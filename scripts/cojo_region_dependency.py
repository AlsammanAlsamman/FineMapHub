#!/usr/bin/env python3
"""Run targeted COJO conditioning experiments to assess regional independence.

This script reuses helper routines from run_cojo_iterative.py to build a region-specific
.ma file and to execute GCTA-COJO while incrementally conditioning on SNPs from one
region and tracking the impact on another region.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from pandas.errors import ParserError

import run_cojo_iterative as cojo  # pylint: disable=import-error


@dataclass
class Region:
    label: str
    chromosome: int
    start: int
    end: int

    def contains(self, chrom: float, position: float) -> bool:
        if math.isnan(chrom) or math.isnan(position):
            return False
        return int(chrom) == self.chromosome and self.start <= position <= self.end


def parse_region(spec: str) -> Region:
    """Parse region specification label:chr:start:end"""
    parts = spec.split(":")
    if len(parts) != 4:
        raise ValueError(f"Region specification '{spec}' must be label:chr:start:end")

    label = parts[0]
    try:
        chromosome = int(parts[1])
    except ValueError as exc:  # pragma: no cover - defensive path
        raise ValueError(f"Invalid chromosome in '{spec}'") from exc

    def parse_coord(value: str) -> int:
        clean = value.replace(",", "").replace("_", "")
        return int(float(clean))

    start = parse_coord(parts[2])
    end = parse_coord(parts[3])
    if end < start:
        start, end = end, start
    return Region(label=label, chromosome=chromosome, start=start, end=end)


def collect_region_snps(gwas_df: pd.DataFrame, region: Region) -> pd.DataFrame:
    """Return GWAS rows within region sorted by ascending p-value."""
    gwas_df = gwas_df.copy()
    gwas_df["CHR_numeric"] = pd.to_numeric(gwas_df["CHR"], errors="coerce")
    gwas_df["POS_numeric"] = pd.to_numeric(gwas_df["POS"], errors="coerce")
    gwas_df["P_numeric"] = pd.to_numeric(gwas_df["P"], errors="coerce")

    subset = gwas_df[
        (gwas_df["CHR_numeric"] == region.chromosome)
        & (gwas_df["POS_numeric"].between(region.start, region.end))
        & (gwas_df["P_numeric"].between(0, 1, inclusive="neither"))
    ].copy()

    if subset.empty:
        return subset

    subset = subset.dropna(subset=["SNPID"])
    subset = subset.sort_values("P_numeric", ascending=True)
    return subset


def top_signal_from_df(
    df: pd.DataFrame,
    region: Region,
    p_column: str,
    chr_column: str = "Chr",
    pos_column: str = "bp",
) -> Optional[Dict[str, float]]:
    """Return the top signal (highest -log10(p)) within a region."""
    if p_column not in df.columns:
        fallback = "p" if p_column != "p" else "pC"
        if fallback in df.columns:
            p_column = fallback
        else:
            return None

    df = df.copy()
    df["chr_numeric"] = pd.to_numeric(df[chr_column], errors="coerce")
    df["pos_numeric"] = pd.to_numeric(df[pos_column], errors="coerce")
    df["p_numeric"] = pd.to_numeric(df[p_column], errors="coerce")

    subset = df[
        (df["chr_numeric"] == region.chromosome)
        & (df["pos_numeric"].between(region.start, region.end))
        & (df["p_numeric"].between(0, 1, inclusive="neither"))
    ].copy()

    if subset.empty:
        return None

    subset.loc[:, "p_numeric"] = subset["p_numeric"].clip(lower=1e-300)
    subset["neg_log_p"] = -1.0 * subset["p_numeric"].map(lambda x: math.log10(x))
    subset = subset.sort_values("neg_log_p", ascending=False)

    top_row = subset.iloc[0]
    return {
        "snp": top_row.get("SNP") or top_row.get("SNPID"),
        "neg_log_p": top_row["neg_log_p"],
        "p": top_row["p_numeric"],
    }


def baseline_signal_from_gwas(gwas_df: pd.DataFrame, region: Region) -> Optional[Dict[str, float]]:
    """Return baseline top signal using raw GWAS p-values."""
    df = gwas_df.copy()
    df["CHR_numeric"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["POS_numeric"] = pd.to_numeric(df["POS"], errors="coerce")
    df["P_numeric"] = pd.to_numeric(df["P"], errors="coerce")

    subset = df[
        (df["CHR_numeric"] == region.chromosome)
        & (df["POS_numeric"].between(region.start, region.end))
        & (df["P_numeric"].between(0, 1, inclusive="neither"))
    ].copy()

    if subset.empty:
        return None

    subset.loc[:, "P_numeric"] = subset["P_numeric"].clip(lower=1e-300)
    subset["neg_log_p"] = -1.0 * subset["P_numeric"].map(lambda x: math.log10(x))
    subset = subset.sort_values("neg_log_p", ascending=False)
    top_row = subset.iloc[0]
    return {
        "snp": top_row["SNPID"],
        "neg_log_p": top_row["neg_log_p"],
        "p": top_row["P_numeric"],
    }


def write_snplist(path: Path, snps: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for snp in snps:
            handle.write(f"{snp}\n")


def ensure_executable(path: Optional[str], candidates: Iterable[str], label: str) -> str:
    if path:
        return path
    found = cojo.find_executable(list(candidates))
    if not found:
        raise RuntimeError(f"{label} executable not found; tried {', '.join(candidates)}")
    return found


def read_cojo_table(file_path: Path) -> pd.DataFrame:
    """Read a GCTA COJO output table allowing space or tab delimiters."""
    try:
        return pd.read_csv(file_path, sep="\t")
    except ParserError:
        return pd.read_csv(file_path, delim_whitespace=True, engine="python")


def run_cojo_selection(
    *,
    region: Region,
    region_snps_df: pd.DataFrame,
    significance_threshold: float,
    selection_dir: Path,
    bfile_prefix: Path,
    ma_file: Path,
    chromosome: int,
    gcta_path: str,
    cojo_collinear: Optional[float],
    log_file: Path,
) -> List[str]:
    """Run GCTA --cojo-slct to preselect independent SNPs within a region."""
    if region_snps_df.empty:
        return []

    if "P_numeric" not in region_snps_df.columns:
        if "P" in region_snps_df.columns:
            region_snps_df = region_snps_df.copy()
            region_snps_df["P_numeric"] = pd.to_numeric(region_snps_df["P"], errors="coerce")
        else:
            return []

    candidate_df = region_snps_df[region_snps_df["P_numeric"] < significance_threshold]
    if candidate_df.empty:
        return []

    selection_dir.mkdir(parents=True, exist_ok=True)
    candidate_snplist = selection_dir / f"{region.label}_candidates.snplist"
    write_snplist(candidate_snplist, candidate_df["SNPID"].astype(str).tolist())

    output_prefix = selection_dir / f"{region.label}_cojo_slct"
    cmd = [
        gcta_path,
        "--bfile",
        str(bfile_prefix),
        "--chr",
        str(chromosome),
        "--maf",
        "0.01",
        "--cojo-file",
        str(ma_file),
        "--cojo-slct",
        "--cojo-p",
        str(significance_threshold),
        "--cojo-snp",
        str(candidate_snplist),
        "--out",
        str(output_prefix),
    ]
    if cojo_collinear is not None and cojo_collinear > 0:
        cmd.extend(["--cojo-collinear", str(cojo_collinear)])

    print(
        f"Running GCTA-COJO selection: {' '.join(cmd)}",
        file=sys.stderr,
    )
    with log_file.open("a", encoding="utf-8") as log_handle:
        subprocess.run(cmd, stdout=log_handle, stderr=log_handle, check=False)

    jma_file = output_prefix.with_suffix(".jma.cojo")
    if not jma_file.exists():
        print(
            f"COJO selection for {region.label} produced no .jma.cojo file; falling back to iterative selection.",
            file=sys.stderr,
        )
        return []

    jma_df = read_cojo_table(jma_file)
    if jma_df.empty or "SNP" not in jma_df.columns:
        print(
            f"COJO selection for {region.label} returned no SNPs; falling back to iterative selection.",
            file=sys.stderr,
        )
        return []

    order_col = "pJ" if "pJ" in jma_df.columns else ("p" if "p" in jma_df.columns else None)
    if order_col:
        jma_df = jma_df.sort_values(order_col, ascending=True)

    return jma_df["SNP"].astype(str).tolist()


def run_conditioning_sequence(
    *,
    region: Region,
    other_region: Region,
    region_snps: Sequence[Dict[str, str]],
    gwas_df: pd.DataFrame,
    ma_file: Path,
    bfile_prefix: Path,
    output_dir: Path,
    chromosome: int,
    gcta_path: str,
    log_file: Path,
    cojo_collinear: Optional[float] = 0.99,
    significance_threshold: float = 5e-8,
    preselected_snps: Optional[Sequence[str]] = None,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    """Condition on SNPs sequentially and track impact on the other region."""
    results: List[Dict[str, object]] = []

    # Baseline entry
    # Record the unconditioned reference point before we begin iterative COJO runs.
    baseline = baseline_signal_from_gwas(gwas_df, other_region)
    if baseline:
        results.append(
            {
                "mode": "baseline",
                "conditioning_region": region.label,
                "other_region": other_region.label,
                "step": 0,
                "total_conditioned": 0,
                "conditioned_snp": "",
                "top_other_snp": baseline["snp"],
                "top_other_p": baseline["p"],
                "top_other_neg_logp": baseline["neg_log_p"],
                "output_prefix": "baseline",
                "attempt_index": 0,
            }
        )

    region_snps_df = pd.DataFrame(region_snps)
    conditioned_snps: List[str] = []
    conditioned_set = set()
    skipped_snps: set[str] = set()
    attempt_index = 0
    last_cma_df: Optional[pd.DataFrame] = None
    last_top_signal: Optional[Dict[str, float]] = None
    final_source = "baseline"
    pending_snps: Optional[List[str]] = None

    if preselected_snps is not None:
        pending_snps = [str(snp) for snp in preselected_snps if snp]

    if not region_snps_df.empty:
        if "P_numeric" not in region_snps_df.columns:
            if "P" in region_snps_df.columns:
                region_snps_df["P_numeric"] = pd.to_numeric(region_snps_df["P"], errors="coerce")
            else:
                region_snps_df["P_numeric"] = float("nan")

    def pick_from_gwas() -> Optional[str]:
        if pending_snps is not None:
            return None
        if region_snps_df.empty:
            return None
        subset = region_snps_df[
            (region_snps_df["P_numeric"] < significance_threshold)
            & (~region_snps_df["SNPID"].isin(conditioned_set))
            & (~region_snps_df["SNPID"].isin(skipped_snps))
        ]
        if subset.empty:
            return None
        best_row = subset.sort_values("P_numeric", ascending=True).iloc[0]
        return str(best_row["SNPID"])

    def pick_from_cma(cma_df: Optional[pd.DataFrame]) -> Optional[str]:
        if pending_snps is not None:
            return None
        if cma_df is None or cma_df.empty:
            return None

        snp_cols = [col for col in ("SNP", "SNPID", "rsid", "RSID", "ID") if col in cma_df.columns]
        p_col = "pC" if "pC" in cma_df.columns else ("p" if "p" in cma_df.columns else None)
        chr_col = next((col for col in ("Chr", "chr", "CHR") if col in cma_df.columns), None)
        pos_col = next((col for col in ("bp", "BP", "Pos", "POS") if col in cma_df.columns), None)

        if not snp_cols or p_col is None or chr_col is None or pos_col is None:
            return None

        work_df = cma_df.copy()
        work_df["chr_numeric"] = pd.to_numeric(work_df[chr_col], errors="coerce")
        work_df["pos_numeric"] = pd.to_numeric(work_df[pos_col], errors="coerce")
        work_df["p_numeric"] = pd.to_numeric(work_df[p_col], errors="coerce")

        snp_col = snp_cols[0]
        work_df["SNPID"] = work_df[snp_col].astype(str)

        subset = work_df[
            (work_df["chr_numeric"] == region.chromosome)
            & (work_df["pos_numeric"].between(region.start, region.end))
            & (work_df["p_numeric"].between(0, 1, inclusive="neither"))
            & (work_df["p_numeric"] < significance_threshold)
            & (~work_df["SNPID"].isin(conditioned_set))
            & (~work_df["SNPID"].isin(skipped_snps))
        ]

        if subset.empty:
            return None

        best_row = subset.sort_values("p_numeric", ascending=True).iloc[0]
        return str(best_row["SNPID"])

    while True:
        if pending_snps is not None:
            if not pending_snps:
                break
            next_snp = pending_snps.pop(0)
            if next_snp in conditioned_set or next_snp in skipped_snps:
                continue
        else:
            if not conditioned_snps:
                next_snp = pick_from_gwas()
            else:
                next_snp = pick_from_cma(last_cma_df)

        if next_snp is None:
            break

        attempt_index += 1
        candidate_snps = conditioned_snps + [next_snp]
        cond_file = output_dir / f"{region.label}_step_{attempt_index}.snplist"
        write_snplist(cond_file, candidate_snps)

        output_prefix = output_dir / f"{region.label}_step_{attempt_index}"
        success = cojo.run_gcta_cojo(
            str(bfile_prefix),
            str(ma_file),
            str(cond_file),
            str(output_prefix),
            chromosome,
            gcta_path,
            str(log_file),
            cojo_collinear=cojo_collinear,
        )
        if not success:
            raise RuntimeError(f"GCTA failed for {output_prefix}")

        cma_file = output_prefix.with_suffix(".cma.cojo")
        if not cma_file.exists():
            skipped_snps.add(next_snp)
            results.append(
                {
                    "mode": "collinearity_skip",
                    "conditioning_region": region.label,
                    "other_region": other_region.label,
                    "step": len(conditioned_snps),
                    "total_conditioned": len(conditioned_snps),
                    "conditioned_snp": next_snp,
                    "top_other_snp": "NA",
                    "top_other_p": float("nan"),
                    "top_other_neg_logp": float("nan"),
                    "output_prefix": str(output_prefix),
                    "attempt_index": attempt_index,
                }
            )
            print(
                f"Skipping {next_snp} in {region.label}: GCTA reported collinearity; keeping previous SNP set.",
                file=sys.stderr,
            )
            continue

        cma_df = read_cojo_table(cma_file)
        top_signal = top_signal_from_df(cma_df, other_region, p_column="pC")
        if top_signal is None:
            top_signal = {
                "snp": "NA",
                "neg_log_p": float("nan"),
                "p": float("nan"),
            }

        conditioned_snps = candidate_snps
        current_step = len(conditioned_snps)
        results.append(
            {
                "mode": "sequential",
                "conditioning_region": region.label,
                "other_region": other_region.label,
                "step": current_step,
                "total_conditioned": current_step,
                "conditioned_snp": next_snp,
                "top_other_snp": top_signal["snp"],
                "top_other_p": top_signal["p"],
                "top_other_neg_logp": top_signal["neg_log_p"],
                "output_prefix": str(output_prefix),
                "attempt_index": attempt_index,
            }
        )

        conditioned_set.add(next_snp)
        last_cma_df = cma_df
        last_top_signal = top_signal
        final_source = "sequential"

    if conditioned_snps:
        cond_file = output_dir / f"{region.label}_all.snplist"
        write_snplist(cond_file, conditioned_snps)
        output_prefix = output_dir / f"{region.label}_all"
        success = cojo.run_gcta_cojo(
            str(bfile_prefix),
            str(ma_file),
            str(cond_file),
            str(output_prefix),
            chromosome,
            gcta_path,
            str(log_file),
            cojo_collinear=cojo_collinear,
        )
        if not success:
            raise RuntimeError(f"GCTA failed for joint conditioning of {region.label}")

        cma_file = output_prefix.with_suffix(".cma.cojo")
        if not cma_file.exists():
            raise RuntimeError(f"Expected result missing: {cma_file}")

        cma_df = read_cojo_table(cma_file)
        top_signal = top_signal_from_df(cma_df, other_region, p_column="pC")
        if top_signal is None:
            top_signal = {
                "snp": "NA",
                "neg_log_p": float("nan"),
                "p": float("nan"),
            }
        results.append(
            {
                "mode": "all_at_once",
                "conditioning_region": region.label,
                "other_region": other_region.label,
                "step": len(conditioned_snps),
                "total_conditioned": len(conditioned_snps),
                "conditioned_snp": "ALL",
                "top_other_snp": top_signal["snp"],
                "top_other_p": top_signal["p"],
                "top_other_neg_logp": top_signal["neg_log_p"],
                "output_prefix": str(output_prefix),
                "attempt_index": attempt_index + 1,
            }
        )
        last_top_signal = top_signal
        final_source = "all_at_once"

    baseline_p = baseline["p"] if baseline else float("nan")
    baseline_neg = baseline["neg_log_p"] if baseline else float("nan")

    if last_top_signal is not None:
        final_p = last_top_signal.get("p", float("nan"))
        final_neg = last_top_signal.get("neg_log_p", float("nan"))
    else:
        final_p = baseline_p
        final_neg = baseline_neg

    is_independent = True
    if final_p is not None and not math.isnan(final_p):
        if final_p < significance_threshold:
            is_independent = False
    elif baseline_p is not None and not math.isnan(baseline_p):
        if baseline_p < significance_threshold:
            is_independent = False

    def safe_float(value: Optional[float]) -> Optional[float]:
        if value is None:
            return None
        if isinstance(value, float) and math.isnan(value):
            return None
        return value

    direction_summary = {
        "conditioning_region": region.label,
        "other_region": other_region.label,
        "baseline_p": safe_float(baseline_p),
        "baseline_neg_logp": safe_float(baseline_neg),
        "final_p": safe_float(final_p),
        "final_neg_logp": safe_float(final_neg),
        "conditioned_snp_count": len(conditioned_snps),
        "final_source": final_source,
        "independence_status": "independent" if is_independent else "dependent",
    }

    return results, direction_summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="COJO regional dependency experiment")
    parser.add_argument("--gwas-file", required=True)
    parser.add_argument("--bfile-prefix", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample-size", type=int, required=True)
    parser.add_argument(
        "--region",
        action="append",
        dest="regions",
        required=True,
        help="Region specification label:chr:start:end (provide twice)",
    )
    parser.add_argument("--chromosome", type=int, default=None)
    parser.add_argument("--gcta-path", default=None)
    parser.add_argument("--plink-path", default=None)
    parser.add_argument("--significance-threshold", type=float, default=5e-8)
    parser.add_argument("--freq-file", default=None)
    parser.add_argument("--generate-freq", action="store_true")
    parser.add_argument(
        "--cojo-collinear",
        type=float,
        default=0.99,
        help="Threshold to pass to GCTA via --cojo-collinear (set <=0 to disable)",
    )
    parser.add_argument(
        "--use-cojo-slct",
        action="store_true",
        help="Use GCTA --cojo-slct to pre-select independent SNPs within each region",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="Size of SNP batches to condition in a single run (use with --chunk-index)",
    )
    parser.add_argument(
        "--chunk-index",
        type=int,
        default=None,
        help="Zero-based index of the SNP batch to process (use with --chunk-size)",
    )

    args = parser.parse_args(argv)

    if len(args.regions) != 2:
        parser.error("Exactly two --region specifications are required")

    if args.chunk_size is not None and args.chunk_size <= 0:
        parser.error("--chunk-size must be a positive integer")

    if args.chunk_index is not None:
        if args.chunk_index < 0:
            parser.error("--chunk-index must be >= 0")
        if args.chunk_size is None:
            parser.error("--chunk-size is required when --chunk-index is provided")

    regions = [parse_region(r) for r in args.regions]
    region_map = {region.label: region for region in regions}

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gwas_file = Path(args.gwas_file)
    bfile_prefix = Path(args.bfile_prefix)

    if not cojo.validate_inputs(str(gwas_file), str(bfile_prefix), str(output_dir)):
        return 1

    chromosome = args.chromosome
    if chromosome is None:
        chromosome = cojo.detect_chromosome(str(gwas_file))
        if chromosome is None:
            sys.stderr.write("ERROR: Could not determine chromosome. Use --chromosome.\n")
            return 1

    gcta_path = ensure_executable(args.gcta_path, ["gcta64", "gcta", "gcta-1.94.1", "gcta_1.94.1"], "GCTA")

    cojo_collinear = None if args.cojo_collinear is not None and args.cojo_collinear <= 0 else args.cojo_collinear

    freq_file = args.freq_file
    if args.generate_freq:
        plink_path = ensure_executable(args.plink_path, ["plink", "plink2", "plink1.9"], "PLINK")
        freq_file = cojo.generate_freq_file(str(bfile_prefix), str(output_dir), plink_path)

    # Prepare .ma file
    ma_file = Path(
        cojo.prepare_ma_file(
            str(gwas_file),
            str(bfile_prefix),
            str(output_dir),
            args.sample_size,
            freq_file,
        )
    )

    gwas_df = pd.read_csv(gwas_file, sep="\t")
    region_snps: Dict[str, pd.DataFrame] = {}
    for region in regions:
        subset = collect_region_snps(gwas_df, region)
        if subset.empty:
            raise RuntimeError(f"Region {region.label} contains no SNPs in GWAS file")
        region_snps[region.label] = subset

    log_file = output_dir / "cojo_region_dependency.log"
    summary_rows: List[Dict[str, object]] = []
    conditioned_summary: Dict[str, List[str]] = {}
    selection_summary: Dict[str, List[str]] = {}
    direction_summaries: List[Dict[str, object]] = []

    def slice_chunk(seq: Optional[Sequence[str]]) -> List[str]:
        if not seq:
            return []
        if args.chunk_size is None or args.chunk_index is None:
            return list(seq)
        start = args.chunk_index * args.chunk_size
        end = start + args.chunk_size
        return list(seq[start:end])

    preselected_by_region: Dict[str, Optional[List[str]]] = {}
    chunk_queue_summary: Dict[str, List[str]] = {}
    selection_root = output_dir / "selection"
    for region in regions:
        selection_list: List[str] = []
        if args.use_cojo_slct:
            selection_list = run_cojo_selection(
                region=region,
                region_snps_df=region_snps[region.label],
                significance_threshold=args.significance_threshold,
                selection_dir=selection_root / region.label,
                bfile_prefix=bfile_prefix,
                ma_file=ma_file,
                chromosome=chromosome,
                gcta_path=gcta_path,
                cojo_collinear=cojo_collinear,
                log_file=log_file,
            )
        selection_summary[region.label] = selection_list

        base_queue: Optional[List[str]] = selection_list if selection_list else None
        if args.chunk_size is not None and base_queue is None:
            region_df = region_snps[region.label].copy()
            region_df = region_df.sort_values("P_numeric", ascending=True)
            base_queue = region_df["SNPID"].astype(str).tolist()

        if args.chunk_size is not None:
            chunk_queue = slice_chunk(base_queue)
            preselected_by_region[region.label] = chunk_queue
            chunk_queue_summary[region.label] = chunk_queue
        else:
            preselected_by_region[region.label] = base_queue
            chunk_queue_summary[region.label] = base_queue if base_queue else []
    if not args.use_cojo_slct:
        for region in regions:
            selection_summary.setdefault(region.label, [])

    order = [(regions[0], regions[1]), (regions[1], regions[0])]
    for conditioning_region, other_region in order:
        cond_dir = output_dir / f"conditioning_{conditioning_region.label}"
        cond_dir.mkdir(exist_ok=True)

        sequential_results, direction_assessment = run_conditioning_sequence(
            region=conditioning_region,
            other_region=other_region,
            region_snps=region_snps[conditioning_region.label].to_dict("records"),
            gwas_df=gwas_df,
            ma_file=ma_file,
            bfile_prefix=bfile_prefix,
            output_dir=cond_dir,
            chromosome=chromosome,
            gcta_path=gcta_path,
            log_file=log_file,
            cojo_collinear=cojo_collinear,
            significance_threshold=args.significance_threshold,
            preselected_snps=preselected_by_region.get(conditioning_region.label),
        )
        summary_rows.extend(sequential_results)
        direction_summaries.append(direction_assessment)
        conditioned_summary[conditioning_region.label] = [
            row["conditioned_snp"]
            for row in sequential_results
            if row.get("mode") == "sequential"
        ]

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / "region_dependency_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    assessment_df = pd.DataFrame(direction_summaries)
    assessment_path = output_dir / "region_dependency_assessment.tsv"
    assessment_df.to_csv(assessment_path, sep="\t", index=False)

    metadata = {
        "gwas_file": str(gwas_file),
        "bfile_prefix": str(bfile_prefix),
        "sample_size": args.sample_size,
        "chromosome": chromosome,
        "regions": {label: region_map[label].__dict__ for label in region_map},
        "summary_file": str(summary_path),
        "assessment_file": str(assessment_path),
        "ma_file": str(ma_file),
        "significance_threshold": args.significance_threshold,
        "freq_file": freq_file,
        "cojo_collinear": cojo_collinear,
        "conditioned_snps_by_region": conditioned_summary,
        "cojo_slct_selected": selection_summary,
        "use_cojo_slct": args.use_cojo_slct,
        "chunk_size": args.chunk_size,
        "chunk_index": args.chunk_index,
        "chunk_snp_queue": chunk_queue_summary,
        "independence_assessment": direction_summaries,
    }
    metadata_path = output_dir / "region_dependency_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    print(f"Results written to {summary_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
