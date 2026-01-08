#!/usr/bin/env python3
"""
Generate Target Analysis Excel Report
Creates a comprehensive Excel report for a single target analysis with COJO results
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Generate Excel report for target analysis COJO results"
    )
    
    parser.add_argument(
        "--target-analysis",
        required=True,
        help="Target analysis name"
    )
    parser.add_argument(
        "--cojo-dir",
        required=True,
        help="Directory containing COJO results"
    )
    parser.add_argument(
        "--loci-dir",
        required=True,
        help="Directory containing loci data"
    )
    parser.add_argument(
        "--loci-info-file",
        required=False,
        default="inputs/LOCI_info.txt",
        help="Path to loci information file (optional)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output Excel file path"
    )
    
    return parser.parse_args()


def load_loci_info(loci_info_file: str) -> Dict[str, Dict]:
    """Load loci information from file if available"""
    loci_info = {}
    
    loci_path = Path(loci_info_file)
    if not loci_path.exists():
        print(f"Note: Loci info file not found: {loci_info_file}")
        return loci_info
    
    try:
        df = pd.read_csv(loci_path, sep='\t')
        print(f"Loci info file columns: {list(df.columns)}")
        
        for _, row in df.iterrows():
            # Try multiple possible column names for locus identifier
            locus_id = None
            for col in ['LOCUS_ID', 'Locus', 'locus', 'LOCUS', 'locus_id', 'Locus_ID']:
                if col in row.index and pd.notna(row[col]):
                    locus_id = str(row[col])
                    break
            
            if locus_id:
                # Store with both original ID and LOC_ prefixed version
                loci_info[locus_id] = row.to_dict()
                if not locus_id.startswith('LOC_'):
                    loci_info[f'LOC_{locus_id}'] = row.to_dict()
        
        print(f"Loaded information for {len(loci_info)} loci")
        if loci_info:
            sample_key = list(loci_info.keys())[0]
            print(f"Sample loci info keys: {sample_key} -> {list(loci_info[sample_key].keys())}")
    except Exception as e:
        print(f"Warning: Could not load loci info: {e}")
    
    return loci_info


def collect_loci_results(cojo_dir: str, loci_dir: str) -> Tuple[List[Dict], List[str]]:
    """Collect COJO results for all loci in the target analysis"""
    cojo_path = Path(cojo_dir)
    loci_path = Path(loci_dir)
    
    results = []
    loci_names = []
    
    # Find all LOC_* directories
    loci_dirs = sorted([d for d in cojo_path.iterdir() if d.is_dir() and d.name.startswith("LOC_")])
    
    print(f"Found {len(loci_dirs)} loci directories")
    
    for locus_dir in loci_dirs:
        locus_name = locus_dir.name
        loci_names.append(locus_name)
        
        cojo_results_dir = locus_dir / "cojo_results"
        if not cojo_results_dir.exists():
            print(f"  Skipping {locus_name}: no cojo_results directory")
            continue
        
        # Read SNP list
        snp_list_file = cojo_results_dir / "final_cond_snplist.txt"
        snps = []
        if snp_list_file.exists():
            with open(snp_list_file, 'r') as f:
                snps = [line.strip() for line in f if line.strip()]
        
        # Read iteration summary
        summary_file = cojo_results_dir / "iteration_summary.json"
        summary = {}
        if summary_file.exists():
            try:
                with open(summary_file, 'r') as f:
                    summary = json.load(f)
            except Exception as e:
                print(f"  Warning: Could not read summary for {locus_name}: {e}")
        
        # Read GWAS data for leading SNP p-value and chromosome/position
        gwas_file = loci_path / locus_name / f"{locus_name}.loci.tsv"
        leading_snp = None
        leading_pval = None
        chromosome = None
        position = None
        
        if gwas_file.exists():
            try:
                # Try reading with automatic compression detection
                gwas_df = pd.read_csv(gwas_file, sep='\t', compression='infer', low_memory=False)
                
                if len(gwas_df) > 0:
                    # Get leading SNP (minimum p-value)
                    if 'P' in gwas_df.columns:
                        min_idx = gwas_df['P'].idxmin()
                        leading_snp = gwas_df.loc[min_idx, 'SNPID'] if 'SNPID' in gwas_df.columns else None
                        leading_pval = gwas_df.loc[min_idx, 'P']
                    
                    # Extract chromosome and position from first row
                    chr_col = None,
            'chromosome': chromosome,
            'position': position
                    pos_col = None
                    for col in ['CHR', 'Chr', 'chr', 'Chromosome', 'chromosome', 'CHROM']:
                        if col in gwas_df.columns:
                            chr_col = col
                            break
                    for col in ['POS', 'Pos', 'pos', 'BP', 'Position', 'position']:
                        if col in gwas_df.columns:
                            pos_col = col
                            break
                    
                    if chr_col and pd.notna(gwas_df[chr_col].iloc[0]):
                        chromosome = - try multiple lookup keys
        info = loci_info.get(locus_name, loci_info.get(locus_num, {}))
        
        # Try multiple possible column names for chromosome, position, gene
        chromosome = None
        position = None
        gene = None
        
        if info:
            # Try chromosome
            for col in ['CHR', 'Chr', 'chr', 'Chromosome', 'chromosome', 'CHROM', 'chrom']:
                if col in info and pd.notna(info[col]):
                    chromosome = info[col]
                    break
            # Try position
            for col in ['BP', 'bp', 'POS', 'Pos', 'pos', 'Position', 'position']:
                if col in info and pd.notna(info[col]):
                    position = info[col]
                    break
            # Try gene
            for col in ['Gene', 'gene', 'GENE', 'Gene_Name', 'gene_name']:
                if col in info and pd.notna(info[col]):
                    gene = info[col]
                    break
        
        # Fall back to data extracted from GWAS file
        if chromosome is None and result.get('chromosome') is not None:
            chromosome = result['chromosome']
        if position is None and result.get('position') is not None:
            position = result['position']
        
        row = {
            'Locus': locus_name,
            'Chromosome': chromosome if chromosome is not None else '',
            'Position': position if position is not None else '',
            'Gene': gene if gene is not None else ''
                                position = gwas_df[pos_col].iloc[0]
                        else:
                            position = gwas_df[pos_col].iloc[0]
            except Exception as e:
                print(f"  Warning: Could not read GWAS data for {locus_name}: {e}")
        
        results.append({
            'locus_name': locus_name,
            'snps': snps,
            'count': len(snps),
            'summary': summary,
            'leading_snp': leading_snp,
            'leading_pval': leading_pval
        })
    
    return results, loci_names


def create_loci_summary_sheet(results: List[Dict], loci_info: Dict[str, Dict]) -> pd.DataFrame:
    """Create summary sheet with one row per locus"""
    rows = []
    
    for result in results:
        locus_name = result['locus_name']
        locus_num = locus_name.replace('LOC_', '')
        
        # Get loci info if available
        info = loci_info.get(locus_name, loci_info.get(locus_num, {}))
        
    # Cache GWAS dataframes to avoid re-reading
    gwas_cache = {}
    
    for result in results:
        locus_name = result['locus_name']
        
        for snp in result['snps']:
            # Try to get SNP details from GWAS file
            gwas_file = loci_path / locus_name / f"{locus_name}.loci.tsv"
            snp_details = {'SNPID': snp}
            
            if gwas_file.exists():
                try:
                    # Use cached dataframe if available
                    if locus_name not in gwas_cache:
                        gwas_cache[locus_name] = pd.read_csv(gwas_file, sep='\t', compression='infer', low_memory=False)
                    
                    gwas_df = gwas_cache[locus_name]
                    snp_row = gwas_df[gwas_df['SNPID'] == snp]
                    
                    if len(snp_row) > 0:
                        snp_details.update(snp_row.iloc[0].to_dict())
                except Exception as e:
                    # Only show warning once per locus, not per SNP
                    if locus_name not in gwas_cache:
                        print(f"  Warning: Could not read GWAS file for {locus_name}: {e}")
                        gwas_cache[locus_name] = None
            
            # Try multiple column names for each field
            row = {
                'Locus': locus_name,
                'SNP': snp,
                'Chromosome': snp_details.get('CHR', snp_details.get('Chr', snp_details.get('chr', ''))),
                'Position': snp_details.get('POS', snp_details.get('Pos', snp_details.get('pos', snp_details.get('BP', '')))),
                'Reference_Allele': snp_details.get('A1', snp_details.get('a1', snp_details.get('Allele1', ''))),
                'Alternative_Allele': snp_details.get('A2', snp_details.get('a2', snp_details.get('Allele2', ''))),
                'P_Value': snp_details.get('P', snp_details.get('p', snp_details.get('PVAL', ''))),
                'Beta': snp_details.get('BETA', snp_details.get('Beta', snp_details.get('beta', snp_details.get('B', '')))),
                'SE': snp_details.get('SE', snp_details.get('se', snp_details.get('StdErr', ''))),
                'Frequency': snp_details.get('FREQ', snp_details.get('Freq', snp_details.get('freq', snp_details.get('MAF', '')))
    loci_path = Path(loci_dir)
    
    for result in results:
        locus_name = result['locus_name']
        
        for snp in result['snps']:
            # Try to get SNP details from GWAS file
            gwas_file = loci_path / locus_name / f"{locus_name}.loci.tsv"
            snp_details = {'SNPID': snp}
            
            if gwas_file.exists():
                try:
                    gwas_df = pd.read_csv(gwas_file, sep='\t')
                    snp_row = gwas_df[gwas_df['SNPID'] == snp]
                    if len(snp_row) > 0:
                        snp_details.update(snp_row.iloc[0].to_dict())
                except Exception as e:
                    print(f"  Warning: Could not read SNP details for {snp} in {locus_name}: {e}")
            
            row = {
                'Locus': locus_name,
                'SNP': snp,
                'Chromosome': snp_details.get('CHR', ''),
                'Position': snp_details.get('POS', ''),
                'Reference_Allele': snp_details.get('A1', ''),
                'Alternative_Allele': snp_details.get('A2', ''),
                'P_Value': snp_details.get('P', ''),
                'Beta': snp_details.get('BETA', ''),
                'SE': snp_details.get('SE', ''),
                'Frequency': snp_details.get('FREQ', '')
            }
            rows.append(row)
    
    return pd.DataFrame(rows)


def create_summary_statistics(results: List[Dict], target_analysis: str) -> pd.DataFrame:
    """Create overall summary statistics"""
    total_loci = len(results)
    total_snps = sum(r['count'] for r in results)
    loci_with_snps = sum(1 for r in results if r['count'] > 0)
    avg_snps_per_locus = total_snps / total_loci if total_loci > 0 else 0
    
    stats = {
        'Metric': [
            'Target Analysis',
            'Total Loci',
            'Loci with Conditioning SNPs',
            'Total Conditioning SNPs',
            'Average SNPs per Locus',
            'Min SNPs in a Locus',
            'Max SNPs in a Locus'
        ],
        'Value': [
            target_analysis,
            total_loci,
            loci_with_snps,
            total_snps,
            f"{avg_snps_per_locus:.2f}",
            min(r['count'] for r in results) if results else 0,
            max(r['count'] for r in results) if results else 0
        ]
    }
    
    return pd.DataFrame(stats)


def write_excel_report(
    output_file: str,
    summary_stats: pd.DataFrame,
    loci_summary: pd.DataFrame,
    snp_details: pd.DataFrame
):
    """Write all data to Excel file with multiple sheets"""
    print(f"Writing Excel report to: {output_file}")
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Sheet 1: Summary Statistics
        summary_stats.to_excel(writer, sheet_name='Summary', index=False)
        
        # Sheet 2: Loci Summary
        loci_summary.to_excel(writer, sheet_name='Loci_Summary', index=False)
        
        # Sheet 3: SNP Details
        if len(snp_details) > 0:
            snp_details.to_excel(writer, sheet_name='SNP_Details', index=False)
        
        # Auto-adjust column widths
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]
            for column in worksheet.columns:
                max_length = 0
                column_letter = column[0].column_letter
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = min(max_length + 2, 50)
                worksheet.column_dimensions[column_letter].width = adjusted_width
    
    print("✓ Excel report created successfully")


def main():
    args = parse_arguments()
    
    print("=" * 60)
    print("TARGET ANALYSIS EXCEL REPORT GENERATOR")
    print("=" * 60)
    print(f"Target analysis: {args.target_analysis}")
    print(f"COJO directory: {args.cojo_dir}")
    print(f"Loci directory: {args.loci_dir}")
    print(f"Output file: {args.output}")
    print("=" * 60)
    
    # Load loci information
    loci_info = load_loci_info(args.loci_info_file)
    
    # Collect results
    print("\nCollecting COJO results...")
    results, loci_names = collect_loci_results(args.cojo_dir, args.loci_dir)
    
    if len(results) == 0:
        print("ERROR: No loci results found")
        sys.exit(1)
    
    print(f"✓ Collected results for {len(results)} loci")
    
    # Create sheets
    print("\nGenerating report sheets...")
    summary_stats = create_summary_statistics(results, args.target_analysis)
    loci_summary = create_loci_summary_sheet(results, loci_info)
    snp_details = create_snp_details_sheet(results, args.loci_dir)
    
    print(f"  Summary statistics: {len(summary_stats)} rows")
    print(f"  Loci summary: {len(loci_summary)} rows")
    print(f"  SNP details: {len(snp_details)} rows")
    
    # Write Excel file
    print("\nWriting Excel file...")
    write_excel_report(args.output, summary_stats, loci_summary, snp_details)
    
    print("\n" + "=" * 60)
    print("REPORT GENERATION COMPLETE")
    print("=" * 60)
    print(f"Output: {args.output}")
    print(f"Total loci: {len(results)}")
    print(f"Total conditioning SNPs: {sum(r['count'] for r in results)}")
    print("=" * 60)


if __name__ == '__main__':
    main()
