#!/usr/bin/env python3
"""
Iterative GCTA-COJO Conditional Analysis - Standalone Script
Performs iterative conditional analysis until no significant SNPs remain
"""

import pandas as pd
import subprocess
import argparse
import sys
import os
import json
from pathlib import Path

def detect_chromosome(gwas_file):
    """Detect chromosome from GWAS file"""
    try:
        gwas_df = pd.read_csv(gwas_file, sep='\t', nrows=100)
        
        # Check for CHR or CHROM column
        chrom_col = None
        for col in ['CHR', 'CHROM', 'chrom', 'chr', 'chromosome', 'Chromosome']:
            if col in gwas_df.columns:
                chrom_col = col
                break
        
        if chrom_col is None:
            print("ERROR: Could not find chromosome column in GWAS file", file=sys.stderr)
            return None
        
        # Get unique chromosomes
        unique_chroms = gwas_df[chrom_col].unique()
        
        if len(unique_chroms) == 0:
            print("ERROR: No chromosome data found", file=sys.stderr)
            return None
        elif len(unique_chroms) > 1:
            print(f"WARNING: Multiple chromosomes found: {unique_chroms}")
            print(f"Using first chromosome: {unique_chroms[0]}")
        
        chrom = int(unique_chroms[0])
        print(f"Auto-detected chromosome: {chrom}")
        return chrom
        
    except Exception as e:
        print(f"ERROR: Failed to detect chromosome: {e}", file=sys.stderr)
        return None

def validate_inputs(gwas_file, bfile_prefix, output_dir):
    """Validate that all required input files exist"""
    errors = []
    
    # Check GWAS file
    if not os.path.exists(gwas_file):
        errors.append(f"GWAS file not found: {gwas_file}")
    
    # Check PLINK files
    for ext in ['.bed', '.bim', '.fam']:
        plink_file = f"{bfile_prefix}{ext}"
        if not os.path.exists(plink_file):
            errors.append(f"PLINK file not found: {plink_file}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return False
    
    return True

def prepare_ma_file(gwas_file, bfile_prefix, output_dir, sample_size, freq_file=None):
    """Convert GWAS data to GCTA-COJO .ma format"""
    print(f"Preparing .ma file from: {gwas_file}")
    
    # Read GWAS data
    gwas_df = pd.read_csv(gwas_file, sep='\t')
    
    # Validate required columns
    required_cols = ['SNPID', 'CHR', 'POS', 'EA', 'NEA', 'BETA', 'SE', 'P', 'N']
    missing_cols = [col for col in required_cols if col not in gwas_df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)
    
    # Read PLINK .bim file to filter SNPs present in reference panel
    bim_file = f"{bfile_prefix}.bim"
    bim_df = pd.read_csv(bim_file, sep='\t', header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    
    print(f"GWAS variants: {len(gwas_df)}")
    print(f"Reference panel variants: {len(bim_df)}")
    
    # Match by SNP ID
    matched_df = gwas_df[gwas_df['SNPID'].isin(bim_df['SNP'])]
    print(f"Matched variants: {len(matched_df)}")
    
    if len(matched_df) == 0:
        print("ERROR: No overlapping SNPs between GWAS and reference panel", file=sys.stderr)
        sys.exit(1)
    
    # Load frequency data if available
    freq_dict = {}
    if freq_file and os.path.exists(freq_file):
        print(f"Loading allele frequencies from: {freq_file}")
        freq_df = pd.read_csv(freq_file, sep=r'\s+')
        for _, row in freq_df.iterrows():
            freq_dict[row['SNP']] = row['MAF']
    
    # Create .ma format: SNP A1 A2 freq b se p N
    ma_data = []
    for _, row in matched_df.iterrows():
        snp = row['SNPID']
        freq = freq_dict.get(snp, 0.5)  # Default to 0.5 if not available
        
        ma_data.append({
            'SNP': snp,
            'A1': row['EA'],
            'A2': row['NEA'],
            'freq': freq,
            'b': row['BETA'],
            'se': row['SE'],
            'p': row['P'],
            'N': sample_size if sample_size else row['N']
        })
    
    ma_df = pd.DataFrame(ma_data)
    
    # Remove rows with missing or invalid data
    ma_df = ma_df.dropna()
    ma_df = ma_df[(ma_df['se'] > 0) & (ma_df['p'] > 0) & (ma_df['p'] <= 1)]
    
    print(f"Final .ma file variants: {len(ma_df)}")
    
    # Save .ma file (no header for GCTA-COJO)
    ma_file = os.path.join(output_dir, "region.ma")
    ma_df.to_csv(ma_file, sep='\t', index=False, header=False)
    
    return ma_file

def create_initial_snplist(gwas_file, output_dir, significance_threshold):
    """Create initial conditioning SNP list with most significant SNP"""
    gwas_df = pd.read_csv(gwas_file, sep='\t')
    
    # Find most significant SNP
    most_sig_snp = gwas_df.loc[gwas_df['P'].idxmin(), 'SNPID']
    
    print(f"Initial conditioning SNP: {most_sig_snp} (p={gwas_df['P'].min():.2e})")
    
    # Save initial SNP list
    snplist_file = os.path.join(output_dir, "initial_cond.snplist")
    with open(snplist_file, 'w') as f:
        f.write(f"{most_sig_snp}\n")
    
    return snplist_file, [most_sig_snp]

def generate_freq_file(bfile_prefix, output_dir, plink_path):
    """Generate allele frequency file using PLINK"""
    freq_output = os.path.join(output_dir, "reference_freq")
    
    cmd = [
        plink_path,
        '--bfile', bfile_prefix,
        '--freq',
        '--out', freq_output
    ]
    
    print(f"Generating allele frequencies: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        freq_file = f"{freq_output}.frq"
        if os.path.exists(freq_file):
            print(f"Frequency file created: {freq_file}")
            return freq_file
        else:
            print("Warning: Frequency file generation failed")
            return None
    except subprocess.CalledProcessError as e:
        print(f"Warning: PLINK frequency generation failed: {e}")
        return None

def run_gcta_cojo(
    bfile_prefix,
    ma_file,
    cond_snplist,
    output_prefix,
    chrom,
    gcta_path,
    log_file,
    *,
    cojo_collinear=0.99,
):
    """Run GCTA-COJO conditional analysis with optional collinearity guard."""
    cmd = [
        gcta_path,
        '--bfile', bfile_prefix,
        '--chr', str(chrom),
        '--maf', '0.01',
        '--cojo-file', ma_file,
        '--cojo-cond', cond_snplist,
        '--out', output_prefix,
    ]

    if cojo_collinear is not None:
        cmd.extend(['--cojo-collinear', str(cojo_collinear)])
    
    print(f"Running GCTA-COJO: {' '.join(cmd)}")
    
    with open(log_file, 'a') as log:
        result = subprocess.run(cmd, stdout=log, stderr=log)
    
    return result.returncode == 0

def parse_cojo_results(cojo_file, significance_threshold):
    """Parse COJO results and find the TOP significant SNP only"""
    if not os.path.exists(cojo_file):
        return None
    
    try:
        df = pd.read_csv(cojo_file, sep='\t')
        
        if 'pC' not in df.columns or len(df) == 0:
            return None
        
        # Filter significant SNPs
        sig_snps = df[df['pC'] < significance_threshold]
        
        if len(sig_snps) == 0:
            return None
        
        # Return ONLY the most significant SNP
        top_snp = sig_snps.loc[sig_snps['pC'].idxmin()]
        
        print(f"  Found top significant SNP: {top_snp['SNP']} (pC={top_snp['pC']:.2e})")
        
        return top_snp['SNP']
    
    except Exception as e:
        print(f"Error parsing COJO results: {e}", file=sys.stderr)
        return None

def save_snp_list(snps, output_file):
    """Save SNP list to file"""
    with open(output_file, 'w') as f:
        for snp in snps:
            f.write(f"{snp}\n")

def merge_iteration_results(output_dir, iteration):
    """Merge .cma.cojo and .given.cojo files"""
    cma_file = os.path.join(output_dir, f"iteration_{iteration}.cma.cojo")
    given_file = os.path.join(output_dir, f"iteration_{iteration}.given.cojo")
    output_file = os.path.join(output_dir, f"iteration_{iteration}_merged.txt")
    
    if not os.path.exists(cma_file):
        return False
    
    cma_df = pd.read_csv(cma_file, sep='\t')
    
    if os.path.exists(given_file):
        given_df = pd.read_csv(given_file, sep='\t')
        
        # Add missing columns to given_df
        for col in ['n', 'freq_geno', 'bC', 'bC_se', 'pC']:
            if col not in given_df.columns:
                given_df[col] = 'NA'
        
        # Ensure same columns
        given_df = given_df.reindex(columns=cma_df.columns, fill_value='NA')
        
        # Combine
        merged_df = pd.concat([cma_df, given_df], ignore_index=True)
    else:
        merged_df = cma_df.copy()
    
    # Sort by chromosome and position
    merged_df['Chr'] = pd.to_numeric(merged_df['Chr'], errors='coerce')
    merged_df['bp'] = pd.to_numeric(merged_df['bp'], errors='coerce')
    merged_df = merged_df.sort_values(['Chr', 'bp'])
    
    # Save
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    return True

def find_executable(names):
    """Find executable from a list of possible names"""
    import shutil
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return None

def main():
    parser = argparse.ArgumentParser(
        description='Iterative GCTA-COJO conditional analysis - ONE SNP per iteration'
    )
    parser.add_argument('--gwas-file', required=True,
                       help='Input GWAS summary statistics file (.tsv) [REQUIRED]')
    parser.add_argument('--bfile-prefix', required=True,
                       help='PLINK binary file prefix (without .bed/.bim/.fam) [REQUIRED]')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for results [REQUIRED]')
    parser.add_argument('--chromosome', type=int, default=None,
                       help='Chromosome number (auto-detected from GWAS file if not provided)')
    parser.add_argument('--gcta-path', default=None,
                       help='Path to GCTA executable (auto-detected if not provided)')
    parser.add_argument('--plink-path', default=None,
                       help='Path to PLINK executable (auto-detected if not provided)')
    parser.add_argument('--gcta-module', default='gcta/1.94.1',
                       help='GCTA module name to load (default: gcta/1.94.1)')
    parser.add_argument('--plink-module', default='plink/1.90',
                       help='PLINK module name to load (default: plink/1.90)')
    parser.add_argument('--python-module', default='python/3.9',
                       help='Python module name to load (default: python/3.9)')
    parser.add_argument('--sample-size', type=int, required=True,
                       help='Sample size [REQUIRED]')
    parser.add_argument('--freq-file', default=None,
                       help='PLINK frequency file (.frq) for allele frequencies')
    parser.add_argument('--generate-freq', action='store_true',
                       help='Generate frequency file using PLINK before analysis')
    parser.add_argument('--max-iterations', type=int, default=20,
                       help='Maximum number of iterations (default: 20)')
    parser.add_argument('--significance-threshold', type=float, default=5e-8,
                       help='Significance threshold for pC (default: 5e-8)')
    parser.add_argument('--load-modules', action='store_true',
                       help='Load environment modules before running (requires module system)')
    
    args = parser.parse_args()
    
    # Auto-detect executables if not provided
    if args.gcta_path is None:
        args.gcta_path = find_executable(['gcta64', 'gcta', 'gcta-1.94.1', 'gcta_1.94.1'])
        if args.gcta_path is None:
            print("ERROR: GCTA executable not found. Please specify with --gcta-path or load gcta module", file=sys.stderr)
            sys.exit(1)
        print(f"Auto-detected GCTA: {args.gcta_path}")
    
    if args.plink_path is None:
        args.plink_path = find_executable(['plink', 'plink2', 'plink1.9'])
        if args.plink_path is None:
            print("ERROR: PLINK executable not found. Please specify with --plink-path or load plink module", file=sys.stderr)
            sys.exit(1)
        print(f"Auto-detected PLINK: {args.plink_path}")
    
    # Load modules if requested
    if args.load_modules:
        print("Loading environment modules...")
        for module in [args.gcta_module, args.plink_module, args.python_module]:
            try:
                subprocess.run(['module', 'load', module], check=False, 
                             capture_output=True, text=True)
                print(f"  Loaded: {module}")
            except Exception as e:
                print(f"  Warning: Could not load module {module}: {e}")
    
    print("="*60)
    print("ITERATIVE GCTA-COJO CONDITIONAL ANALYSIS")
    print("="*60)
    print(f"GWAS file: {args.gwas_file}")
    print(f"Reference panel: {args.bfile_prefix}")
    print(f"Output directory: {args.output_dir}")
    print(f"Chromosome: {args.chromosome}")
    print(f"Sample size: {args.sample_size}")
    print(f"Max iterations: {args.max_iterations}")
    print(f"Significance threshold: {args.significance_threshold}")
    print(f"GCTA path: {args.gcta_path}")
    print(f"PLINK path: {args.plink_path}")
    if args.load_modules:
        print(f"Modules: {args.gcta_module}, {args.plink_module}, {args.python_module}")
    if args.generate_freq:
        print("Will generate frequency file using PLINK")
    print("="*60)
    
    # Validate inputs
    if not validate_inputs(args.gwas_file, args.bfile_prefix, args.output_dir):
        sys.exit(1)
    
    # Auto-detect chromosome if not provided
    if args.chromosome is None:
        print("\nAuto-detecting chromosome from GWAS file...")
        args.chromosome = detect_chromosome(args.gwas_file)
        if args.chromosome is None:
            print("ERROR: Could not detect chromosome. Please specify with --chromosome", file=sys.stderr)
            sys.exit(1)
    
    # Generate frequency file if requested
    freq_file = args.freq_file
    if args.generate_freq and not freq_file:
        print("\nGenerating allele frequency file...")
        freq_file = generate_freq_file(args.bfile_prefix, args.output_dir, args.plink_path)
    
    # Prepare .ma file
    ma_file = prepare_ma_file(
        args.gwas_file, 
        args.bfile_prefix, 
        args.output_dir,
        args.sample_size,
        freq_file
    )
    
    # Create initial conditioning SNP list
    initial_snplist, current_snps = create_initial_snplist(
        args.gwas_file,
        args.output_dir,
        args.significance_threshold
    )
    
    # Initialize tracking
    iteration = 1
    iteration_summary = {
        'iterations': [],
        'final_snps': [],
        'total_iterations': 0,
        'converged': False
    }
    
    log_file = os.path.join(args.output_dir, "cojo_analysis.log")
    
    # Iterative COJO analysis
    while iteration <= args.max_iterations:
        print(f"\n{'='*60}")
        print(f"ITERATION {iteration}")
        print(f"{'='*60}")
        print(f"Current conditioning SNPs: {len(current_snps)}")
        print(f"SNPs: {', '.join(current_snps[:5])}{'...' if len(current_snps) > 5 else ''}")
        
        # Create current SNP list
        current_snplist = os.path.join(args.output_dir, f"iteration_{iteration}.snplist")
        save_snp_list(current_snps, current_snplist)
        
        # Run COJO
        output_prefix = os.path.join(args.output_dir, f"iteration_{iteration}")
        
        success = run_gcta_cojo(
            args.bfile_prefix,
            ma_file,
            current_snplist,
            output_prefix,
            args.chromosome,
            args.gcta_path,
            log_file
        )
        
        if not success:
            print(f"WARNING: GCTA-COJO failed at iteration {iteration}")
            print("Continuing to finalization...")
            break
        
        # Parse results
        cojo_file = f"{output_prefix}.cma.cojo"
        new_snp = parse_cojo_results(cojo_file, args.significance_threshold)
        
        # Merge results
        merge_iteration_results(args.output_dir, iteration)
        
        # Record iteration
        iteration_summary['iterations'].append({
            'iteration': iteration,
            'conditioning_snps': len(current_snps),
            'new_snp_found': new_snp is not None,
            'new_snp': new_snp if new_snp else 'None'
        })
        
        # Check convergence
        if new_snp is None:
            print(f"\nCONVERGENCE: No significant SNPs found at iteration {iteration}")
            iteration_summary['converged'] = True
            break
        
        # Add new SNP to conditioning list
        if new_snp not in current_snps:
            current_snps.append(new_snp)
            print(f"  Added SNP to conditioning list: {new_snp}")
        else:
            print(f"  SNP already in list (convergence): {new_snp}")
            iteration_summary['converged'] = True
            break
        
        iteration += 1
    
    # Finalize
    iteration_summary['total_iterations'] = iteration
    iteration_summary['final_snps'] = current_snps
    
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"Total iterations: {iteration_summary['total_iterations']}")
    print(f"Final conditioning SNPs: {len(current_snps)}")
    print(f"Converged: {iteration_summary['converged']}")
    
    # Save final SNP list
    final_snplist = os.path.join(args.output_dir, "final_cond_snplist.txt")
    save_snp_list(current_snps, final_snplist)
    print(f"\nFinal SNP list saved: {final_snplist}")
    
    # Save iteration summary
    summary_file = os.path.join(args.output_dir, "iteration_summary.json")
    with open(summary_file, 'w') as f:
        json.dump(iteration_summary, f, indent=2)
    print(f"Iteration summary saved: {summary_file}")
    
    # Merge all results
    print("\nMerging all iteration results...")
    all_merged = []
    for i in range(1, iteration + 1):
        merged_file = os.path.join(args.output_dir, f"iteration_{i}_merged.txt")
        if os.path.exists(merged_file):
            df = pd.read_csv(merged_file, sep='\t')
            df['iteration'] = i
            all_merged.append(df)
    
    if all_merged:
        final_merged = pd.concat(all_merged, ignore_index=True)
        final_merged_file = os.path.join(args.output_dir, "final_cojo_results_merged.txt")
        final_merged.to_csv(final_merged_file, sep='\t', index=False)
        print(f"Final merged results saved: {final_merged_file}")
    
    print(f"\n{'='*60}")
    print("SUCCESS: Analysis completed")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
