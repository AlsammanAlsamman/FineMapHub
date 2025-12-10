#!/usr/bin/env python3
"""
Generate COJO Excel Report
Creates comprehensive Excel report with multiple sheets summarizing COJO conditional analysis results
"""

import argparse
import json
import os
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate Excel report from COJO conditional analysis results'
    )
    parser.add_argument('--cojo-dir', required=True,
                       help='Directory containing COJO results (results/05_cojo)')
    parser.add_argument('--loci-info-file', required=True,
                       help='File containing loci information')
    parser.add_argument('--target-analyses', nargs='+', required=True,
                       help='List of target analysis names')
    parser.add_argument('--populations', type=str, required=True,
                       help='JSON string mapping target analysis to population')
    parser.add_argument('--output', required=True,
                       help='Output Excel file path')
    
    return parser.parse_args()

def load_loci_info(loci_info_file):
    """Load loci information from file"""
    print(f"Loading loci information from: {loci_info_file}")
    
    if not os.path.exists(loci_info_file):
        print(f"WARNING: Loci info file not found: {loci_info_file}")
        return {}
    
    loci_info = {}
    try:
        with open(loci_info_file, 'r') as f:
            header = f.readline().strip().split('\t')
            print(f"Header columns: {header}")
            
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= len(header):
                        row_dict = dict(zip(header, parts))
                        # Use the 'locus' column value as the key
                        locus_id = row_dict.get('locus', parts[0])
                        loci_info[locus_id] = row_dict
                        
            # Debug: show first locus
            if loci_info:
                first_locus = next(iter(loci_info.keys()))
                print(f"Example locus: {first_locus}")
                print(f"Example data: {loci_info[first_locus]}")
                
    except Exception as e:
        print(f"ERROR loading loci info: {e}")
        import traceback
        traceback.print_exc()
        return {}
    
    print(f"Loaded info for {len(loci_info)} loci")
    return loci_info

def collect_cojo_results(cojo_dir, target_analyses):
    """Collect COJO results from all target analyses and loci"""
    print(f"\nCollecting COJO results from: {cojo_dir}")
    print(f"Target analyses: {', '.join(target_analyses)}")
    
    results = defaultdict(lambda: defaultdict(dict))
    all_loci = set()
    
    for target_analysis in target_analyses:
        analysis_dir = os.path.join(cojo_dir, target_analysis)
        
        if not os.path.exists(analysis_dir):
            print(f"WARNING: Analysis directory not found: {analysis_dir}")
            continue
        
        # Find all LOC_* directories
        for locus_dir in Path(analysis_dir).glob("LOC_*"):
            locus_name = locus_dir.name
            
            # Skip merged loci
            if "-" in locus_name or "newLOC" in locus_name:
                continue
            
            all_loci.add(locus_name)
            
            cojo_results_dir = locus_dir / "cojo_results"
            
            # Read final SNP list
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
                except:
                    pass
            
            results[locus_name][target_analysis] = {
                'snps': snps,
                'count': len(snps),
                'summary': summary
            }
    
    print(f"Found results for {len(all_loci)} loci across {len(target_analyses)} target analyses")
    return results, sorted(all_loci)

def create_locus_summary_sheet(results, all_loci, target_analyses, populations, loci_info):
    """Create Sheet 1: Locus-level summary"""
    print("\nGenerating Sheet 1: Locus Summary")
    
    rows = []
    
    for locus in all_loci:
        row = {'Locus': locus}
        
        # Add locus information immediately after Locus column
        if locus in loci_info:
            info_dict = loci_info[locus]
            print(f"Adding info for {locus}: {list(info_dict.keys())}")
            for key, value in info_dict.items():
                if key.lower() != 'locus':  # Skip locus ID column (already added)
                    row[key] = value
        else:
            print(f"WARNING: No loci info found for {locus}")
        
        # Add count and SNPs for each target analysis
        for target_analysis in target_analyses:
            if target_analysis in results[locus]:
                row[f'{target_analysis}_count'] = results[locus][target_analysis]['count']
                row[f'{target_analysis}_snps'] = ', '.join(results[locus][target_analysis]['snps'])
            else:
                row[f'{target_analysis}_count'] = 0
                row[f'{target_analysis}_snps'] = ''
        
        # Calculate unique SNPs per reference panel
        ref_panel_snps = defaultdict(set)
        all_unique_snps = set()
        
        for target_analysis in target_analyses:
            if target_analysis in results[locus] and target_analysis in populations:
                pop = populations[target_analysis]
                snps = results[locus][target_analysis]['snps']
                ref_panel_snps[pop].update(snps)
                all_unique_snps.update(snps)
        
        # Add unique SNP counts and IDs per reference panel
        for pop in sorted(ref_panel_snps.keys()):
            row[f'{pop}_unique_snps'] = len(ref_panel_snps[pop])
            row[f'{pop}_unique_snp_ids'] = ', '.join(sorted(ref_panel_snps[pop]))
        
        # Add total unique SNPs across all reference panels
        row['total_unique_snps'] = len(all_unique_snps)
        row['total_unique_snp_ids'] = ', '.join(sorted(all_unique_snps))
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    print(f"Created locus summary with {len(df)} rows")
    return df

def create_snp_detail_sheet(results, all_loci, target_analyses, populations, loci_info):
    """Create Sheet 2: SNP-level details (one row per SNP-per-locus combination)"""
    print("\nGenerating Sheet 2: SNP Details")
    
    rows = []
    
    for locus in all_loci:
        # Collect all unique SNPs for this locus across all target analyses
        snp_to_analyses = defaultdict(list)
        
        for target_analysis in target_analyses:
            if target_analysis in results[locus]:
                for snp in results[locus][target_analysis]['snps']:
                    snp_to_analyses[snp].append(target_analysis)
        
        # Create one row per SNP for this locus
        for snp in sorted(snp_to_analyses.keys()):
            row = {
                'Locus': locus,
                'SNP': snp,
                'Target_Analyses': ', '.join(snp_to_analyses[snp]),
                'Num_Analyses': len(snp_to_analyses[snp])
            }
            
            # Add locus information
            if locus in loci_info:
                for key, value in loci_info[locus].items():
                    if key.lower() != 'locus' and key not in row:
                        row[f'Locus_{key}'] = value
            
            # Add reference panels where this SNP appears
            ref_panels = set()
            for target_analysis in snp_to_analyses[snp]:
                if target_analysis in populations:
                    ref_panels.add(populations[target_analysis])
            row['Reference_Panels'] = ', '.join(sorted(ref_panels))
            
            rows.append(row)
    
    df = pd.DataFrame(rows)
    print(f"Created SNP details with {len(df)} rows")
    return df

def create_distribution_summary_sheet(results, all_loci, target_analyses):
    """Create Sheet 3: Distribution summary (loci by number of conditioning SNPs)"""
    print("\nGenerating Sheet 3: Distribution Summary")
    
    # Count loci by number of SNPs for each target analysis
    distribution = defaultdict(lambda: defaultdict(int))
    
    for locus in all_loci:
        for target_analysis in target_analyses:
            count = 0
            if target_analysis in results[locus]:
                count = results[locus][target_analysis]['count']
            
            if count == 0:
                distribution[target_analysis]['0 SNPs'] += 1
            elif count == 1:
                distribution[target_analysis]['1 SNP'] += 1
            elif count == 2:
                distribution[target_analysis]['2 SNPs'] += 1
            elif count == 3:
                distribution[target_analysis]['3 SNPs'] += 1
            elif count >= 4:
                distribution[target_analysis]['4+ SNPs'] += 1
    
    # Create summary rows
    rows = []
    categories = ['0 SNPs', '1 SNP', '2 SNPs', '3 SNPs', '4+ SNPs']
    
    for category in categories:
        row = {'Category': category}
        for target_analysis in target_analyses:
            row[target_analysis] = distribution[target_analysis][category]
        rows.append(row)
    
    # Add total row
    total_row = {'Category': 'Total Loci'}
    for target_analysis in target_analyses:
        total_row[target_analysis] = len(all_loci)
    rows.append(total_row)
    
    df = pd.DataFrame(rows)
    print(f"Created distribution summary with {len(df)} rows")
    return df

def main():
    args = parse_arguments()
    
    print("="*60)
    print("COJO EXCEL REPORT GENERATOR")
    print("="*60)
    print(f"COJO directory: {args.cojo_dir}")
    print(f"Loci info file: {args.loci_info_file}")
    print(f"Target analyses: {', '.join(args.target_analyses)}")
    print(f"Output file: {args.output}")
    print("="*60)
    
    # Parse populations JSON
    try:
        populations = json.loads(args.populations.replace("'", '"'))
        print(f"Populations mapping: {populations}")
    except Exception as e:
        print(f"ERROR parsing populations: {e}")
        sys.exit(1)
    
    # Load loci information
    loci_info = load_loci_info(args.loci_info_file)
    
    # Collect COJO results
    results, all_loci = collect_cojo_results(args.cojo_dir, args.target_analyses)
    
    if len(all_loci) == 0:
        print("ERROR: No loci found in COJO results")
        sys.exit(1)
    
    # Generate sheets
    print("\n" + "="*60)
    print("GENERATING EXCEL SHEETS")
    print("="*60)
    
    sheet1 = create_locus_summary_sheet(results, all_loci, args.target_analyses, populations, loci_info)
    sheet2 = create_snp_detail_sheet(results, all_loci, args.target_analyses, populations, loci_info)
    sheet3 = create_distribution_summary_sheet(results, all_loci, args.target_analyses)
    
    # Write to Excel
    print(f"\nWriting Excel file: {args.output}")
    try:
        with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
            sheet1.to_excel(writer, sheet_name='Locus_Summary', index=False)
            sheet2.to_excel(writer, sheet_name='SNP_Details', index=False)
            sheet3.to_excel(writer, sheet_name='Distribution_Summary', index=False)
        
        print("✓ Excel file created successfully")
    except Exception as e:
        print(f"ERROR writing Excel file: {e}")
        sys.exit(1)
    
    print("\n" + "="*60)
    print("REPORT GENERATION COMPLETE")
    print("="*60)
    print(f"Output: {args.output}")
    print(f"Total loci: {len(all_loci)}")
    print(f"Total target analyses: {len(args.target_analyses)}")
    print(f"Total SNP entries: {len(sheet2)}")
    print("="*60)

if __name__ == '__main__':
    main()
