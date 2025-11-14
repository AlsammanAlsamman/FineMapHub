#!/usr/bin/env python3
"""
Extract loci from population-harmonized GWAS data using regions file
Creates individual locus folders with extracted GWAS data
"""

import argparse
import pandas as pd
import logging
import sys
from pathlib import Path

def setup_logging():
    """Setup basic logging for the script"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

logger = setup_logging()

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Extract loci from population-harmonized GWAS data')
    
    parser.add_argument('--gwas-file', required=True,
                       help='Path to population-harmonized GWAS file (.popharmonized.tsv)')
    parser.add_argument('--loci', required=True,
                       help='Path to loci regions file (chr, start, end, locus format)')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for target analysis (e.g., results/04_loci/analysis_name)')
    parser.add_argument('--analysis-name', required=True,
                       help='Target analysis name for logging and organization')
    
    return parser.parse_args()

def read_loci_file(loci_file):
    """
    Read loci regions file with chr, start, end, locus format
    
    Parameters
    ----------
    loci_file : str
        Path to loci regions file
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with chr, start, end, locus columns
    """
    try:
        logger.info(f"Reading loci regions file: {loci_file}")
        loci_df = pd.read_csv(loci_file, sep='\t')
        
        # Check required columns
        required_cols = ['chr', 'start', 'end', 'locus']
        missing_cols = [col for col in required_cols if col not in loci_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in loci file: {missing_cols}")
        
        # Convert chr to string and ensure positions are integers
        loci_df['chr'] = loci_df['chr'].astype(str)
        loci_df['start'] = pd.to_numeric(loci_df['start'], errors='coerce')
        loci_df['end'] = pd.to_numeric(loci_df['end'], errors='coerce')
        
        # Remove rows with invalid positions
        loci_df = loci_df.dropna(subset=['start', 'end']).reset_index(drop=True)
        
        logger.info(f"Successfully read {len(loci_df)} loci regions")
        return loci_df
        
    except Exception as e:
        logger.error(f"Error reading loci file: {e}")
        return None

def read_gwas_file(gwas_file):
    """
    Read population-harmonized GWAS file
    
    Parameters
    ---------- 
    gwas_file : str
        Path to population-harmonized GWAS file
        
    Returns
    -------
    pandas.DataFrame
        GWAS DataFrame
    """
    try:
        logger.info(f"Reading GWAS file: {gwas_file}")
        
        # Check file size
        file_size_mb = Path(gwas_file).stat().st_size / (1024 * 1024)
        logger.info(f"GWAS file size: {file_size_mb:.2f} MB")
        
        # Read with progress monitoring for large files
        if file_size_mb > 100:
            logger.info("Large file detected - reading with progress monitoring...")
            gwas_df = pd.read_csv(gwas_file, sep='\t', low_memory=False)
        else:
            gwas_df = pd.read_csv(gwas_file, sep='\t')
        
        logger.info(f"Successfully read GWAS data: {len(gwas_df)} variants")
        
        # Ensure chromosome column is string for consistent matching
        if 'CHR' in gwas_df.columns:
            gwas_df['CHR'] = gwas_df['CHR'].astype(str)
        elif 'chrom' in gwas_df.columns:
            gwas_df['chrom'] = gwas_df['chrom'].astype(str)
        
        return gwas_df
        
    except Exception as e:
        logger.error(f"Error reading GWAS file: {e}")
        return None

def extract_locus_data(gwas_df, chrom, start, end, locus_name):
    """
    Extract GWAS data for a specific locus region
    
    Parameters
    ----------
    gwas_df : pandas.DataFrame
        GWAS DataFrame
    chrom : str
        Chromosome
    start : int
        Start position
    end : int  
        End position
    locus_name : str
        Locus name for logging
        
    Returns
    -------
    pandas.DataFrame
        Filtered GWAS DataFrame for the locus
    """
    logger.info(f"Extracting locus {locus_name}: chr{chrom}:{start}-{end}")
    
    # Debug: Show GWAS columns
    logger.info(f"Available GWAS columns: {list(gwas_df.columns)}")
    
    # Determine chromosome and position columns (case insensitive)
    chr_col = None
    pos_col = None
    
    for col in gwas_df.columns:
        col_upper = col.upper()
        if col_upper in ['CHR', 'CHROM', 'CHROMOSOME']:
            chr_col = col
        elif col_upper in ['POS', 'POSITION', 'BP', 'BASE_PAIR']:
            pos_col = col
    
    if not chr_col or not pos_col:
        logger.error(f"Could not find chromosome and position columns. Available: {list(gwas_df.columns)}")
        return None
    
    logger.info(f"Using columns - Chr: {chr_col}, Pos: {pos_col}")
    
    # Normalize chromosome format (remove 'chr' prefix if present)
    chrom_normalized = str(chrom).replace('chr', '').replace('CHR', '')
    
    # Filter by chromosome and position range
    logger.info(f"Filtering chr {chrom_normalized} between {start}-{end}")
    
    # Convert chromosome column to string and normalize
    gwas_df_temp = gwas_df.copy()
    gwas_df_temp['temp_chr'] = gwas_df_temp[chr_col].astype(str).str.replace('chr', '').str.replace('CHR', '')
    
    locus_df = gwas_df_temp[
        (gwas_df_temp['temp_chr'] == chrom_normalized) &
        (gwas_df_temp[pos_col] >= start) &
        (gwas_df_temp[pos_col] <= end)
    ].copy()
    
    # Remove temporary column
    locus_df = locus_df.drop('temp_chr', axis=1)
    
    logger.info(f"Found {len(locus_df)} variants in locus {locus_name}")
    
    # Debug: Show sample of filtered data
    if len(locus_df) > 0:
        logger.info(f"Sample variant: chr{locus_df.iloc[0][chr_col]}:{locus_df.iloc[0][pos_col]}")
    else:
        logger.warning(f"No variants found for chr{chrom_normalized}:{start}-{end}")
    
    return locus_df

def process_loci_extraction(gwas_file, loci_file, output_dir, analysis_name):
    """
    Process loci extraction for all regions in the loci file
    
    Parameters
    ----------
    gwas_file : str
        Path to population-harmonized GWAS file
    loci_file : str
        Path to loci regions file
    output_dir : str
        Output directory for target analysis
    analysis_name : str
        Target analysis name
    """
    
    logger.info(f"Starting loci extraction for analysis: {analysis_name}")
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Read input files
    logger.info("Reading input files...")
    gwas_df = read_gwas_file(gwas_file)
    if gwas_df is None:
        logger.error("Failed to read GWAS file")
        return False
        
    loci_df = read_loci_file(loci_file)
    if loci_df is None:
        logger.error("Failed to read loci file")
        return False
    
    # Process each locus
    extracted_count = 0
    failed_count = 0
    
    for _, locus_row in loci_df.iterrows():
        try:
            chrom = locus_row['chr']
            start = int(locus_row['start'])
            end = int(locus_row['end'])
            locus_name = locus_row['locus']
            
            # Create locus-specific output directory
            locus_dir = Path(output_dir) / locus_name
            locus_dir.mkdir(parents=True, exist_ok=True)
            
            # Extract locus data
            locus_data = extract_locus_data(gwas_df, chrom, start, end, locus_name)
            
            if locus_data is not None and len(locus_data) > 0:
                # Save locus data
                output_file = locus_dir / f"{locus_name}.loci.tsv"
                locus_data.to_csv(output_file, sep='\t', index=False)
                logger.info(f"✅ Saved locus {locus_name}: {output_file} ({len(locus_data)} variants)")
                extracted_count += 1
            else:
                logger.warning(f"❌ No variants found for locus {locus_name}")
                failed_count += 1
                
        except Exception as e:
            logger.error(f"❌ Error processing locus {locus_row.get('locus', 'unknown')}: {e}")
            failed_count += 1
    
    # Summary
    logger.info(f"\n=== LOCI EXTRACTION SUMMARY ===")
    logger.info(f"Analysis: {analysis_name}")
    logger.info(f"Total loci processed: {len(loci_df)}")
    logger.info(f"Successfully extracted: {extracted_count}")
    logger.info(f"Failed extractions: {failed_count}")
    logger.info(f"Output directory: {output_dir}")
    
    return extracted_count > 0

def main():
    """Main function"""
    args = parse_arguments()
    
    logger.info("=== LOCI EXTRACTION ===")
    logger.info(f"GWAS file: {args.gwas_file}")
    logger.info(f"Loci file: {args.loci}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Analysis name: {args.analysis_name}")
    logger.info("=======================")
    
    # Check pandas import
    try:
        import pandas as pd
        logger.info(f"Pandas version: {pd.__version__}")
    except ImportError as e:
        logger.error(f"Failed to import pandas: {e}")
        return 1
    
    # Validate input files exist
    if not Path(args.gwas_file).exists():
        logger.error(f"GWAS file not found: {args.gwas_file}")
        return 1
    
    if not Path(args.loci).exists():
        logger.error(f"Loci file not found: {args.loci}")
        return 1
    
    # Process loci extraction
    success = process_loci_extraction(
        args.gwas_file,
        args.loci,
        args.output_dir,
        args.analysis_name
    )
    
    if success:
        logger.info("✅ Loci extraction completed successfully!")
        return 0
    else:
        logger.error("❌ Loci extraction failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())