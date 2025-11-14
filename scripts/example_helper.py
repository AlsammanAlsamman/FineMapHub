#!/usr/bin/env python3
"""
Example helper script for FineMapHub pipeline.
Demonstrates proper CLI argument handling - never reads YAML files directly.
"""

import argparse
import os
import sys
from pathlib import Path


def setup_logging():
    """Setup basic logging to stdout/stderr."""
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def process_sample(input_file: str, sample_name: str, output_dir: str) -> bool:
    """
    Process sample data - placeholder implementation.
    
    Args:
        input_file: Path to input sample file
        sample_name: Name of the sample
        output_dir: Output directory for results
        
    Returns:
        bool: True if processing successful
    """
    logger = setup_logging()
    
    logger.info(f"Processing sample: {sample_name}")
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output directory: {output_dir}")
    
    # Validate input file exists
    if not Path(input_file).exists():
        logger.error(f"Input file does not exist: {input_file}")
        return False
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Placeholder processing - read input file and create output
    try:
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Create sample output file
        output_file = Path(output_dir) / f"{sample_name}_processed.txt"
        with open(output_file, 'w') as f:
            f.write(f"Processed sample: {sample_name}\n")
            f.write(f"Input content length: {len(content)} characters\n")
            f.write(f"Processing completed successfully\n")
        
        logger.info(f"Created output file: {output_file}")
        return True
        
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        return False


def main():
    """Main function - handles CLI arguments and orchestrates processing."""
    parser = argparse.ArgumentParser(
        description="Example helper script for FineMapHub pipeline"
    )
    
    parser.add_argument(
        "--input", 
        required=True, 
        help="Input sample file path"
    )
    parser.add_argument(
        "--sample", 
        required=True, 
        help="Sample name/identifier"
    )
    parser.add_argument(
        "--output-dir", 
        required=True, 
        help="Output directory for results"
    )
    parser.add_argument(
        "--done-marker", 
        required=True, 
        help="Path for completion marker file"
    )
    
    args = parser.parse_args()
    logger = setup_logging()
    
    logger.info("Starting example processing")
    logger.info(f"Arguments: {vars(args)}")
    
    # Process the sample
    success = process_sample(args.input, args.sample, args.output_dir)
    
    if success:
        # Create completion marker
        Path(args.done_marker).parent.mkdir(parents=True, exist_ok=True)
        with open(args.done_marker, 'w') as f:
            f.write(f"Example processing completed for sample: {args.sample}\n")
            f.write(f"Timestamp: $(date)\n")
        
        logger.info(f"Created completion marker: {args.done_marker}")
        logger.info("Processing completed successfully")
        sys.exit(0)
    else:
        logger.error("Processing failed")
        sys.exit(1)


if __name__ == "__main__":
    main()