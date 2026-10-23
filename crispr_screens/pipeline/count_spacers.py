#!/usr/bin/env python3
"""
sgRNA Spacer Counter with UMI Support
=====================================

Counts sgRNA sequences from FASTQ files with optional UMI extraction.

Original version by Julia Joung
Refactored and extended by Claude

Features:
    - Single-end or paired-end FASTQ processing
    - Configurable key sequences for different library designs
    - UMI extraction and counting (paired-end only)
    - Support for gzipped FASTQ files
    - Detailed statistics output
    - Progress reporting

Usage:
    # Paired-end with UMI
    python count_spacers.py -f read1.fastq -r read2.fastq -i library.csv -o output_prefix
    
    # Single-end (no UMI)
    python count_spacers.py -f read1.fastq -i library.csv -o output_prefix
    
    # Without leading G in guide
    python count_spacers.py -f read1.fastq -i library.csv -o output -no-g
"""

import csv
import gzip
import argparse
import sys
from pathlib import Path
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Iterator, TextIO

import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class LibraryConfig:
    """Configuration for library-specific sequence keys."""
    
    # Forward read (R1) - contains sgRNA sequence
    fkey_start: int = 0          # Start index of key search region
    fkey_end: int = 45           # End index of key search region
    fkey: str = "CACCTTGTTG"     # Key sequence before guide (mU6 promoter)
    guide_length: int = 20       # Length of sgRNA sequence
    
    # Reverse read (R2) - contains UMI
    rkey_start: int = 5          # Start index of key search region
    rkey_end: int = 35           # End index of key search region
    rkey1: str = "ATCCTGCAGG"    # First key flanking UMI
    rkey2: str = "GGATCCAAAA"    # Second key flanking UMI
    
    # Options
    guide_g: bool = True         # Whether library has leading G before spacer


# Preset configurations for common library types
LIBRARY_PRESETS = {
    'mU6': LibraryConfig(
        fkey_start=0,
        fkey_end=45,
        fkey="CACCTTGTTG",
    ),
    'hU6': LibraryConfig(
        fkey_start=20,
        fkey_end=45,
        fkey="CGAAACACC",
    ),
}


@dataclass
class CountStats:
    """Statistics from counting run."""
    num_reads: int = 0
    perfect_matches: int = 0
    non_perfect_matches: int = 0
    key_not_found: int = 0
    reads_with_n: int = 0
    pair_mismatches: int = 0
    
    @property
    def percent_matched(self) -> float:
        total = self.perfect_matches + self.non_perfect_matches
        if total == 0:
            return 0.0
        return round(self.perfect_matches / total * 100, 1)
    
    def to_dict(self) -> dict:
        return {
            'num_reads': self.num_reads,
            'perfect_matches': self.perfect_matches,
            'non_perfect_matches': self.non_perfect_matches,
            'key_not_found': self.key_not_found,
            'reads_with_n': self.reads_with_n,
            'pair_mismatches': self.pair_mismatches,
            'percent_matched': self.percent_matched,
        }


# =============================================================================
# File I/O
# =============================================================================

def open_fastq(filepath: str) -> TextIO:
    """Open FASTQ file, handling gzip compression."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def load_guide_library(input_file: str) -> List[str]:
    """
    Load guide sequences from CSV file.
    
    Handles multiple formats:
    - Single column of sequences (no header)
    - Multi-column CSV with 'sgRNA Sequence' or similar column
    - First column contains sequences if no recognized header
    """
    guides = []
    
    with open(input_file, 'r') as f:
        reader = csv.reader(f)
        rows = list(reader)
    
    if not rows:
        return guides
    
    # Check if first row is a header
    header = rows[0]
    seq_col_idx = 0  # Default to first column
    has_header = False
    
    # Look for sequence column by name
    sequence_col_names = ['sgRNA Sequence', 'sgrna sequence', 'sequence', 'Sequence', 
                          'sgRNA', 'sgrna', 'guide', 'Guide', 'spacer', 'Spacer']
    
    for i, col_name in enumerate(header):
        col_clean = col_name.strip()
        if col_clean in sequence_col_names or col_clean.lower() in [s.lower() for s in sequence_col_names]:
            seq_col_idx = i
            has_header = True
            break
    
    # If no recognized column name, check if first column looks like a header
    if not has_header:
        first_val = header[0].strip().upper()
        # If first value doesn't look like a sequence (too short, contains non-ACGT), skip it
        if len(first_val) < 15 or not all(c in 'ACGTN' for c in first_val):
            has_header = True
    
    # Process rows
    start_row = 1 if has_header else 0
    for row in rows[start_row:]:
        if row and len(row) > seq_col_idx:
            seq = row[seq_col_idx].strip().upper()
            if seq and not seq.startswith('#') and len(seq) >= 15:
                # Validate it looks like a sequence
                if all(c in 'ACGTN' for c in seq):
                    guides.append(seq)
    
    return guides


# =============================================================================
# Core Counting Logic
# =============================================================================

def extract_guide_from_read(
    sequence: str,
    config: LibraryConfig
) -> Optional[str]:
    """
    Extract sgRNA sequence from forward read.
    
    Returns None if key sequence not found.
    """
    key_region = sequence[config.fkey_start:config.fkey_end]
    key_index = key_region.find(config.fkey)
    
    if key_index < 0:
        return None
    
    start = key_index + config.fkey_start + len(config.fkey)
    guide = sequence[start:start + config.guide_length]
    
    return guide if len(guide) == config.guide_length else None


def extract_umi_from_read(
    sequence: str,
    config: LibraryConfig
) -> Optional[str]:
    """
    Extract UMI sequence from reverse read.
    
    Returns None if flanking keys not found.
    """
    # Find first key
    region1 = sequence[config.rkey_start:config.rkey_end]
    index1 = region1.find(config.rkey1)
    
    # Find second key (search in extended region)
    region2 = sequence[config.rkey_start + 20:config.rkey_end + 20]
    index2 = region2.find(config.rkey2)
    
    if index1 < 0 or index2 < 0:
        return None
    
    # Extract UMI between the two keys
    umi_start = index1 + config.rkey_start + len(config.rkey1)
    umi_end = index2 + config.rkey_start + 20
    
    return sequence[umi_start:umi_end]


def count_spacers(
    r1_fastq_file: str,
    r2_fastq_file: Optional[str],
    guide_library: List[str],
    config: LibraryConfig,
    progress_interval: int = 1_000_000
) -> Tuple[Dict[str, int], Dict[str, Dict[str, int]], CountStats]:
    """
    Count sgRNA sequences from FASTQ files.
    
    Parameters
    ----------
    r1_fastq_file : str
        Path to forward read FASTQ
    r2_fastq_file : str or None
        Path to reverse read FASTQ (optional, for UMI extraction)
    guide_library : list
        List of expected guide sequences
    config : LibraryConfig
        Configuration for sequence extraction
    progress_interval : int
        Print progress every N reads
    
    Returns
    -------
    tuple
        (count_dict, umi_dict, stats)
        - count_dict: guide -> total count
        - umi_dict: guide -> {umi: count} (empty if no R2)
        - stats: CountStats object
    """
    # Initialize dictionaries
    count_dict = {guide: 0 for guide in guide_library}
    umi_dict = {guide: {} for guide in guide_library}
    stats = CountStats()
    
    # Prepare key sequence (add G if needed)
    effective_fkey = config.fkey
    if config.guide_g:
        effective_fkey += "G"
    
    # Create modified config with effective key
    effective_config = LibraryConfig(
        fkey_start=config.fkey_start,
        fkey_end=config.fkey_end,
        fkey=effective_fkey,
        guide_length=config.guide_length,
        rkey_start=config.rkey_start,
        rkey_end=config.rkey_end,
        rkey1=config.rkey1,
        rkey2=config.rkey2,
        guide_g=config.guide_g
    )
    
    # Process reads - paired or single-end
    if r2_fastq_file:
        # Paired-end with UMI extraction
        with open_fastq(r1_fastq_file) as f1, open_fastq(r2_fastq_file) as f2:
            f_iter = FastqGeneralIterator(f1)
            r_iter = FastqGeneralIterator(f2)
            
            for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in zip(f_iter, r_iter):
                stats.num_reads += 1
                
                # Progress reporting
                if stats.num_reads % progress_interval == 0:
                    print(f"  Processed {stats.num_reads:,} reads...", file=sys.stderr)
                
                # Verify read pair matching
                if f_id.split()[0] != r_id.split()[0]:
                    stats.pair_mismatches += 1
                    continue
                
                # Skip reads with N
                f_seq_upper = f_seq.upper()
                if 'N' in f_seq_upper:
                    stats.reads_with_n += 1
                    continue
                
                # Extract guide sequence
                guide = extract_guide_from_read(f_seq_upper, effective_config)
                
                if guide is None:
                    stats.key_not_found += 1
                    continue
                
                # Check if guide is in library
                if guide in count_dict:
                    count_dict[guide] += 1
                    stats.perfect_matches += 1
                    
                    # Extract and count UMI
                    umi = extract_umi_from_read(r_seq.upper(), effective_config)
                    if umi:
                        umi_dict[guide][umi] = umi_dict[guide].get(umi, 0) + 1
                else:
                    stats.non_perfect_matches += 1
    else:
        # Single-end, no UMI
        with open_fastq(r1_fastq_file) as f1:
            f_iter = FastqGeneralIterator(f1)
            
            for f_id, f_seq, f_q in f_iter:
                stats.num_reads += 1
                
                # Progress reporting
                if stats.num_reads % progress_interval == 0:
                    print(f"  Processed {stats.num_reads:,} reads...", file=sys.stderr)
                
                # Skip reads with N
                f_seq_upper = f_seq.upper()
                if 'N' in f_seq_upper:
                    stats.reads_with_n += 1
                    continue
                
                # Extract guide sequence
                guide = extract_guide_from_read(f_seq_upper, effective_config)
                
                if guide is None:
                    stats.key_not_found += 1
                    continue
                
                # Check if guide is in library
                if guide in count_dict:
                    count_dict[guide] += 1
                    stats.perfect_matches += 1
                else:
                    stats.non_perfect_matches += 1
    
    return count_dict, umi_dict, stats


# =============================================================================
# Output Functions
# =============================================================================

def write_count_file(
    count_dict: Dict[str, int],
    output_file: str
) -> None:
    """Write guide counts to CSV."""
    sorted_guides = OrderedDict(sorted(count_dict.items()))
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        for guide, count in sorted_guides.items():
            writer.writerow([guide, count])


def write_umi_file(
    umi_dict: Dict[str, Dict[str, int]],
    output_file: str
) -> None:
    """Write UMI counts to CSV."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        for guide, umis in sorted(umi_dict.items()):
            # Sort UMIs by count (descending)
            sorted_umis = dict(sorted(umis.items(), key=lambda x: x[1], reverse=True))
            writer.writerow([guide, str(sorted_umis)])


def write_stats_file(
    stats: CountStats,
    count_dict: Dict[str, int],
    output_file: str
) -> None:
    """Write statistics to text file."""
    counts = np.array(list(count_dict.values()))
    guides_with_reads = np.count_nonzero(counts)
    guides_no_reads = len(counts) - guides_with_reads
    percent_no_reads = round(guides_no_reads / len(counts) * 100, 1) if len(counts) > 0 else 0
    
    # Calculate skew ratio
    top_10 = np.percentile(counts, 90)
    bottom_10 = np.percentile(counts, 10)
    if bottom_10 > 0:
        skew_ratio = f"{top_10 / bottom_10:.2f}"
    else:
        skew_ratio = "N/A (bottom 10% is zero)"
    
    with open(output_file, 'w') as f:
        f.write(f"Number of reads processed: {stats.num_reads:,}\n")
        f.write(f"Number of perfect guide matches: {stats.perfect_matches:,}\n")
        f.write(f"Number of non-perfect guide matches: {stats.non_perfect_matches:,}\n")
        f.write(f"Number of reads where key was not found: {stats.key_not_found:,}\n")
        f.write(f"Number of reads with N: {stats.reads_with_n:,}\n")
        f.write(f"Percentage of guides that matched perfectly: {stats.percent_matched}%\n")
        f.write(f"Guides with zero reads: {guides_no_reads:,} ({percent_no_reads}%)\n")
        f.write(f"Skew ratio (90th/10th percentile): {skew_ratio}\n")


def run_counting(
    input_file: str,
    r1_fastq_file: str,
    r2_fastq_file: Optional[str],
    output_prefix: str,
    config: Optional[LibraryConfig] = None,
    guide_g: bool = True
) -> Tuple[Dict[str, int], CountStats]:
    """
    Run full counting pipeline.
    
    Parameters
    ----------
    input_file : str
        Path to library CSV (guide sequences)
    r1_fastq_file : str
        Path to forward read FASTQ
    r2_fastq_file : str or None
        Path to reverse read FASTQ (optional, for UMI extraction)
    output_prefix : str
        Prefix for output files
    config : LibraryConfig, optional
        Library configuration (uses default if None)
    guide_g : bool
        Whether guides have leading G
    
    Returns
    -------
    tuple
        (count_dict, stats)
    """
    if config is None:
        config = LibraryConfig()
    config.guide_g = guide_g
    
    print(f"Loading guide library from {input_file}...")
    guides = load_guide_library(input_file)
    print(f"  Loaded {len(guides):,} guides")
    
    print(f"Processing FASTQ files...")
    print(f"  R1: {r1_fastq_file}")
    if r2_fastq_file:
        print(f"  R2: {r2_fastq_file}")
    else:
        print(f"  R2: None (single-end mode, no UMI)")
    
    count_dict, umi_dict, stats = count_spacers(
        r1_fastq_file=r1_fastq_file,
        r2_fastq_file=r2_fastq_file,
        guide_library=guides,
        config=config
    )
    
    print(f"\nWriting output files...")
    write_count_file(count_dict, f"{output_prefix}_count.csv")
    if r2_fastq_file:
        write_umi_file(umi_dict, f"{output_prefix}_umi.csv")
    write_stats_file(stats, count_dict, f"{output_prefix}_stat.txt")
    
    print(f"\nResults:")
    print(f"  Total reads: {stats.num_reads:,}")
    print(f"  Perfect matches: {stats.perfect_matches:,} ({stats.percent_matched}%)")
    print(f"  Key not found: {stats.key_not_found:,}")
    
    return count_dict, stats


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Count sgRNA sequences from paired-end FASTQ files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python count_spacers.py -f read1.fastq -r read2.fastq -i library.csv -o sample1

    # Without leading G
    python count_spacers.py -f read1.fastq -r read2.fastq -i library.csv -o sample1 -no-g

    # With gzipped files
    python count_spacers.py -f read1.fastq.gz -r read2.fastq.gz -i library.csv -o sample1

    # Using hU6 preset
    python count_spacers.py -f read1.fastq -r read2.fastq -i library.csv -o sample1 --preset hU6

Output files:
    {output}_count.csv  - Guide counts
    {output}_umi.csv    - UMI counts per guide (only with paired-end)
    {output}_stat.txt   - Run statistics
        """
    )
    
    parser.add_argument('-f', '--fastq1', type=str, dest='r1_fastq_file',
                        required=True, help='Read 1 (forward) FASTQ file')
    parser.add_argument('-r', '--fastq2', type=str, dest='r2_fastq_file',
                        required=False, default=None, 
                        help='Read 2 (reverse) FASTQ file (optional, for UMI extraction)')
    parser.add_argument('-i', '--input', type=str, dest='input_file',
                        required=True, help='Guide library CSV file')
    parser.add_argument('-o', '--output', type=str, dest='output_file',
                        required=True, help='Output file prefix')
    parser.add_argument('-no-g', dest='guide_g', action='store_false',
                        help='Guides do NOT have leading G before spacer')
    parser.add_argument('--preset', type=str, choices=list(LIBRARY_PRESETS.keys()),
                        help='Use preset configuration for library type')
    parser.add_argument('--fkey', type=str,
                        help='Custom forward key sequence (overrides preset)')
    parser.set_defaults(guide_g=True)
    
    args = parser.parse_args()
    
    # Get configuration
    if args.preset:
        config = LIBRARY_PRESETS[args.preset]
    else:
        config = LibraryConfig()
    
    # Override with custom key if provided
    if args.fkey:
        config.fkey = args.fkey.upper()
    
    # Run counting
    run_counting(
        input_file=args.input_file,
        r1_fastq_file=args.r1_fastq_file,
        r2_fastq_file=args.r2_fastq_file,
        output_prefix=args.output_file,
        config=config,
        guide_g=args.guide_g
    )


if __name__ == '__main__':
    main()
