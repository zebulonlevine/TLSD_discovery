#!/usr/bin/env python3
"""
CRISPRi Screen Pipeline with SLURM Support
==========================================

End-to-end pipeline for CRISPRi screen analysis with cluster job submission.

Workflow:
    1. submit_counting_jobs() - Submit SLURM jobs to count sgRNAs from FASTQs
    2. check_job_status()     - Monitor job completion
    3. aggregate_counts()     - Recombine counts from all jobs
    4. submit_mageck_jobs()   - Submit SLURM jobs for MAGeCK analysis
    5. collect_results()      - Generate plots and summary

Stages can be run independently or as a full pipeline.

Original sgRNA counting by Julia Joung
Pipeline by Zeb Levine with assistance from Claude Opus 4.1
"""

import os
import sys
import re
import argparse
import subprocess
import shlex
import time
import yaml
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
from datetime import datetime
import random
from random import shuffle

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

# Repo path configuration: exports TLSD_DATA / TLSD_OUT / TLSD_REPO so configs
# and library tables can use ${TLSD_DATA}-style paths (expanded on load).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from common import paths as _paths  # noqa: E402


# =============================================================================
# Negative-control (NTC) identification
# =============================================================================

NTC_LABELS = frozenset({
    'negative_control', 'negative control', 'non-targeting', 'nontargeting',
    'ntc', 'neg_control', 'negative', 'safe_harbor',
})


def is_ntc(names: pd.Series) -> pd.Series:
    """Exact-label NTC detection on the leading token of a gene id.

    Substring matching false-flags real genes -- KNTC1 contains "NTC" and
    NEGR1 contains "NEG" -- which silently pools their sgRNAs into the null
    distribution and deletes the genes from the results. See CLAUDE.md section 3
    ("Control-label precision ... so real genes are never false-flagged as
    controls").

    Gene ids at every stage downstream of prepare_mageck_input() have the form
    'SYMBOL_suffix' (A1CF_P1P2, non-targeting_7_negative_control, NTC_412), so
    comparing the token before the first underscore is both sufficient and
    unambiguous. Matching is exact against NTC_LABELS, so it fails *closed*
    (an unrecognised control spelling becomes a fake gene, which is loud)
    rather than open (a real gene silently becomes a control).
    """
    head = names.astype(str).str.strip().str.split('_').str[0].str.lower()
    return head.isin(NTC_LABELS)


def is_ntc_pseudogene(gene_ids: pd.Series) -> pd.Series:
    """True for the synthetic NTC_{j} null genes built during scoring.

    Distinct from is_ntc(): these are the pseudogenes used to estimate the FDR,
    not the raw control sgRNAs.
    """
    return gene_ids.astype(str).str.startswith('NTC_')


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class SlurmConfig:
    """SLURM job submission settings."""
    partition: str = "20"
    mem: str = "8gb"
    cpus: int = 4
    time: str = "4:00:00"
    email: str = ""
    mail_type: str = "FAIL"  # NONE, BEGIN, END, FAIL, ALL
    # Conda env to activate inside the job. Empty = emit nothing, so jobs
    # inherit the submitting shell (the old behaviour; existing configs that
    # don't set this keep working). Setting it makes jobs self-contained.
    conda_env: str = ""
    # Path to conda.sh. Empty = derive from $CONDA_EXE at submit time.
    conda_sh: str = ""



@dataclass 
class LibraryReference:
    """A library reference file for counting."""
    name: str                    # e.g., "LibA", "extras"
    sequence_file: str           # Path to CSV with sgRNA sequences
    annotation_file: str = ""    # Path to full annotation CSV (optional)
    no_g: bool = True            # Whether guides lack leading G
    preset: str = "mU6"          # Key sequence preset


@dataclass
class PipelineConfig:
    """Full pipeline configuration."""
    
    # Project info
    project_name: str = "screen"
    
    # Input files
    sample_table: str = ""              # CSV: sample_name, r1_fastq, [r2_fastq], [lane]
    library_table: str = ""             # CSV: library_name, sequence_file, annotation_file
    comparisons_table: str = ""         # CSV: comparison_name, treatment, control
    
    # Directories
    fastq_dir: str = ""                 # Base directory for FASTQ files
    output_dir: str = "./results"       # Output directory
    counts_dir: str = "./counts"        # Directory for count files
    jobs_dir: str = "./jobs"            # Directory for SLURM scripts
    logs_dir: str = "./logs"            # Directory for job logs
    
    # Scripts
    count_script: str = "count_spacers.py"  # Path to counting script
    
    # Analysis parameters
    fdr: float = 0.05
    counts_threshold: int = 50           # Min combined counts per replicate pair
    min_treat_mean: int = 20
    
    # SLURM settings
    slurm: SlurmConfig = field(default_factory=SlurmConfig)
    
    # Execution control
    dry_run: bool = False               # Generate scripts but don't submit


def load_config(config_file: str) -> PipelineConfig:
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        data = yaml.safe_load(f)
    data = {k: _paths.expand(v) for k, v in data.items()}

    # Handle nested slurm config
    slurm_data = {k: _paths.expand(v) for k, v in data.pop('slurm', {}).items()}
    slurm_config = SlurmConfig(**slurm_data)
    
    config = PipelineConfig(**data)
    config.slurm = slurm_config
    
    return config


# =============================================================================
# Input Table Handling
# =============================================================================

def load_sample_table(filepath: str) -> pd.DataFrame:
    """
    Load sample table CSV.
    
    Required columns:
        - sample_name: Logical sample name (for aggregation)
        - r1_fastq: R1 FASTQ filename
    
    Optional columns:
        - r2_fastq: R2 FASTQ filename (for UMI extraction)
        - lane: Lane identifier (L001, L002, etc.)
        - fastq_id: Unique ID for this FASTQ pair (auto-generated if missing)
    """
    df = pd.read_csv(filepath)
    required = ['sample_name', 'r1_fastq']
    
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Sample table missing columns: {missing}")
    
    # Make r2_fastq optional
    if 'r2_fastq' not in df.columns:
        df['r2_fastq'] = None
    
    # Auto-generate fastq_id if not present
    if 'fastq_id' not in df.columns:
        if 'lane' in df.columns:
            df['fastq_id'] = df['sample_name'] + '_' + df['lane']
        else:
            df['fastq_id'] = df['sample_name'] + '_' + df.index.astype(str)
    
    return df


def load_library_table(filepath: str) -> pd.DataFrame:
    """
    Load library reference table CSV.
    
    Required columns:
        - library_name: Identifier for this library subset
        - sequence_file: CSV with sgRNA sequences (for counting)
    
    Optional columns:
        - annotation_file: CSV with full annotations (for aggregation)
        - no_g: Whether guides lack leading G (default: True)
        - preset: Key sequence preset - mU6 or hU6 (default: mU6)
        - library_filter: Value to filter Library column by (e.g., "Weissman")
        - gene_column: Column name to use for gene aggregation (e.g., "Gene_Single")
        - sequence_column: Column name for sequence matching (default: "sgRNA Sequence")
        - guide_length: Override guide length for counting (default: 20)
    """
    df = pd.read_csv(filepath)
    
    if 'library_name' not in df.columns or 'sequence_file' not in df.columns:
        raise ValueError("Library table requires 'library_name' and 'sequence_file'")
    
    # Set defaults
    if 'no_g' not in df.columns:
        df['no_g'] = True
    if 'preset' not in df.columns:
        df['preset'] = 'mU6'
    if 'annotation_file' not in df.columns:
        df['annotation_file'] = df['sequence_file']
    for c in ('sequence_file', 'annotation_file'):
        df[c] = df[c].map(_paths.expand)   # ${TLSD_DATA}/... -> absolute
    if 'library_filter' not in df.columns:
        df['library_filter'] = ''
    if 'gene_column' not in df.columns:
        df['gene_column'] = ''
    if 'sequence_column' not in df.columns:
        df['sequence_column'] = 'sgRNA Sequence'
    if 'guide_length' not in df.columns:
        df['guide_length'] = 20
    
    return df


def load_comparisons_table(filepath: str) -> pd.DataFrame:
    """
    Load comparisons table CSV.
    
    Required columns:
        - comparison_name: Unique name for this comparison
        - treatment: Treatment sample name(s), comma-separated
        - control: Control sample name(s), comma-separated
    
    Replicate pairing:
        Treatment and control lists are position-matched for replicate pairs.
        Example: treatment="T1,T2", control="C1,C2" means:
            - Replicate 1: T1 vs C1
            - Replicate 2: T2 vs C2
        
        For count thresholding, each replicate pair must have >= threshold
        combined counts. Guides failing ANY pair are excluded.
    """
    df = pd.read_csv(filepath)
    required = ['comparison_name', 'treatment', 'control']
    
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Comparisons table missing columns: {missing}")
    
    return df


def parse_replicate_pairs(treatment: str, control: str) -> List[Tuple[str, str]]:
    """
    Parse treatment and control strings into replicate pairs.
    
    Returns list of (treatment_sample, control_sample) tuples.
    Position-matched: T1,T2 vs C1,C2 -> [(T1,C1), (T2,C2)]
    """
    treat_samples = [s.strip() for s in treatment.split(',')]
    ctrl_samples = [s.strip() for s in control.split(',')]
    
    if len(treat_samples) != len(ctrl_samples):
        raise ValueError(
            f"Treatment and control must have same number of samples for replicate pairing. "
            f"Got {len(treat_samples)} treatments, {len(ctrl_samples)} controls."
        )
    
    return list(zip(treat_samples, ctrl_samples))


# =============================================================================
# SLURM Job Generation and Submission
# =============================================================================

def generate_slurm_header(
    job_name: str,
    output_file: str,
    config: SlurmConfig
) -> str:
    """Generate SLURM script header."""
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={output_file}",
        f"#SBATCH --error={output_file.replace('.out', '.err')}",
        f"#SBATCH --partition={config.partition}",
        f"#SBATCH --mem={config.mem}",
        f"#SBATCH --cpus-per-task={config.cpus}",
        f"#SBATCH --time={config.time}",
    ]
    
    if config.email:
        lines.append(f"#SBATCH --mail-user={config.email}")
        lines.append(f"#SBATCH --mail-type={config.mail_type}")

    lines.append("")

    # Activate the analysis env inside the job. Without this the job inherits
    # whatever environment the submitting shell had, which has silently failed
    # before (ModuleNotFoundError: pandas._libs.window.aggregations) and left
    # comparison directories containing only a log file.
    if config.conda_env:
        conda_sh = config.conda_sh
        if not conda_sh:
            conda_exe = os.environ.get('CONDA_EXE', '')
            if conda_exe:
                conda_sh = os.path.join(
                    os.path.dirname(os.path.dirname(conda_exe)),
                    'etc', 'profile.d', 'conda.sh')
        if not conda_sh:
            raise ValueError(
                "slurm.conda_env is set but conda.sh could not be located. "
                "Set slurm.conda_sh explicitly, or submit from an active conda shell."
            )
        lines.append(f"source {conda_sh}")
        lines.append(f"conda activate {config.conda_env}")
        lines.append("")

    return "\n".join(lines)


def submit_job(script_path: str, dry_run: bool = False) -> Optional[str]:
    """Submit SLURM job and return job ID."""
    if dry_run:
        print(f"  [DRY RUN] Would submit: {script_path}")
        return None
    
    result = subprocess.run(["sbatch", script_path], capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  ERROR submitting {script_path}: {result.stderr}")
        return None
    
    match = re.search(r'Submitted batch job (\d+)', result.stdout)
    if match:
        job_id = match.group(1)
        print(f"  Submitted {Path(script_path).name} -> Job {job_id}")
        return job_id
    return None


def check_jobs_complete(job_ids: List[str]) -> Tuple[bool, Dict[str, str]]:
    """Check if all SLURM jobs are complete."""
    if not job_ids:
        return True, {}
    
    result = subprocess.run(
        ["sacct", "-j", ",".join(job_ids), "--format=JobID,State", "--noheader", "-P"],
        capture_output=True, text=True
    )
    
    status = {}
    for line in result.stdout.strip().split('\n'):
        if line and '|' in line:
            parts = line.split('|')
            jid = parts[0].split('.')[0]
            state = parts[1] if len(parts) > 1 else "UNKNOWN"
            if jid in job_ids:
                status[jid] = state
    
    complete_states = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT'}
    all_complete = all(status.get(jid, 'PENDING') in complete_states for jid in job_ids)
    
    return all_complete, status


def wait_for_jobs(job_ids: List[str], poll_interval: int = 60) -> Dict[str, str]:
    """Wait for all jobs to complete."""
    if not job_ids:
        return {}
    
    print(f"\nWaiting for {len(job_ids)} jobs...")
    
    while True:
        complete, status = check_jobs_complete(job_ids)
        
        n_complete = sum(1 for s in status.values() if s in {'COMPLETED', 'FAILED', 'CANCELLED'})
        n_running = sum(1 for s in status.values() if s == 'RUNNING')
        n_pending = len(job_ids) - n_complete - n_running
        
        print(f"  [{datetime.now().strftime('%H:%M:%S')}] "
              f"Done: {n_complete}, Running: {n_running}, Pending: {n_pending}")
        
        if complete:
            break
        time.sleep(poll_interval)
    
    failed = [jid for jid, s in status.items() if s != 'COMPLETED']
    if failed:
        print(f"\n  WARNING: {len(failed)} jobs failed: {failed}")
    
    return status


# =============================================================================
# Stage 1: Counting Jobs
# =============================================================================

def generate_counting_jobs(
    sample_table: pd.DataFrame,
    library_table: pd.DataFrame,
    config: PipelineConfig
) -> List[Tuple[str, str, str]]:
    """Generate SLURM scripts for counting jobs."""
    
    os.makedirs(config.jobs_dir, exist_ok=True)
    os.makedirs(config.logs_dir, exist_ok=True)
    os.makedirs(config.counts_dir, exist_ok=True)
    
    jobs = []
    
    for _, sample in sample_table.iterrows():
        fastq_id = sample['fastq_id']
        r1 = os.path.join(config.fastq_dir, sample['r1_fastq']) if config.fastq_dir else sample['r1_fastq']
        
        # Handle optional R2
        r2_val = sample.get('r2_fastq')
        has_r2 = pd.notna(r2_val) and str(r2_val).strip() != ''
        if has_r2:
            r2 = os.path.join(config.fastq_dir, r2_val) if config.fastq_dir else r2_val
            r2_arg = f"-r {r2} "
        else:
            r2_arg = ""
        
        for _, lib in library_table.iterrows():
            lib_name = lib['library_name']
            seq_file = lib['sequence_file']
            no_g_val = lib['no_g']
            # Default to True (use -no-g flag) if not specified
            no_g_flag = "-no-g" if (pd.isna(no_g_val) or no_g_val) else ""
            preset = lib.get('preset', 'mU6')
            if pd.isna(preset) or not preset:
                preset = 'mU6'
            
            job_name = f"count_{fastq_id}_{lib_name}"
            output_prefix = os.path.join(config.counts_dir, f"{fastq_id}_{lib_name}")
            count_file = f"{output_prefix}_count.csv"
            log_file = os.path.join(config.logs_dir, f"{job_name}.out")
            
            cmd = (
                f"python {config.count_script} "
                f"-f {r1} {r2_arg}"
                f"-i {seq_file} "
                f"-o {output_prefix} "
                f"--preset {preset} {no_g_flag}"
            ).strip()
            
            script = generate_slurm_header(job_name, log_file, config.slurm)
            script += f"\n{cmd}\n"
            
            script_path = os.path.join(config.jobs_dir, f"{job_name}.sh")
            with open(script_path, 'w') as f:
                f.write(script)
            
            jobs.append((script_path, count_file, job_name))
    
    return jobs


def submit_counting_jobs(
    sample_table: pd.DataFrame,
    library_table: pd.DataFrame,
    config: PipelineConfig
) -> List[str]:
    """Generate and submit all counting jobs."""
    
    print(f"\n{'='*60}")
    print("STAGE 1: Submitting counting jobs")
    print(f"{'='*60}")
    
    jobs = generate_counting_jobs(sample_table, library_table, config)
    
    n_samples = len(sample_table)
    n_libs = len(library_table)
    print(f"\n  {n_samples} FASTQ pairs × {n_libs} libraries = {len(jobs)} jobs")
    
    job_ids = []
    skipped = 0
    
    for script_path, count_file, job_name in jobs:
        if os.path.exists(count_file):
            skipped += 1
            continue
        job_id = submit_job(script_path, config.dry_run)
        if job_id:
            job_ids.append(job_id)
    
    if skipped:
        print(f"  Skipped {skipped} jobs (count files exist)")
    
    # Save job tracking
    os.makedirs(config.output_dir, exist_ok=True)
    tracking = os.path.join(config.output_dir, 'counting_jobs.txt')
    with open(tracking, 'w') as f:
        f.write('\n'.join(job_ids))
    
    print(f"\n  Submitted {len(job_ids)} jobs")
    return job_ids


# =============================================================================
# Stage 2: Count Aggregation
# =============================================================================

def aggregate_counts(
    sample_table: pd.DataFrame,
    library_table: pd.DataFrame,
    config: PipelineConfig
) -> pd.DataFrame:
    """Aggregate counts from all counting jobs."""
    
    print(f"\n{'='*60}")
    print("STAGE 2: Aggregating counts")
    print(f"{'='*60}")
    
    # Load all library annotations
    all_annotations = []
    for _, lib in library_table.iterrows():
        ann = pd.read_csv(lib['annotation_file'])
        
        # Apply library filter if specified
        lib_filter = lib.get('library_filter', '')
        if pd.notna(lib_filter) and lib_filter and 'Library' in ann.columns:
            ann = ann[ann['Library'] == lib_filter].copy()
            print(f"  Filtered {lib['library_name']} to Library='{lib_filter}': {len(ann)} sgRNAs")
        
        # Apply custom gene column if specified
        gene_col = lib.get('gene_column', '')
        if pd.notna(gene_col) and gene_col and gene_col in ann.columns:
            ann['_gene_override'] = ann[gene_col]
            print(f"  Using '{gene_col}' for gene aggregation in {lib['library_name']}")
        
        # Determine sequence column for this library
        seq_col_name = lib.get('sequence_column', 'sgRNA Sequence')
        if pd.isna(seq_col_name) or not seq_col_name:
            seq_col_name = 'sgRNA Sequence'
        if seq_col_name in ann.columns:
            ann['_match_sequence'] = ann[seq_col_name].str.upper().str.strip()
            print(f"  Using '{seq_col_name}' for sequence matching in {lib['library_name']}")
        else:
            # Fallback to sgRNA Sequence
            for col in ['sgRNA Sequence', 'sequence', 'Sequence']:
                if col in ann.columns:
                    ann['_match_sequence'] = ann[col].str.upper().str.strip()
                    break
        
        # Store guide length for this library
        guide_len_val = lib.get('guide_length', 20)
        guide_len = int(guide_len_val) if pd.notna(guide_len_val) else 20
        ann['_guide_length'] = guide_len
        
        ann['_library'] = lib['library_name']
        all_annotations.append(ann)
    
    library = pd.concat(all_annotations, ignore_index=True)
    print(f"  Total: {len(library)} sgRNAs from {len(library_table)} library files")
    
    # Create sequence lookup dict with prefix support
    # For 21bp guides counted as 20bp, we'll match on first 20 chars
    seq_to_idx = {}  # sequence -> list of row indices
    seq_prefixes = {}  # 20bp prefix -> list of (full_seq, row_idx) for longer guides
    
    for idx, row in library.iterrows():
        seq = row['_match_sequence']
        if pd.isna(seq):
            continue
        seq = str(seq).upper()
        
        # Exact match lookup
        if seq not in seq_to_idx:
            seq_to_idx[seq] = []
        seq_to_idx[seq].append(idx)
        
        # Prefix lookup for sequences longer than 20bp
        if len(seq) > 20:
            prefix = seq[:20]
            if prefix not in seq_prefixes:
                seq_prefixes[prefix] = []
            seq_prefixes[prefix].append((seq, idx))
    
    print(f"  Annotated sequences: {len(seq_to_idx)}")
    print(f"  Sequences with >20bp (prefix matching enabled): {len(seq_prefixes)}")
    
    # Track all sequences and their counts across all samples
    all_seq_counts = {}  # seq -> {sample: count}
    
    # Initialize count columns
    sample_names = list(sample_table['sample_name'].unique())
    for sample_name in sample_names:
        library[sample_name] = 0.0
    
    # Aggregate counts per sample
    for sample_name in sample_names:
        sample_rows = sample_table[sample_table['sample_name'] == sample_name]
        
        for _, row in sample_rows.iterrows():
            fastq_id = row['fastq_id']
            
            for _, lib in library_table.iterrows():
                lib_name = lib['library_name']
                count_file = os.path.join(config.counts_dir, f"{fastq_id}_{lib_name}_count.csv")
                
                if not os.path.exists(count_file):
                    print(f"  WARNING: Missing {count_file}")
                    continue
                
                counts = pd.read_csv(count_file, names=['sequence', 'counts'])
                counts['sequence'] = counts['sequence'].str.upper().str.strip()
                
                for _, cnt_row in counts.iterrows():
                    seq = cnt_row['sequence']
                    cnt = cnt_row['counts']
                    
                    # Track for unmatched reporting
                    if seq not in all_seq_counts:
                        all_seq_counts[seq] = {}
                    all_seq_counts[seq][sample_name] = all_seq_counts[seq].get(sample_name, 0) + cnt
                    
                    # Try exact match first
                    if seq in seq_to_idx:
                        for idx in seq_to_idx[seq]:
                            library.at[idx, sample_name] += cnt
                    # Try prefix match for 20bp counted sequences
                    elif len(seq) == 20 and seq in seq_prefixes:
                        for full_seq, idx in seq_prefixes[seq]:
                            library.at[idx, sample_name] += cnt
    
    print(f"  Aggregated {len(sample_names)} samples")
    
    # Find unmatched sequences (not matched by exact or prefix)
    matched_seqs = set(seq_to_idx.keys()) | set(seq_prefixes.keys())
    all_counted_seqs = set(all_seq_counts.keys())
    unmatched_seqs = all_counted_seqs - matched_seqs
    
    # Also remove sequences that matched via prefix
    truly_unmatched = set()
    for seq in unmatched_seqs:
        if len(seq) == 20 and seq in seq_prefixes:
            continue  # This matched via prefix
        truly_unmatched.add(seq)
    
    if truly_unmatched:
        print(f"\n  WARNING: {len(truly_unmatched)} sequences in count files have no annotation!")
        print(f"           These will be EXCLUDED from analysis.")
        
        # Build unmatched dataframe
        unmatched_data = []
        for seq in truly_unmatched:
            row_data = {'sequence': seq}
            for sn in sample_names:
                row_data[sn] = all_seq_counts[seq].get(sn, 0)
            unmatched_data.append(row_data)
        
        unmatched_df = pd.DataFrame(unmatched_data)
        unmatched_df['total_counts'] = unmatched_df[sample_names].sum(axis=1)
        unmatched_df = unmatched_df.sort_values('total_counts', ascending=False)
        
        unmatched_file = os.path.join(config.output_dir, 'unmatched_sequences.csv')
        unmatched_df.to_csv(unmatched_file, index=False)
        
        print(f"           Saved to: {unmatched_file}")
        print(f"           Top unmatched by counts:")
        for _, row in unmatched_df.head(5).iterrows():
            print(f"             {row['sequence']}: {int(row['total_counts'])} total counts")
    else:
        print(f"\n  ✓ All counted sequences have annotations (including prefix matches)")
    
    # Also copy sgRNA Sequence to standard column for downstream
    if '_match_sequence' in library.columns and 'sgRNA Sequence' not in library.columns:
        library['sgRNA Sequence'] = library['_match_sequence']
    
    # Save aggregated counts (only annotated sequences)
    agg_file = os.path.join(config.output_dir, 'aggregated_counts.csv')
    library.to_csv(agg_file, index=False)
    print(f"  Saved: {agg_file}")
    
    return library


def prepare_mageck_input(
    library: pd.DataFrame,
    sample_names: List[str],
    output_file: str
) -> pd.DataFrame:
    """Format for MAGeCK input, excluding unannotated sequences."""
    
    library = library.copy()
    seq_col = 'sgRNA Sequence'
    gene_col = 'Gene'
    
    # Debug: show what columns exist
    print(f"  Library columns: {list(library.columns)}")
    
    # Drop existing Gene column if present to avoid conflicts
    if gene_col in library.columns:
        print(f"  Dropping existing '{gene_col}' column")
        library = library.drop(columns=[gene_col])
    
    # Initialize Gene column as empty
    library[gene_col] = ''
    
    # Priority 1: Use _gene_override if specified (from library_table gene_column)
    if '_gene_override' in library.columns:
        # Debug: show _gene_override values
        print(f"  _gene_override column found, sample values: {library['_gene_override'].head(3).tolist()}")
        
        override_mask = library['_gene_override'].notna() & \
                       (library['_gene_override'].astype(str).str.strip() != '') & \
                       (library['_gene_override'].astype(str).str.upper() != 'NA')
        library.loc[override_mask, gene_col] = library.loc[override_mask, '_gene_override'].astype(str)
        n_override = override_mask.sum()
        if n_override > 0:
            print(f"  Applied gene_column override to {n_override} sgRNAs")
            # Show sample of overridden values
            sample_genes = library.loc[override_mask, gene_col].head(3).tolist()
            print(f"    Sample overridden genes: {sample_genes}")
        else:
            print(f"  WARNING: _gene_override column exists but no valid values found!")
            print(f"    Non-null count: {library['_gene_override'].notna().sum()}")
            print(f"    Sample raw values: {library['_gene_override'].head(5).tolist()}")
    else:
        print(f"  No _gene_override column found")
    
    # Priority 2: Fill remaining from standard columns
    gene_columns = ['NCBI', 'Ensembl', 'Original', 'gene', 'Gene_Symbol', 'Symbol']
    available_cols = [col for col in gene_columns if col in library.columns]
    print(f"  Available gene columns for fallback: {available_cols}")
    
    for col in available_cols:
        # Fill empty/NA values from this column, but only if source has valid value
        fill_mask = (library[gene_col].astype(str).str.strip() == '') | \
                    (library[gene_col].isna())
        source_valid = library[col].notna() & \
                      (library[col].astype(str).str.strip() != '') & \
                      (library[col].astype(str).str.upper() != 'NA')
        fill_mask = fill_mask & source_valid
        if fill_mask.any():
            n_filled = fill_mask.sum()
            library.loc[fill_mask, gene_col] = library.loc[fill_mask, col].astype(str)
            print(f"    Filled {n_filled} from '{col}'")
    
    # Debug: show final Gene column
    print(f"  Final Gene column sample: {library[gene_col].head(5).tolist()}")
    print(f"  Unique genes: {library[gene_col].nunique()}")
    
    # Identify rows with no valid gene annotation
    unknown_mask = (library[gene_col].isna() | 
                   (library[gene_col].astype(str).str.upper() == 'NA') |
                   (library[gene_col].astype(str).str.strip() == ''))
    
    n_unknown = unknown_mask.sum()
    if n_unknown > 0:
        print(f"\n  WARNING: {n_unknown} sgRNAs have no gene annotation!")
        print(f"           These will be EXCLUDED from MAGeCK analysis.")
        
        # Save unknown sequences for inspection
        output_dir = os.path.dirname(output_file)
        unknown_file = os.path.join(output_dir, 'unknown_gene_sequences.csv')
        unknown_df = library[unknown_mask].copy()
        unknown_df['total_counts'] = unknown_df[list(sample_names)].sum(axis=1)
        unknown_df = unknown_df.sort_values('total_counts', ascending=False)
        unknown_df.to_csv(unknown_file, index=False)
        print(f"           Saved to: {unknown_file}")
        
        # Show top offenders
        if len(unknown_df) > 0:
            print(f"           Top unknown by counts:")
            for _, row in unknown_df.head(5).iterrows():
                print(f"             {row[seq_col]}: {int(row['total_counts'])} counts")
        
        # Exclude unknown from analysis
        library = library[~unknown_mask].copy()
        print(f"           Proceeding with {len(library)} annotated sgRNAs")
    
    # Rename NTCs (after removing unknowns)
    # Match the NTC label exactly, via the shared is_ntc() helper. The old
    # substring regex also caught real genes: KNTC1 (contains "NTC") and
    # NEGR1 (contains "NEG").
    ntc_mask = is_ntc(library[gene_col])

    # report anything the old regex would have swallowed, so silent
    # mislabelling can't come back
    _gene_norm = library[gene_col].astype(str).str.strip().str.lower()
    _legacy = _gene_norm.str.contains('neg|non-targeting|ntc|negative_control', na=False)
    _fp = _legacy & ~ntc_mask
    if _fp.any():
        print(f"  NTC regex false positives avoided ({_fp.sum()} rows): "
              f"{sorted(library.loc[_fp, gene_col].unique())}")

    # The NTC_{j} pseudogene namespace used during scoring collides with any
    # real symbol literally named NTC. None exists in the libraries in use;
    # fail loudly rather than silently mis-scoring if that ever changes.
    _real = library.loc[~ntc_mask, gene_col].astype(str).str.split('_').str[0]
    assert not (_real.str.lower() == 'ntc').any(), (
        "A real gene is named 'NTC', which collides with the NTC_{j} pseudogene "
        "labels used by the scorers. Rename the pseudogene prefix before proceeding."
    )

    n_ntc = int(ntc_mask.sum())
    print(f"  NTCs identified: {n_ntc}")
    library.loc[ntc_mask, gene_col] = [f'non-targeting_{i}' for i in range(n_ntc)]
    
    # Add isoform suffix - but NOT for rows with custom gene_column override
    # (custom gene columns like Gene_Single typically already have the desired naming)
    if 'Isoforms' in library.columns:
        if '_gene_override' in library.columns:
            # Only add isoform suffix to rows WITHOUT override
            no_override_mask = library['_gene_override'].isna() | \
                              (library['_gene_override'].astype(str).str.strip() == '') | \
                              (library['_gene_override'].astype(str).str.upper() == 'NA')
            if no_override_mask.any():
                library.loc[no_override_mask, gene_col] = \
                    library.loc[no_override_mask, gene_col].astype(str) + '_' + \
                    library.loc[no_override_mask, 'Isoforms'].astype(str).apply(lambda x: str(x).split(',')[0])
                print(f"  Added isoform suffix to {no_override_mask.sum()} rows (skipped {(~no_override_mask).sum()} with gene_column override)")
            else:
                print(f"  Skipping isoform suffix (all rows have custom gene_column)")
        else:
            # No override column - add isoform suffix to all
            library[gene_col] = library[gene_col].astype(str) + '_' + \
                library['Isoforms'].astype(str).apply(lambda x: str(x).split(',')[0])
            print(f"  Added isoform suffix from 'Isoforms' column")
    
    # Create sgRNA names
    library['sgRNA'] = ''
    for gene in library[gene_col].unique():
        mask = library[gene_col] == gene
        library.loc[mask, 'sgRNA'] = [f'{gene}_sg{i}' for i in range(mask.sum())]
    
    # Final output
    output_df = library[['sgRNA', gene_col] + list(sample_names)]
    output_df.columns = ['sgRNA', 'Gene'] + list(sample_names)
    
    # Debug: verify Gene column before saving
    print(f"  Final output Gene column sample:")
    for _, row in output_df.head(5).iterrows():
        print(f"    {row['sgRNA']} -> {row['Gene']}")
    
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved MAGeCK input: {output_file}")
    print(f"    {len(output_df)} sgRNAs, {output_df['Gene'].nunique()} genes")

    # Persist the sgRNA -> Gene -> sequence map. counts_for_mageck.txt carries
    # the final sgRNA/Gene names but not the spacer sequence, while
    # aggregated_counts.csv carries the sequence but not the names (they are
    # built here). Supplementary Table S1 needs all three, and joining those two
    # files on row order breaks silently as soon as a guide is dropped as
    # unannotated -- so write the mapping while all three are in one frame.
    if seq_col in library.columns:
        ann_file = os.path.join(os.path.dirname(output_file), 'sgrna_annotation.csv')
        library[['sgRNA', gene_col, seq_col]].rename(
            columns={gene_col: 'Gene'}).to_csv(ann_file, index=False)
        print(f"  Saved sgRNA annotation: {ann_file}")

    return library


# =============================================================================
# Stage 3: MAGeCK Jobs
# =============================================================================

def generate_mageck_jobs(
    comparisons: pd.DataFrame,
    counts_file: str,
    config: PipelineConfig
) -> List[Tuple[str, str, str]]:
    """Generate SLURM scripts for MAGeCK + U-test analysis."""
    
    jobs = []
    
    # Get absolute path to this script for U-test
    pipeline_script = os.path.abspath(__file__)
    
    for _, comp in comparisons.iterrows():
        comp_name = comp['comparison_name']
        treatment = comp['treatment']
        control = comp['control']
        
        job_name = f"mageck_{comp_name}"
        comp_dir = os.path.join(config.output_dir, comp_name)
        os.makedirs(comp_dir, exist_ok=True)
        
        output_prefix = os.path.join(comp_dir, comp_name)
        log_file = os.path.join(config.logs_dir, f"{job_name}.out")
        sgrna_file = f"{output_prefix}.sgrna_summary.txt"
        thresholded_file = f"{output_prefix}_thresholded_counts.txt"
        
        # Threshold command (calls this script's threshold subcommand)
        threshold_cmd = (
            f"python {pipeline_script} threshold "
            f"--counts {counts_file} "
            f"--treatment {treatment} "
            f"--control {control} "
            f"--min-counts {config.counts_threshold} "
            f"--output {thresholded_file}"
        )
        
        # MAGeCK command (uses thresholded counts)
        mageck_cmd = f"mageck test -k {thresholded_file} -t {treatment} -c {control} -n {output_prefix}"
        
        # U-test command (calls this script's utest subcommand)
        utest_cmd = (
            f"python {pipeline_script} utest "
            f"--sgrna-summary {sgrna_file} "
            f"--output-prefix {output_prefix} "
            f"--fdr {config.fdr} "
            f"--min-treat-mean {config.min_treat_mean}"
        )
        
        script = generate_slurm_header(job_name, log_file, config.slurm)
        script += f"\n# Threshold counts (require >= {config.counts_threshold} per replicate pair)\n{threshold_cmd}\n"
        script += f"\n# Run MAGeCK\n{mageck_cmd}\n"
        script += f"\n# Run U-test analysis\n{utest_cmd}\n"
        
        script_path = os.path.join(config.jobs_dir, f"{job_name}.sh")
        with open(script_path, 'w') as f:
            f.write(script)
        
        jobs.append((script_path, comp_dir, job_name))
    
    return jobs


def threshold_counts(
    counts_file: str,
    treatment: str,
    control: str,
    min_counts: int,
    output_file: str
) -> pd.DataFrame:
    """
    Threshold counts based on replicate pairs.
    
    For each replicate pair (treatment_i, control_i), requires:
        treatment_i + control_i >= min_counts
    
    Guides failing ANY replicate pair are excluded.
    """
    df = pd.read_table(counts_file)
    
    # Parse replicate pairs
    pairs = parse_replicate_pairs(treatment, control)
    
    print(f"Thresholding counts:")
    print(f"  Input: {len(df)} sgRNAs")
    print(f"  Replicate pairs: {len(pairs)}")
    print(f"  Min counts per pair: {min_counts}")
    
    # Build mask: guide passes if ALL replicate pairs meet threshold
    pass_mask = pd.Series(True, index=df.index)
    
    for i, (treat_sample, ctrl_sample) in enumerate(pairs):
        pair_sum = df[treat_sample] + df[ctrl_sample]
        pair_pass = pair_sum >= min_counts
        pass_mask = pass_mask & pair_pass
        
        n_pass = pair_pass.sum()
        print(f"    Pair {i+1} ({treat_sample} + {ctrl_sample}): {n_pass} pass")
    
    df_filtered = df[pass_mask].copy()
    
    print(f"  Output: {len(df_filtered)} sgRNAs ({len(df) - len(df_filtered)} removed)")
    
    df_filtered.to_csv(output_file, sep='\t', index=False)
    return df_filtered


def run_threshold_standalone(
    counts_file: str,
    treatment: str,
    control: str,
    min_counts: int,
    output_file: str
) -> None:
    """Run count thresholding as standalone command (called by SLURM jobs)."""
    threshold_counts(counts_file, treatment, control, min_counts, output_file)


def submit_mageck_jobs(
    comparisons: pd.DataFrame,
    counts_file: str,
    config: PipelineConfig
) -> List[str]:
    """Submit MAGeCK jobs."""
    
    print(f"\n{'='*60}")
    print("STAGE 3: Submitting MAGeCK jobs")
    print(f"{'='*60}")
    
    jobs = generate_mageck_jobs(comparisons, counts_file, config)
    print(f"\n  {len(jobs)} comparisons")
    
    job_ids = []
    for script_path, comp_dir, job_name in jobs:
        comp_name = job_name.replace('mageck_', '')
        sgrna_file = os.path.join(comp_dir, f"{comp_name}.sgrna_summary.txt")
        
        if os.path.exists(sgrna_file):
            print(f"  [SKIP] {job_name} - results exist")
            continue
        
        job_id = submit_job(script_path, config.dry_run)
        if job_id:
            job_ids.append(job_id)
    
    tracking = os.path.join(config.output_dir, 'mageck_jobs.txt')
    with open(tracking, 'w') as f:
        f.write('\n'.join(job_ids))
    
    print(f"\n  Submitted {len(job_ids)} jobs")
    return job_ids


# =============================================================================
# Stage 4: Results Processing
# =============================================================================

def mann_whitney_analysis(
    sgrna_summary: str,
    output_prefix: str,
    fdr: float = 0.05,
    min_treat_mean: int = 20
) -> Tuple[pd.DataFrame, pd.DataFrame, float]:
    """Mann-Whitney U test analysis."""
    
    df_mageck = pd.read_table(sgrna_summary)
    df = df_mageck[df_mageck['treat_mean'] > min_treat_mean].copy()
    # Exact leading-token match, not substring: KNTC1 contains "NTC" and would
    # otherwise be dropped from the gene set and pooled into the null.
    ntc_where = is_ntc(df['Gene'])
    df_ntc = df[ntc_where]
    df_targeting = df[~ntc_where]
    
    ntc_pvalues = list(df_ntc['p.twosided'])
    ntc_p_lfc = list(zip(df_ntc['p.twosided'], df_ntc['LFC']))
    
    genes = df_targeting['Gene'].unique()
    gene_results = {}
    
    for gene in genes:
        df_gene = df_targeting[df_targeting['Gene'] == gene].sort_values('p.twosided')
        lfc = df_gene.iloc[:3]['LFC'].mean()
        if len(ntc_pvalues) > 0:
            _, pvalue = mannwhitneyu(list(df_gene['p.twosided'])[:3], ntc_pvalues, alternative='two-sided')
        else:
            pvalue = 1.0
        gene_results[gene] = [lfc, pvalue]
    
    # Null distribution
    random.seed(42)
    for j in range(len(genes)):
        if len(ntc_p_lfc) >= 5:
            shuffle(ntc_p_lfc)
            selected = ntc_p_lfc[:5]
            selected_p = [x[0] for x in selected]
            selected_lfc = np.mean([x[1] for x in sorted(selected, key=lambda x: x[0])][:3])
            _, pval = mannwhitneyu(selected_p, ntc_pvalues, alternative='two-sided')
            gene_results[f'NTC_{j}'] = [selected_lfc, pval]
    
    results = pd.DataFrame(gene_results).T
    results.columns = ['epsilon', 'pvalue']
    results = results.reset_index().rename(columns={'index': 'gene_id'})
    results['gene'] = results['gene_id'].apply(lambda x: x.split('_')[0])
    
    ntc_med = results[results['gene'] == 'NTC']['epsilon'].median()
    if pd.notna(ntc_med):
        results['epsilon'] -= ntc_med
    
    results['product'] = results['epsilon'] * (-np.log10(results['pvalue'].clip(lower=1e-300)))
    
    # FDR threshold
    threshold = 0
    for t in np.arange(0, abs(results['product']).max(), 0.1):
        above = results[abs(results['product']) > t]
        if len(above) > 0 and (is_ntc_pseudogene(above['gene_id']).sum() / len(above)) < fdr:
            threshold = t
            break
    
    hits = results[abs(results['product']) > threshold]
    
    results.sort_values('product', ascending=False).to_csv(f'{output_prefix}_all_genes.csv', index=False)
    hits.sort_values('product', ascending=False).to_csv(f'{output_prefix}_hits_fdr{fdr}.csv', index=False)
    
    return results, hits, threshold


def _line_score(parallel_dist, perp_dist, mode: str = 'perp_weighted', eps: float = 1e-10):
    """Per-sgRNA line score from signed along-line (parallel) and off-line
    (perpendicular) distances to the reference line.

    mode='perp_weighted' (default, original behavior):
        parallel * |parallel| / (|parallel| + perp). The signed projection is
        scaled by |parallel| / (|parallel| + perp), so guides on the line keep
        their full projection while guides far off the line shrink toward zero.
    mode='parallel':
        the raw signed distance along the line, no perpendicular normalization.
    """
    parallel_dist = np.asarray(parallel_dist, dtype=float)
    perp_dist = np.asarray(perp_dist, dtype=float)
    if mode == 'parallel':
        return parallel_dist
    if mode == 'perp_weighted':
        abs_parallel = np.abs(parallel_dist)
        return parallel_dist * abs_parallel / (abs_parallel + perp_dist + eps)
    raise ValueError(f"score_mode must be 'perp_weighted' or 'parallel', got {mode!r}")


def sgrna_line_score_utest(
    sgrna_file1: str,
    sgrna_file2: str,
    slope: float,
    intercept: float,
    output_prefix: str,
    comp1_name: str = None,
    comp2_name: str = None,
    top_n: int = 3,
    fdr: float = 0.05,
    min_count: int = 20,
    score_mode: str = 'perp_weighted'
) -> Tuple[pd.DataFrame, pd.DataFrame, float]:
    """
    Calculate line score at sgRNA level, then U-test top N per gene against NTCs.
    
    This performs a two-comparison analysis:
    1. For each sgRNA, calculate line score based on LFC from both comparisons
    2. For each gene, select top N sgRNAs by |score|
    3. U-test those scores against NTC score distribution
    
    Parameters
    ----------
    sgrna_file1, sgrna_file2 : str
        MAGeCK sgrna_summary.txt files from two comparisons
    slope, intercept : float
        Reference line parameters (y = slope*x + intercept)
    output_prefix : str
        Prefix for output files
    comp1_name, comp2_name : str, optional
        Names for the two comparisons (default: parent directory names)
    top_n : int
        Number of top sgRNAs per gene for U-test (default: 3)
    fdr : float
        FDR threshold for hit calling
    min_count : int
        Minimum treat_mean for sgRNA inclusion
    
    Returns
    -------
    results : DataFrame
        All genes with epsilon (mean score), pvalue, product
    hits : DataFrame
        Genes passing FDR threshold
    threshold : float
        Product threshold for FDR
    """
    # Load sgRNA summaries
    df1 = pd.read_table(sgrna_file1)
    df2 = pd.read_table(sgrna_file2)
    
    # Filter by count threshold
    df1 = df1[df1['treat_mean'] > min_count].copy()
    df2 = df2[df2['treat_mean'] > min_count].copy()
    
    # Get names
    if comp1_name is None:
        comp1_name = os.path.basename(os.path.dirname(sgrna_file1))
    if comp2_name is None:
        comp2_name = os.path.basename(os.path.dirname(sgrna_file2))
    
    print(f"sgRNA line score U-test analysis:")
    print(f"  Comparison 1: {comp1_name} ({len(df1)} sgRNAs)")
    print(f"  Comparison 2: {comp2_name} ({len(df2)} sgRNAs)")
    print(f"  Reference line: y = {slope}x + {intercept}")
    print(f"  Top N sgRNAs per gene: {top_n}")
    
    # Merge on sgRNA
    merged = df1[['sgrna', 'Gene', 'LFC']].merge(
        df2[['sgrna', 'LFC']], 
        on='sgrna', 
        how='inner',
        suffixes=('_1', '_2')
    )
    merged.columns = ['sgrna', 'gene', 'lfc1', 'lfc2']
    
    print(f"  Merged sgRNAs: {len(merged)}")
    
    # Calculate line score for each sgRNA using raw LFC values
    # Convention: negative LFC = depletion (matches MAGeCK)
    x = merged['lfc1'].values
    y = merged['lfc2'].values
    m = slope
    b = intercept
    
    norm = np.sqrt(1 + m**2)
    
    # Perpendicular distance (absolute)
    perp_dist = np.abs(y - m*x - b) / norm
    
    # Parallel distance (signed)
    parallel_dist = (x + m * (y - b)) / norm
    
    # Score from along-line (parallel) and off-line (perpendicular) distances;
    # mode selects perpendicular-weighted projection vs. raw distance along line.
    merged['score'] = _line_score(parallel_dist, perp_dist, mode=score_mode)
    merged['parallel_dist'] = parallel_dist
    merged['perp_dist'] = perp_dist
    
    # Identify NTCs -- exact leading-token match, not substring. A case-insensitive
    # 'NTC' substring also matches the real gene KNTC1 (see is_ntc docstring).
    ntc_mask = is_ntc(merged['gene'])
    
    df_ntc = merged[ntc_mask]
    df_targeting = merged[~ntc_mask]
    
    print(f"  NTC sgRNAs: {len(df_ntc)}")
    print(f"  Targeting sgRNAs: {len(df_targeting)}")
    
    # NTC score distribution for U-test
    ntc_scores = df_ntc['score'].tolist()
    
    # Gene-level analysis: top N sgRNAs per gene by absolute score
    genes = df_targeting['gene'].unique()
    gene_results = {}
    
    for gene in genes:
        df_gene = df_targeting[df_targeting['gene'] == gene].copy()
        
        # Sort by absolute score (most extreme first), then take top N
        df_gene['abs_score'] = df_gene['score'].abs()
        df_gene = df_gene.sort_values('abs_score', ascending=False)
        top_sgrnas = df_gene.head(top_n)
        
        # Mean score of top N (this is "epsilon" for the gene)
        mean_score = top_sgrnas['score'].mean()
        
        # U-test: compare top N scores to NTC distribution
        if len(ntc_scores) > 0 and len(top_sgrnas) > 0:
            _, pvalue = mannwhitneyu(
                top_sgrnas['score'].tolist(), 
                ntc_scores, 
                alternative='two-sided'
            )
        else:
            pvalue = 1.0
        
        gene_results[gene] = [mean_score, pvalue, len(df_gene)]
    
    # Null distribution from NTCs
    random.seed(42)
    ntc_scores_for_null = list(zip(df_ntc['score'], df_ntc['gene']))
    
    for j in range(len(genes)):
        if len(ntc_scores_for_null) >= 5:
            shuffle(ntc_scores_for_null)
            selected = ntc_scores_for_null[:5]
            selected_scores = [x[0] for x in selected]
            # Sort by absolute score, take top N
            selected_sorted = sorted(selected_scores, key=lambda x: abs(x), reverse=True)[:top_n]
            mean_score = np.mean(selected_sorted)
            _, pval = mannwhitneyu(selected_sorted, ntc_scores, alternative='two-sided')
            gene_results[f'NTC_{j}'] = [mean_score, pval, 5]
    
    # Build results DataFrame
    results = pd.DataFrame(gene_results).T
    results.columns = ['epsilon', 'pvalue', 'n_sgrnas']
    results = results.reset_index().rename(columns={'index': 'gene_id'})
    results['gene'] = results['gene_id'].apply(lambda x: x.split('_')[0] if not x.startswith('NTC_') else 'NTC')
    
    # Center on NTC median
    ntc_med = results[results['gene'] == 'NTC']['epsilon'].median()
    if pd.notna(ntc_med):
        results['epsilon'] -= ntc_med
        print(f"  NTC median score: {ntc_med:.4f} (subtracted)")
    
    # Product score
    results['product'] = results['epsilon'] * (-np.log10(results['pvalue'].clip(lower=1e-300)))
    
    # FDR threshold
    threshold = 0
    for t in np.arange(0, abs(results['product']).max(), 0.1):
        above = results[abs(results['product']) > t]
        if len(above) > 0:
            ntc_frac = is_ntc_pseudogene(above['gene_id']).sum() / len(above)
            if ntc_frac < fdr:
                threshold = t
                break
    
    hits = results[abs(results['product']) > threshold].copy()
    hits = hits[~is_ntc_pseudogene(hits['gene_id'])]
    
    print(f"  Product threshold (FDR={fdr}): {threshold:.2f}")
    print(f"  Hits: {len(hits)} genes")
    
    # Save results
    results.sort_values('product', ascending=False).to_csv(
        f'{output_prefix}_all_genes.csv', index=False
    )
    hits.sort_values('product', ascending=False).to_csv(
        f'{output_prefix}_hits_fdr{fdr}.csv', index=False
    )
    
    # Also save sgRNA-level scores
    merged.to_csv(f'{output_prefix}_sgrna_scores.csv', index=False)
    
    print(f"  Saved: {output_prefix}_all_genes.csv")
    print(f"  Saved: {output_prefix}_hits_fdr{fdr}.csv")
    print(f"  Saved: {output_prefix}_sgrna_scores.csv")
    
    return results, hits, threshold


COLORS = {'red': '#DC0000FF', 'blue': '#3C5488FF', 'salmon': '#F39B7FFF', 'green': '#00A087FF', 'grey': '#7F7F7F'}


def plot_volcano(df, hits, threshold, output_file, title="", fdr=0.05, label_genes=None, label_n_top=10):
    """Generate volcano plot."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    df_ntc = df[is_ntc_pseudogene(df['gene_id'])]
    hits_pos = hits[(hits['epsilon'] > 0) & (~is_ntc_pseudogene(hits['gene_id']))]
    hits_neg = hits[(hits['epsilon'] < 0) & (~is_ntc_pseudogene(hits['gene_id']))]
    
    ax.scatter(df['epsilon'], -np.log10(df['pvalue']), c=COLORS['salmon'], s=5, alpha=0.7, zorder=1)
    ax.scatter(hits_pos['epsilon'], -np.log10(hits_pos['pvalue']), c=COLORS['red'], s=8, label=f'Pos (n={len(hits_pos)})', zorder=2)
    ax.scatter(hits_neg['epsilon'], -np.log10(hits_neg['pvalue']), c=COLORS['blue'], s=8, label=f'Neg (n={len(hits_neg)})', zorder=2)
    ax.scatter(df_ntc['epsilon'], -np.log10(df_ntc['pvalue']), c=COLORS['grey'], s=3, alpha=0.5, zorder=1)
    
    x = np.linspace(0.01, 10, 100)
    ax.plot(x, threshold/x, '--k', lw=0.8)
    ax.plot(-x, threshold/x, '--k', lw=0.8, label=f'FDR={fdr}')
    
    if label_genes:
        for gene in label_genes:
            sub = df[df['gene'].str.upper() == gene.upper()]
            if len(sub):
                r = sub.iloc[0]
                ax.annotate(gene, (r['epsilon'], -np.log10(r['pvalue'])), fontsize=9)
                ax.scatter(r['epsilon'], -np.log10(r['pvalue']), c=COLORS['green'], s=30, zorder=3, edgecolors='k', lw=0.5)
    
    if label_n_top:
        for _, r in hits[~is_ntc_pseudogene(hits['gene_id'])].head(label_n_top).iterrows():
            ax.annotate(r['gene'], (r['epsilon'], -np.log10(r['pvalue'])), fontsize=8)
    
    xlim = max(abs(df['epsilon'].min()), df['epsilon'].max()) + 0.5
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(0, -np.log10(df['pvalue'].min()) + 0.5)
    ax.set_xlabel('Phenotype (ε)')
    ax.set_ylabel('-log₁₀(p)')
    ax.set_title(title)
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def plot_scatter(df, x_col, y_col, output_file, highlight_genes=None, log=True):
    """Scatter plot comparing conditions."""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    x, y = df[x_col].values, df[y_col].values
    if log:
        x, y = np.log10(x + 1), np.log10(y + 1)
    
    r, _ = pearsonr(x, y)
    ax.scatter(x, y, c=COLORS['blue'], s=3, alpha=0.5)
    
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, '--', c='grey', lw=0.8)
    
    if highlight_genes and 'gene' in df.columns:
        for gene in highlight_genes:
            sub = df[df['gene'].str.upper() == gene.upper()]
            if len(sub):
                r_ = sub.iloc[0]
                xi, yi = (np.log10(r_[x_col]+1), np.log10(r_[y_col]+1)) if log else (r_[x_col], r_[y_col])
                ax.scatter(xi, yi, c=COLORS['red'], s=50, zorder=3, edgecolors='k')
                ax.annotate(gene, (xi, yi), fontsize=9)
    
    ax.set_xlabel(f'log₁₀({x_col}+1)' if log else x_col)
    ax.set_ylabel(f'log₁₀({y_col}+1)' if log else y_col)
    ax.set_title(f'r = {r:.3f}')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def process_results(comparisons: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Collect and summarize results from MAGeCK + U-test jobs."""
    
    print(f"\n{'='*60}")
    print("STAGE 4: Collecting results")
    print(f"{'='*60}")
    
    summary_rows = []
    
    for _, comp in comparisons.iterrows():
        comp_name = comp['comparison_name']
        comp_dir = os.path.join(config.output_dir, comp_name)
        
        # Check for U-test results
        hits_file = os.path.join(comp_dir, f"{comp_name}_hits_fdr{config.fdr}.csv")
        all_genes_file = os.path.join(comp_dir, f"{comp_name}_all_genes.csv")
        
        if not os.path.exists(hits_file):
            print(f"  WARNING: Missing results for {comp_name}")
            continue
        
        df_hits = pd.read_csv(hits_file)
        df_all = pd.read_csv(all_genes_file)
        
        n_pos = len(df_hits[(df_hits['epsilon'] > 0) & (~is_ntc_pseudogene(df_hits['gene_id']))])
        n_neg = len(df_hits[(df_hits['epsilon'] < 0) & (~is_ntc_pseudogene(df_hits['gene_id']))])
        n_genes = len(df_all[~is_ntc_pseudogene(df_all['gene_id'])])
        
        summary_rows.append({
            'comparison': comp_name,
            'total_genes': n_genes,
            'positive_hits': n_pos,
            'negative_hits': n_neg,
            'total_hits': n_pos + n_neg,
            'fdr': config.fdr
        })
        
        print(f"  {comp_name}: {n_pos} pos, {n_neg} neg hits")
    
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(config.output_dir, 'analysis_summary.csv'), index=False)
    return summary


def run_utest_standalone(
    sgrna_summary: str,
    output_prefix: str,
    fdr: float = 0.05,
    min_treat_mean: int = 20
) -> None:
    """Run U-test analysis as standalone command (called by SLURM jobs)."""
    
    print(f"Running U-test analysis...")
    print(f"  Input: {sgrna_summary}")
    print(f"  Output prefix: {output_prefix}")
    print(f"  FDR: {fdr}")
    
    df_all, df_hits, threshold = mann_whitney_analysis(
        sgrna_summary, output_prefix, fdr, min_treat_mean
    )
    
    n_pos = len(df_hits[(df_hits['epsilon'] > 0) & (~is_ntc_pseudogene(df_hits['gene_id']))])
    n_neg = len(df_hits[(df_hits['epsilon'] < 0) & (~is_ntc_pseudogene(df_hits['gene_id']))])
    
    print(f"  Positive hits: {n_pos}")
    print(f"  Negative hits: {n_neg}")
    print(f"  Product threshold: {threshold:.2f}")
    print(f"  Results saved to {output_prefix}_*.csv")


# =============================================================================
# Full Pipeline
# =============================================================================

def run_pipeline(config: PipelineConfig, wait: bool = False):
    """Run complete pipeline."""
    os.makedirs(config.output_dir, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"CRISPRi Screen Pipeline: {config.project_name}")
    print(f"{'='*60}")
    
    sample_table = load_sample_table(config.sample_table)
    library_table = load_library_table(config.library_table)
    comparisons = load_comparisons_table(config.comparisons_table)
    
    print(f"\n  Samples: {len(sample_table)}")
    print(f"  Libraries: {len(library_table)}")
    print(f"  Comparisons: {len(comparisons)}")
    
    # Stage 1
    job_ids = submit_counting_jobs(sample_table, library_table, config)
    
    if wait and job_ids:
        wait_for_jobs(job_ids)
    elif job_ids:
        print("\n  Jobs submitted. Re-run with --wait or after completion.")
        return
    
    # Stage 2
    library = aggregate_counts(sample_table, library_table, config)
    sample_names = list(sample_table['sample_name'].unique())
    counts_file = os.path.join(config.output_dir, 'counts_for_mageck.txt')
    prepare_mageck_input(library, sample_names, counts_file)
    
    # Stage 3
    job_ids = submit_mageck_jobs(comparisons, counts_file, config)
    
    if wait and job_ids:
        wait_for_jobs(job_ids)
    elif job_ids:
        print("\n  Jobs submitted. Re-run with --wait or after completion.")
        return
    
    # Stage 4
    summary = process_results(comparisons, config)
    
    print(f"\n{'='*60}")
    print("Complete!")
    print(f"{'='*60}")
    print(f"\n{summary.to_string(index=False)}")


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='CRISPRi Screen Pipeline')
    sub = parser.add_subparsers(dest='command')
    
    run_p = sub.add_parser('run', help='Run full pipeline')
    run_p.add_argument('--config', required=True)
    run_p.add_argument('--wait', action='store_true')
    run_p.add_argument('--dry-run', action='store_true')
    
    count_p = sub.add_parser('count', help='Submit counting jobs')
    count_p.add_argument('--config', required=True)
    count_p.add_argument('--dry-run', action='store_true')
    
    agg_p = sub.add_parser('aggregate', help='Aggregate counts')
    agg_p.add_argument('--config', required=True)
    
    mageck_p = sub.add_parser('mageck', help='Submit MAGeCK + U-test jobs')
    mageck_p.add_argument('--config', required=True)
    mageck_p.add_argument('--dry-run', action='store_true')
    
    results_p = sub.add_parser('results', help='Collect results summary')
    results_p.add_argument('--config', required=True)
    
    status_p = sub.add_parser('status', help='Check job status')
    status_p.add_argument('--config', required=True)
    
    # U-test subcommand (called by SLURM jobs)
    utest_p = sub.add_parser('utest', help='Run U-test analysis (used by SLURM jobs)')
    utest_p.add_argument('--sgrna-summary', required=True, help='MAGeCK sgrna_summary.txt file')
    utest_p.add_argument('--output-prefix', required=True, help='Output file prefix')
    utest_p.add_argument('--fdr', type=float, default=0.05, help='FDR threshold')
    utest_p.add_argument('--min-treat-mean', type=int, default=20, help='Min treatment mean')
    
    # Threshold subcommand (called by SLURM jobs)
    thresh_p = sub.add_parser('threshold', help='Threshold counts by replicate pairs (used by SLURM jobs)')
    thresh_p.add_argument('--counts', required=True, help='Counts file (MAGeCK format)')
    thresh_p.add_argument('--treatment', required=True, help='Treatment sample(s), comma-separated')
    thresh_p.add_argument('--control', required=True, help='Control sample(s), comma-separated')
    thresh_p.add_argument('--min-counts', type=int, default=50, help='Min combined counts per replicate pair')
    thresh_p.add_argument('--output', required=True, help='Output thresholded counts file')
    
    # Line-score U-test subcommand (two-comparison analysis)
    line_p = sub.add_parser('line-score-utest', 
        help='Two-comparison analysis: sgRNA line score with U-test on top N')
    line_p.add_argument('--sgrna1', required=True, help='First MAGeCK sgrna_summary.txt')
    line_p.add_argument('--sgrna2', required=True, help='Second MAGeCK sgrna_summary.txt')
    line_p.add_argument('--line', required=True, 
        help='Reference line as "slope,intercept" (e.g., "1,0" for y=x)')
    line_p.add_argument('--output', required=True, help='Output prefix for result files')
    line_p.add_argument('--names', help='Comma-separated comparison names')
    line_p.add_argument('--top-n', type=int, default=3, 
        help='Top N sgRNAs per gene for U-test (default: 3)')
    line_p.add_argument('--fdr', type=float, default=0.05, help='FDR threshold')
    line_p.add_argument('--min-treat-mean', type=int, default=20, help='Min treatment mean')
    line_p.add_argument('--score-mode', choices=['perp_weighted', 'parallel'],
        default='perp_weighted',
        help="Scoring: 'perp_weighted' (default) downweights guides off the line; "
             "'parallel' uses raw distance along the line")
    
    # Volcano plot subcommand (manual use)
    volcano_p = sub.add_parser('volcano', help='Generate volcano plot from results')
    volcano_p.add_argument('--results', required=True, help='Path to _all_genes.csv file')
    volcano_p.add_argument('--hits', required=True, help='Path to _hits_fdr*.csv file')
    volcano_p.add_argument('--output', required=True, help='Output PDF path')
    volcano_p.add_argument('--title', default='', help='Plot title')
    volcano_p.add_argument('--fdr', type=float, default=0.05)
    volcano_p.add_argument('--label-genes', default='', help='Comma-separated genes to label')
    volcano_p.add_argument('--label-n-top', type=int, default=10, help='Label top N hits')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    # Commands that don't need config
    if args.command == 'utest':
        run_utest_standalone(
            args.sgrna_summary,
            args.output_prefix,
            args.fdr,
            args.min_treat_mean
        )
        return
    
    if args.command == 'threshold':
        run_threshold_standalone(
            args.counts,
            args.treatment,
            args.control,
            args.min_counts,
            args.output
        )
        return
    
    if args.command == 'line-score-utest':
        # Parse line parameters
        line_parts = args.line.split(',')
        slope = float(line_parts[0])
        intercept = float(line_parts[1])
        
        # Parse names
        names = [None, None]
        if args.names:
            names = [s.strip() for s in args.names.split(',')]
        
        sgrna_line_score_utest(
            args.sgrna1,
            args.sgrna2,
            slope,
            intercept,
            args.output,
            comp1_name=names[0],
            comp2_name=names[1] if len(names) > 1 else None,
            top_n=args.top_n,
            fdr=args.fdr,
            min_count=args.min_treat_mean,
            score_mode=args.score_mode
        )
        return
    
    if args.command == 'volcano':
        df_all = pd.read_csv(args.results)
        df_hits = pd.read_csv(args.hits)
        # Estimate threshold from hits
        if len(df_hits) > 0:
            threshold = abs(df_hits['product']).min()
        else:
            threshold = 0
        label_genes = [g.strip() for g in args.label_genes.split(',') if g.strip()] if args.label_genes else None
        plot_volcano(df_all, df_hits, threshold, args.output, 
                    title=args.title, fdr=args.fdr, 
                    label_genes=label_genes, label_n_top=args.label_n_top)
        print(f"Saved volcano plot to {args.output}")
        return
    
    # Commands that need config
    config = load_config(args.config)
    if hasattr(args, 'dry_run') and args.dry_run:
        config.dry_run = True
    
    if args.command == 'run':
        run_pipeline(config, args.wait)
    elif args.command == 'count':
        sample_table = load_sample_table(config.sample_table)
        library_table = load_library_table(config.library_table)
        submit_counting_jobs(sample_table, library_table, config)
    elif args.command == 'aggregate':
        sample_table = load_sample_table(config.sample_table)
        library_table = load_library_table(config.library_table)
        library = aggregate_counts(sample_table, library_table, config)
        # Also write counts_for_mageck.txt, as run_pipeline() does. Without this
        # a staged run (count -> aggregate -> mageck) scores against whatever
        # counts_for_mageck.txt happens to be on disk, silently skipping the
        # NTC relabelling in prepare_mageck_input().
        sample_names = list(sample_table['sample_name'].unique())
        counts_file = os.path.join(config.output_dir, 'counts_for_mageck.txt')
        prepare_mageck_input(library, sample_names, counts_file)
    elif args.command == 'mageck':
        comparisons = load_comparisons_table(config.comparisons_table)
        counts_file = os.path.join(config.output_dir, 'counts_for_mageck.txt')
        submit_mageck_jobs(comparisons, counts_file, config)
    elif args.command == 'results':
        comparisons = load_comparisons_table(config.comparisons_table)
        process_results(comparisons, config)
    elif args.command == 'status':
        for tf in ['counting_jobs.txt', 'mageck_jobs.txt']:
            p = os.path.join(config.output_dir, tf)
            if os.path.exists(p):
                with open(p) as f:
                    jids = [l.strip() for l in f if l.strip()]
                if jids:
                    print(f"\n{tf}:")
                    _, status = check_jobs_complete(jids)
                    for j, s in status.items():
                        print(f"  {j}: {s}")


if __name__ == '__main__':
    main()
