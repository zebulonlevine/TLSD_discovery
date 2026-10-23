#!/usr/bin/env bash
# Run MAGeCK test inside a named conda env.
#
# Args:
#   $1  counts.txt        MAGeCK-format count table (sgRNA, gene, samples...)
#   $2  output_prefix     full prefix (incl. directory) for MAGeCK outputs
#   $3  ctrl_ids          comma-separated control sample names matching the count header
#   $4  trt_ids           comma-separated treatment sample names matching the count header
#   $5  conda_env         conda env name where mageck is installed
#
# MAGeCK normalization uses non-targeting sgRNAs as the size-factor anchors.
# LFCs in sgrna_summary.txt come from median-normalized counts and are what
# the gene-level MWU pipeline (02) consumes.

set -euo pipefail

COUNTS="$1"
PREFIX="$2"
CTRL_IDS="$3"
TRT_IDS="$4"
CONDA_ENV="${5:-mageck}"

if [[ ! -f "$COUNTS" ]]; then
    echo "ERROR: count file not found: $COUNTS" >&2
    exit 1
fi

# Build NTC sgRNA list (control sgRNAs for MAGeCK)
NTC_LIST="$(dirname "$PREFIX")/ntc_sgrnas.txt"
awk -F'\t' 'NR>1 && $2=="non-targeting"{print $1}' "$COUNTS" > "$NTC_LIST"
N_NTC=$(wc -l < "$NTC_LIST")
echo "[mageck] using $N_NTC non-targeting sgRNAs as normalization controls"

# Activate conda env
# (works whether or not the user already has conda initialized)
if [[ -z "${CONDA_EXE:-}" ]]; then
    if command -v conda >/dev/null 2>&1; then
        CONDA_EXE="$(command -v conda)"
    else
        echo "ERROR: conda not on PATH. Please initialize conda or pass an env in PATH." >&2
        exit 1
    fi
fi
CONDA_BASE="$("$CONDA_EXE" info --base)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# Verify mageck is available
if ! command -v mageck >/dev/null 2>&1; then
    echo "ERROR: mageck not found in env '$CONDA_ENV'. Install with:" >&2
    echo "       conda install -c bioconda mageck" >&2
    exit 1
fi
echo "[mageck] using $(which mageck) ($(mageck --version 2>&1 | head -1))"

mkdir -p "$(dirname "$PREFIX")"

# Run MAGeCK test
#   --norm-method median        : median ratio normalization (matches our pipeline)
#   --gene-lfc-method median    : gene-level LFC aggregation (not used downstream but harmless)
#   --control-sgrna             : use NTC list for null model fitting
#   --adjust-method fdr         : BH FDR
mageck test \
    --count-table "$COUNTS" \
    --control-id "$CTRL_IDS" \
    --treatment-id "$TRT_IDS" \
    --output-prefix "$PREFIX" \
    --norm-method median \
    --gene-lfc-method median \
    --control-sgrna "$NTC_LIST" \
    --adjust-method fdr

echo
echo "[mageck] outputs:"
echo "  ${PREFIX}.sgrna_summary.txt    <- consumed by 02_gene_scores.py"
echo "  ${PREFIX}.gene_summary.txt     (RRA gene scores, not used downstream)"
echo "  ${PREFIX}.normalized.txt       (median-normalized counts)"
