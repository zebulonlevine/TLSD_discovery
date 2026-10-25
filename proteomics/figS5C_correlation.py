"""Fig S5C: correlation between the SILAC and total-protein comparisons (sgDDRGK1 vs sgNT).

Proteins quantified in both experiments (matched on accession and gene) with both turnover deltas;
Pearson correlation of nine comparisons, clustered.

  SILAC:  sgDDRGK1 / sgNT ratio of total-normalized signal, day 3 and day 8, heavy (old) and light
          (new); delta_old_log2, delta_new_log2 (silac_turnover.py)
  total:  log2FC columns supplied with all_ttest_results.csv (growing; quiescent; RPL29 sgDDRGK1
          quiescent vs RPS3 sgNT quiescent)

Note: the four SILAC sgDDRGK1 / sgNT columns are the linear ratios (named log2_* in the notebook);
the published panel was drawn from these.

Outputs (output root, proteomics/): figS5C_correlation_heatmap.pdf, figS5C_correlation.csv

    python proteomics/figS5C_correlation.py
"""
import os
import sys
import warnings

import pandas as pd
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402
from settings import OUT_DIR  # noqa: E402
from silac_turnover import build_silac, load_localization  # noqa: E402
from total_abundance import build_total  # noqa: E402

SILAC_COLS = ['log2_sgDDRGK1_vs_NT_8d_L', 'log2_sgDDRGK1_vs_NT_3d_L',
              'log2_sgDDRGK1_vs_NT_8d_H', 'log2_sgDDRGK1_vs_NT_3d_H',
              'delta_old_log2', 'delta_new_log2']
TOTAL_COLS = ['log2FC_S3_sgDDRGK1_growing-S3_sgNT_growing',
              'log2FC_S3_sgDDRGK1_quiescent-S3_sgNT_quiescent',
              'log2FC_L29_sgDDRGK1_quiescent-S3_sgNT_quiescent']
KEY = ['Reference', 'Gene']


def correlation(silac, total):
    combined = pd.merge(total[KEY + TOTAL_COLS], silac[KEY + SILAC_COLS], on=KEY, how='inner')
    combined = combined.dropna(subset=['delta_old_log2', 'delta_new_log2'])
    return combined[SILAC_COLS + TOTAL_COLS].corr(), len(combined)


def main():
    warnings.filterwarnings('ignore', message='After omitting NaNs')   # SEM of proteins seen in < 2 replicates = NaN
    loc_df = load_localization()
    corr, n = correlation(build_silac(loc_df), build_total(loc_df, topology=False))
    corr.to_csv(paths.outp(OUT_DIR, 'figS5C_correlation.csv', mkdir=True))
    sns.clustermap(corr, cmap='vlag', vmin=-1, vmax=1, center=0)
    plt.savefig(paths.outp(OUT_DIR, 'figS5C_correlation_heatmap.pdf'), dpi=300, bbox_inches='tight')
    plt.close('all')
    print(f'{n} proteins in both experiments; wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
