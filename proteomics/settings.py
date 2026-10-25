"""Shared inputs and constants for the proteomics analyses.

Input files are read from the data root (common/paths.py) under proteomics/:

  protein_quant_30972.tsv        SILAC-TMT protein quantification (RPS3-HaloTag line; 12-plex x heavy/light)
  all_ttest_results.csv          total-abundance TMT protein table (15 samples)
  mmc4_published_OpenCell.xlsx   organelle-IP localization table (sheet 'subcell_annotation_summary'), ref. 170
  silac_topology.tsv,            UniProt membrane topology + MitoCarta 3.0 flags for each table's proteins,
  total_topology.tsv             written by annotation/uniprot_topology_pipeline.py
"""

# ---------------------------------------------------------------------------
# Input files (relative to the data root)
# ---------------------------------------------------------------------------
SILAC_TSV = 'proteomics/protein_quant_30972.tsv'
TOTAL_CSV = 'proteomics/all_ttest_results.csv'
LOCALIZATION_XLSX = 'proteomics/mmc4_published_OpenCell.xlsx'
LOCALIZATION_SHEET = 'subcell_annotation_summary'
SILAC_TOPOLOGY = 'proteomics/silac_topology.tsv'
TOTAL_TOPOLOGY = 'proteomics/total_topology.tsv'

OUT_DIR = 'proteomics'          # under the output root

# ---------------------------------------------------------------------------
# SILAC-TMT plex (RPS3-HaloTag cells; heavy = old protein, light = new protein)
# ---------------------------------------------------------------------------
BAIT = 'RPS3'
CHANNELS = ['126', '127n', '127c', '128n', '128c', '129n', '129c', '130n', '130c', '131n', '131c', '132n']
SAMPLE_NAMES = [
    'NT_3d_1', 'NT_3d_2', 'NT_3d_3',
    'sgDDRGK1_3d_1', 'sgDDRGK1_3d_2', 'sgDDRGK1_3d_3',
    'NT_8d_1', 'NT_8d_2', 'NT_8d_3',
    'sgDDRGK1_8d_1', 'sgDDRGK1_8d_2', 'sgDDRGK1_8d_3',
]
LABELS = ('H', 'L')
RAW_LABEL_TOKENS = {'H': 'heavy', 'L': 'light'}
RAW_TEMPLATE = '{label_token}~rq_{channel}_sn sum'
SAMPLE_TEMPLATE = BAIT + '_{sample}_{label}'

# day-8 / day-3 ratio per condition and label: (out_name, condition, late, early, label)
TURNOVER_COMPARISONS = [
    ('turnover_old_NT', 'NT', '8d', '3d', 'H'),
    ('turnover_old_sgDDRGK1', 'sgDDRGK1', '8d', '3d', 'H'),
    ('turnover_new_NT', 'NT', '8d', '3d', 'L'),
    ('turnover_new_sgDDRGK1', 'sgDDRGK1', '8d', '3d', 'L'),
]
# sgDDRGK1 - sgNT on the log2 turnover ratios
DELTA_COMPARISONS = [
    ('delta_old_log2', 'turnover_old_sgDDRGK1_log2', 'turnover_old_NT_log2'),
    ('delta_new_log2', 'turnover_new_sgDDRGK1_log2', 'turnover_new_NT_log2'),
]
T_TEST_DF = 4                   # 3 + 3 replicates - 2

# ---------------------------------------------------------------------------
# Total-abundance TMT plex
# ---------------------------------------------------------------------------
# Column names given to all_ttest_results.csv on load. The table arrives with its own
# fold changes and tests (columns 3-50, unused except the three log2FC columns in Fig S5C);
# the analysis recomputes the comparisons from the replicate columns.
TOTAL_COLUMNS = [
    'Reference', 'Gene', 'Annotation',
    'p.val_S3_sgDDRGK1_growing-S3_sgNT_growing', 'log2FC_S3_sgDDRGK1_growing-S3_sgNT_growing', 'S3_sgNT_growing',
    'q.val_S3_sgDDRGK1_growing-S3_sgNT_growing', 'sig_S3_sgDDRGK1_growing-S3_sgNT_growing',
    'p.val_S3_sgNT_quiescent-S3_sgNT_growing', 'log2FC_S3_sgNT_quiescent-S3_sgNT_growing', 'S3_sgDDRGK1_growing',
    'q.val_S3_sgNT_quiescent-S3_sgNT_growing', 'sig_S3_sgNT_quiescent-S3_sgNT_growing', 'S3_sgNT_quiescent',
    'p.val_S3_sgNT_quiescent-S3_sgDDRGK1_growing', 'log2FC_S3_sgNT_quiescent-S3_sgDDRGK1_growing',
    'q.val_S3_sgNT_quiescent-S3_sgDDRGK1_growing', 'sig_S3_sgNT_quiescent-S3_sgDDRGK1_growing',
    'p.val_anova', 'q.val_anova', 'anova_sig',
    'p.val_S3_sgDDRGK1_quiescent-S3_sgNT_growing', 'log2FC_S3_sgDDRGK1_quiescent-S3_sgNT_growing', 'S3_sgDDRGK1_quiescent',
    'q.val_S3_sgDDRGK1_quiescent-S3_sgNT_growing', 'sig_S3_sgDDRGK1_quiescent-S3_sgNT_growing',
    'p.val_S3_sgDDRGK1_quiescent-S3_sgDDRGK1_growing', 'log2FC_S3_sgDDRGK1_quiescent-S3_sgDDRGK1_growing',
    'q.val_S3_sgDDRGK1_quiescent-S3_sgDDRGK1_growing', 'sig_S3_sgDDRGK1_quiescent-S3_sgDDRGK1_growing',
    'p.val_S3_sgDDRGK1_quiescent-S3_sgNT_quiescent', 'log2FC_S3_sgDDRGK1_quiescent-S3_sgNT_quiescent',
    'q.val_S3_sgDDRGK1_quiescent-S3_sgNT_quiescent', 'sig_S3_sgDDRGK1_quiescent-S3_sgNT_quiescent',
    'p.val_L29_sgDDRGK1_quiescent-S3_sgNT_growing', 'log2FC_L29_sgDDRGK1_quiescent-S3_sgNT_growing', 'L29_sgDDRGK1_quiescent',
    'q.val_L29_sgDDRGK1_quiescent-S3_sgNT_growing', 'sig_L29_sgDDRGK1_quiescent-S3_sgNT_growing',
    'p.val_L29_sgDDRGK1_quiescent-S3_sgDDRGK1_growing', 'log2FC_L29_sgDDRGK1_quiescent-S3_sgDDRGK1_growing',
    'q.val_L29_sgDDRGK1_quiescent-S3_sgDDRGK1_growing', 'sig_L29_sgDDRGK1_quiescent-S3_sgDDRGK1_growing',
    'p.val_L29_sgDDRGK1_quiescent-S3_sgNT_quiescent', 'log2FC_L29_sgDDRGK1_quiescent-S3_sgNT_quiescent',
    'q.val_L29_sgDDRGK1_quiescent-S3_sgNT_quiescent', 'sig_L29_sgDDRGK1_quiescent-S3_sgNT_quiescent',
    'p.val_L29_sgDDRGK1_quiescent-S3_sgDDRGK1_quiescent', 'log2FC_L29_sgDDRGK1_quiescent-S3_sgDDRGK1_quiescent',
    'q.val_L29_sgDDRGK1_quiescent-S3_sgDDRGK1_quiescent', 'sig_L29_sgDDRGK1_quiescent-S3_sgDDRGK1_quiescent', 'cluster',
    'S3_sgNT_growing_1', 'S3_sgNT_growing_2', 'S3_sgNT_growing_3',
    'S3_sgDDRGK1_growing_1', 'S3_sgDDRGK1_growing_2', 'S3_sgDDRGK1_growing_3',
    'S3_sgNT_quiescent_1', 'S3_sgNT_quiescent_2', 'S3_sgNT_quiescent_3',
    'S3_sgDDRGK1_quiescent_1', 'S3_sgDDRGK1_quiescent_2', 'S3_sgDDRGK1_quiescent_3',
    'L29_sgDDRGK1_quiescent_1', 'L29_sgDDRGK1_quiescent_2', 'L29_sgDDRGK1_quiescent_3',
]
TOTAL_REPLICATES = TOTAL_COLUMNS[-15:]          # log2 TMT abundances, one column per sample
# (numerator, denominator): sgDDRGK1 vs sgNT in growing and quiescent RPS3 cells, and quiescent
# RPL29 sgDDRGK1 vs RPS3 sgNT (no RPL29 sgNT sample was collected)
TOTAL_COMPARISONS = [
    ('S3_sgDDRGK1_growing', 'S3_sgNT_growing'),
    ('S3_sgDDRGK1_quiescent', 'S3_sgNT_quiescent'),
    ('L29_sgDDRGK1_quiescent', 'S3_sgNT_quiescent'),
]
QUIESCENT_FC = 'S3_sgDDRGK1_quiescent_vs_S3_sgNT_quiescent_log2FC'
GROWING_FC = 'S3_sgDDRGK1_growing_vs_S3_sgNT_growing_log2FC'

# ---------------------------------------------------------------------------
# Localization categories
# ---------------------------------------------------------------------------
KEEP_CATEGORIES = ['RPL', 'RPS', 'proteasome', 'nucleus', 'mitochondrion', 'ER', 'cytosol',
                   'lysosome', 'Golgi', 'plasma_membrane']

RPL_GENES = '''RPL10A RPL8 RPL3 RPL4 RPL11 RPL9 RPL6 RPL7A RPLP0 RPL12 RPL13A RPL13
RPL23 RPL14 RPL27A RPL15 RPL10 RPL5 RPL18 RPL19 RPL18A RPL21 RPL17 RPL22 RPL23A
RPL26 RPL24 RPL27 RPL28 RPL35 RPL29 RPL7 RPL30 RPL31 RPL32 RPL35A RPL34 RPL36
RPL37 RPL38 RPL39 RPL40 RPL41 RPL36A RPL37A RPLP1 RPLP2'''.split()

RPS_GENES = '''RPS3A RPSA RPS3 RPS9 RPS4X RPS2 RPS6 RPS5 RPS7 RPS15A RPS8 RPS16
RPS20 RPS10 RPS14 RPS23 RPS12 RPS18 RPS29 RPS13 RPS11 RPS17 RPS15 RPS19 RPS21
RPS24 RPS25 RPS26 RPS27 RPS28 RPS30 RPS27A RACK1'''.split()

CATEGORY_OVERRIDES = {'RPL': RPL_GENES, 'RPS': RPS_GENES}     # applied after the localization map
DESCRIPTION_RULES = {'collagen': 'ollagen'}                   # case-insensitive description match

# Adjusted categories (the ones in the figures and tables): nucleolus counts as nucleus;
# stress granule, P-body and 'translation' count as cytosol; no collagen description rule.
CYTO_CATEGORIES = ['cytosol', 'stress_granule', 'p-body', 'translation']
NUC_CATEGORIES = ['nucleus', 'nucleolus']
ADJ_CATEGORY = 'category_locAdj'

# subsampled Mann-Whitney U: category vs rest of the proteome, n <= 50 per draw, 1,000 draws,
# median p, BH across categories (seeded inside silac_proteomics.mann_whitney_subsampled)
MWU = dict(alternative='two-sided', fdr=True, max_n=50, n_resamples=1000)

# ---------------------------------------------------------------------------
# Membrane topology (total-abundance data; Fig 4E, S6C)
# ---------------------------------------------------------------------------
TOPOLOGY_CATEGORIES = ['soluble (cytosolic)', 'soluble (secreted/lumenal)', 'single-pass TM', 'multi-pass TM',
                       'TA', 'soluble (mito)', 'multi-pass TM (mito)', 'single-pass TM (mito)', 'TA (mito)']
TOPOLOGY_MAP = {   # uniprot_topology_pipeline class (+ '_mito' if in MitoCarta 3.0) -> simplified label
    'soluble': 'soluble (cytosolic)',
    'secreted_or_lumenal': 'soluble (secreted/lumenal)',
    'soluble_mito': 'soluble (mito)', 'secreted_or_lumenal_mito': 'soluble (mito)',
    'single_TM_unclear': 'single-pass TM', 'TypeIII_or_unusual': 'single-pass TM',
    'Type_I_or_SP_single_TM': 'single-pass TM', 'signal_anchored_TypeII': 'single-pass TM',
    'single_TM_unclear_mito': 'single-pass TM (mito)', 'TypeIII_or_unusual_mito': 'single-pass TM (mito)',
    'Type_I_or_SP_single_TM_mito': 'single-pass TM (mito)', 'signal_anchored_TypeII_mito': 'single-pass TM (mito)',
    'SP_plus_multipass': 'multi-pass TM', 'multipass_no_SP': 'multi-pass TM',
    'SP_plus_multipass_mito': 'multi-pass TM (mito)', 'multipass_no_SP_mito': 'multi-pass TM (mito)',
    'tail_anchored': 'TA', 'tail_anchored_mito': 'TA (mito)',
    'beta_barrel': 'beta_barrel', 'beta_barrel_mito': 'beta_barrel (mito)',
}

# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------
# highlighted in the Fig S6D/E volcanoes (Table S9 column isTransloconOrUPR)
TRANSLOCON = ['SEC61B', 'SEC61A1', 'SEC61G', 'SEC61A2']
UPR_ER = ['HSPA5', 'DNAJB9', 'HERPUD1', 'HYOU1', 'PDIA4', 'PDIA6', 'SDF2L1', 'CALR', 'HSP90B1']

# Fig S5F strip plot (Table S8 column ER_GeneSets_delta)
GENE_SETS = {
    'Sec61': ['SEC61A1', 'SEC61B', 'SEC61G', 'SEC61A2'],
    'TRAP': ['SSR1', 'SSR2', 'SSR3', 'SSR4'],
    'OST (core)': ['RPN1', 'RPN2', 'DDOST', 'DAD1', 'TMEM258', 'OST4', 'KRTCAP2'],
    'OST (STT3A/B)': ['STT3A', 'OSTC', 'STT3B', 'MAGT1', 'TUSC3'],
    'Sec62/63': ['SEC62', 'SEC63'],
    'PAT': ['CCDC47', 'WDR83OS'],                                   # Asterix = WDR83OS
    'BOS': ['NCLN', 'NOMO1', 'NOMO2', 'NOMO3', 'TMEM147'],
    'GEL': ['TMCO1', 'RAB5IF'],                                     # RAB5IF = C5orf15
    'EMC': ['EMC1', 'EMC2', 'EMC3', 'EMC4', 'MMGT1', 'EMC6', 'EMC7', 'EMC8', 'EMC9', 'EMC10'],   # MMGT1 = EMC5
}
DELTA_SETS = {
    'Sec61 (core)': GENE_SETS['Sec61'],
    'translocon-associated': (GENE_SETS['TRAP'] + GENE_SETS['Sec62/63'] + GENE_SETS['PAT']
                              + GENE_SETS['BOS'] + GENE_SETS['GEL']),
    'OST': GENE_SETS['OST (core)'] + GENE_SETS['OST (STT3A/B)'],
    'EMC': GENE_SETS['EMC'],
    '40S (RPS)': RPS_GENES,
    '60S (RPL)': RPL_GENES,
}
BACKGROUND = {'All proteins': None, 'Cytosol': 'cytosol', 'ER': 'ER'}
SET_COLORS = {'40S (RPS)': '#3C7DA6', '60S (RPL)': '#3C7DA6'}

# Fig 7C
FOCUS_GENE = 'R3HDM4'           # TLSD1
