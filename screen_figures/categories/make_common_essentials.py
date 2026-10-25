#!/usr/bin/env python3
"""Build CommonEssentials.csv (Fig S1E/F highlight, Table S2 CommonEssential column)
from the DepMap Achilles common-essential gene list.

Input (not committed; download from DepMap, https://depmap.org/portal/download/):
    ${TLSD_DATA}/DepMap/Achilles_common_essentials.csv   one column 'gene', e.g. "AAAS (8086)"
Output: screen_figures/categories/CommonEssentials.csv  (gene, category)

Converted from Screen_Visualization/Categories/depMapParse.ipynb.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from common import paths  # noqa: E402

df = pd.read_csv(paths.inp('DepMap', 'Achilles_common_essentials.csv'))
essentials = df.gene.apply(lambda x: x.split(' (')[0])
pd.DataFrame({'gene': essentials, 'category': ['DepMap Common Essential'] * len(essentials)}).to_csv(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CommonEssentials.csv'), index=False)
