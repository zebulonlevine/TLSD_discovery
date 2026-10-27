"""Write the UniProt accession lists that uniprot_topology_pipeline.py takes as input.

One accession per protein, in table order ('sp|P12931-2|SRC_HUMAN' -> 'P12931-2'):
  silac_uniprotIDs.csv   SILAC table (contaminants and decoys removed)
  total_uniprotIDS.csv   total-protein table

    python proteomics/annotation/export_accessions.py --outdir <dir>
"""
import argparse
import os
import sys

import pandas as pd

_PROT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PROT)
sys.path.insert(0, os.path.dirname(_PROT))
from common import paths  # noqa: E402
import silac_proteomics as sp  # noqa: E402
from settings import SILAC_TSV, TOTAL_CSV  # noqa: E402
from total_abundance import uniprot_accession  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    silac = sp.drop_contaminants_and_decoys(sp.load_protein_quant(paths.repo(SILAC_TSV)))
    total = pd.read_csv(paths.repo(TOTAL_CSV))
    for name, refs in (('silac_uniprotIDs.csv', silac['Protein Id']), ('total_uniprotIDS.csv', total.iloc[:, 0])):
        uniprot_accession(refs).to_csv(os.path.join(a.outdir, name), index=False, header=False)
        print(f'{name}: {len(refs)} accessions')


if __name__ == '__main__':
    main()
