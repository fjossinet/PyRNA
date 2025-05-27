#!/usr/bin/env python
"""
This script takes as arguments the data source and a chain name. 
It prints the helices from the 2D infered from the 3D coordinates
"""

import sys
from pyrna.db import PDBProvider, PDBProviderParams
from pyrna.parsers import parse_pdb, PDBFileProvider, PDBFileProviderParams
from pyrna.model import SecondaryStructure

def find_helices(pdb_id, chain_name):
    ss = None
    if "/" in pdb_id:
        params = PDBFileProviderParams()
        params.set_input_file(pdb_id)
        params.set_chain_name(chain_name)
        ss = PDBFileProvider().get_secondary_structure(params)
    elif len(pdb_id) == 4 :
        params = PDBProviderParams()
        params.set_pdb_id(pdb_id)
        params.set_chain_name(chain_name)
        ss = PDBProvider().get_secondary_structure(params)
    else:
        print("Data source unknown")
    if ss:
        for h in ss.helices:
            print(f"{h.name}")
            print(f"strand1 : {h.location.blocks[0].start}-{h.location.blocks[0].end}")
            print(f"strand2 : {h.location.blocks[1].start}-{h.location.blocks[1].end}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: find_helices.py data_source chain_name (try for example: find_helices.py 1HR2 A)")
        sys.exit()
    find_helices(sys.argv[1], sys.argv[2])