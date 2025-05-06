#!/usr/bin/env python
"""
This script takes as first argument a PDB ID, extracts the RNA sequences and prints them using the FASTA format
"""

import sys
from pyrna.db import PDB
from pyrna.parsers import parse_pdb, to_fasta

def fetch(pdb_id):
    pdb = PDB()
    content = pdb.get_entry(pdb_id)
    print(to_fasta([ts.rna for ts in parse_pdb(content)]))  

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: pdb2_fasta.py PDB_ID (try for example: pdb2_fasta.py 1HR2)")
        sys.exit()
    fetch(sys.argv[1])