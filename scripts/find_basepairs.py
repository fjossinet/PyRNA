#!/usr/bin/env python
"""
This script takes as first argument a PDB ID and searches for basepairs.
"""

import sys
from pyrna.db import PDB
from pyrna.parsers import parse_pdb

def find_bps(pdb_id, chain_name):
    for ts in parse_pdb(PDB().get_entry(pdb_id)):
        for bp in ts.find_canonical_basepairs(chain_name):
            print(f"{bp.location.start()} - {bp.location.end()}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: find_basepairs.py PDB_ID chain_name (try for example: find_basepairs.py 1HR2 A)")
        sys.exit()
    find_bps(sys.argv[1], sys.argv[2])