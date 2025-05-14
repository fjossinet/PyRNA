#!/usr/bin/env python
"""
This script takes as first argument a PDB ID and searches for basepairs.
"""

import sys
from pyrna.db import PDB
from pyrna.parsers import parse_pdb

def find_bps(pdb_id):
    for ts in parse_pdb(PDB().get_entry(pdb_id)):
        ts.find_canonical_basepairs()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: find_basepairs.py PDB_ID (try for example: find_basepairs.py 1HR2)")
        sys.exit()
    find_bps(sys.argv[1])