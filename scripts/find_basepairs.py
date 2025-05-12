#!/usr/bin/env python
"""
This script takes as first argument a PDB ID and searches for basepairs.
"""

import sys
from itertools import product
from pyrna.db import PDB
from pyrna.parsers import parse_pdb
from pyrna.model import DonorEndoA, DonorExoA, AcceptorEndoA, AcceptorExoA

def find_bps(pdb_id):
    pdb = PDB()
    content = pdb.get_entry(pdb_id)
    for ts in parse_pdb(content):
        for i in range(0, len(ts.residues)-1):
            r = ts.residues[i]
            donors = [a for a in r.atoms if isinstance(a, DonorEndoA) or isinstance(a, DonorExoA)]
            acceptors = [a for a in r.atoms if isinstance(a, AcceptorEndoA) or isinstance(a, AcceptorExoA)]
            for j in range(i+1, len(ts.residues)):
                next_r = ts.residues[j]
                print(f"{r}{i} - {next_r}{j}")
                next_acceptors = [a for a in next_r.atoms if isinstance(a, AcceptorEndoA) or isinstance(a, AcceptorExoA)]
                hbonds = list(product(donors, next_acceptors))
                for hbond in hbonds:
                    print(f"{hbond[0]} ~ {hbond[1]}")
                next_donors = [a for a in next_r.atoms if isinstance(a, DonorEndoA) or isinstance(a, DonorExoA)]
                hbonds = list(product(acceptors, next_donors))
                for hbond in hbonds:
                    print(f"{hbond[0]} ~ {hbond[1]}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: find_basepairs.py PDB_ID (try for example: find_basepairs.py 1HR2)")
        sys.exit()
    find_bps(sys.argv[1])