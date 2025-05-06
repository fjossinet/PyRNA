#!/usr/bin/env python
"""
This script takes as first argument a 2D structure described in a VIENNA file and prints it using the RNAML format
"""

import sys, os
from pyrna.parsers import parse_vienna, to_rnaml

def convert(vienna_file):
    with open(os.path.join(vienna_file)) as f:
        rna, bps = parse_vienna(f.read())
        print(to_rnaml(rna, bps)) 
            
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: vienna_to_rnaml.py vienna_file (try for example: ./vienna_to_rnaml.py ./data/sample.vienna)")
        sys.exit()
    convert(sys.argv[1])