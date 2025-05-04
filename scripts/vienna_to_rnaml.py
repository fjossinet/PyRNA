#!/usr/bin/env python
"""
This script takes as first argument a 2D structure described in a VIENNA file and prints it using the RNAML format
"""

import sys, os
from pyrna.parsers import parse_vienna

def to_rnaml(vienna_file):
    with open(os.path.join(working_dir,f)) as f:
        print(parse_vienna(f.read()))
    
if __name__ == "__main__":
    to_rnaml(sys.argv[1])