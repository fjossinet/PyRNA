#!/usr/bin/env python
"""
This script takes as first argument a 2D structure described in a VIENNA file and prints it using the RNAML format
"""

import sys, os
from pyrna.parsers import parse_vienna

def to_rnaml(vienna_file):
    with open(os.path.join(vienna_file)) as f:
        rna, bps = parse_vienna(f.read())
        for i in range(0,len(bps)):
            print("""\
<rnaml version="1.1">
  <molecule id="{rna_name}">
    <identity>
      <name>{rna_name}</name>
    </identity>
    <sequence length="">
        <seq-data>
        {rna_sequence}
        </seq-data>
    </sequence>
    <structure>
        <model id="{model_id}">		        				    
        <str-annotation>""".format(rna_name = rna[i].name, rna_sequence = rna[i].sequence, model_id = i+1))
            for bp in bps[i]:
                print("""\
            <base-pair">
                <base-id-5p>
                    <base-id>
                        <position>{bp_start}</position>
                    </base-id>
                </base-id-5p>
                <base-id-3p>
                    <base-id>
                        <position>{bp_end}</position>
                    </base-id>
                </base-id-3p>
                <edge-5p>W</edge-5p>
                <edge-3p>W</edge-3p>
                <bond-orientation>cis</bond-orientation>
            </base-pair>""".format(bp_start = bp.location.start(), bp_end = bp.location.end()))
            print("""\
        </str-annotation>	
        </model>
    </structure>		
  </molecule>			
</rnaml>""")
            
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: vienna_to_rnaml.py vienna_file (try for example: ./vienna_to_rnaml.py ./data/sample.vienna)")
        sys.exit()
    to_rnaml(sys.argv[1])