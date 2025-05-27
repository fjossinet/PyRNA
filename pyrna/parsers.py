import re
from pyrna.model import RNA, BasePair, TertiaryStructure

def to_pdb(tertiary_structure, location = None):
    """
    Convert a TertiaryStructure object into PDB data

    Parameters:
    -----------
    - tertiary_structure: a TertiaryStructure object (see pyrna.model)
    - location (default: None): a Location object (see pyrna.model). Restrict the export to the atoms of the residues enclosed by this location.

    Returns:
    ------
    the PDB data as a String
    """
    lines= []
    i = 1
    keys = []

    for k in tertiary_structure.residues:
        if location and location.has_position(k) or not location:
            keys.append(k)

    keys.sort() #the absolute position are sorted

    for key in keys:
        atoms = tertiary_structure.residues[key]['atoms']
        for atom in atoms:
            lines.append("%-6s%5u  %-4s%3s %s%4u    %8.3f%8.3f%8.3f"%("ATOM", i, atom['name'], tertiary_structure.rna.sequence[key-1], tertiary_structure.rna.name[0], key, atom['coords'][0], atom['coords'][1], atom['coords'][2]))
            i += 1

    lines.append("END")

    return '\n'.join(lines)

def to_fasta(molecules, single_line=False):
    """
    Convert a list of Molecule objects into FASTA data

    Parameters:
    ---------
    - molecules: a list of Molecule objects (see pyrna.model)
    - single_line (default: False): if True, each molecular sequence will we exported into a single line

    Returns:
    ------
    the FASTA data as a String
    """
    outputs = []
    for molecule in molecules:
        outputs.append(molecule.to_fasta(single_line))
    return '\n'.join(outputs)

def to_rnaml(rna, bps):
    """
    Convert an RNA object and its list of base-pairs into RNAML data

    Parameters:
    -----------
    - rna: an RNA object (see pyrna.model)
    - bps: a list of BasePair objects

    Returns:
    ------
    the RNAML data as a String
    """
    for i in range(0,len(bps)):
        data = """\
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
            <str-annotation>""".format(rna_name = rna[i].name, rna_sequence = rna[i].sequence, model_id = i+1)
        for bp in bps[i]:
            data += """
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
                </base-pair>""".format(bp_start = bp.location.start(), bp_end = bp.location.end())
            
        data += """
            </str-annotation>	
        </model>
    </structure>		
  </molecule>			
</rnaml>"""
    return data


def parse_fasta(fasta_data):
    """
    Parse FASTA data

    Parameters:
    ---------
    - fasta_data: the Fasta data as a String

    Returns:
    ------
    a list of RNA objects (according to the value of the parameter type) (see pyrna.model)
    """
    molecules = []
    pieces = []
    molecule_name = None
    for line in fasta_data.split('\n'):
        if re.match('>',line):
            if molecule_name and len(pieces) > 0:
                m = RNA(sequence = ''.join(pieces), name = molecule_name.strip())
                if m != None:
                    molecules.append(m)
            molecule_name = line[1:]
            pieces = []
        else:
            pieces.append(line.strip().upper())
    #last molecule
    if molecule_name and len(pieces) > 0:
        m = RNA(sequence = ''.join(pieces), name = molecule_name.strip())
        if m != None:
            molecules.append(m)
    return molecules

def parse_vienna(vienna_data):
    """
    Parse Vienna data

    Parameters:
    ---------
     - vienna_data: the Vienna data as a String

    Returns:
    ------
    tuple containg a list of RNA objects and a list of secondary structures, each 2D described as a list of BasePair objects
    """
    name = None
    secondary_structures = []
    rnas = []
    current_bn = []
    current_sequence = []
    for line in vienna_data.split('\n'):
        if re.match('^[\.()\{\}\[\]]+$', line):
            current_bn.append(line)
        elif re.match('^>', line):
            if len(current_sequence):
                rnas.append(RNA(name = name, sequence = ''.join(current_sequence)))
                if len(current_bn):
                    secondary_structures.append(parse_bn(''.join(current_bn)))
            name = line[1:]
            current_bn = []
            current_sequence = []
        elif len(line.strip()):
            current_sequence.append(line.strip())

    #last one
    if len(current_sequence):
        rnas.append(RNA(name = name, sequence = ''.join(current_sequence)))
        if len(current_bn):
            secondary_structures.append(parse_bn(''.join(current_bn)))
    return rnas, secondary_structures

def parse_bn(bn):
    """
    Parse a bracket notation. The function supports characters like '(', ')', '[', ']', '{' and '}'

    Parameters:
    ---------
     - bn: the bracket notation as a String

    Returns:
    ------
    a list of base pairs
    """

    i = 0
    lastPairedPos = []
    lastPairedSymbol = []
    basePairs = []

    for s in list(bn):
        i+=1
        if s in ['(','{','[']:
            lastPairedPos.append(i)
            lastPairedSymbol.append(s)
        elif s in [')','}',']']:
            basePairs.append(BasePair(lastPairedPos.pop(),i))

    return basePairs

def parse_pdb(pdb_data):
    """
    Parse PDB data.

    Parameters:
    ---------
     - pdb_data: the PDB data as a list of lines

    Returns:
    ------
    a list of TertiaryStructure objects (see pyrna.model). if the PDB data describes a tertiary structure made with several molecular chains, this method will return one TertiaryStructure object per chain.
    """
    molecules = []
    chains = []
    tertiary_structures = []
    current_chain = None
    current_residue_name = None
    current_residue_pos = None
    absolute_position = -1
    current_molecule = None
    residues = []
    current_3D = None
    title = "N.A."

    for line in pdb_data:
        header = line[0:6].strip()
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip().upper()
        chain_name = line[21:22].strip()
        residue_pos = line[22:27].strip()

        if (header == "ATOM" or header == "HETATM") and not residue_name in ["FMN","PRF","HOH","MG","OHX","MN","ZN", "SO4", "CA", "UNK", "AMO"] and not atom_name in ["MG","K", "NA", "SR", "CL", "CD", "ACA"] and len(chain_name):
            if chain_name != current_chain: #new chain
                current_residue_name = residue_name
                current_residue_pos = residue_pos
                current_chain = chain_name
                absolute_position = 1
                residues = []
                current_molecule = RNA(sequence="", name = current_chain)
                residues.append(current_residue_name)
                current_3D = TertiaryStructure(current_molecule)
                current_3D.title = re.sub(' +', ' ', title)
                tertiary_structures.append(current_3D)

            elif current_residue_pos != residue_pos: # new residue
                current_residue_name = residue_name
                current_residue_pos = residue_pos
                absolute_position += 1

            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            current_3D.add_atom(atom_name, current_residue_name, current_chain, absolute_position, [x,y,z])

        elif header == 'TITLE':
            title += line[10:]

        elif header == "TER":
            current_chain = None
            current_residue_pos = None
            current_molecule = None
            residues = []

    return tertiary_structures
