import re
from itertools import groupby, product
from operator import itemgetter

class Block:
    """
    A continuous range of molecular positions, with a single start and end point.
    """
    def __init__(self, start, end):
        if start < end:
            self.start = start
            self.end = end
        else:
            self.start = end
            self.end = start

class Location:
    """
    A Location defines a range of molecular positions, continuous or not. A location is made with Block objects.
    """
    def __init__(self, start = None, end = None, single_positions = None, nested_lists = None):
        """
        To instantiate a Location, you can:
        - set a start and end position: Location(start=34, end=69). The location will contain all the positions between the start and end.
        - list all the single positions (sorted or not) to be contained in the location: Location(single_positions=[34, 56, 57, 58, 67, 68, 69])
        - list the ranges of continuous positions as nested lists: Location(nested_lists=[[34,34], [56,58], [67,69]])
        """
        self.blocks = []
        if start and end:
            self.add_block(Block(start, end))
        elif single_positions:
            single_positions.sort()
            for k, g in groupby(enumerate(single_positions), lambda i,x:i-x):
                _range = map(itemgetter(1), g)
                self.blocks.append(Block(min(_range), max(_range)))
        elif nested_lists:
            for nested_list in nested_lists:
                self.blocks.append(Block(min(nested_list), max(nested_list)))

    def add_block(self, block):
        blocks_to_remove = []

        for _block in self.blocks:
            if block.is_before(_block) and not block.is_beside(_block):
                break
            elif block.intersects(_block) or block.is_beside(_block):
                block.merge(_block)
                blocks_to_remove.append(_block)
                #its necessary to continue to see if the new Block can merge with other blocks
                continue
            elif len(blocks_to_remove):
                break

        for block_to_remove in blocks_to_remove:
            self.blocks.remove(block_to_remove)

        self.blocks.append(block)
        self.blocks = sorted(self.blocks, key=lambda block: block.start)

    def remove_location(self, location):
        """
        Return a new Location object from the difference between the current Location and the Location given as argument.
        Difference means all the positions not found in the Location given as argument
        """
        single_positions_1 = self.get_single_positions()
        single_positions_2 = location.get_single_positions()

        diff = list(set(single_positions_1) - set(single_positions_2))

        return Location(single_positions = diff)

    def remove_locations(self, locations):
        """
        Return a new Location object from the difference between the current Location with all the Locations given in a list as argument.
        Difference means all the positions not found in the Locations given as argument
        """
        single_positions_1 = self.get_single_positions()
        single_positions_2 = []

        for location in locations:
            single_positions_2 += location.get_single_positions()

        diff = list(set(single_positions_1) - set(single_positions_2))

        return Location(single_positions = diff)


    def get_single_positions(self):
        """
        Returns:
        ------
        all the single positions making this Location as a list.
        """
        single_positions = []
        for block in self.blocks:
            single_positions += range(block.start, block.end+1)
        return single_positions

    def has_position(self, position):
        """
        Test if the location encloses a single position.
        Parameters:
        ---------
        position: an integer
        """
        return position in self.get_single_positions()

    def start(self):
        return self.blocks[0].start

    def end(self):
        return self.blocks[-1].end

class Molecule:
    def __init__(self, name):
        self.modified_residues = []
        self.name = name
        self.organism = None
        self.sequence = ""

    def to_fasta(self, single_line=False):
        lines = []
        lines.append(">" + self.name)
        if single_line:
            lines.append(self.sequence)
        else:
            c = 0
            while c < len(self.sequence):
                d = min(len(self.sequence), c + 79)
                lines.append(self.sequence[c:d])
                c += 79
        return '\n'.join(lines)

    def __add__(self, seq):
        if seq.__class__ == str:
            self.sequence = ''.join([self.sequence, seq])

    def __sub__(self, length):
        if length.__class__ == int and length <= len(self.sequence):
            self.sequence = self.sequence[0: len(self.sequence)-length]

    def __len__(self):
        return len(self.sequence)

    def __iter__(self):
        return iter(self.sequence)

    def __getslice__(self, i, j):
        return self.sequence.__getslice__(i, j)

    def __getitem__(self, i):
        return self.sequence.__getitem__(i)


class RNA(Molecule):
    def __init__(self, sequence, name = 'rna'):
        Molecule.__init__(self, name)

        for residue in list(sequence):
            self.add_residue(residue)

    def add_residue(self, residue):
        if residue in modified_ribonucleotides:
            self.modified_residues.append((residue, len(self.sequence)+1))
            residue = modified_ribonucleotides[residue]
        if residue in ['A', 'U', 'G', 'C']:
            self.sequence = ''.join([self.sequence, residue])
        elif residue in ['.', '_', '-']:
            self.sequence = ''.join([self.sequence, '-'])
        else:
            #print "Unknown residue "+residue
            self.sequence = ''.join([self.sequence, residue])

    def get_complement(self):
        """
        Returns:
        ------
        the complement sequence as a string.
        """
        basecomplement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
        letters = list(self.sequence)
        letters = [basecomplement[base] if base in basecomplement else base for base in letters]
        return ''.join(letters)
    
    
class BasePair:
    def __init__(self, pos1, pos2, edge1 = "WC", edge2 = "WC", orientation = "cis"):
        """
        Parameters:
        ---------
        edge1: string. Can be WC, H or SE
        edge2: string. Can be WC, H or SE
        orientation: string. Can be cis or trans
        """
        self.location = Location(nested_lists=[[pos1,pos1],[pos2,pos2]])
        self.edge1 = edge1
        self.edge2 = edge2
        self.orientation = orientation

class Helix:
    def __init__(self, start, end, length):
        """
        Parameters:
        ---------
        start: the position of the first residue in the first strand 
        end: the position of the last residue in the second strand 
        length: helix size, meaning the number of basepairs
        """
        self.location = Location(nested_lists=[[start,start+length-1],[end-length+1,end]])

class SecondaryStructure:
    """
    This class stores everything describing a secondary structure: helices, junctions, tertiary interactions...
    """
    def __init__(self, rna):
        self.name = "2D"
        self.rna = rna
        self.tertiary_interactions = []
        self.helices = []
        self.junctions = []
        self.base_pairs = []

class Atom:
    def __init__(self, name, x, y, z):
        self.name = name
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return self.name

class EndocyclicA(Atom):
    def __init__(self, name, x, y, z):
        Atom.__init__(self,name, x, y, z)

    def __str__(self):
        return self.name

class DonorEndoA(EndocyclicA):
    def __init__(self, name, x, y, z):
        EndocyclicA.__init__(self,name, x, y, z)

    def __str__(self):
        return self.name

class AcceptorEndoA(EndocyclicA):
    def __init__(self, name, x, y, z):
        EndocyclicA.__init__(self,name, x, y, z)

    def __str__(self):
        return self.name

class ExocyclicA(Atom):
    def __init__(self, name, x, y, z):
        Atom.__init__(self,name, x, y, z)

    def __str__(self):
        return self.name

class DonorExoA(ExocyclicA):
    def __init__(self, name, x, y, z):
        ExocyclicA.__init__(self,name, x, y, z)

    def __str__(self):
        return self.name

class AcceptorExoA:
    def __init__(self, name, x, y, z):
        ExocyclicA.__init__(self,name, x, y, z)

    def __str__(self):
        return self.name

class Residue3D:
    def __init__(self, chain_name, position):
        self.atoms = [] #a list of Atom objects
        self.position = position
        self.chain_name = chain_name

    def add_atom(self, atom_name, coords):
        self.atoms.append(Atom(atom_name, coords[0], coords[1], coords[2]))
    
    def get_WC_atoms(self):
        return []

class Guanine3D(Residue3D):
    def __init__(self, chain_name, position):
        Residue3D.__init__(self, chain_name, position)

    def add_atom(self, atom_name, coords):
        match atom_name:
            case "N1": self.atoms.append(DonorEndoA(atom_name, coords[0], coords[1], coords[2]))
            case "N2": self.atoms.append(DonorExoA(atom_name, coords[0], coords[1], coords[2]))
            case "N3": self.atoms.append(AcceptorEndoA(atom_name, coords[0], coords[1], coords[2]))
            case "O6": self.atoms.append(AcceptorExoA(atom_name, coords[0], coords[1], coords[2]))
            case "N7": self.atoms.append(AcceptorEndoA(atom_name, coords[0], coords[1], coords[2]))
            case _: self.atoms.append(Atom(atom_name, coords[0], coords[1], coords[2]))

    def get_WC_atoms(self):
        target_atom_names = ["N1", "N2", "O6"] 
        return [a for a in self.atoms if a.name in target_atom_names]
    
    def __str__(self):
        return "G"+str(self.position)

class Adenine3D(Residue3D):
    def __init__(self, chain_name, position):
        Residue3D.__init__(self, chain_name, position)

    def add_atom(self, atom_name, coords):
        match atom_name:
            case "N1": self.atoms.append(AcceptorEndoA(atom_name, coords[0], coords[1], coords[2]))
            case "N3": self.atoms.append(AcceptorEndoA(atom_name, coords[0], coords[1], coords[2]))
            case "N6": self.atoms.append(DonorExoA(atom_name, coords[0], coords[1], coords[2]))
            case "N7": self.atoms.append(AcceptorEndoA(atom_name, coords[0], coords[1], coords[2]))
            case _: self.atoms.append(Atom(atom_name, coords[0], coords[1], coords[2]))

    def get_WC_atoms(self):
        target_atom_names = ["N1", "N6"] 
        return [a for a in self.atoms if a.name in target_atom_names]

    def __str__(self):
        return "A"+str(self.position)

class Uracil3D(Residue3D):
    def __init__(self, chain_name, position):
        Residue3D.__init__(self, chain_name, position)

    def add_atom(self, atom_name, coords):
        match atom_name:
            case "O2": self.atoms.append(AcceptorExoA(atom_name, coords[0], coords[1], coords[2]))
            case "N3": self.atoms.append(DonorEndoA(atom_name, coords[0], coords[1], coords[2]))
            case "O4": self.atoms.append(AcceptorExoA(atom_name, coords[0], coords[1], coords[2]))
            case _: self.atoms.append(Atom(atom_name, coords[0], coords[1], coords[2]))

    def get_WC_atoms(self):
        target_atom_names = ["O2", "N3", "O4"] 
        return [a for a in self.atoms if a.name in target_atom_names]

    def __str__(self):
        return "U"+str(self.position)

class Cytosine3D(Residue3D):
    def __init__(self, chain_name, position):
        Residue3D.__init__(self, chain_name, position)

    def add_atom(self, atom_name, coords):
        match atom_name:
            case "O2": self.atoms.append(AcceptorExoA(atom_name, coords[0], coords[1], coords[2]))
            case "N3": self.atoms.append(AcceptorEndoA(atom_name, coords[0], coords[1], coords[2]))
            case "N4": self.atoms.append(DonorExoA(atom_name, coords[0], coords[1], coords[2]))
            case _: self.atoms.append(Atom(atom_name, coords[0], coords[1], coords[2]))

    def get_WC_atoms(self):
        target_atom_names = ["O2", "N3", "N4"]
        return [a for a in self.atoms if a.name in target_atom_names]

    def __str__(self):
        return "C"+str(self.position)

class TertiaryStructure:

    def __init__(self, rna):
        self.rna = rna
        self.name = "N.A."
        self.residues = [] #a list of Residue3D objects

    def add_atom(self, atom_name, residue_name, chain_name, absolute_position, coords):
        atom_name = re.sub("\*", "'", atom_name)
        if atom_name == 'OP1':
            atom_name = 'O1P'
        elif atom_name == 'OP2':
            atom_name = 'O2P'
        elif atom_name == 'OP3':
            atom_name = 'O3P'
        if absolute_position-1 >= len(self.residues):
            self.rna.add_residue(residue_name)
            match(residue_name):
                case "A": self.residues.append(Adenine3D(chain_name, absolute_position))
                case "G": self.residues.append(Guanine3D(chain_name, absolute_position))
                case "U": self.residues.append(Uracil3D(chain_name, absolute_position))
                case "C": self.residues.append(Cytosine3D(chain_name, absolute_position))
                case _: raise RuntimeError("Unknown residue "+residue_name) 
        self.residues[absolute_position-1].add_atom(atom_name,coords)

    def find_canonical_basepairs(self):
        for i in range(0, len(self.residues)-1):
            r = self.residues[i]
            for j in range(i+1, len(self.residues)):
                next_r = self.residues[j]
                match r :
                    case Adenine3D():
                        match next_r:
                            case Uracil3D(): compute_pairing(r, next_r)
                    case Guanine3D():
                        match next_r:
                            case Cytosine3D(): compute_pairing(r,next_r)
                            case Uracil3D(): compute_pairing(r,next_r)
                    case Uracil3D():
                        match next_r:
                            case Guanine3D(): compute_pairing(r,next_r)
                            case Adenine3D(): compute_pairing(r,next_r)
                    case Cytosine3D():
                        match next_r:
                            case Guanine3D(): compute_pairing(r,next_r)

"""
r1: a first instance of Residue3D
r2: a second instance of Residue3D 
"""
def compute_pairing(r1, r2):
    print(f"{r1}-{r2}")
    donors = [a for a in r1.atoms if isinstance(a, DonorEndoA) or isinstance(a, DonorExoA)]
    acceptors = [a for a in r1.atoms if isinstance(a, AcceptorEndoA) or isinstance(a, AcceptorExoA)]
    
    acceptors_2 = [a for a in r2.atoms if isinstance(a, AcceptorEndoA) or isinstance(a, AcceptorExoA)]
    hbonds = list(product(donors, acceptors_2))
    for hbond in hbonds:
        print(f"{hbond[0]} ~ {hbond[1]}")
    donors_2 = [a for a in r2.atoms if isinstance(a, DonorEndoA) or isinstance(a, DonorExoA)]
    hbonds = list(product(acceptors, donors_2))
    for hbond in hbonds:
        print(f"{hbond[0]} ~ {hbond[1]}")
            
modified_ribonucleotides = {
    "T": "U",
    "PSU": "U",
    "I": "A",
    "N": "U",
    "S": "U",
    "+A": "A",
    "+C": "C",
    "+G": "G",
    "+I": "I",
    "+T": "U",
    "+U": "U",
    "PU": "A",
    "YG": "G",
    "1AP": "G",
    "1MA": "A",
    "1MG": "G",
    "2DA": "A",
    "2DT": "U",
    "2MA": "A",
    "2MG": "G",
    "4SC": "C",
    "4SU": "U",
    "5IU": "U",
    "5MC": "C",
    "5MU": "U",
    "5NC": "C",
    "6MP": "A",
    "7MG": "G",
    "A23": "A",
    "AD2": "A",
    "AET": "A",
    "AMD": "A",
    "AMP": "A",
    "APN": "A",
    "ATP": "A",
    "AZT": "U",
    "CCC": "C",
    "CMP": "A",
    "CPN": "C",
    "DAD": "A",
    "DCT": "C",
    "DDG": "G",
    "DG3": "G",
    "DHU": "U",
    "DOC": "C",
    "EDA": "A",
    "G7M": "G",
    "GDP": "G",
    "GNP": "G",
    "GPN": "G",
    "GTP": "G",
    "GUN": "G",
    "H2U": "U",
    "HPA": "A",
    "IPN": "U",
    "M2G": "G",
    "MGT": "G",
    "MIA": "A",
    "OMC": "C",
    "OMG": "G",
    "OMU": "U",
    "ONE": "U",
    "P2U": "P",
    "PGP": "G",
    "PPU": "A",
    "PRN": "A",
    "PST": "U",
    "QSI": "A",
    "QUO": "G",
    "RIA": "A",
    "SAH": "A",
    "SAM": "A",
    "T23": "U",
    "T6A": "A",
    "TAF": "U",
    "TLC": "U",
    "TPN": "U",
    "TSP": "U",
    "TTP": "U",
    "UCP": "U",
    "VAA": "A",
    "YYG": "G",
    "70U": "U",
    "12A": "A",
    "2MU": "U",
    "127": "U",
    "125": "U",
    "126": "U",
    "MEP": "U",
    "TLN": "U",
    "ADP": "A",
    "TTE": "U",
    "PYO": "U",
    "SUR": "U",
    "PSD": "A",
    "S4U": "U",
    "CP1": "C",
    "TP1": "U",
    "NEA": "A",
    "GCK": "C",
    "CH": "C",
    "EDC": "G",
    "DFC": "C",
    "DFG": "G",
    "DRT": "U",
    "2AR": "A",
    "8OG": "G",
    "IG": "G",
    "IC": "C",
    "IGU": "G",
    "IMC": "C",
    "GAO": "G",
    "UAR": "U",
    "CAR": "C",
    "PPZ": "A",
    "M1G": "G",
    "ABR": "A",
    "ABS": "A",
    "S6G": "G",
    "HEU": "U",
    "P": "G",
    "DNR": "C",
    "MCY": "C",
    "TCP": "U",
    "LGP": "G",
    "GSR": "G",
    "X": "G",
    "R": "A",
    "Y": "A",
    "E": "A",
    "GSS": "G",
    "THX": "U",
    "6CT": "U",
    "TEP": "G",
    "GN7": "G",
    "FAG": "G",
    "PDU": "U",
    "MA6": "A",
    "UMP": "U",
    "SC": "C",
    "GS": "G",
    "TS": "U",
    "AS": "A",
    "ATD": "U",
    "T3P": "U",
    "5AT": "U",
    "MMT": "U",
    "SRA": "A",
    "6HG": "G",
    "6HC": "C",
    "6HT": "U",
    "6HA": "A",
    "55C": "C",
    "U8U": "U",
    "BRO": "U",
    "BRU": "U",
    "5IT": "U",
    "ADI": "A",
    "5CM": "C",
    "IMP": "G",
    "THM": "U",
    "URI": "U",
    "AMO": "A",
    "FHU": "P",
    "TSB": "A",
    "CMR": "C",
    "RMP": "A",
    "SMP": "A",
    "5HT": "U",
    "RT": "U",
    "MAD": "A",
    "OXG": "G",
    "UDP": "U",
    "6MA": "A",
    "5IC": "C",
    "SPT": "U",
    "TGP": "G",
    "BLS": "A",
    "64T": "U",
    "CB2": "C",
    "DCP": "C",
    "ANG": "G",
    "BRG": "G",
    "Z": "A",
    "AVC": "A",
    "5CG": "G",
    "UDP": "U",
    "UMS": "U",
    "BGM": "G",
    "SMT": "U",
    "DU": "U",
    "CH1": "C",
    "GH3": "G",
    "GNG": "G",
    "TFT": "U",
    "U3H": "U",
    "MRG": "G",
    "ATM": "U",
    "GOM": "A",
    "UBB": "U",
    "A66": "A",
    "T66": "U",
    "C66": "C",
    "3ME": "A",
    "A3P": "A",
    "ANP": "A",
    "FA2": "A",
    "9DG": "G",
    "GMU": "U",
    "UTP": "U",
    "5BU": "U",
    "APC": "A",
    "DI": "I",
    "UR3": "U",
    "3DA": "A",
    "DDY": "C",
    "TTD": "U",
    "TFO": "U",
    "TNV": "U",
    "MTU": "U",
    "6OG": "G",
    "E1X": "A",
    "FOX": "A",
    "CTP": "C",
    "D3T": "U",
    "TPC": "C",
    "7DA": "A",
    "7GU": "U",
    "2PR": "A",
    "CBR": "C",
    "I5C": "C",
    "5FC": "C",
    "GMS": "G",
    "2BT": "U",
    "8FG": "G",
    "MNU": "U",
    "AGS": "A",
    "NMT": "U",
    "NMS": "U",
    "UPG": "U",
    "G2P": "G",
    "2NT": "U",
    "EIT": "U",
    "TFE": "U",
    "P2T": "U",
    "2AT": "U",
    "2GT": "U",
    "2OT": "U",
    "BOE": "U",
    "SFG": "G",
    "CSL": "I",
    "PPW": "G",
    "IU": "U",
    "D5M": "A",
    "ZDU": "U",
    "DGT": "U",
    "UD5": "U",
    "S4C": "C",
    "DTP": "A",
    "5AA": "A",
    "2OP": "A",
    "PO2": "A",
    "DC": "C",
    "DA": "A",
    "LOF": "A",
    "ACA": "A",
    "BTN": "A",
    "PAE": "A",
    "SPS": "A",
    "TSE": "A",
    "A2M": "A",
    "NCO": "A",
    "A5M": "C",
    "M5M": "C",
    "S2M": "U",
    "MSP": "A",
    "P1P": "A",
    "N6G": "G",
    "MA7": "A",
    "FE2": "G",
    "AKG": "G",
    "SIN": "G",
    "PR5": "G",
    "GOL": "G",
    "XCY": "G",
    "5HU": "U",
    "CME": "C",
    "EGL": "G",
    "LC": "C",
    "LHU": "U",
    "LG": "G",
    "PUY": "U",
    "PO4": "U",
    "PQ1": "U",
    "ROB": "U",
    "O2C": "C",
    "C30": "C",
    "C31": "C",
    "C32": "C",
    "C33": "C",
    "C34": "C",
    "C35": "C",
    "C36": "C",
    "C37": "C",
    "C38": "C",
    "C39": "C",
    "C40": "C",
    "C41": "C",
    "C42": "C",
    "C43": "C",
    "C44": "C",
    "C45": "C",
    "C46": "C",
    "C47": "C",
    "C48": "C",
    "C49": "C",
    "C50": "C",
    "A30": "A",
    "A31": "A",
    "A32": "A",
    "A33": "A",
    "A34": "A",
    "A35": "A",
    "A36": "A",
    "A37": "A",
    "A38": "A",
    "A39": "A",
    "A40": "A",
    "A41": "A",
    "A42": "A",
    "A43": "A",
    "A44": "A",
    "A45": "A",
    "A46": "A",
    "A47": "A",
    "A48": "A",
    "A49": "A",
    "A50": "A",
    "G30": "G",
    "G31": "G",
    "G32": "G",
    "G33": "G",
    "G34": "G",
    "G35": "G",
    "G36": "G",
    "G37": "G",
    "G38": "G",
    "G39": "G",
    "G40": "G",
    "G41": "G",
    "G42": "G",
    "G43": "G",
    "G44": "G",
    "G45": "G",
    "G46": "G",
    "G47": "G",
    "G48": "G",
    "G49": "G",
    "G50": "G",
    "T30": "U",
    "T31": "U",
    "T32": "U",
    "T33": "U",
    "T34": "U",
    "T35": "U",
    "T36": "U",
    "T37": "U",
    "T38": "U",
    "T39": "U",
    "T40": "U",
    "T41": "U",
    "T42": "U",
    "T43": "U",
    "T44": "U",
    "T45": "U",
    "T46": "U",
    "T47": "U",
    "T48": "U",
    "T49": "U",
    "T50": "U",
    "U30": "U",
    "U31": "U",
    "U32": "U",
    "U33": "U",
    "U34": "U",
    "U35": "U",
    "U36": "U",
    "U37": "U",
    "U38": "U",
    "U39": "U",
    "U40": "U",
    "U41": "U",
    "U42": "U",
    "U43": "U",
    "U44": "U",
    "U45": "U",
    "U46": "U",
    "U47": "U",
    "U48": "U",
    "U49": "U",
    "U50": "U",
    "UFP": "U",
    "UFR": "U",
    "UCL": "U",
    "3DR": "U",
    "CBV": "C",
    "HFA": "A",
    "MMA": "A",
    "DCZ": "C",
    "GNE": "C",
    "A1P": "A",
    "6IA": "A",
    "CTG": "G",
    "5FU": "U",
    "2AD": "A",
    "T2T": "U",
    "XUG": "G",
    "2ST": "U",
    "5PY": "U",
    "4PC": "C",
    "US1": "U",
    "M5C": "C",
    "DG": "G",
    "DA": "A",
    "DT": "U",
    "DC": "C",
    "P5P": "A",
    "FMU": "U"
}
