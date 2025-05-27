from urllib.request import urlopen
from pyrna.model import SecondaryStructure, SecondaryStructureProvider, SecondaryStructureProviderParams
from pyrna.parsers import parse_pdb


class PDBProviderParams(SecondaryStructureProviderParams):

    """
    This class stores the parameters needed by a SecondaryStructureProvider to connect to its data source and construct a SecondaryStructure object
    """
    
    def __init__(self):
        SecondaryStructureProviderParams.__init__(self) 

    def set_pdb_id(self, pdb_id):
        self.params["pdb_id"] =  pdb_id

    def get_pdb_id(self):
        return self.params["pdb_id"]
    
    def set_chain_name(self, chain_name):
        self.params["chain_name"] = chain_name

    def get_chain_name(self):
        return self.params["chain_name"]

class PDBProvider(SecondaryStructureProvider):   
    """
    A SecondaryStructureProvider using the data stored in the Protein Database http://www.rcsb.org/
    """

    def __init__(self):
        SecondaryStructureProvider.__init__(self) 

    def get_entry(self, pdb_id):
        """
        Return the content of a PDB entry as a list of lines
        """
        response = urlopen("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s"%pdb_id)
        content = str(response.read())
        return content.split("\\n")
    
    def get_secondary_structure(self, params):
        for ts in parse_pdb(self.get_entry(params.get_pdb_id())):
            if ts.rna.name == params.get_chain_name():
                ss = SecondaryStructure(rna = ts.rna, base_pairs=ts.find_canonical_basepairs(chain_name = ts.rna.name))
                return ss

