from urllib.request import urlopen

class PDB:
    """
    Wrapper for the PDB database http://www.rcsb.org/
    """

    def __init__(self):
        pass

    def get_entry(self, pdb_id):
        """
        Return the content of a PDB entry as a list of lines
        """
        response = urlopen("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s"%pdb_id)
        content = str(response.read())
        return content.split("\\n")