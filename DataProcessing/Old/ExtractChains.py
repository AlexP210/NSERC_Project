import os
from Bio import PDB


class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, chain_letters, overwrite=False, struct=None):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """
        chain_letters = [chain.upper() for chain in chain_letters]

        # Input/output files
        pdb_fn = os.path.split(pdb_path)[1]
        pdb_id = pdb_fn[3:7]
        out_name = "pdb%s_%s.ent" % (pdb_id, "".join(chain_letters))
        out_path = os.path.join(self.out_dir, out_name)
        print("OUT PATH:",out_path)
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(chain_letters), pdb_id, out_name))
            return out_path

        print("Extracting chain%s %s from %s..." % (plural,
                ", ".join(chain_letters), pdb_fn))

        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters))

        return out_path


class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


if __name__ == "__main__":
    """ Parses PDB id's desired chains, and creates new PDB structures. """
    import sys
    if not len(sys.argv) == 4:
        print(f"Usage: python {__file__} <chains.txt> <chains location> <pdb location>" % __file__)
        sys.exit()

    chains_to_get_path = sys.argv[1]
    chains_save_path = sys.argv[2]
    pdb_download_path = sys.argv[3]

    pdbList = PDB.PDBList()
    splitter = ChainSplitter(chains_save_path)

    with open(chains_to_get_path) as pdb_textfile:
        for line in pdb_textfile:
            pdb_id = line[:4].lower()
            chain = line[4].upper()
            # if we don't have the pdb, download it
            if not os.path.isfile(os.path.join(pdb_download_path, f"pdb{pdb_id}.ent")):
                pdb_fn = pdbList.retrieve_pdb_file(pdb_id, file_format = "pdb", pdir=pdb_download_path)
            # I we have it and haven't split it, do so
            if os.path.isfile(pdb_fn) and not os.path.isfile( os.path.join(chains_save_path, f"pdb{pdb_id}_{chain}") ): 
                splitter.make_pdb(pdb_fn, chain)