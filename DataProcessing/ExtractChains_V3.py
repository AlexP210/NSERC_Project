# Changes from V2:
#   1. Works with CIF files, and does not require a chains.txt

import os
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO
from Bio import PDB
from Bio import SeqIO

from My_Library import *

class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = MMCIFParser(QUIET=True)
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
        pdb_id = pdb_fn.split(".")[0]
        out_name = "%s_%s.ent" % (pdb_id, "".join(chain_letters))
        out_path = os.path.join(self.out_dir, out_name)
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("        Chain%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(chain_letters), pdb_id, out_name))
            return out_path

        # print("Extracting chain%s %s from %s..." % (plural,
        #         ", ".join(chain_letters), pdb_fn))

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
    if len(sys.argv) != 2:
        print(f"\nUsage: python ExtractChains_V3.py <Root of Species Tree>\n")
        sys.exit()

    root_dir = sys.argv[1]

    print("1. EXTRACTING CHAINS")
    for species_fol in os.listdir(root_dir):
        print(species_fol)
        # Set the folder with the selected assemblies, and the folder to save the extracted chains
        assembly_dir = os.path.join(root_dir, species_fol, "Selected_Biological_Assemblies")
        chains_save_dir = os.path.join(root_dir, species_fol, "Extracted_Chains")
        if not os.path.exists(chains_save_dir): os.mkdir(chains_save_dir)

        pdbList = PDB.PDBList()
        splitter = ChainSplitter(chains_save_dir)

        for assembly_filename in os.listdir(assembly_dir):
            path = os.path.join(assembly_dir, assembly_filename)
            pdb_id = assembly_filename.split(".")[0]
            for chain in get_chain_ids(path):
                print(f"    Extracting {pdb_id}:{chain}")
                # If we have it and haven't split it, do so
                if not os.path.isfile( os.path.join(chains_save_dir, f"{pdb_id}_{chain}.ent") ): 
                    splitter.make_pdb(os.path.join(assembly_dir, assembly_filename), chain)
            