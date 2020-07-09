# Differences over V1:
# 1. Works with CIF format instead of PDB
# 2. Also deletes all Biological assemblies not corresponding to desired number of monomers

from os.path import join, exists
from os import listdir, remove, rename, mkdir
from Bio import SeqIO
import sys
from Bio.PDB.MMCIFParser import MMCIFParser
import shutil

def filename_to_id(filename):
    return filename.split("/")[-1].split("-")[0]

def num_chains(filename):
    global assembly_directory
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure(filename.split("-")[0], join(assembly_directory, filename))
    return len(list(struct.get_chains()))

if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) != 3:
        print(f"\nUsage: $ python PickRepresentativeAssembly.py <Root of Species Tree> <Number of Monomers>\n")
        sys.exit()
    root_directory = sys.argv[1]
    n_monomers = int(sys.argv[2])

    print("1. SELECTING ASSEMBLIES")
    for species_fol in listdir(root_directory):
        print(species_fol)
        assembly_directory = join(root_directory, species_fol, "All_Biological_Assemblies")
        id_to_rep_assembly = {}

        # Go through all the PDB IDs, and see which of the biological assembly .cif files to keep
        for filename in listdir(assembly_directory):
            pdb_id = filename_to_id(filename)
            assembly_nr = filename.split("/")[-1].split("-")[-1].split(".")[0]
            try:
                if not pdb_id in id_to_rep_assembly and num_chains(filename) == n_monomers:
                    id_to_rep_assembly[pdb_id] = assembly_nr
                    print(f"    Keeping {pdb_id} as {assembly_nr}")
                else:
                    print(f"    Removing {filename.split('.')[0]}")
                    remove(join(assembly_directory, f"{pdb_id}-{assembly_nr}.cif"))
            except:
                print(f"    ERROR with {filename.split('.')[0]}, removing")
                remove(join(assembly_directory, f"{pdb_id}-{assembly_nr}.cif"))


        # Rename all the biological assembly files to remove the assembly number
        selected_assemblies_folder = "Selected_Biological_Assemblies"
        selected_assemblies_dir = join(root_directory, species_fol, selected_assemblies_folder)
        if not exists(selected_assemblies_dir): mkdir(selected_assemblies_dir)
        print("    Saving selected assemblies ...")
        for pdb_id, assembly_nr in id_to_rep_assembly.items():
            shutil.copy(join(assembly_directory, f"{pdb_id}-{assembly_nr}.cif"), join(selected_assemblies_dir, f"{pdb_id}.cif"))
    

            
