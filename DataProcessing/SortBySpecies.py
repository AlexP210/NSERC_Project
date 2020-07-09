import os
import pypdb
import sys

if __name__ == "__main__":
    # Check Arguments
    if len(sys.argv) < 2:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <Root Folder for Species Tree> <CIF Folders>\n")
        sys.exit()
    root_folder = sys.argv[1]
    cif_folders = sys.argv[2:]

    for folder in cif_folders:
        for filename in os.listdir(folder):
            cif_filepath = os.path.join(folder, filename)
            pdb_id = filename.split(".")[0]
            info = pypdb.get_all_info(pdb_id)
            species_for_chains = set()
            for entity in info["polymer"]:
                if type(entity["Taxonomy"]) == type({}):
                    species_for_chains.add(entity["Taxonomy"]["@name"])
            if len(species_for_chains) == 1:
                if os.path.exists(os.path.join(root_folder, ))
