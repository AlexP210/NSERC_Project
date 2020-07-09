import os
from Bio import PDB
import pypdb

if __name__ == "__main__":
    """ Downloads CIF files listed for PDB IDs listed in chains.txt """
    import sys
    if not len(sys.argv) == 3:
        print(f"\nUsage: $ python DownloadCIFs.py <IdList.txt> <Root Folder of CIF Tree>\n")
        sys.exit()

    ids_to_get_path = sys.argv[1]
    root_folder = sys.argv[2]
    if not os.path.exists(root_folder): os.mkdir(root_folder)

    pdbList = PDB.PDBList()
    print("1. DOWNLOADING .CIF FILES")

    with open(ids_to_get_path) as pdb_textfile:
        for line in pdb_textfile:
            pdb_id = line[:4]

            pdb_id = line[:4]
            info = pypdb.get_all_info(pdb_id)
            species_for_chains = []
            if type(info["polymer"]) == type([]):
                # Why in some cases is there no "Taxonomy"?
                for entity in info["polymer"]:
                    if "Taxonomy" in entity and type(entity["Taxonomy"]) == type({}): # Expand this to include cases where there is more than 1 species
                        if entity["Taxonomy"]["@name"] not in species_for_chains: species_for_chains.append(entity["Taxonomy"]["@name"])
            else:
                entity = info["polymer"]
                if "Taxonomy" in entity and type(entity["Taxonomy"]) == type({}): # Expand this to include cases where there is more than 1 species
                    if entity["Taxonomy"]["@name"] not in species_for_chains: species_for_chains.append(entity["Taxonomy"]["@name"])
            if len(species_for_chains) == 1:
                full_species_name = species_for_chains[0]
                species_folder_name = "_".join(full_species_name.split(" "))
                species_folder_path = os.path.join(root_folder, species_folder_name)
                if not os.path.exists(species_folder_path):
                    os.mkdir(species_folder_path)
                download_folder_path = os.path.join(species_folder_path, "Structure_Downloads")
                if not os.path.exists(download_folder_path):
                    os.mkdir(download_folder_path)
                # if we don't have the pdb, download it
                if not os.path.isfile(os.path.join(download_folder_path, f"{pdb_id}.cif")):
                    print(f"Downloading {pdb_id}.cif")
                    pdb_fn = pdbList.retrieve_pdb_file(pdb_id, pdir=download_folder_path)
                else:
                    print(f"Already have {pdb_id}.cif")