import os
from Bio import PDB
import pypdb
import My_Library as ml

if __name__ == "__main__":
    """ Downloads CIF files listed for PDB IDs listed in chains.txt """
    import sys
    if not len(sys.argv) == 4:
        print(f"\nUsage: $ python DownloadCIFs.py <IdList.txt> <Root Folder of CIF Tree> <Number of Monomers>\n")
        sys.exit()

    ids_to_get_path = sys.argv[1]
    root_folder = sys.argv[2]
    if not os.path.exists(root_folder): os.mkdir(root_folder)
    n_monomers = int(sys.argv[3])

    pdbList = PDB.PDBList(verbose=True)
    print("1. DOWNLOADING .CIF FILES")

    with open(ids_to_get_path) as pdb_textfile:
        for line in pdb_textfile:
            download = True
            # Get the info for this pdb
            pdb_id = line[:4]
            print(pdb_id)
            info = pypdb.get_all_info(pdb_id)
            # Initialize our list of lists for the species from which each chain comes
            species_for_chains = []
            # If "polymer" is not a list, make it one. If it has different number of 
            # entities than the number of monomers, do not continue with it
            if type(info["polymer"]) != type([]):
                info["polymer"] = [info["polymer"],]
            if len(info["polymer"]) != n_monomers:
                download = False
            # If we have the right number of entities, then get the species for each.
            for entity in info["polymer"]:
                # If there are is no "Taxonomy", we can't classify it; flag this pdb for discarding
                # AKA: not downloading
                if "Taxonomy" not in entity:
                    download = False
                    break
                # If there is only one species, turn it into a list
                if type(entity["Taxonomy"]) != type([]):
                    entity["Taxonomy"] = [entity["Taxonomy"],]
                # Initialize the list of species for this chain
                species_for_entity = []
                # Go through the species
                for species in entity["Taxonomy"]:
                    if species["@name"] not in species_for_entity: species_for_entity.append(species["@name"])
                species_for_chains.append(species_for_entity)
            # Get the species name, if there is one
            # NOTE: Right now, the chain-species mapping {"A":[H Sapiens, E Coli], "B":[H Sapiens, E Coli]}
            # Only returns H Sapiens as the species, to avoid duplicating data
            # Alternate approaches: Have it duplicated, OR add it to the species that has more/less structures
            if download:
                exists_path, species_name = ml.find_species(species_for_chains)
                if exists_path:
                    full_species_name = species_name
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
                        pdb_fn = pdbList.retrieve_pdb_file(pdb_id, pdir=download_folder_path, file_format="mmCif")
                    else:
                        print(f"Already have {pdb_id}.cif")
                # If we tripped download = False, or there is no path through the species_for_chains list, then
                # exclude that structure
            else:
                print(f"Excluding {pdb_id}.cif")
            print()