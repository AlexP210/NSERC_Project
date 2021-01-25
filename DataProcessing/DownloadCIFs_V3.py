### Changes Over V2 ###
# 1. Modified to account for new PDB API format. May change as PyPDB library develops in response to the changes.

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
            info = pypdb.get_all_info(pdb_id)
            # Initialize our container; is a list of list, showing the source species for each chain
            species_for_chains = []
            # If "polymer" is not a list, make it one. If it has different number of 
            # entities than the number of monomers, do not continue with it
            if int(info["rcsb_entry_info"]["polymer_entity_count_protein"]) != n_monomers:
                download = False
            # Get the number of assemblies for this structure
            n_assemblies = info["rcsb_entry_container_identifiers"]["assembly_ids"]
            # Get the names of the entities in this structure (These are NOT NECESSARILY the chains- could be licands, nucleic acids, etc)
            entities = info["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
            # Go through each entity, and check the species of the protein entities
            for entity in entities:
                entity_info = pypdb.get_all_info(f"{pdb_id}/{entity}", url_root="http://data.rcsb.org/rest/v1/core/polymer_entity/")
                # Verify the current entity is a protein before we go on
                if entity_info["entity_poly"]["rcsb_entity_polymer_type"] != "Protein": continue
                # Go through each source organism for this entity, and if it's a natural source, save it
                species_for_entity = []
                number_of_natural_chains = 0
                for source in entity_info["rcsb_entity_source_organism"]:
                    if source["source_type"] == "natural":
                        species_for_entity.append(source["ncbi_scientific_name"])
                        number_of_natural_chains += 1
                species_for_chains.append(species_for_entity)
            # Check if the number of natural chains matches the number of monomers of the symmetry group
            if number_of_natural_chains != n_monomers: download = False
            # Check if there is one organism in which all chains are found
            exists_path, species_name = ml.find_species(species_for_chains)
            if not exists_path: download = False

            # Get the species name, if there is one
            # NOTE: Right now, the chain-species mapping {"A":[H Sapiens, E Coli], "B":[H Sapiens, E Coli]}
            # Only returns H Sapiens as the species, to avoid duplicating data
            # Alternate approaches: Have it duplicated, OR add it to the species that has more/less structures
            if download:
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
                # If we tripped download = False, then exclude that structure
            else:
                print(f"Excluding {pdb_id}.cif")
            print()