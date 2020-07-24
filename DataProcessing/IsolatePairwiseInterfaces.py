import My_Library as ml
import os
import sys

if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) < 3:
        print(f"$ Usage: python {os.path.basename(__file__)} <Species Root Directory> <Distance Cutoff>")
        sys.exit()
    root_directory = sys.argv[1]
    distance_cutoff = float(sys.argv[2])

    # Go through each species folder
    for species_folder in os.listdir(root_directory):
        print(species_folder)

        # Set the directories
        species_directory = os.path.join(root_directory, species_folder)
        clean_structures_directory = os.path.join(species_directory, "Clean_Structures")
        pairwise_interfaces_directory = os.path.join(species_directory, "Pairwise_Interfaces")
        if not os.path.exists(pairwise_interfaces_directory): os.mkdir(pairwise_interfaces_directory)

        # for each structure, get the chain ids, and take all pairs
        for clean_structure_filename in os.listdir(clean_structures_directory):
            extension = os.path.basename(clean_structure_filename).split(".")[-1]
            pdb_id = clean_structure_filename.split(".")[0]
            print(f"    {pdb_id}")
            clean_structure_path = os.path.join(clean_structures_directory, clean_structure_filename)
            chain_letters = ml.get_chain_ids(clean_structure_path)
            for a in range(1, len(chain_letters)):
                for b in range(a):
                    chain_letter_a = chain_letters[a]
                    chain_letter_b = chain_letters[b]
                    print(f"        {chain_letter_a}{chain_letter_b}")
                    pairwise_interfaces_filename = f"{pdb_id}{chain_letter_a}{chain_letter_b}.{extension}"
                    pairwise_interface_path = os.path.join( pairwise_interfaces_directory, pairwise_interfaces_filename )
                    ml.get_interfaces(clean_structure_path, distance_cutoff, pairwise_interface_path, [chain_letter_a, chain_letter_b])
