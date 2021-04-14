import Bio.Blast.NCBIWWW as blast
import os
import My_Library as ml

if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) != 2:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <Species Root Directory>")
        sys.exit()
    root_directory = sys.argv[1]

    print(" ")
    for species_fol in os.listdir(root_directory):
        print(species_fol)
        structure_directory = os.path.join(root_directory, species_fol, "Duplicates_Removed")
        chains_directory = os.path.join(root_directory, species_fol, "Extracted_Chains")
        ids_to_chains = {}
        sequences = {}
        for filename in os.listdir(chains_directory):
            pdb_id, chain_id = filename.split(".")[0].split("_")
            if pdb_id not in ids_to_chains: ids_to_chains[pdb_id] = [chain_id,]
            else: ids_to_chains[pdb_id].append(chain_id)
            filepath = os.path.join(chains_directory, filename)
            sequences[f"{pdb_id}_{chain_id}"] = ml.load_sequence(filepath)
        for name, sequence in sequences.items():
            result = blast.qblast("blastp", "swissprot", sequence)
            print(result)


