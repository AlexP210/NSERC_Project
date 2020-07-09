import csv
import pypdb
import sys
from Bio import PDB
import time
import os

if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) != 6:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <Paralogs.txt> <conversions.txt> <output.txt> <paralogs complexes dir> <chains.txt>\n") 
        sys.exit()
    paralogs_list_path = sys.argv[1]
    conversions_table_path = sys.argv[2]
    output_table_path = sys.argv[3]
    paralogs_saveplace = sys.argv[4]
    chains_textfile_path = sys.argv[5]


    # Find which yeast genes have 3D structure, and store their OLN to Uniprot accession code mapping
    print("--------------------- Loading Yeast Genes ---------------------")
    with open(conversions_table_path, "r", encoding="UTF-8-sig") as conversion_table:
        accession_codes = {}
        rows = csv.DictReader(conversion_table)
        for row in rows:
            if row["OLN"] not in accession_codes and int(row["3D"]) == 3:
                accession_codes[row["OLN"]] = row["S-P Num"]
        
    
    # Go through the Marchant paralogs, and for each pair find the complexes that contain both
    print("--------------------- Identifying Complexes Containing Paralogous Chains ---------------------")
    with open(paralogs_list_path, "r", encoding="UTF-8-sig") as paralogous_heteromers_list:
        rows = csv.DictReader(paralogous_heteromers_list)
        HET_complexes = {}
        for row in rows:
            if row["Both_3D"] == "TRUE":
                P1_query = pypdb.make_query(accession_codes[row["P1"]])
                P1_list = pypdb.do_search(P1_query)
                P2_query = pypdb.make_query(accession_codes[row["P2"]])
                P2_list = pypdb.do_search(P2_query)
                HET_complex_list = list(set(P1_list) & set(P2_list))
                if len(HET_complex_list) > 0:
                    pair_name = f"{accession_codes[row['P1']]}&{accession_codes[row['P2']]}"
                    HET_complexes[pair_name] = HET_complex_list

    # Write down the PDBs that are associated with each pair into output.txt
    with open(output_table_path, "w", encoding="UTF-8-sig") as output:
        writer = csv.DictWriter(output, ["P1_Accession", "P2_Accession", "PDB_IDs"])
        for pair, complex_list in HET_complexes.items():
            p1_acc, p2_acc = pair.split("&")
            writer.writerow({"P1_Accession": p1_acc, "P2_Accession": p2_acc, "PDB_IDs": ";".join(complex_list)})
    
    # Go through each HET complex, and count how many have only 2 chains => not many; need to use all,
    # and see which have an interface between P1 and P2
    print("--------------------- Downloading Complexes Containing Paralogous Chains ---------------------")
    with open(chains_textfile_path, "w", encoding="UTF-8-sig") as chains_textfile:
        downloaded = {}
        for pair, complex_list in HET_complexes.items():
            print()
            print(f"Downloading {pair}")
            dl = PDB.PDBList()
            for pdb in complex_list:
                print(f"     Trying {pdb}")
                dl.retrieve_pdb_file(pdb, pdir=paralogs_saveplace, file_format="pdb")
                time.sleep(3)
                if os.path.isfile(os.path.join(paralogs_saveplace, f"pdb{pdb.lower()}.ent")):
                    downloaded[pair] = pdb
                    break
        print("Recording labels for paralogous chains...")
        for pair, pdb in downloaded.items():
            info = pypdb.get_all_info(pdb)
            for entity in info["polymer"]:
                if type(entity) == type({}) and entity["@type"] == "protein" and "macroMolecule" in entity.keys() and "accession" in entity["macroMolecule"].keys() and "@id" in entity["macroMolecule"]["accession"].keys():
                    if entity["macroMolecule"]["accession"]["@id"].lower() in [i.lower() for i in pair.split("&")]:
                        if type(entity["chain"]) == type({}):
                            chains_textfile.write(f"{pdb.upper()}{entity['chain']['@id'].upper()}\n")
                        # If same chain appears multiple times, we take only first occurrence (i.e. chain letter A)
                        # Might want to take all; we'll see
                        elif type(entity["chain"]) == type([]):
                            chains_textfile.write(f"{pdb.upper()}{entity['chain'][0]['@id'][0].upper()}\n")



                
