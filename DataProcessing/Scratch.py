# import os
# import sys
# import pandas as pd
# import My_Library as ml
# import matplotlib.pyplot as plt
# import numpy as np
# import sklearn.decomposition as sk
# from Bio.PDB import MMCIFParser

# if __name__ == "__main__":

#     c = 0
#     for species_folder in os.listdir(r"E:\NSERC_Data\Heterotrimers_07_12\Heterotrimer_07-12\Heterotrimers\Data"):
#         species_directory = os.path.join(r"E:\NSERC_Data\Heterotrimers_07_12\Heterotrimer_07-12\Heterotrimers\Data", species_folder)
#         downloaded_structures_directory = os.path.join(species_directory, "Selected_Biological_Assemblies")
#         for downloaded_structure_filename in os.listdir(downloaded_structures_directory):
#             print(species_folder, downloaded_structure_filename)
#             downloaded_structure_path = os.path.join(downloaded_structures_directory, downloaded_structure_filename)
#             parser = MMCIFParser(QUIET=True)
#             structure = parser.get_structure(f"{downloaded_structure_filename.split('.')[0]}", downloaded_structure_path)
#             if len(list(structure.get_chains())) == 3:
#                 c += 1
#                 print(len(ml.load_sequence(downloaded_structure_path)))
#                 for chain in structure.get_chains():
#                     print(chain.get_id())
#                     #print(ml.load_chain_sequence(structure, chain.get_id()))
#                     print()
#             if c == 3:
#                 assert False
            
            

x = 3
def test():
    global x
    x += 1
print(x)

