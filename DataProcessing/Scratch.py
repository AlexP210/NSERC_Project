import os
import sys
import pandas as pd
import My_Library as ml
import matplotlib.pyplot as plt
import numpy as np
import sklearn.decomposition as sk
from Bio.PDB import MMCIFParser, NeighborSearch, MMCIFIO, Select, PDBIO, PDBParser
import csv

if __name__ == "__main__":
    first_results = set()
    second_results = set()
    t1_hist = {}
    t2_hist = {}
    t3_hist = {}
    for species_folder in os.listdir(r"E:\NSERC_Data\Heterotrimers_07_12\Heterotrimer_07-12\Heterotrimers\Data"):
        species_directory = os.path.join(r"E:\NSERC_Data\Heterotrimers_07_12\Heterotrimer_07-12\Heterotrimers\Data", species_folder)
        chains_csv = os.path.join(species_directory, "ChainsWithSignificances2.csv")
        chains = pd.read_csv(chains_csv)
        t1_count = 0
        t2_count = 0
        t3_count = 0
        chains_hist = {}
        for index, row in chains.iterrows():
            pdb_id = row["PDB_ID"]
            if pdb_id not in chains_hist:
                chains_hist[pdb_id] = [(row["AlignmentScoreAdjusted_pvalue"] < 0.01) and (row["TM_pvalue"] < 0.01),]
            elif pdb_id in chains_hist:
                chains_hist[pdb_id].append((row["AlignmentScoreAdjusted_pvalue"] < 0.01) and (row["TM_pvalue"] < 0.01))
        for pdb_id, l in chains_hist.items():
            if sum(l) == 1: t1_count += 1
            if sum(l) == 2: t2_count += 1
            if sum(l) == 3: t3_count += 1
        t1_hist[species_folder] = t1_count
        t2_hist[species_folder] = t2_count
        t3_hist[species_folder] = t3_count

    print(sum(t1_hist.values()))
    print(sum(t2_hist.values()))
    print(sum(t3_hist.values()))


        



    #     try:
    #         print(species_folder)
    #         species_directory = os.path.join(r"E:\NSERC_Data\Heterotrimers_07_12\Heterotrimer_07-12\Heterotrimers\Data", species_folder)
    #         first_path = os.path.join(species_directory, "ChainsWithSignificances.csv")
    #         second_path = os.path.join(species_directory, "ChainsWithSignificances2.csv")
    #         first = pd.read_csv(first_path)
    #         second = pd.read_csv(second_path)
    #         for row in first[first["p_value"] <= 0.01].iterrows():
    #             i = row[1]
    #             name1 = f"{species_folder}_{i['PDB_ID']}_{i['Chain_1']}{i['Chain_2']}"
    #             name2 = f"{species_folder}_{i['PDB_ID']}_{i['Chain_2']}{i['Chain_1']}"
    #             if not name1 in first_results or name2 in first_results:
    #                 first_results.add(name1)
    #         for row in second[(second["AlignmentScoreAdjusted_pvalue"] <= 0.05) & (second["TM_pvalue"] <= 0.05)].iterrows():
    #             i = row[1]
    #             name1 = f"{species_folder}_{i['PDB_ID']}_{i['Chain_1']}{i['Chain_2']}"
    #             name2 = f"{species_folder}_{i['PDB_ID']}_{i['Chain_2']}{i['Chain_1']}"
    #             if not name1 in second_results or name2 in second_results:
    #                 second_results.add(name1)
    #     except KeyError:
    #         print("     Empty")
    # print(f"Significant with PCA, but not 'and':")
    # print(first_results - second_results)
    # print()
    # print(f"Significant with 'and', but not PCA:")
    # print(second_results - first_results)

            
            




    


