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
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    structure = parser.get_structure("1aif", r"E:\NSERC_Data\Test\Test_Species\Complex_Interfaces\1aif.cif")
    io.set_structure(structure)
    res_list = []
    for residue in structure.get_residues():
        residue_dict =  {
                "Resname":residue.get_resname(),
                "Seq":residue.get_id()[1],
                "Chain":residue.get_parent().get_id()
            }
        print("My interface: ", residue_dict)
        res_list.append(residue_dict)
    with open(r"C:\Users\alexp\Desktop\School\3rd_Year\NSERC\Scratch\1aif_cocomaps_residues.csv", "r", encoding="UTF-8-sig") as cocomaps:
        rows = csv.DictReader(cocomaps)
        unique_cocomaps = []
        for row in rows:
            residue = {"Resname":row["Resname"],"Seq":int(row["Seq"]),"Chain":row["Chain"]}
            print("Cocomaps interface: ", residue)
            if residue not in res_list:
                assert False
            if residue not in unique_cocomaps:
                unique_cocomaps.append(residue)

    print("Length of my interface: ", len(res_list))
    print("Length of cocomaps interface: ", len(unique_cocomaps))

    
    
    # for residue in structure.get_residues():
    #     if residue.get_id()[1] == 39 and residue.get_resname() == "GLN" and residue.get_parent().get_id() == "H":
    #         res1 = residue
    #         print(res1, res1.get_parent().get_id())
    #     if residue.get_id()[1] == 45 and residue.get_resname() == "PRO" and residue.get_parent().get_id() == "L":
    #         res2 = residue
    #         print(res2, res2.get_parent().get_id())
    # min_dist = float("inf")
    # a1 = None
    # a2 = None
    # for atom1 in res1.get_atoms():
    #     for atom2 in res2.get_atoms():
    #         if np.linalg.norm(atom1.get_coord() - atom2.get_coord()) < min_dist:
    #             min_dist = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
    #             a1 = atom1
    #             a2 = atom2
    # print(min_dist, a1, a2)

    #io.save(r"E:\NSERC_Data\Test\Test_Species\Complex_Interfaces\1a0qH.ent", select=SelectIntefrace(structure))
                


            
            




    


