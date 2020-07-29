import os
import sys
import pandas as pd
import My_Library as ml
import matplotlib.pyplot as plt
import numpy as np
import sklearn.decomposition as sk
from Bio.PDB import MMCIFParser

if __name__ == "__main__":
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("", r"E:\NSERC_Data\2C2\Data\Mus_musculus\Clean_Structures\1a0q.cif")
    for residue in structure.get_residues():
        print(residue)
        for atom in residue.get_atoms():
            print(atom)
        break


