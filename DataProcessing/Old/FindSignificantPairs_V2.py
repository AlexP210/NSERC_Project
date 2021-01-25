# Changes over V1:
# 1. Uses PCA

import os
import sys
import pandas as pd
import My_Library as ml
import matplotlib.pyplot as plt
import numpy as np
import sklearn.decomposition as sk
from mpl_toolkits.mplot3d import Axes3D
import math as m

if __name__ == "__main__":

    
    # Check the arguments
    usage = f"$ python {__file__.split('/')[-1]} <Species Root Directory>"
    if len(sys.argv) != 2:
        print(usage)
        sys.exit()
    root_directory = sys.argv[1]


    # Read in the random comparisons
    # random_comparisons_path = os.path.join(root_directory, "..", "RandomComparisons.csv")
    random_comparisons_path = r"C:\Users\alexp\NSERC_Project\DataProcessing\RandomComparisons.csv"
    random_comparisons = pd.read_csv(random_comparisons_path, encoding="UTF-8-sig")

    # Get the features from the random comparisons
    print("     Doing PCA")
    features = ["PID", "TM", "Alignment_Score", "Length_Difference"]
    # # Log'd
    # log_features = pd.DataFrame(columns=features)
    # for i in log_features:
    #     if i in ("Alignment_Score"):
    #         log_features[i] = np.log(random_comparisons[i] - min(random_comparisons[i]) + 1)
    #     if i in ("Length_Difference"):
    #         log_features[i] = np.log(random_comparisons[i] + 1)
    #     elif i in ("PID", "TM"):
    #         log_features[i] = np.log(max(random_comparisons[i]) - random_comparisons[i] + 1)
    

    # Now, the actual data
    for species_folder in os.listdir(root_directory):
        print(species_folder)
        # Set paths for this species folder
        species_directory = os.path.join(root_directory, species_folder)
        
        heteromer_comparisons_path = os.path.join(species_directory, "Chains.csv")

        # Import the data
        heteromer_comparisons = pd.read_csv(heteromer_comparisons_path)

        # If we actually have data to analyze
        if len(heteromer_comparisons.index) > 0:

            # heteromer_comparisons = ml.preprocess_data(heteromer_comparisons, features, pre_pca_means, pre_pca_stdevs)
            heteromer_comparisons = ml.statisticalize(pca, heteromer_comparisons, features, pre_pca_means, pre_pca_stdevs, post_pca_means, post_pca_stdevs)
            print(heteromer_comparisons)
        print("     Saving CSV")
        save_path = os.path.join(species_directory, "ChainsWithSignificances.csv")
        heteromer_comparisons.to_csv(save_path, encoding="UTF-8")


