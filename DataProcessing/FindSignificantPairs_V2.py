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

if __name__ == "__main__":

    
    # Check the arguments
    usage = f"$ python {__file__.split('/')[-1]} <Species Root Directory>"
    if len(sys.argv) != 2:
        print(usage)
        sys.exit()
    root_directory = sys.argv[1]


    # Read in the random comparisons
    random_comparisons_path = os.path.join(root_directory, "..", "RandomComparisons.csv")
    random_comparisons = pd.read_csv(random_comparisons_path, encoding="UTF-8-sig")

    # Get the CDFs from the random comparisons
    print("     Doing PCA")
    features = random_comparisons[["PID", "TM", "Alignment_Score"]]
    pca = sk.PCA()
    transformed_features = pca.fit_transform(features)

    # plt.scatter(transformed_features[:,0], transformed_features[:,1])
    # #plt.show()
    # plt.scatter(transformed_features[:, 0], transformed_features[:,2])
    # #plt.show()
    # plt.scatter(transformed_features[:,1], transformed_features[:,2])
    # plt.show()
    # plt.hist(transformed_features[:,0])
    # #plt.show()
    # plt.hist(transformed_features[:,1])
    # #plt.show()
    # plt.hist(transformed_features[:,2])
    # plt.show()

    fig = plt.figure()
    ax = Axes3D(fig)

    ax.scatter(transformed_features[:,0], transformed_features[:,1], transformed_features[:,2])
    plt.show()

    fig = plt.figure()
    ax = Axes3D(fig)

    ax.scatter(random_comparisons["PID"], random_comparisons["Alignment_Score"], random_comparisons["TM"])
    plt.show()
    
    # for species_folder in os.listdir(root_directory):
    #     print(species_folder)
    #     # Set paths for this species folder
    #     species_directory = os.path.join(root_directory, species_folder)
        
    #     heteromer_comparisons_path = os.path.join(species_directory, "Chains.csv")

    #     # Import the data
    #     heteromer_comparisons = pd.read_csv(heteromer_comparisons_path)

    #     if len(heteromer_comparisons.index) > 0:

    #         # Go through each of the heteromer alignments, and find the probability of getting something as significant, or more
    #         print("     Calculating P-Values")
    #         p_AlignmentScores = 1 - cdf_AlignmentScore(heteromer_comparisons["Alignment_Score"])
    #         p_PercentIdentities = 1 - cdf_PercentIdentity(heteromer_comparisons["PID"])
    #         p_TMScores = 1 - cdf_TMScore(heteromer_comparisons["TM"])

    #         # Add the significances to the dataframe
    #         heteromer_comparisons["AlignmentScore_pvalue"] = p_AlignmentScores
    #         heteromer_comparisons["PID_pvalue"] = p_PercentIdentities
    #         heteromer_comparisons["TM_pvalue"] = p_TMScores   
    
    #     # Save the new csv file
    #     print("     Saving CSV")
    #     save_path = os.path.join(species_directory, "ChainsWithSignificances.csv")
    #     heteromer_comparisons.to_csv(save_path, encoding="UTF-8")


