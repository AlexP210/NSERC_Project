import os
import sys
import pandas as pd
import My_Library as ml
import matplotlib.pyplot as plt
import numpy as np
import sklearn.decomposition as sk

if __name__ == "__main__":

    
    # Check the arguments
    usage = f"$ python {__file__.split('/')[-1]} <Species Root Directory>"
    if len(sys.argv) != 2:
        print(usage)
        sys.exit()
    root_directory = sys.argv[1]


    # Read in the random comparisons
    random_comparisons_path = os.path.join(root_directory, "..", "RandomComparisons.csv")
    random_comparisons = pd.read_csv(random_comparisons_path)

    # plt.scatter(random_comparisons["Alignment_Score"], random_comparisons["PID"])
    # plt.show()
    # plt.scatter(random_comparisons["Alignment_Score"], random_comparisons["TM"])
    # plt.show()
    # plt.scatter(random_comparisons["TM"], random_comparisons["PID"])
    # plt.show()

    features = ["Alignment_Score", "PID", "TM"]
    feature_frame = random_comparisons[features]
    pca = sk.PCA(n_components=3)
    pca.fit(feature_frame)
    fitted_data = pca.fit_transform(feature_frame)
    print(fitted_data)

    plt.scatter(fitted_data[0,:], fitted_data[1,:])
    plt.show()
    plt.scatter(fitted_data[0,:], fitted_data[2,:])
    plt.show()
    plt.scatter(fitted_data[2,:], fitted_data[1,:])
    plt.show()


