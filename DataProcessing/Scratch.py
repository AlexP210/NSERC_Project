import matplotlib.pyplot as plt
import os
import pandas as pd

for species_folder in os.listdir(r"E:\NSERC_Data\Heterotrimers\Data"):
    species_directory = os.path.join(r"E:\NSERC_Data\Heterotrimers\Data", species_folder)
    random_comparisons_path = os.path.join(species_directory, "RandomComparisons.csv")
    random_comparisons = pd.read_csv(random_comparisons_path)
    plt.hist(random_comparisons["Alignment_Score"])
    plt.show()

