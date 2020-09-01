import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc

homo_sapiens_directory = r"E:\NSERC_Data\Test\Homo_Sapiens"

method_1 = pd.read_csv(os.path.join(homo_sapiens_directory, "ChainsWithSignificances.csv"))
method_2 = pd.read_csv(os.path.join(homo_sapiens_directory, "ChainsWithSignificances2.csv"))
labels = pd.read_csv(os.path.join(homo_sapiens_directory, "Labels.csv"))
labels_dict = {}
all_positives = 0
all_negatives = 0
for index, row in labels.iterrows():
    pdb_id = row["PDB_ID"]
    homologous = row["Homologous"]
    labels_dict[pdb_id] = homologous
    if homologous:
        all_positives += 1
    else:
        all_negatives += 1

thresholds = np.linspace(0, 1, 200)
method_1_TPRs = np.zeros(shape=thresholds.shape)
method_1_FPRs = np.zeros(shape=thresholds.shape)
method_2_TPRs = np.zeros(shape=thresholds.shape)
method_2_FPRs = np.zeros(shape=thresholds.shape)


for threshold_idx in range(thresholds.size):
    threshold = thresholds[threshold_idx]
    print(threshold)
    method_1["Significance"] = method_1["p_value"] < threshold
    method_1_dict = {}
    for index, row in method_1.iterrows():
        method_1_dict[row["PDB_ID"]] = row["Significance"]
    # method_2["Significance"] = (method_2["TM_pvalue"] < threshold) & (method_2["AlignmentScoreAdjusted_pvalue"] < threshold)
    # method_2["Significance"] = method_2["AlignmentScoreAdjusted_pvalue"] < threshold   
    # method_2["Significance"] = method_2["TM_pvalue"] < threshold

    method_2_dict = {}
    for index, row in method_2.iterrows():
        method_2_dict[row["PDB_ID"]] = row["Significance"]

    method_1_TPs = 0
    method_1_FPs = 0
    method_2_TPs = 0
    method_2_FPs = 0
    for pdb_id, label in labels_dict.items():
        if label == method_1_dict[pdb_id] and label: method_1_TPs += 1
        elif label != method_1_dict[pdb_id] and not label: method_1_FPs += 1
        if label == method_2_dict[pdb_id] and label: method_2_TPs += 1
        elif label != method_2_dict[pdb_id] and not label: method_2_FPs += 1
    
    method_1_TPRs[threshold_idx] = method_1_TPs / all_positives
    method_1_FPRs[threshold_idx] = method_1_FPs / all_negatives

    method_2_TPRs[threshold_idx] = method_2_TPs / all_positives
    method_2_FPRs[threshold_idx] = method_2_FPs / all_negatives

plt.plot(method_1_FPRs, method_1_TPRs, label=f"PCA Method (AUC={round(auc(method_1_FPRs, method_1_TPRs), 3)}")
plt.plot(method_2_FPRs, method_2_TPRs, label=f"AND Method (AUC={round(auc(method_2_FPRs, method_2_TPRs), 3)}")
plt.xlabel("False Positive Rates")
plt.ylabel("True Positive Rates")
plt.title("ROC Curves for PCA & AND Methods")
plt.legend()
plt.show()

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('TPR and FPR vs Significance Threshold')

ax1.plot(thresholds[:10], method_1_TPRs[:10], label="True Positive Rate")
ax1.plot(thresholds[:10], method_1_FPRs[:10], label="False Positive Rate")
ax1.set(xlabel="Significance Threshold", ylabel="True Positive Rate / False Positive Rate")
ax1.set_title("PCA Method")
ax1.legend()

ax2.plot(thresholds[:10], method_2_TPRs[:10], label="True Positive Rate")
ax2.plot(thresholds[:10], method_2_FPRs[:10], label="False Positive Rate")
ax2.set(xlabel="Significance Threshold", ylabel="True Positive Rate / False Positive Rate")
ax2.set_title("AND Method")
ax2.legend()

plt.show()



    

