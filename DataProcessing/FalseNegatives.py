import pandas as pd

labels = pd.read_csv(r"E:\NSERC_Data\Test\Homo_Sapiens\Labels.csv")
data = pd.read_csv(r"E:\NSERC_Data\Test\Homo_Sapiens\ChainsWithSignificances2.csv")
# data["Significant"] = (data["AlignmentScoreAdjusted_pvalue"] < 0.01) & (data["TM_pvalue"] < 0.01)
# data["Significant"] = data["AlignmentScoreAdjusted_pvalue"] < 0.01
data["Significant"] = data["TM_pvalue"] < 0.01


print(data[data["PDB_ID"] == "1jkg"]["AlignmentScoreAdjusted_pvalue"])
print(data[data["PDB_ID"] == "1jkg"]["TM_pvalue"])



print(f"False Negatives are:")

for index, row in labels.iterrows():
    pdb_id = row["PDB_ID"]
    homologous = row["Homologous"]
    if homologous and not data[data["PDB_ID"] == pdb_id]["Significant"].all():
        print(pdb_id)

print()

print(f"False Positives are:")
for index, row in labels.iterrows():
    pdb_id = row["PDB_ID"]
    homologous = row["Homologous"]
    if not homologous and data[data["PDB_ID"] == pdb_id]["Significant"].all():
        print(pdb_id)
