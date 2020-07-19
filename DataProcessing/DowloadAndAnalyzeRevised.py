import os
import sys

# Script Folder
DataProcessing_dir = f"{os.path.dirname(__file__)}"

# Processing Scripts
ScrapePeriodicTable = os.path.join(DataProcessing_dir, "ScrapePeriodicTable_V2.py")
DownloadCIFs = os.path.join(DataProcessing_dir, "DownloadCIFs.py")
BiologicalAssemblies = os.path.join(DataProcessing_dir, "BiologicalAssemblies_V2.py")
PickRepresentativeAssembly = os.path.join(DataProcessing_dir, "PickRepresentativeAssembly.py")
RemoveDuplicates = os.path.join(DataProcessing_dir, "RemoveDuplicates_V2.py")
CleanStructures = os.path.join(DataProcessing_dir, "CleanStructures_V2.py")
ExtractChains = os.path.join(DataProcessing_dir, "ExtractChains_V3.py")
CompareChains = os.path.join(DataProcessing_dir, "CompareChains.py")
RemoveRedundantChains = os.path.join(DataProcessing_dir, "RemoveRedundantChains.py")
RandomComparisons = os.path.join(DataProcessing_dir, "RandomComparisons_V4.py")
FindHomologs = os.path.join(DataProcessing_dir, "FindSignificantPairs.py")
ExtractInterfaces = os.path.join(DataProcessing_dir, "ExtractInterfaces.py")
CompareInterfaces = os.path.join(DataProcessing_dir, "CompareInterfaces.py")

def nohupify(call):
    if sys.platform in ("linux", "linux2") and False:
        return "nohup " + call + " &"
    else: return call

# Main
if __name__ == "__main__":
    usage = f"\n$ python {os.path.basename(__file__)} <Periodic Table URL> <Current Working Directory> <Number of Monomers> <Max Sequence Similarity> <Max Structural Similarity> <Distance Cutoff> <Symmetries> <Testing>\n"
    if len(sys.argv) < 8:
        print(usage)
        sys.exit()
    url = sys.argv[1]
    cwd = sys.argv[2]
    if not os.path.exists(cwd): os.mkdir(cwd)
    n_monomers = sys.argv[3]
    max_sequence_similarity = sys.argv[4]
    max_structural_similarity = sys.argv[5]
    distance_cutoff = float(sys.argv[6])
    symmetry_groups = " ".join(sys.argv[7:-1])
    testing = bool(sys.argv[-1])

    # Activate the environment, get the cwd, set up folder for logs
    print("Initializing ...")
    os.system("activate Bioinformatics")
    logs_dir = os.path.join(cwd, "Logs")
    if not os.path.exists(logs_dir): os.mkdir(logs_dir)

    # Scrape the periodic table
    print("Scraping Periodic Table ...")
    call = f"python -u {ScrapePeriodicTable} '{url}' {os.path.join(cwd, 'IDs.txt')} true > {os.path.join(logs_dir, 'ScrapePeriodicTable_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    # if testing: input("Waiting ...")

    # Download the PDBs
    print("Downloading CIF Files ...")
    call = f"python -u {DownloadCIFs} {os.path.join(cwd, 'IDs.txt')} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'DownloadCIFs_Log.txt')}"
    call = nohupify(call)
    os.system(call)

#    if testing: input("Waiting ...")

    # Make the biological assemblies
    print("Forming Biological Assemblies ...")
    call = f"python -u {BiologicalAssemblies} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'BiologicalAssemblies_Log.txt')}"
    call = nohupify(call)
    os.system(call)

#     if testing: input("Waiting ...")

    # Pick the representative biological assembly for each PDB
    print("Picking Representative Assemblies ...")
    call = f"python -u {PickRepresentativeAssembly} {os.path.join(cwd, 'Data')} {n_monomers} > {os.path.join(logs_dir, 'PickRepresentativeAssembly_Log.txt')}"
    call = nohupify(call)
    os.system(call)

#     if testing: input("Waiting ...")

    # Remove the duplicated PDBs
    print("Removing Duplicates ...")
    call = f"python -u {RemoveDuplicates} {os.path.join(cwd, 'Data')} 1 1 > {os.path.join(logs_dir, 'RemoveDuplicates_Log.txt')}"
    call = nohupify(call)
    os.system(call)

#     if testing: input("Waiting ...")


    # Clean the PDBs
    print("Cleaning structures ...")
    call = f"python -u {CleanStructures} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'CleanStructures_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    # # if testing: input("Waiting ...")

    # # Extract the chains
    print("Extracting Chains ...")
    call = f"python -u {ExtractChains} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'ExtractChains_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    # # if testing: input("Waiting ...")

    # # Compare the chains
    print("Comparing Chains ...")
    call = f"python -u {CompareChains} {os.path.join(cwd, 'Data')} {symmetry_groups} > {os.path.join(logs_dir, 'CompareChains_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    # # if testing: input("Waiting ...")

    # # Filter Redundant Chains
    print("Filtering Redundant Chains ...")
    call = f"python -u {RemoveRedundantChains} {os.path.join(cwd, 'Data')} {max_sequence_similarity} {max_structural_similarity} > {os.path.join(logs_dir, 'FilterRedundantChains_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    # # if testing: input("Waiting ...")

    # # Random Comparisons
    print("Comparing random chains ...")
    call = f"python -u {RandomComparisons} {os.path.join(cwd, 'Data')} 10000 > {os.path.join(logs_dir, 'RandomComparisons_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    # # if testing: input("Waiting ...")

    # # Find the significant pairs
    print("Finding significant heteromers ...")
    call = f"python -u {FindHomologs} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'FindingSignificantHeteromers_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    # # if testing: input("Waiting ...")

    # # Extract interface
    print("Extracting interfaces ...")
    call = f"python -u {ExtractInterfaces} {os.path.join(cwd, 'Data')} {distance_cutoff} > {os.path.join(logs_dir, 'ExtractInterfaces_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    # #if testing: input("Waiting ...")

    # # Compare interface
    print("Comparing interfaces ...")
    call = f"python -u {CompareInterfaces} {os.path.join(cwd, 'Data')} {symmetry_groups} > {os.path.join(logs_dir, 'CompareInterfaces_Log.txt')}"
    call = nohupify(call)
    os.system(call)

    print("Done")
