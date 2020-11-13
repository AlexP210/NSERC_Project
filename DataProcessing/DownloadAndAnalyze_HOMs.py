import os
import sys

# Script Folder
DataProcessing_dir = f"{os.path.dirname(__file__)}"

# Processing Scripts
# Set-Up
ScrapePeriodicTable = os.path.join(DataProcessing_dir, "ScrapePeriodicTable_V2.py")
DownloadCIFs = os.path.join(DataProcessing_dir, "DownloadCIFs_V2.py") # CIF
RenameFolders = os.path.join(DataProcessing_dir, "RenameFolders.py")
BiologicalAssemblies = os.path.join(DataProcessing_dir, "BiologicalAssemblies_V2.py") # CIF
PickRepresentativeAssembly = os.path.join(DataProcessing_dir, "PickRepresentativeAssembly.py") # CIF
RemoveDuplicates = os.path.join(DataProcessing_dir, "RemoveDuplicates_V2.py") # CIF
CleanStructures = os.path.join(DataProcessing_dir, "CleanStructures_HOMs.py") # CIF
# Compare the Structures
CompareHOMs = os.path.join(DataProcessing_dir, "CompareHOMs.py")
RemoveRedundantHOMs = os.path.join(DataProcessing_dir, "RemoveRedundantHOMs.py")
RandomComparisons_HOMs = os.path.join(DataProcessing_dir, "RandomComparisons_HOMs.py")
FindSignificantPairs_HOMS = os.path.join(DataProcessing_dir, "FindSignificantPairs_HOMs.py")


# Compare the Protein Chains
SeparateChains = os.path.join(DataProcessing_dir, "SeparateChains.py") # PDB
CompareChains = os.path.join(DataProcessing_dir, "CompareChains.py")
RemoveRedundantChains = os.path.join(DataProcessing_dir, "RemoveRedundantChains_V2.py") # PDB
RandomComparisons = os.path.join(DataProcessing_dir, "RandomComparisons_V4.py")
FindSignificantPairs = os.path.join(DataProcessing_dir, "FindSignificantPairs.py")
# Compare the whole-complex interfaces
IsolateComplexInterfaces = os.path.join(DataProcessing_dir, "IsolateComplexInterfaces.py")
RemoveEmptyComplexInterfaces = os.path.join(DataProcessing_dir, "RemoveEmptyInterfaces.py")
SeparateComplexInterfaces = os.path.join(DataProcessing_dir, "SeparateComplexInterfaces.py")
CompareComplexInterfaces = os.path.join(DataProcessing_dir, "CompareComplexInterfaces.py")
RemoveRedundantComplexInterfaces = os.path.join(DataProcessing_dir, "RemoveRedundantComplexInterfaces.py")
RandomComplexInterfaceComparisonsAndStatistics = os.path.join(DataProcessing_dir, "RandomComplexInterfaceComparisonsAndStatistics.py")
# # Compare the pairwise-interfaces for each pair of chains
# NOTE: Not yet ready
# IsolatePairwiseInterfaces = os.path.join(DataProcessing_dir, "IsolatePairwiseInterfaces.py")
# RemoveEmptyPairwiseInterfaces = None
# SeparatePairwiseInterfaces = os.path.join(DataProcessing_dir, "SeparatePairwiseInterfaces.py")

# Main
if __name__ == "__main__":
    usage = f"\n$ python {os.path.basename(__file__)} <Periodic Table URL> <Current Working Directory> <Number of Monomers> <Max Sequence Similarity> <Max Structural Similarity> <Distance Cutoff> <Symmetries>\n"
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
    if len(sys.argv) > 8: symmetry_groups = " ".join(sys.argv[7:-1])
    elif len(sys.argv) == 8: symmetry_groups = sys.argv[7]

    # Activate the environment, get the cwd, set up folder for logs
    print("Initializing ...")
    os.system("activate Bioinformatics")
    logs_dir = os.path.join(cwd, "Logs")
    if not os.path.exists(logs_dir): os.mkdir(logs_dir)

    # # Scrape the periodic table
    # print("Scraping Periodic Table ...")
    # call = f"python -u {ScrapePeriodicTable} '{url}' {os.path.join(cwd, 'IDs.txt')} true > {os.path.join(logs_dir, 'ScrapePeriodicTable_Log.txt')}"
    # os.system(call)

    # # Download the PDBs
    # print("Downloading CIF Files ...")
    # call = f"python -u {DownloadCIFs} {os.path.join(cwd, 'IDs.txt')} {os.path.join(cwd, 'Data')} {n_monomers} > {os.path.join(logs_dir, 'DownloadCIFs_Log.txt')}"
    # os.system(call)

    # Sanitize the Folder names
    print("Sanitizing Folder Names ...")
    call = f"python -u {RenameFolders} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'RenameFolders_Log.txt')}"
    os.system(call)

    # Make the biological assemblies
    print("Forming Biological Assemblies ...")
    call = f"python -u {BiologicalAssemblies} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'BiologicalAssemblies_Log.txt')}"
    os.system(call)

    # Pick the representative biological assembly for each PDB
    print("Picking Representative Assemblies ...")
    call = f"python -u {PickRepresentativeAssembly} {os.path.join(cwd, 'Data')} {n_monomers} > {os.path.join(logs_dir, 'PickRepresentativeAssembly_Log.txt')}"
    os.system(call)

    # Remove the duplicated PDBs
    print("Removing Duplicates ...")
    call = f"python -u {RemoveDuplicates} {os.path.join(cwd, 'Data')} 1 1 > {os.path.join(logs_dir, 'RemoveDuplicates_Log.txt')}"
    os.system(call)

    # Clean the PDBs
    print("Cleaning structures ...")
    call = f"python -u {CleanStructures} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'CleanStructures_Log.txt')}"
    os.system(call)

    # # Compare the HOMs
    print("Comparing HOMs ...")
    call = f"python -u {CompareHOMs} {os.path.join(cwd, 'Data')} {symmetry_groups} > {os.path.join(logs_dir, 'CompareHOMs_Log.txt')}"
    os.system(call)

    # # Filter Redundant Chains
    print("Filtering Redundant HOMs ...")
    call = f"python -u {RemoveRedundantHOMs} {os.path.join(cwd, 'Data')} {max_sequence_similarity} {max_structural_similarity} > {os.path.join(logs_dir, 'FilterRedundantHOMs_Log.txt')}"
    os.system(call)

    # # Random Comparisons
    print("Comparing random HOMs ...")
    call = f"python -u {RandomComparisons_HOMs} {os.path.join(cwd, 'Data')} 10000 > {os.path.join(logs_dir, 'RandomComparisonsHOMs_Log.txt')}"
    os.system(call)

    # # Find Significant Pairs
    print("Finding significant pairs of HOMs ...")
    call = f"python -u {FindSignificantPairs_HOMS} {os.path.join(cwd, 'Data')} > {os.path.join(logs_dir, 'FindSignificantPairsHOMs_Log.txt')}"
    os.system(call)


    # # # Isolate complex interfaces
    # print("Extracting interfaces ...")
    # call = f"python -u {IsolateComplexInterfaces} {os.path.join(cwd, 'Data')} {distance_cutoff} > {os.path.join(logs_dir, 'IsolateComplexInterfaces_Log.txt')}"
    # os.system(call)

    # # # Separate complex interfaces
    # print("Separating interfaces ...")
    # call = f"python -u {SeparateComplexInterfaces} {os.path.join(cwd, 'Data')} {distance_cutoff} > {os.path.join(logs_dir, 'SeparateComplexInterfaces_Log.txt')}"
    # os.system(call)


    # # # Compare interfaces
    # print("Comparing interfaces ...")
    # call = f"python -u {CompareComplexInterfaces} {os.path.join(cwd, 'Data')} {symmetry_groups} > {os.path.join(logs_dir, 'CompareComplexInterfaces_Log.txt')}"
    # os.system(call)

    print("Done")
