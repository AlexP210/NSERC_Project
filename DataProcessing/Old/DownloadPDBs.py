import os
from Bio import PDB

if __name__ == "__main__":
    """ Downloads PDB files listed for PDB IDs listed in chains.txt """
    import sys
    if not len(sys.argv) == 3:
        print(f"Usage: python DownloadPDBs.py <chainsList.txt> <PDB Directory>")
        sys.exit()

    chains_to_get_path = sys.argv[1]
    pdb_download_path = sys.argv[2]
    if not os.path.exists(pdb_download_path): os.mkdir(pdb_download_path)

    pdbList = PDB.PDBList()

    with open(chains_to_get_path) as pdb_textfile:
        for line in pdb_textfile:
            pdb_id = line[:4].lower()
            chain = line[4].upper()
            # if we don't have the pdb, download it
            if not os.path.isfile(os.path.join(pdb_download_path, f"pdb{pdb_id}.ent")):
                pdb_fn = pdbList.retrieve_pdb_file(pdb_id, file_format = "pdb", pdir=pdb_download_path)