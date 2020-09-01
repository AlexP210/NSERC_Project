# Changes from V1:
# 1. Only saves the PDB IDs of complexes in the periodic table. Then, BiologicalAssemblies will tell us which chains we need.

from bs4 import BeautifulSoup
import requests
import pypdb
import os
import sys

# Change this for use with non-dimers
n_monomers = 2

def get_pdbIDs(url, strict):
    """Gets the PDB IDs for structures at periodicproteincomplexes.org, at the given URL.

    :param url: URL to scrape from at periodicproteincomplexes.org
    :type url: string (url)
    :yield: PDB ID tag
    :rtype: string
    """
    html = requests.get(url).text
    soup = BeautifulSoup(html, "lxml")
    proteins = soup.findAll("td")
    ids = []
    for protein in proteins:
        ids.append(protein.text.lower())
    return ids

def find_chains(pdbID, strict):
    """Finds the PDB chain labels for a PDB structure. Uses PyPDB to make a information query to PDB.

    :param pdbID: PDB ID for a structure
    :type pdbID: string
    :param strict: Whether or not to omit PDBs that have more subunits than the specified n_monomers, defaults to False
    :type strict: bool, optional
    :return: A list of chain tags
    :rtype: list of strings
    """
    global n_monomers

    polymers = pypdb.get_all_info(pdbID)["polymer"]
    # Next line checks if the PDB file listed in the periodic table is really what we need. Can be invalid in a few ways:
    # 1. Has too few monomers: Definitely invalid because can't do all the comparisons needed
    # 2. Has too many monomers: Can be invalid depending on how strict we are about dimers. If we are ok with a small extra chain, take the longest 2 chains listed. Control which is done through "strict"
    # 3. Has only 1 chain: In this case, it's definitely invalid. Due to the way the results are returned, polymers will be a dictionary in this case rather than a singleton list
    if len(polymers) < n_monomers or (len(polymers) > n_monomers and strict) or type(polymers) == type({}):
        return None
    chains = []
    for entity_dict in sorted(polymers, key=lambda d: int(d["@length"]), reverse=True)[0:2]:
        if len(entity_dict["chain"]) > 1: chains.append(entity_dict["chain"][0]["@id"])
        else: chains.append(entity_dict["chain"]["@id"])
    return chains
    
if __name__ == "__main__":

    if len(sys.argv) != 4: 
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <url> <chains.txt> <strict>\n")
        sys.exit()
    url_raw = sys.argv[1]
    # This is new - powershell does not know how to pass in the "url" parameter without the double quotes
    url = ""
    for c in url_raw:
        if c != "'": url += c

    outfile = sys.argv[2]
    strict = bool(sys.argv[3].lower())

    print("1. SCRAPING IDs ...")
    ids = get_pdbIDs(url, strict)
    print("2. WRITING IDs ... ")
    with open(outfile, "w") as out:
        for Id in ids:
            print(Id)
            out.write(Id + "\n")
                




