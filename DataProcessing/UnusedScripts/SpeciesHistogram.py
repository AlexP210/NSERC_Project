import pypdb
import seaborn as sns
import matplotlib.pyplot as plt

heterotrimer_list_path = "E:\\NSERC_Data\\Heterotrimers_C3\\Heterotrimers.txt"
with open(heterotrimer_list_path, "r") as f:
    heterotrimer_list = f.readlines()

species_dict = {}
n_lists = 0
n_diff = 0
frequencies = {}
for pdb_id in heterotrimer_list:
    print(pdb_id[:4])
    try:
        info = pypdb.get_all_info(pdb_id[:4])
        species_names = set()
        for entity in info["polymer"]:
            if type(entity["Taxonomy"]) == type({}):
                species_names.add(entity["Taxonomy"]["@name"].lower())
            else:
                print("Type is list...")
                n_lists += 1
        
        if len(species_names) == 1:
            for species_name in species_names:
                print(f"     {species_name}")
                species_dict[pdb_id] = species_name
                if species_name in frequencies:
                    frequencies[species_name] += 1
                else:
                    frequencies[species_name] = 1
                
        else:
            print(f"Different species: {species_names}")
            n_diff += 1
    except KeyError:
        pass

print(frequencies)
print(species_names)
plt.bar(x=list(frequencies.keys()), height=list(frequencies.values()))
plt.show()
print(f"Number of different: {n_diff}")
print(f"Number of lists: {n_lists}")
    


