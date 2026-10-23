import pandas as pd
import itertools


# Load each source and reduce it to a set of uppercase PDB codes
def pdb_set(filename):
    df = pd.read_csv(filename)
    return set(df["pdb"].str.upper())


files = {
    "psiblast": "Results/psiblast_hits_annotated.csv",
    "pfam": "Results/Interpro_PFam_FN3_homologs.csv",
    "cpdf": "Results/FN3_homologs_CPDF.csv",
    "interpro": "Results/Interpro_FN3_homologs.csv",
    "prosite": "Results/Interpro_Prosite_FN3_homologs.csv",
}

sets = {}
for name, filename in files.items():
    sets[name] = pdb_set(filename)
    print("{name}", name, " had :", len(sets[name]))

# Every PDB code seen anywhere.
# union(*sets.values()) combines all the sets at once.
all_pdbs = set().union(*sets.values())

# One row per PDB code, with a True/False column per source
rows = []
for pdb in sorted(all_pdbs):
    row = {"pdb": pdb}
    for name in sets:
        row["in_" + name] = pdb in sets[name]
    rows.append(row)

overlap = pd.DataFrame(rows)

# The source columns, in the order they were listed
source_columns = ["in_" + name for name in sets]

# True counts as 1, so this adds up how many sources found each entry
overlap["n_sources"] = overlap[source_columns].sum(axis=1)

overlap = overlap.sort_values("n_sources", ascending=False)
overlap.to_csv("pdb_overlap.csv", index=False)

print(len(overlap), "unique PDB codes across", len(sets), "sources\n")
# print(overlap["n_sources"].value_counts())

# How many entries each source found alone
for name in sets:
    only = overlap[(overlap["n_sources"] == 1) & overlap["in_" + name]]
    print(name, "only:", len(only))


# combinations gives every unordered pair from the list of names
for a, b in itertools.combinations(sets, 2):
    print(a, "&", b, ":", len(sets[a] & sets[b]))

# And the entries every source agrees on
print("in all sources:", len(set.intersection(*sets.values())))
