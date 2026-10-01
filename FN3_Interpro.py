# ---------- Libraries ----------

import requests
import json

# Lets me pause between requests
import time
import os
import pandas as pd

# ---------- InterPro Data to a CSV ----------
with open("Results/structure-matching-IPR003961.json") as file:
    interPro_results = json.load(file)

rows = []

for record in interPro_results:
    meta = record["metadata"]

    for entry in record["entries"]:
        # PFam has e-values for each copy of the protein domain
        locations = entry.get("entry_protein_locations")
        if locations is None:
            locations = []

        n_domains_found = len(locations)

        rows.append(
            {
                "pdb": meta["accession"].upper(),
                "chain": entry.get("chain"),
                "description": meta["name"],
                "method": meta["experiment_type"],
                "resolution": meta["resolution"],
                "source": "interpro_IPR003961",
            }
        )
interpro = pd.DataFrame(rows)
interpro = interpro.drop_duplicates(["pdb", "chain"])
interpro.to_csv("Interpro_FN3_homologs.csv", index=False)

print(len(interpro), "unique chains saved")
print(interpro["method"].value_counts())
# Check if I have any rows with no protein chains
print(interpro["chain"].isna().sum(), "rows with no chain")


# ---------- InterPro PROSITE Data to a CSV ----------
with open("Results/structure-matching-PS50853.json") as file:
    interPro_results = json.load(file)

rows = []

for record in interPro_results:
    meta = record["metadata"]

    for entry in record["entries"]:
        # PFam has e-values for each copy of the protein domain
        locations = entry.get("entry_protein_locations")
        if locations is None:
            locations = []

        n_domains_found = len(locations)

        rows.append(
            {
                "pdb": meta["accession"].upper(),
                "chain": entry.get("chain"),
                "description": meta["name"],
                "method": meta["experiment_type"],
                "resolution": meta["resolution"],
                "source": "interpro_IPR003961",
            }
        )
interpro = pd.DataFrame(rows)
interpro = interpro.drop_duplicates(["pdb", "chain"])
interpro.to_csv("Interpro_Prosite_FN3_homologs.csv", index=False)

print(len(interpro), "unique chains saved")
print(interpro["method"].value_counts())
# Check if I have any rows with no protein chains
print(interpro["chain"].isna().sum(), "rows with no chain")
