# ---------- Libraries ----------

import requests
import json

# Lets me pause between requests
import time
import os
import pandas as pd

# ---------- InterPro Data to a CSV ----------
with open("Results/structure-matching-PF00041.json") as file:
    interPro_results = json.load(file)

rows = []

for record in interPro_results:
    meta = record["metadata"]

    for entry in record["entries"]:
        # PFam has e-values for each copy of the protein domain
        locations = entry.get("entry_protein_locations")
        if locations is None:
            locations = []

        scores = []
        for location in locations:
            if location.get("score") is not None:
                scores.append(location["score"])

        if scores:
            best_e_value = min(scores)
        else:
            best_e_value = None

        n_domains_found = len(locations)

        rows.append(
            {
                "pdb": meta["accession"].upper(),
                "chain": entry.get("chain"),
                "description": meta["name"],
                "method": meta["experiment_type"],
                "resolution": meta["resolution"],
                "pfam_best_evalue": best_e_value,
                "n_fn3_domains": n_domains_found,
                "uniprot": entry.get("protein"),
                "source": "interpro_PF00041",
            }
        )
interpro = pd.DataFrame(rows)
interpro = interpro.drop_duplicates(["pdb", "chain"])
interpro.to_csv("Interpro_PFam_FN3_homologs.csv", index=False)

print(len(interpro), "unique chains saved")
print(interpro["method"].value_counts())
print(interpro["n_fn3_domains"].value_counts().head())
# Check if I have any rows with no protein chains
print(interpro["chain"].isna().sum(), "rows with no chain")
