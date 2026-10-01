# the for loop was written by AI
#
# # ---------- Libraries ----------

import requests
import re

# Lets me pause between requests
import time
import os
import pandas as pd

# ---------- InterPro Data to a CSV ----------
with open("FN3_homologs_CPDF.txt") as file:
    cpdf_results = file.read()
rows = []

# finditer walks through every match in the file.
# (\d+) captures the MMDB number, (\w{4}) the PDB code.
# re.finditer yields every match in the text, unlike re.search which stops at the first.
# \s* allows any spacing after the colon; \s+ covers the gap between the two fields.
# .group(1) and .group(2) are the first and second sets of parentheses, in order

for match in re.finditer(r"MMDB ID:\s*(\d+)\s+PDB ID:\s*(\w{4})", cpdf_results):
    rows.append(
        {
            "mmdb_id": match.group(1),
            "pdb": match.group(2).upper(),
        }
    )

mmdb = pd.DataFrame(rows)
mmdb = mmdb.drop_duplicates("pdb")
mmdb.to_csv("FN3_homologs_CPDF.csv", index=False)

print(len(mmdb), "entries saved")

