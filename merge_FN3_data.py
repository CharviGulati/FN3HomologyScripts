# ---------- Merge Interpro and PSI-Blast results ----------
# the following merge was written with the help of Claude.ai

import pandas as pd

psiblast = pd.read_csv("Results/psiblast_final_hits_annotated.csv")
interpro = pd.read_csv("Results/Interpro_FN3_homologs.csv")

# Same case on all, or the merge silently misses matches
psiblast["pdb"] = psiblast["pdb"].str.upper()
interpro["pdb"] = interpro["pdb"].str.upper()


# ---- PSI-BLAST: one row per entry, keeping the best-scoring chain ----
psiblast_entries = psiblast.sort_values("e-value").drop_duplicates("pdb")
psiblast_entries = psiblast_entries[
    ["pdb", "chain", "description", "e-value", "pct identity", "method", "resolution"]
]
psiblast_entries = psiblast_entries.rename(columns={"chain": "psiblast_chain"})
psiblast_entries["in_psiblast"] = True


# ---- InterPro: one row per entry ----
interpro_entries = interpro[
    ["pdb", "chain", "description", "method", "resolution", "source"]
]
interpro_entries = interpro_entries.drop_duplicates("pdb")
interpro_entries = interpro_entries.rename(columns={"chain": "interpro_chain"})
interpro_entries["in_interpro"] = True


# ---- Merge ----
# how="outer" keeps entries that appear in only one table.
# suffixes label the columns that exist in both (description, method, resolution).
merged = psiblast_entries.merge(
    interpro_entries, on="pdb", how="outer", suffixes=("_psiblast", "_interpro")
)

# Entries missing from a source have NaN in its flag column; make those False
for col in ["in_psiblast", "in_interpro"]:
    merged[col] = merged[col].fillna(False)

# How many sources found each entry
merged["n_sources"] = merged["in_psiblast"].astype(int) + merged["in_interpro"].astype(
    int
)

# One chain column: PSI-BLAST's where it exists, otherwise InterPro's
merged["chain"] = merged["psiblast_chain"].fillna(merged["interpro_chain"])

merged = merged.sort_values(["n_sources", "e-value"], ascending=[False, True])
merged.to_csv("merged_fn3_candidates.csv", index=False)

print(len(merged), "unique PDB entries")
print(merged["n_sources"].value_counts())
print(merged[["in_psiblast", "in_interpro"]].sum())


# ---------- One method and resolution column ----------

# Both sources report these; take PSI-BLAST's, fall back to InterPro's
merged["method"] = merged["method_psiblast"].fillna(merged["method_interpro"])
merged["resolution"] = merged["resolution_psiblast"].fillna(
    merged["resolution_interpro"]
)
merged["description"] = merged["description_psiblast"].fillna(
    merged["description_interpro"]
)


# ---------- Filtering for good methods ----------

method_nmr = merged["method"].str.contains("NMR", case=False, na=False)
good_xray = merged["resolution"] <= 2.0
final_dataset = merged[method_nmr | good_xray]

print(
    len(merged),
    "--->",
    len(final_dataset),
    "proteins left after quality method filtering",
)


# ---------- Filtering top 100 E-values for alignment ----------

top_evals = final_dataset.sort_values("e-value").head(100)
print("Aligning", len(top_evals), "chains")
print("E-value range:", top_evals["e-value"].min(), "to", top_evals["e-value"].max())


# ---------- Checking for the known proteins ----------

knowns = ["1TEN", "1QR4", "3TES", "4U3H"]

print("--- knowns in the merged set ---")
print(
    merged[merged["pdb"].isin(knowns)][
        ["pdb", "chain", "e-value", "method", "resolution", "n_sources"]
    ]
)

print("--- knowns that survived the quality method filter ---")
print(
    final_dataset[final_dataset["pdb"].isin(knowns)][
        ["pdb", "chain", "e-value", "method", "resolution"]
    ]
)

print("--- knowns in the top 100 of E-values ---")
print(top_evals[top_evals["pdb"].isin(knowns)][["pdb", "chain", "e-value"]])


# ---------- Part 5: run alignments on every chain ----------
# align_csv = "rcsb_fn3_alignments.csv"

# # If alignments were already run, load them instead of redoing ~400 jobs
# if os.path.exists(align_csv):
#     aln = pd.read_csv(align_csv)
#     print("Loaded saved alignments:", len(aln))

# else:
#     # One dictionary per chain goes in here
#     aln_rows = []

#     # Go through the annotated table row by row
#     for i in range(len(top_evals)):
#         # This row's PDB ID and chain
#         pdb_id = top_evals["pdb"].iloc[i]
#         chain_id = top_evals["chain"].iloc[i]

#         # Rigid alignment: gives the TM-score and RMSD used for the cutoff
#         tm = run_alignment(pdb_id, chain_id, "tm-align")

#         # Flexible alignment: gives the number of blocks (twists)
#         fc = run_alignment(pdb_id, chain_id, "fatcat-flexible")

#         # One row per chain. If an alignment failed, run_alignment returned {},
#         # so .get() gives None instead of crashing.
#         row = {
#             "pdb": pdb_id,
#             "chain": chain_id,
#             # From TM-align
#             "tm_score": tm.get("tm_score"),
#             "tm_score_rmsd": tm.get("rmsd"),
#             "tm_score_n_aligned": tm.get("n_aligned"),
#             "tm_score_seq_identity": tm.get("seq_identity"),
#             # From FATCAT-flexible
#             "fatcat_tm": fc.get("tm_score"),
#             "fatcat_rmsd": fc.get("rmsd"),
#             "fatcat_seq_identity": fc.get("seq_identity"),
#             "fatcat_blocks": fc.get("n_blocks"),
#             "fatcat_score": fc.get("similarity_score"),
#         }

#         aln_rows.append(row)

#         # Progress: which chain, how far through
#         print(
#             i + 1,
#             "/",
#             len(top_evals),
#             pdb_id,
#             chain_id,
#             "TM:",
#             row["tm_score"],
#             "blocks:",
#             row["fatcat_blocks"],
#         )

#         # Pause between chains to stay under RCSB's rate limit
#         time.sleep(0.5)

#     # List of dictionaries -> table, then save it
#     aln = pd.DataFrame(aln_rows)
#     aln.to_csv(align_csv, index=False)


# # Join annotations and alignment scores, matching on pdb AND chain
# full = top_evals.merge(aln, on=["pdb", "chain"])
# full.to_csv("fn3_homologs_full.csv", index=False)


# # ---------- Part 6: checks and filter ----------

# # How many alignments failed
# print(full["tm_score"].isna().sum(), "chains with no TM-align result")
# print(full["fatcat_blocks"].isna().sum(), "chains with no FATCAT result")

# # Best matches first
# print(
#     full.sort_values("tm_score", ascending=False)[
#         ["pdb", "chain", "desc", "tm_score", "tm_score_rmsd", "fatcat_blocks"]
#     ].head(30)
# )

# # The cutoff: similar enough (TM), closely superimposed (RMSD)
# homologs = full[
#     (full["tm_score"] >= TM_CUTOFF) & (full["tm_score_rmsd"] <= RMSD_CUTOFF)
# ]

# print(len(full), "->", len(homologs), "homologs")
# homologs.to_csv("fn3_homologs_final.csv", index=False)

# print(
#     homologs[homologs["fatcat_blocks"] > 1][
#         ["pdb", "chain", "desc", "fatcat_tm", "fatcat_rmsd", "fatcat_blocks"]
#     ]
# )
