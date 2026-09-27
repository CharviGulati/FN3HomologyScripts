# ---------- Merge Interpro and PSI-Blast results ----------


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
