# This script is for looking at the RSCB PDB database to find homologous structures to 1FNA (tenth domain of fibronection 3)

# I find search up all the possible proteins that would match the 1FNA_A structure then put a condition on it that it should be a protein and not
# DNA or RNA molecule.

# Then I need to filter for NMR structures only and then I look at structures resolved by x-ray crystallography with a resolution < 2.0A

# After that, I need to filter based on similarity percentage. I made the cutoff 85% similarity because I believe everything over that will just be other monobodies and we
# are hoping to find dimers and oligomers

# https://rcsbsearchapi.readthedocs.io/en/latest/api.html <-- to build the structure node
# https://search.rcsb.org/#search-api <-- all SEARCH RSCB API documentation
# https://alignment.rcsb.org/#api-intro <-- all alignment API
# https://www.rcsb.org/docs/search-and-browse/advanced-search/overview-advanced-search <-- great resource
# https://data.rcsb.org/data-attributes.html
# https://data.rcsb.org/index.html#data-api


# ---------- Libraries ----------

import requests
import json

# Lets me pause between requests
import time
import os
import pandas as pd

# ---------- Part 1: search ----------

# The address of RCSB's search service
ALIGNMENT_URL = "https://alignment.rcsb.org/api/v1/structures"
SUBMIT_URL = "/submit"
RESULTS_URL = "/results"
MONITOR_ALIGNMENT = (
    "https://alignment.rcsb.org/api/v1/structures/results?uuid={job_ticket}"
)


SEQUENCE_IDENTITY_CUTTOFF = 0.70
TM_CUTOFF = 0.8
RMSD_CUTOFF = 3


# ---------- Part 1: annotation helpers ----------


# Entry -> (method, resolution)
def get_method_and_resolution(pdb_id):
    url = "https://data.rcsb.org/rest/v1/core/entry/" + pdb_id
    response = requests.get(url)
    if response.status_code != 200:
        return None, None
    info = response.json()["rcsb_entry_info"]
    method = info["experimental_method"]
    resolutions = info.get("resolution_combined")
    if resolutions is None:
        resolution = None
    else:
        resolution = resolutions[0]
    return method, resolution


# ---------- Part 2: ALIGNMENT CALCULATIONS AND HELPERS ----------


# def alignment_query(pdb_id, chain_id, alignment_method):
#     method = {"name": alignment_method}

#     if alignment_method == "fatcat-flexible":
#         method["parameters"] = {
#             "rmsd_cutoff": 3,
#             "afp_dist_cutoff": 1600,
#             "fragment_length": 8,
#             "max_num_twists": 5,
#         }

#     # Outside the if, so it runs for every method
#     return {
#         "options": {"return_sequence_data": False},
#         "context": {
#             "mode": "pairwise",
#             "method": method,
#             "structures": [
#                 {"entry_id": "1FNA", "selection": {"asym_id": "A"}},
#                 {"entry_id": pdb_id, "selection": {"asym_id": chain_id}},
#             ],
#         },
#     }


# def submit(query):
#     data = {"query": json.dumps(query)}
#     response = requests.post(url=(ALIGNMENT_URL + SUBMIT_URL), data=data)
#     response.raise_for_status()
#     return response.text


# def monitor_job_status_of_submit_query(ticket, max_attempts=50):
#     for attempt in range(max_attempts):
#         response = requests.get(
#             url=(ALIGNMENT_URL + RESULTS_URL), params={"uuid": ticket}
#         )

#         if response.status_code == 404:
#             time.sleep(5)
#             continue

#         response.raise_for_status()
#         data = response.json()
#         job_status = data["info"]["status"]

#         if job_status == "COMPLETE":
#             return data

#         if job_status == "ERROR":
#             print("Alignment failed:", data["info"].get("message", data["info"]))
#             return None

#         time.sleep(5)

#     print("Gave up after", max_attempts, "attempts")
#     return None


# def parse_result(result):
#     alignment = result["results"][0]
#     summary = alignment["summary"]

#     scores = {}
#     for s in summary["scores"]:
#         scores[s["type"]] = s["value"]

#     return {
#         "tm_score": scores.get("TM-score"),
#         "rmsd": scores.get("RMSD"),
#         "seq_identity": scores.get("sequence-identity"),
#         "similarity_score": scores.get("similarity-score"),
#         "n_aligned": summary.get("n_aln_residue_pairs"),
#         "n_blocks": len(alignment["structure_alignment"]),
#     }


# def run_alignment(pdb_id, chain_id, method_name):
#     try:
#         ticket = submit(alignment_query(pdb_id, chain_id, method_name))
#         result = monitor_job_status_of_submit_query(ticket)
#     except requests.exceptions.RequestException as e:
#         print(pdb_id, chain_id, method_name, "request failed:", e)
#         return {}

#     if result is None:
#         return {}

#     return parse_result(result)


# ---------- Part 3: build the table, or load it if already built ----------
annotated_csv = "Results/psiblast_final_hits_annotated.csv"
filtered_annotated_csv = "Results/psiblast_filtered_hits_annotated.csv"

if os.path.exists(annotated_csv):
    dataFrame = pd.read_csv(annotated_csv)
    print("Loaded", len(dataFrame), "annotated hits from pre-saved file")
else:
    dataFrame = pd.read_csv("Results/psiblast_final_hits.csv")
    print("Loaded", len(dataFrame), "protein chains from PSI-Blast results")

    # Since PSI-Blast does not report the methods I have to individually look each up so I can filter based on that

    cache = {}
    method_list = []
    resolution_list = []

    for i in range(len(dataFrame)):
        pdb_id = dataFrame["pdb"].iloc[i]

        if pdb_id not in cache:
            cache[pdb_id] = get_method_and_resolution(pdb_id)
            time.sleep(0.1)

        method, resolution = cache[pdb_id]
        method_list.append(method)
        resolution_list.append(resolution)

        # print only every 100 lines, not every line
        if i % 100 == 0:
            print(i, "/", len(dataFrame), flush=True)

    dataFrame["method"] = method_list
    dataFrame["resolution"] = resolution_list
    # dataFrame.to_csv(annotated_csv, index=False)

    print("Methods: _________________________________")
    print(dataFrame["method"].value_counts())


# ---------- Filtering for good methods ----------

method_nmr = dataFrame["method"].str.contains("NMR", case=False, na=False)
method_xray_resolution = dataFrame["resolution"]
good_xray = method_xray_resolution <= 2.55
final_dataset = dataFrame[method_nmr | good_xray]

print(
    len(dataFrame),
    "--->",
    len(final_dataset),
    "protiens left after quality method filtering",
)

final_dataset.to_csv(filtered_annotated_csv, index=False)

# ---------- Checking for the known proteins ----------

knowns = ["1TEN", "1QR4", "3TES", "4U3H"]

print("--- knowns in the full PSI-BLAST set ---")
print(
    dataFrame[dataFrame["pdb"].isin(knowns)][
        ["pdb", "chain", "e-value", "pct identity", "method", "resolution"]
    ]
)

print("--- knowns surviving the quality filter ---")
print(
    final_dataset[final_dataset["pdb"].isin(knowns)][
        ["pdb", "chain", "e-value", "method", "resolution"]
    ]
)
