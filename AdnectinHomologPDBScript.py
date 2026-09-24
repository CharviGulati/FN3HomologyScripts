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


# ---------- Libraries ----------

import requests
import json

# Lets me pause between requests
import time
import os
import pandas as pd

# ---------- Part 1: search ----------

# The address of RCSB's search service
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ALIGNMENT_URL = "https://alignment.rcsb.org/api/v1/structures"
SUBMIT_URL = "/submit"
RESULTS_URL = "/results"
MONITOR_ALIGNMENT = (
    "https://alignment.rcsb.org/api/v1/structures/results?uuid={job_ticket}"
)

SEQUENCE_IDENTITY_CUTTOFF = 0.70
TM_CUTOFF = 0.8
RMSD_CUTOFF = 3

# Condition: chains shaped like 1FNA chain A
structure_node = {
    "type": "terminal",
    "service": "structure",
    "parameters": {
        "value": {"entry_id": "1FNA", "asym_id": "A"},
        "operator": "relaxed_shape_match",
        "target_search_space": "polymer_entity_instance",  # looks up the polymer I described and not an assembly
    },
}

# Condition: entry contains only protein
protein_only = {
    "type": "terminal",
    "service": "text",
    "parameters": {
        "attribute": "rcsb_entry_info.selected_polymer_entity_types",
        "operator": "exact_match",
        "value": "Protein (only)",  # this way we can discard DNA and RNA molecules
    },
}

# Condition: solution NMR
nmr = {
    "type": "terminal",
    "service": "text",
    "parameters": {
        "attribute": "exptl.method",
        "operator": "exact_match",
        "value": "SOLUTION NMR",
    },
}

# Condition: X-ray
xray = {
    "type": "terminal",
    "service": "text",
    "parameters": {
        "attribute": "exptl.method",
        "operator": "exact_match",
        "value": "X-RAY DIFFRACTION",
    },
}

# Condition: resolution below or equal to 2.0A
under_2A = {
    "type": "terminal",
    "service": "text",
    "parameters": {
        "attribute": "rcsb_entry_info.resolution_combined",
        "operator": "less_or_equal",
        "value": 2.0,
    },
}


# X-ray AND under 2.0A
good_xray = {"type": "group", "logical_operator": "and", "nodes": [xray, under_2A]}

# NMR OR good X-ray
method_filter = {"type": "group", "logical_operator": "or", "nodes": [nmr, good_xray]}


# Full query, asking for every hit
query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [structure_node, protein_only, method_filter],
    },
    "return_type": "polymer_instance",
    "request_options": {"return_all_hits": True},  # gets rid of the pagination headache
}


# ---------- Part 2: annotation helpers ----------


# Chain == entity number
def get_entity_id(pdb_id, chain_id):
    url = (
        "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/"
        + pdb_id
        + "/"
        + chain_id
    )
    response = requests.get(url)
    if response.status_code != 200:
        return None

    data = response.json()
    return data["rcsb_polymer_entity_instance_container_identifiers"]["entity_id"]


# Entity -> (description, 85% cluster ID)
# Entity -> description (the molecule's name)
def get_description(pdb_id, entity_id):
    url = (
        "https://data.rcsb.org/rest/v1/core/polymer_entity/" + pdb_id + "/" + entity_id
    )
    response = requests.get(url)
    if response.status_code != 200:
        return None
    return response.json()["rcsb_polymer_entity"]["pdbx_description"]

    # cluster = None  # going to build cluster here with sequence identify of 85%

    # for g in data.get("rcsb_polymer_entity_group_membership", []):
    #     if (
    #         g["aggregation_method"] == "sequence_identity"
    #         and g["similarity_cutoff"] == SEQUENCE_IDENTITY_CUTTOFF
    #     ):
    #         cluster = g["group_id"]  # sequence similarity is the cluster

    return description  # cluster


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


# ---------- Part 3: ALIGNMENT CALCULATIONS AND HELPERS ----------


def alignment_query(pdb_id, chain_id, alignment_method):
    method = {"name": alignment_method}

    if alignment_method == "fatcat-flexible":
        method["parameters"] = {
            "rmsd_cutoff": 3,
            "afp_dist_cutoff": 1600,
            "fragment_length": 8,
            "max_num_twists": 5,
        }

    # Outside the if, so it runs for every method
    return {
        "options": {"return_sequence_data": False},
        "context": {
            "mode": "pairwise",
            "method": method,
            "structures": [
                {"entry_id": "1FNA", "selection": {"asym_id": "A"}},
                {"entry_id": pdb_id, "selection": {"asym_id": chain_id}},
            ],
        },
    }


def submit(query):
    data = {"query": json.dumps(query)}
    response = requests.post(url=(ALIGNMENT_URL + SUBMIT_URL), data=data)
    response.raise_for_status()
    return response.text


def monitor_job_status_of_submit_query(ticket, max_attempts=50):
    for attempt in range(max_attempts):
        response = requests.get(
            url=(ALIGNMENT_URL + RESULTS_URL), params={"uuid": ticket}
        )

        if response.status_code == 404:
            time.sleep(5)
            continue

        response.raise_for_status()
        data = response.json()
        job_status = data["info"]["status"]

        if job_status == "COMPLETE":
            return data

        if job_status == "ERROR":
            print("Alignment failed:", data["info"].get("message", data["info"]))
            return None

        time.sleep(5)

    print("Gave up after", max_attempts, "attempts")
    return None


def parse_result(result):
    alignment = result["results"][0]
    summary = alignment["summary"]

    scores = {}
    for s in summary["scores"]:
        scores[s["type"]] = s["value"]

    return {
        "tm_score": scores.get("TM-score"),
        "rmsd": scores.get("RMSD"),
        "seq_identity": scores.get("sequence-identity"),
        "similarity_score": scores.get("similarity-score"),
        "n_aligned": summary.get("n_aln_residue_pairs"),
        "n_blocks": len(alignment["structure_alignment"]),
    }


def run_alignment(pdb_id, chain_id, method_name):
    try:
        ticket = submit(alignment_query(pdb_id, chain_id, method_name))
        result = monitor_job_status_of_submit_query(ticket)
    except requests.exceptions.RequestException as e:
        print(pdb_id, chain_id, method_name, "request failed:", e)
        return {}

    if result is None:
        return {}

    return parse_result(result)


print(run_alignment("8GZW", "F", "fatcat-flexible"))
print(run_alignment("8GZW", "F", "tm-align"))


# ---------- Part 4: build the table, or load it if already built ----------

csv_name = "rcsb_fn3_annotated.csv"

# If the CSV exists from a previous run, just load it
if os.path.exists(csv_name):
    df = pd.read_csv(csv_name)
    print("Loaded saved table:", len(df), "chains")

else:
    # Run the search
    response = requests.post(SEARCH_URL, json=query)

    # Collect chain IDs like "1QR4.B"
    ids = []
    if response.status_code != 204:
        for hit in response.json()["result_set"]:
            ids.append(hit["identifier"])
    print(len(ids), "chains found")

    rows = []

    for n, full_id in enumerate(ids):
        pdb_id, chain_id = full_id.split(".")

        entity_id = get_entity_id(pdb_id, chain_id)
        if entity_id is None:
            description = None
        else:
            description = get_description(pdb_id, entity_id)

        method, resolution = get_method_and_resolution(pdb_id)

        rows.append(
            {
                "pdb": pdb_id,
                "chain": chain_id,
                "desc": description,
                "method": method,
                "resolution": resolution,
            }
        )

        time.sleep(0.1)

    df = pd.DataFrame(rows)
    df.to_csv(csv_name, index=False)


# ---------- Part 5: run alignments on every chain ----------
align_csv = "rcsb_fn3_alignments.csv"

# If alignments were already run, load them instead of redoing ~400 jobs
if os.path.exists(align_csv):
    aln = pd.read_csv(align_csv)
    print("Loaded saved alignments:", len(aln))

else:
    # One dictionary per chain goes in here
    aln_rows = []

    # Go through the annotated table row by row
    for i in range(len(df)):
        # This row's PDB ID and chain
        pdb_id = df["pdb"].iloc[i]
        chain_id = df["chain"].iloc[i]

        # Rigid alignment: gives the TM-score and RMSD used for the cutoff
        tm = run_alignment(pdb_id, chain_id, "tm-align")

        # Flexible alignment: gives the number of blocks (twists)
        fc = run_alignment(pdb_id, chain_id, "fatcat-flexible")

        # One row per chain. If an alignment failed, run_alignment returned {},
        # so .get() gives None instead of crashing.
        row = {
            "pdb": pdb_id,
            "chain": chain_id,
            # From TM-align
            "tm_score": tm.get("tm_score"),
            "tm_score_rmsd": tm.get("rmsd"),
            "tm_score_n_aligned": tm.get("n_aligned"),
            "tm_score_seq_identity": tm.get("seq_identity"),
            # From FATCAT-flexible
            "fatcat_tm": fc.get("tm_score"),
            "fatcat_rmsd": fc.get("rmsd"),
            "fatcat_seq_identity": fc.get("seq_identity"),
            "fatcat_blocks": fc.get("n_blocks"),
            "fatcat_score": fc.get("similarity_score"),
        }

        aln_rows.append(row)

        # Progress: which chain, how far through
        print(
            i + 1,
            "/",
            len(df),
            pdb_id,
            chain_id,
            "TM:",
            row["tm_score"],
            "blocks:",
            row["fatcat_blocks"],
        )

        # Pause between chains to stay under RCSB's rate limit
        time.sleep(0.5)

    # List of dictionaries -> table, then save it
    aln = pd.DataFrame(aln_rows)
    aln.to_csv(align_csv, index=False)


# Join annotations and alignment scores, matching on pdb AND chain
full = df.merge(aln, on=["pdb", "chain"])
full.to_csv("rcsb_fn3_full.csv", index=False)


# ---------- Part 6: checks and filter ----------

# How many alignments failed
print(full["tm_score"].isna().sum(), "chains with no TM-align result")
print(full["fatcat_blocks"].isna().sum(), "chains with no FATCAT result")

# Best matches first
print(
    full.sort_values("tm_score", ascending=False)[
        ["pdb", "chain", "desc", "tm_score", "rmsd", "fatcat_blocks"]
    ].head(30)
)

# The cutoff: similar enough (TM), closely superimposed (RMSD)
homologs = full[(full["tm_score"] >= TM_CUTOFF) & (full["rmsd"] <= RMSD_CUTOFF)]

# Drop engineered FNfn10 variants by name
homologs = homologs[
    ~homologs["desc"].str.contains("monobody|adnectin", case=False, na=False)
]

print(len(full), "->", len(homologs), "homologs")
homologs.to_csv("rcsb_fn3_homologs.csv", index=False)

print(
    homologs[homologs["fatcat_blocks"] > 1][
        ["pdb", "chain", "desc", "tm_score", "fatcat_blocks"]
    ]
)
