# ---------- Merge Interpro and PSI-Blast results ----------
# the following merge was written with the help of Claude.ai

# my_FN3_merge_align.py


# ---------- Libraries ----------

import requests
import json

# Lets me pause between requests
import time
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import pandas as pd

# ---------- Logging ----------


# Sends everything printed to the terminal AND to the open log files.
# Several workers print at once, so each worker's text is held until its line
# is complete, then written whole under a lock - lines never get mixed together.
class Tee:
    def __init__(self):
        # The real terminal, saved before sys.stdout is replaced below
        self.terminal = sys.stdout
        self.files = []
        self.lock = threading.Lock()
        # Each worker thread gets its own unfinished-line buffer
        self.local = threading.local()

    def _emit(self, text):
        with self.lock:
            self.terminal.write(text)
            for f in self.files:
                f.write(text)

    def write(self, text):
        buffer = getattr(self.local, "buffer", "") + text
        if "\n" in buffer:
            # Write every complete line, keep the unfinished part
            complete, buffer = buffer.rsplit("\n", 1)
            self._emit(complete + "\n")
        self.local.buffer = buffer

    def flush(self):
        # Write out anything left that has no newline yet
        leftover = getattr(self.local, "buffer", "")
        if leftover:
            self._emit(leftover)
            self.local.buffer = ""
        with self.lock:
            self.terminal.flush()
            for f in self.files:
                f.flush()


tee = Tee()
sys.stdout = tee
sys.stderr = tee

psiblast = pd.read_csv("Results/psiblast_filtered_hits_annotated.csv")
interpro = pd.read_csv("Results/Interpro_FN3_homologs.csv")

# The address of RCSB's search service
ALIGNMENT_URL = "https://alignment.rcsb.org/api/v1/structures"
SUBMIT_URL = "/submit"
RESULTS_URL = "/results"
MONITOR_ALIGNMENT = (
    "https://alignment.rcsb.org/api/v1/structures/results?uuid={job_ticket}"
)


SEQUENCE_IDENTITY_CUTTOFF = 0.70
TM_CUTOFF_MIN = 0.5
TM_CUTOFF_MAX = 1.0
# One full run (log, failed entries, homologs file) per RMSD cutoff
RMSD_CUTOFFS = [2, 2.5, 3]
EVAL_MIN = 0.0001

# Failed chains timed out at 50 attempts, so give them longer on the re-run
RETRY_MAX_ATTEMPTS = 150

# How many chains are aligned at the same time.
# Raise it to go faster; lower it if many "too many requests" messages show up.
MAX_WORKERS = 4

# One log per RMSD cutoff. The shared steps (merge, alignment) go in all of them,
# so each log is a complete record of its run.
log_files = {}
for rmsd in RMSD_CUTOFFS:
    log_files[rmsd] = open(f"LogFiles/FN3_merge_align_rmsd{rmsd}.log", "w")
    log_files[rmsd].write(
        f"===== FN3_merge_align | RMSD cutoff {rmsd} | started {datetime.now():%Y-%m-%d %H:%M:%S} =====\n"
    )
tee.files = list(log_files.values())


# Same case on all, or the merge silently misses matches
psiblast["pdb"] = psiblast["pdb"].str.upper()
interpro["pdb"] = interpro["pdb"].str.upper()


# ---- PSI-BLAST: one row per entry, keeping all unique chain ----
# I think this should be changed back to what
# it originally was wher eI took the highest E-value chain
psiblast_entries = psiblast.sort_values("e-value").drop_duplicates(["pdb", "chain"])
psiblast_entries = psiblast_entries[
    ["pdb", "chain", "description", "e-value", "pct identity", "method", "resolution"]
]
# psiblast_entries = psiblast_entries.rename(columns={"chain": "psiblast_chain"})
psiblast_entries["in_psiblast"] = True
print(psiblast_entries)


# ---- InterPro: one row per entry ----
interpro_entries = interpro[
    ["pdb", "chain", "description", "method", "resolution", "source"]
]
interpro_entries = interpro_entries.drop_duplicates(["pdb", "chain"])
# interpro_entries = interpro_entries.rename(columns={"chain": "interpro_chain"})
interpro_entries["in_interpro"] = True
print(interpro_entries)


# ---- Merge ----
# how="outer" keeps entries that appear in only one table.
# suffixes label the columns that exist in both (description, method, resolution).
merged = psiblast_entries.merge(
    interpro_entries,
    on=["pdb", "chain"],
    how="outer",
    suffixes=("_psiblast", "_interpro"),
)

# Entries missing from a source have NaN in its flag column; make those False
for col in ["in_psiblast", "in_interpro"]:
    merged[col] = merged[col].fillna(False)

# How many sources found each entry
merged["n_sources"] = merged["in_psiblast"].astype(int) + merged["in_interpro"].astype(
    int
)

# One chain column: PSI-BLAST's where it exists, otherwise InterPro's


merged = merged.sort_values(["n_sources", "e-value"], ascending=[False, True])
# Doesn't depend on the RMSD cutoff, so one file for all of them
merged.to_csv("Results/merged_fn3_candidates.csv", index=False)

print(len(merged), "unique PDB entries\n")
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


# ---------- Filtering by E-value for alignment ----------

# Strong sequence hits, plus InterPro-only entries that have no E-value at all
strong_hits = merged["e-value"] <= EVAL_MIN
no_evalue = merged["e-value"].isna()

top_evals = merged[strong_hits | no_evalue].sort_values("e-value")

print("\nAligning", len(top_evals), "chains")
print("  from PSI-BLAST:", strong_hits.sum())
print("  from InterPro only (no E-value):", no_evalue.sum())
print("E-value range:", top_evals["e-value"].min(), "to", top_evals["e-value"].max())


# ---------- Top chain per PDB entry ----------

# Lowest E-val first. Rows with no E-value go last (na_position),
# and chain breaks ties so the choice is reproducible.
ranked = top_evals.sort_values(["e-value", "chain"], na_position="last")

# The single best chain per entry
to_align = ranked.drop_duplicates("pdb")

# All other chains not aligned
not_aligned = ranked.drop(to_align.index)

print("\n", len(ranked), "chains --->", len(to_align), "entries to align")
print(len(not_aligned), "chains set aside")

# Doesn't depend on the RMSD cutoff, so one file for all of them
not_aligned.to_csv("Results/FN3_not_aligned.csv", index=False)


# ---------- Submit jobs for alignment ----------


# If RCSB answers 429 ("too many requests"), wait and ask again
# instead of counting the chain as failed
def request_with_backoff(send, max_tries=10):
    for attempt in range(max_tries):
        response = send()
        if response.status_code != 429:
            return response

        wait = 10 * (attempt + 1)
        print("RCSB says too many requests - waiting", wait, "s", flush=True)
        time.sleep(wait)

    return response


def submit(query):
    data = {"query": json.dumps(query)}
    response = request_with_backoff(
        lambda: requests.post(url=(ALIGNMENT_URL + SUBMIT_URL), data=data)
    )
    response.raise_for_status()
    return response.text


# ---------- Author chain -> label chain ----------

# The alignment API's asym_id is the mmCIF label_asym_id assigned by the PDB.
# Our chain column holds AUTHOR chain IDs, which often differ.
chain_map_cache = {}


def get_label_chain_map(pdb_id):
    # One request set per entry, cached
    if pdb_id in chain_map_cache:
        return chain_map_cache[pdb_id]

    mapping = {}

    response = requests.get("https://data.rcsb.org/rest/v1/core/entry/" + pdb_id)
    if response.status_code != 200:
        chain_map_cache[pdb_id] = mapping
        return mapping

    entity_ids = response.json()["rcsb_entry_container_identifiers"][
        "polymer_entity_ids"
    ]

    for entity_id in entity_ids:
        url = (
            "https://data.rcsb.org/rest/v1/core/polymer_entity/"
            + pdb_id
            + "/"
            + entity_id
        )
        entity_response = requests.get(url)
        if entity_response.status_code != 200:
            continue

        ids = entity_response.json()["rcsb_polymer_entity_container_identifiers"]

        # Parallel lists: position i in one matches position i in the other
        for auth, label in zip(ids["auth_asym_ids"], ids["asym_ids"]):
            mapping[auth] = label

        time.sleep(0.1)

    chain_map_cache[pdb_id] = mapping
    return mapping


# ---------- Alignment Query ----------


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
                {
                    "entry_id": "1FNF",
                    "selection": {"asym_id": "A", "beg_seq_id": 275, "end_seq_id": 368},
                },
                {"entry_id": pdb_id, "selection": {"asym_id": chain_id}},
            ],
        },
    }


# ---------- Monitor Alignment Jobs ----------


def monitor_job_status_of_submit_query(ticket, max_attempts=50):
    failed_jobs = []
    url = ""
    for attempt in range(max_attempts):
        response = request_with_backoff(
            lambda: requests.get(
                url=(ALIGNMENT_URL + RESULTS_URL), params={"uuid": ticket}
            )
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

        if max_attempts == 50:
            failed_jobs.append(url)

        time.sleep(5)

    print("Gave up after", max_attempts, "attempts")
    print("All the files that failed:", failed_jobs)
    return None


# ---------- Parse Alignment Job ----------


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


def run_alignment(pdb_id, chain_id, method_name, max_attempts=50):
    try:
        ticket = submit(alignment_query(pdb_id, chain_id, method_name))
        result = monitor_job_status_of_submit_query(ticket, max_attempts)
    except requests.exceptions.RequestException as e:
        print(pdb_id, chain_id, method_name, "request failed:", e)
        return {}

    if result is None:
        return {}

    return parse_result(result)


label_chains = []

for i in range(len(to_align)):
    pdb_id = to_align["pdb"].iloc[i]
    auth_chain = to_align["chain"].iloc[i]
    label_chains.append(get_label_chain_map(pdb_id).get(auth_chain))

    if i % 50 == 0:
        print("mapping chains", i, "/", len(to_align), flush=True)

to_align = to_align.copy()
to_align["label_chain"] = label_chains

print(to_align["label_chain"].isna().sum(), "chains could not be mapped")
print(
    (to_align["chain"] != to_align["label_chain"]).sum(),
    "chains where auth differs from label",
)


# ---------- Part 5: run alignments on every chain ----------


# Counts finished chains across all workers, for the progress lines
progress_lock = threading.Lock()
progress = {"done": 0}


def align_one(pdb_id, chain_id, max_attempts, total):
    # Skip anything that could not be mapped
    if pd.isna(chain_id):
        print("skipping", pdb_id, "- no label chain", flush=True)
        return None

    # Rigid alignment: gives the TM-score and RMSD used for the cutoff
    tm = run_alignment(pdb_id, chain_id, "tm-align", max_attempts)

    # Flexible alignment: gives the number of blocks (twists)
    fc = run_alignment(pdb_id, chain_id, "fatcat-flexible", max_attempts)

    # One row per chain. If an alignment failed, run_alignment returned {},
    # so .get() gives None instead of crashing.
    # Keyed on label_chain: the "chain" column holds AUTHOR chain IDs
    row = {
        "pdb": pdb_id,
        "label_chain": chain_id,
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

    # Progress: how many chains are finished so far (in finishing order)
    with progress_lock:
        progress["done"] += 1
        done = progress["done"]

    print(
        done,
        "/",
        total,
        pdb_id,
        chain_id,
        "TM:",
        row["tm_score"],
        "FATCAT (Flexible):",
        row["fatcat_tm"],
        "blocks:",
        row["fatcat_blocks"],
        flush=True,
    )

    # Pause between chains to stay under RCSB's rate limit
    time.sleep(0.3)
    return row


def align_chains(entries, max_attempts=50):
    pdb_ids = list(entries["pdb"])
    chain_ids = list(entries["label_chain"])
    total = len(entries)
    progress["done"] = 0

    print("Aligning", total, "chains,", MAX_WORKERS, "at a time", flush=True)

    # MAX_WORKERS chains run at once. executor.map hands the results back in
    # the same order as pdb_ids, even though the chains finish out of order.
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        results = executor.map(
            align_one,
            pdb_ids,
            chain_ids,
            [max_attempts] * total,
            [total] * total,
        )
        # Drop the chains that were skipped (None)
        aln_rows = [row for row in results if row is not None]

    # List of dictionaries -> table, joined back onto the entry info
    aln = pd.DataFrame(aln_rows)
    return entries.merge(aln, on=["pdb", "label_chain"])


# The alignments don't depend on the RMSD cutoff (it is only used in the
# filter below), so they are run once and shared by every cutoff
aligned_csv = "Results/FN3_aligned.csv"

if os.path.exists(aligned_csv):
    full = pd.read_csv(aligned_csv)
    print("Loaded saved alignments:", len(full))

else:
    full = align_chains(to_align)
    full.to_csv(aligned_csv, index=False)


# ---------- Part 6: per RMSD cutoff: re-run failed chains, then filter ----------

for rmsd in RMSD_CUTOFFS:
    # From here on, only this cutoff's log gets the output
    tee.files = [log_files[rmsd]]
    print("\n==========", "RMSD cutoff", rmsd, "==========", flush=True)

    # ---- Failed chains: no TM-align or no FATCAT result ----
    failed_csv = f"Results/failed_entries_{rmsd}.csv"
    no_result = full["tm_score"].isna() | full["fatcat_tm"].isna()
    full.loc[no_result, ["pdb", "chain", "label_chain"]].to_csv(failed_csv, index=False)
    print(no_result.sum(), "failed chains saved to", failed_csv)

    # ---- Re-run the chains listed in the failed file ----
    failed = pd.read_csv(failed_csv)
    if len(failed) > 0:
        print("\nRe-aligning", len(failed), "failed chains (RMSD", rmsd, "run)")

        failed_keys = full.set_index(["pdb", "label_chain"]).index.isin(
            failed.set_index(["pdb", "label_chain"]).index
        )
        retry = full.loc[failed_keys, to_align.columns]
        retried = align_chains(retry, RETRY_MAX_ATTEMPTS)

        # Swap the new results in for the old failed rows
        full = pd.concat([full[~failed_keys], retried], ignore_index=True)
        # Save, so the next cutoff only re-runs what is still failing
        full.to_csv(aligned_csv, index=False)

        still_failed = retried["tm_score"].isna() | retried["fatcat_tm"].isna()
        print(
            len(retried) - still_failed.sum(),
            "recovered,",
            still_failed.sum(),
            "still failing",
        )
        if still_failed.sum() > 0:
            print(
                retried.loc[still_failed, ["pdb", "chain", "label_chain"]].to_string(
                    index=False
                )
            )

    # ---- Checks ----
    # How many alignments failed
    print(
        "\n", full["tm_score"].isna().sum(), "chains in TOTAL with no TM-align result"
    )
    print(full["fatcat_blocks"].isna().sum(), "chains in TOTAL with no FATCAT result")

    # Best matches first
    print("\n")
    print(
        full.sort_values("tm_score", ascending=False)[
            ["pdb", "chain", "tm_score", "tm_score_rmsd", "fatcat_tm", "fatcat_rmsd"]
        ].head(30)
    )

    # ---- Filter ----
    # The cutoff: similar enough (TM), closely superimposed (RMSD)
    homologs = full[
        (full["tm_score"] >= TM_CUTOFF_MIN)
        & (full["tm_score"] <= TM_CUTOFF_MAX)
        & (full["tm_score_rmsd"] <= rmsd)
        & (full["fatcat_tm"] >= TM_CUTOFF_MIN)
        & (full["fatcat_tm"] <= TM_CUTOFF_MAX)
        & (full["fatcat_rmsd"] <= rmsd)
    ]

    homologs_csv = f"Results/FN3_homologs_final_{rmsd}.csv"
    print("\nRMSD cutoff", rmsd, ":", len(full), "chains ->", len(homologs), "homologs")
    homologs.to_csv(homologs_csv, index=False)
    print("Saved", homologs_csv)
    print(f"===== finished {datetime.now():%Y-%m-%d %H:%M:%S} =====", flush=True)

for f in log_files.values():
    f.close()
tee.files = []
