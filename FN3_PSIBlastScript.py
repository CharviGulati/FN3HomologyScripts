# ---------- Libraries ----------

import requests
import json

# Lets me pause between requests
import time
import os
import pandas as pd

# ____________________________Global Variables____________________________________
BLAST_URL = "https://www.ebi.ac.uk/Tools/services/rest/psiblast"  # https://www.ebi.ac.uk/Tools/services/rest/psiblast/run
EMAIL = " "  # EBI requires one
MAX_ITERATIONS = 10
FILENAME = "psiblast_1fna.csv"


# ____________________________Housekeeping work____________________________________
# get FASTA sequence for reference structure which is 1FNA
fasta_query = requests.get("https://www.rcsb.org/fasta/entry/1FNA").text
print(fasta_query)

# The following are the parameters needed for the PSI-Blast to run
# https://www.ebi.ac.uk/Tools/services/rest/psiblast/parameters <-- this is where all the params are listed
# https://www.ebi.ac.uk/Tools/services/rest/psiblast/parameterdetails/{PARAM_NAME} <-- this is where all the params are explained

params = {
    "email": EMAIL,
    "sequence": fasta_query,
    "database": "pdb",
    # Maximum number of match alignments reported in the result output
    # 5, 10, 20, 50, 100, 150, 200, 250, 500, 750, 1000, 5000
    "alignments": 5000,
    # Maximum number of match score summaries reported in the result output
    # 5, 10, 20, 50, 100, 150, 200, 250, 500, 750, 1000, 5000
    "scores": 5000,
    # 1.0e-6, 1.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 0.1, 0.3 ... use this exact format
    "psithr": "0.1",
}


# ____________________________PSI-Blast Result Helpers____________________________________


def blast_result_types(job_id):
    response = requests.get(BLAST_URL + "/resulttypes/" + job_id)
    if response.status_code != 200:
        print("Run failed with response code: ", response.status_code)
        print("EBI error message: ", response.text)
    response.raise_for_status()
    return response.text


def get_results(job_id, result_type):
    response = requests.get(BLAST_URL + "/result/" + job_id + "/" + result_type)
    if response.status_code != 200:
        print("Run failed with response code: ", response.status_code)
        print("EBI error message: ", response.text)
    response.raise_for_status()
    return response.text


# ____________________________Running PSI-Blast____________________________________


# Submit a PSI-Blast job and get a JOB ID back
# if response code == 200 then the job was submitted successfully
def psiblast_run(params, selected_hits=None):

    files = None
    if selected_hits != None:
        files = {"selectedHits": ("selectedHits.txt", selected_hits)}

    response = requests.post(BLAST_URL + "/run", data=params, files=files)

    if response.status_code != 200:
        print("Reason for failure: ", response.text)

    response.raise_for_status()  # if the error code is a 400 the execution is stopped by this
    return response.text


job_id = psiblast_run(params)
print("Job with ID: ", job_id, " queued")

# ________________________________Waiting for PSI-Blast Results________________________________


# Use Job ID to watch for results
# Returns True if job finished
def psiblast_check_job_completed(job_id, max_attemps=100):
    for attempt in range(max_attemps):
        job_status = requests.get(BLAST_URL + "/status/" + job_id).text

        if job_status == "FINISHED":
            return True

        if job_status in ("ERROR", "FAILURE", "NOT_FOUND"):
            print("Job with ID: ", job_id, "ended with status ", job_status)
            print("Try again")
            return False

        # if none of the statements are true then the job is still running and
        # we just have to wait
        time.sleep(20)

    # if we still don't have results now then that means we maxed out our attempts
    print("Maxed out all attempts for Job ID = ", job_id)
    return False


# ________________________________Save PSI-Blast Results____________________________________

# job_ids = [job_id]

# # Check if the job is completed, and if it is we need to look at the results
# if psiblast_check_job_completed(job_id):
#     job_record = []
#     # need enumerate here because we need an index to label each json
#     # so we can keep track of every iteration and follow it logically
#     for phi_iteration, j_id in enumerate(job_ids, start=1):
#         plain_text_results = get_results(j_id, "json")

#         with open(f"psiblast_iteration_{phi_iteration}.json", "w") as f:
#             f.write(plain_text_results)

#         filename = f"psiblast_iteration_{phi_iteration}.json"
#         job_record.append(
#             {
#                 "PSI Iteration: ": phi_iteration,
#                 "job_id: ": job_id,
#                 "file: ": filename,
#             }
#         )

#         print("Done and saved iteration: ", phi_iteration)

#     pd.DataFrame(job_record).to_csv("psiblast_1fna.csv", index=False)
#     data_layout = json.loads(plain_text_results)
#     print(data_layout.keys())


if os.path.exists(FILENAME):
    job_record = pd.read_csv(FILENAME).to_dict("records")
    print("Loaded", len(job_record), " saved jobs from file")
else:
    job_record = []
    previous_job_id = None
    selected_entries_for_next_iteration = None
    selected_entries_for_previous_iteration = None

    for iteration in range(0, MAX_ITERATIONS):
        this_params = dict(params)

        if previous_job_id is None:
            job_id = psiblast_run(this_params)
        else:
            this_params["previousjobid"] = previous_job_id
            job_id = psiblast_run(
                this_params, selected_hits=selected_entries_for_next_iteration
            )

        print(
            "Iteration: ",
            iteration,
            "submitted to EBI with Job ID: ",
            job_id,
            flush=True,
        )

        # now I need to wait for the previous round to finish so the next round can proceed
        if not psiblast_check_job_completed(job_id):
            print("iteration number: ", iteration, "failed. Process stopped")
            break

        filename = f"psiblast_iteration_{iteration}.json"
        result_text = get_results(job_id, "json")
        with open(filename, "w") as file:
            file.write(get_results(job_id, "json"))

        job_record.append(
            {"iteration": iteration, "job_id": job_id, "filename": filename}
        )

        results = json.loads(result_text)
        new_hits = []

        for single_hit in results["hits"]:
            if single_hit.get("hit_iter") == "new":
                new_hits.append(single_hit)

        print(
            " ",
            len(results["hits"]),
            "hits total. ",
            len(new_hits),
            "are the new hits found in this iteration",
            flush=True,
        )

        selected_entries_for_next_iteration = get_results(job_id, "preselected_ids")

        print(
            " ",
            len(selected_entries_for_next_iteration.split()),
            "hits selected for next round",
            flush=True,
        )

        if iteration > 0 and len(new_hits) == 0:
            print(
                "No new hits found in this iteration. Process stopping. Iteration = ",
                iteration,
            )
            break

        previous_job_id = job_id
        time.sleep(10)

    pd.DataFrame(job_record).to_csv(FILENAME, index=False)


# # The top-level keys: tells you where the list of hits lives
# print(data_layout.keys())
# # How many hits this iteration returned
# print(len(data_layout["hits"]), "hits")
# # The first hit, indented so the nesting is readable
# print(json.dumps(data_layout["hits"][0], indent=2))


# ________________________________Parse PSI-Blast Results____________________________________
# use this if you need to save every single hit from every iteration


def parse_psiBlast(data, iteration):
    rows = []

    query_length = data["query_len"]

    for single_hit in data["hits"]:
        # this will give the top return so instead of keeping track of every single region that matched the query
        # I am only keeping the top hit since it matches the most.
        best = single_hit["hit_hsps"][0]

        # 1FNA_A == 1FNA and A
        pdb_id, chain_id = single_hit["hit_acc"].split("_", 1)

        covered_residues = best["hsp_query_to"] - best["hsp_query_from"] + 1

        rows.append(
            {
                "iteration": iteration,
                "pdb": pdb_id,
                "chain": chain_id,
                "description": single_hit["hit_desc"],
                "organism": single_hit.get("hit_os"),
                # will say "new if the hit appeared in this iteration"
                "hit iteration appearance": single_hit.get("hit_iter"),
                "e-value": best["hsp_expect"],
                # normalized alignment score
                "bit score": best["hsp_bit_score"],
                # percentage of aligned portions with the AA
                "pct identity": best["hsp_identity"],
                # percent of aligned positions that are identical in both or chemically similar
                "pct positive": best["hsp_positive"],
                # alignment legth includes gaps
                "align length": best["hsp_align_len"],
                "query cover": (covered_residues / query_length) * 100,
                # how many regions in total this hit had. a number > 1 means the protein has repeated FN3 domains
                "number of hsps": len(single_hit["hit_hsps"]),
            }
        )
    return rows


# ________________________________Save PSI-Blast Results____________________________________

all_rows = []

for job in job_record:
    with open(job["filename"]) as file:
        data = json.load(file)

    all_rows.extend(parse_psiBlast(data, job["iteration"]))

hits = pd.DataFrame(all_rows)
# hits.to_csv("psiblast_all_iterations.csv", index=False)

final_round = hits["iteration"].max()
final_hits = hits[hits["iteration"] == final_round]

combining_hits = hits.sort_values("e-value").drop_duplicates(["pdb", "chain"])

print("final round:", len(final_hits), "chains")
print("union of all rounds:", len(combining_hits), "chains")

# Use the union: it keeps hits that scored well early but dropped out later
combining_hits.to_csv("psiblast_final_hits.csv", index=False)


# ________________________________Checks____________________________________

# Hits per round
print(hits.groupby("iteration").size())

# New hits per round: a sudden jump means the profile drifted
print(hits[hits["hit iteration appearance"] == "new"].groupby("iteration").size())

# Did the knowns come through, and in which round?
print(
    hits[hits["pdb"].isin(["1TEN", "1QR4", "3TES", "4U3H"])][
        ["iteration", "pdb", "chain", "pct identity", "e-value"]
    ]
)

print(hits[hits["pdb"].isin(["2KIU"])])
