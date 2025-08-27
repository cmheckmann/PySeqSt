"""Functions for running and processing BLAST"""
import re
import sys
import time
import requests

from fp import save_blast
from seq import Seqs

# Define minimum identity, cover and maximum number of gaps to be considered a match
IDENT = 1
COV = 0.9
GAPS = 1

# Define URL for BLAST API
B_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

# Define URl for UniProt API
UP_URL = "https://rest.uniprot.org/idmapping/"
# Define databases to use for UniProt conversion
DBS = ["EMBL-GenBank-DDBJ_CDS", "RefSeq_Protein"]

#Define long messages
E1 = "Provided BLAST.json file has unexpected structure. " \
     "Verify provided file was produced by this program, and hasn't been tampered with."
E2 = "Sequence\n'{dscrptr}':\n'{seq}'\ncould not be added from provided BLAST result.\n" \
     "Verify provided BLAST.json was produced by this program, and hasn't been tampered with.\n"
M1 = "BLAST (ID: {}) is running, initial estimated time to completion is {} seconds..."


def extractseq(blast: dict) -> "Seqs":
    """Process user provided BLAST file to extract all query sequences.
    Return a new Seqs object."""

    seqs = Seqs()
    for report in blast["BlastOutput2"]:
        try:
            dscrptr = report["report"]["results"]["search"]["query_title"]
            seq = report["report"]["results"]["search"]["query_seq"]
        except KeyError:
            sys.exit(E1)
        if seqs.add_seq(dscrptr, seq) != 0:
            sys.exit(E2.format(dscrptr=dscrptr, seq=seq))

    return seqs


def process(blast: dict, seqs: "Seqs") -> int:
    """Go through a BLAST json dictionary and find pdb entries for each query.
    Ignore descriptors not in Seqs object.
    Add a list of those pdb entries to the Seqs object.
    The first entry of the list indicates the type of structure.
    Return the number of sequences that do not match the Seqs object."""

    def _checkhit() -> bool:
        """Determines if a hit sequence matches the query sequence"""

        hsp = hit["hsps"][0]
        # This allows for gaps and tags
        if (hsp["identity"] + hsp["gaps"] >= hsp["align_len"] * IDENT
            and hsp["align_len"] > q_len * COV):

            # Check any gaps are consequtive (likely caused by a tag)
            if hsp["gaps"] == 0 or (len(re.findall(r"-+", hsp["hseq"])) <= GAPS
                                    and len(re.findall(r"-+", hsp["qseq"])) <= GAPS):
                return True
        return False

    def _acc_pdb_from_hit():
        """Appends pdb number and accession number to the respective lists for that query"""

        flag = True
        for description in hit["description"]:

            # PDB number
            if description["id"].startswith("pdb|"):
                # Add only pdbs, removing chain ids
                if description["accession"].split("_", 1)[0] not in pdbs:
                    pdbs.append(description["accession"].split("_", 1)[0])

            # Accession number
            elif flag:
                # We only need one accession per hit, and no duplicates
                if description["accession"] not in accessions:
                    accessions.append(description["accession"])
                    flag = False

    # Loop through the queries
    c = 0
    for report in blast["BlastOutput2"]:
        try:
            d = report["report"]["results"]["search"]["query_title"]
            s = report["report"]["results"]["search"]["query_seq"]
            # Skip sequences not in Seqs object
            if s not in seqs.seqs:
                print(f"BLAST-query '{d}': query-sequence not in user-provided sequences.")
                c += 1
                continue
            q_len = report["report"]["results"]["search"]["query_len"]

            pdbs = ["pdb"]
            accessions: list = []
            for hit in report["report"]["results"]["search"]["hits"]:

                # Skip to next hit if no match
                if not _checkhit():
                    continue

                # Grab any pdb numbers & accession numbers
                _acc_pdb_from_hit()

        except KeyError:
            sys.exit(E1)

        # Convert accession numbers to UniProt
        if len(accessions) > 0:
            print(f"Converting Accession Numbers for BLAST-query: '{d}' to UniProt...")
            accessions = _convert_to_uniprot(accessions)
        # Try to add list of accession numbers to the seq object
        if not seqs.add_accessions(d, accessions, s):
            print(f"BLAST-query '{d}' does not match any user-provided descriptor+sequence pair.")
            c += 1
            continue
        # Add list of pdbs (skip if add_accessions has failed)
        seqs.add_structures(d, pdbs, s)

    # Number of sequences not matching Seqs object
    return c


def run(seqs: "Seqs") -> dict:
    """Format and submit request to ncib BLAST. Check for usage limits at:
    https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo

    Dumps BLAST results to a file, and returns the json dictionary.

    Adapted from https://blast.ncbi.nlm.nih.gov/docs/web_blast.pl and
    https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html#urlapi"""

    query = seqs.get_fasta()

    if not query:
        sys.exit("No sequences to run BLAST on.")

    # Submit query
    response = _b_request({"CMD": "Put",
                         "PROGRAM": "blastp",
                         "DATABASE": "nr",
                         "QUERY": query,})

    # Extract the Request ID
    if _rid := re.search(r"RID = ([A-Z0-9]+)\n", response.text):
        rid = _rid.group(1)
    else:
        sys.exit("No response ID in response")

    # Extract the estimated time to completion
    if _rtoe := re.search(r"RTOE = ([0-9]+)\n", response.text):
        rtoe = int(_rtoe.group(1))
    else:
        sys.exit("No estimated time to completion in response")

    # Wait for estimated time to completion + delay
    delay = 10
    print(M1.format(rid, rtoe + delay))
    time.sleep(rtoe)

    # Query for BLAST response
    while True:
        time.sleep(delay)
        response = _b_request({"CMD": "Get",
                             "RID": rid,})

        match _outcome(response.text):
            # Not yet finished
            case 1:
                delay += 10
                print(f"Job not yet finished, trying again in {delay} seconds...")
            # Failed
            case 2:
                sys.exit(f"BLAST (ID: '{rid}') failed")
            # RID not found
            case 3:
                sys.exit(f"BLAST (ID: '{rid}') expired")
            # It is done
            case 4:
                print(f"BLAST (ID: '{rid}') finished.")
                break
            # Unknown error
            case -1:
                sys.exit("Unknown error retrieving BLAST output")

    # Process BLAST response
    response = _b_request({"CMD": "Get",
                         "RID": rid,
                         "FORMAT_TYPE": "JSON2_S",})
    try:
        blast = response.json()
    except requests.exceptions.JSONDecodeError:
        sys.exit("Server BLAST output could not be decoded. Invalid JSON.")

    # Include query sequence in each output
    try:
        for report in blast["BlastOutput2"]:
            dscrptr = report["report"]["results"]["search"]["query_title"]
            report["report"]["results"]["search"]["query_seq"] = seqs.get_seq(dscrptr)
    except KeyError:
        sys.exit("Server BLAST JSON output has unexpected structure.")

    # Print to file
    save_blast(blast)
    return blast


def _convert_to_uniprot(a: list) -> list:
    """Convert accession numbers to UniProt
    Adapted from https://www.uniprot.org/help/id_mapping_prog
    and https://www.uniprot.org/api-documentation/idmapping"""

    # Prepare the list to be returned
    accessions = []
    # Submit one job for each databank
    job_ids = []
    for db in DBS:
        job_ids.append(_up_post({"from": db, "to": "UniProtKB", "ids": ",".join(a)}))

    # Check job status for each job
    delay = 15
    for job_id in job_ids:
        while True:
            response = _up_get(f"status/{job_id}")
            if response.get("jobStatus") in ("NEW", "RUNNING"):
                print("Running...")
                time.sleep(delay)
            elif response.get("jobStatus") is None:
                break
            else:
                print(f"UniProt job status: {response.get("jobStatus")}")
                sys.exit("Unexpected error from UniProt.")

    # Get results
    for job_id in job_ids:
        result = _up_get(f"uniprotkb/results/stream/{job_id}")["results"]

        # Grab unique accession numbers from the results
        for entry in result:
            try:
                if entry["to"]["primaryAccession"] not in accessions:
                    accessions.append(entry["to"]["primaryAccession"])
            except KeyError:
                pass

    return accessions


def _b_request(payload :dict) -> "requests.Response":
    """Send a request to the BLAST API"""

    try:
        response = requests.get(B_URL, params=payload, timeout=10)
        response.raise_for_status()
    except (requests.HTTPError, requests.exceptions.Timeout):
        sys.exit("HTTP-error while querying BLAST server!")
    return response


def _outcome(s: str) -> int:
    """Process the BLAST status."""
    if re.search(r"\s+Status=WAITING\n", s):
        return 1
    if re.search(r"\s+Status=FAILED\n", s):
        return 2
    if re.search(r"\s+Status=UNKNOWN\n", s):
        return 3
    if re.search(r"\s+Status=READY\n", s):
        return 4
    return -1


def _up_post(payload :dict) -> str:
    """Send a POST request to the UniProt ID Mapping API
    https://www.uniprot.org/help/id_mapping"""

    try:
        response = requests.post(UP_URL + "run", data=payload, timeout=10)
        response.raise_for_status()
    except (requests.HTTPError, requests.exceptions.Timeout):
        sys.exit("HTTP-error while querying UniProt server!")
    try:
        return response.json()["jobId"]
    except requests.exceptions.JSONDecodeError:
        sys.exit("Unexpected response format from UniProt.")


def _up_get(path: str) -> dict:
    """Send a GET request to the UniProt ID Mapping API
    https://www.uniprot.org/help/id_mapping"""

    try:
        response = requests.get(UP_URL + path, timeout=10)
        response.raise_for_status()
    except (requests.HTTPError, requests.exceptions.Timeout):
        sys.exit("HTTP-error while querying UniProt server!")
    try:
        return response.json()
    except requests.exceptions.JSONDecodeError:
        sys.exit("Unexpected response format from UniProt.")
