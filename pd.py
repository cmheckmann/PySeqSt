"""Functions for interacting with the rscb pdb"""
import json
import sys
import requests

from seq import Seqs

# Define URL for Search API
SURL = "https://search.rcsb.org/rcsbsearch/v2/query"
# Define URL for download API
PDB_URL = "https://files.rcsb.org/download/"


def download(pdb: str) -> list:
    """Download the PDBx/mmCIF file corresponsing to a pdb entry"""

    try:
        response = requests.get(PDB_URL + pdb + ".cif", timeout=10)
        response.raise_for_status()
    except (requests.HTTPError, requests.exceptions.Timeout):
        sys.exit("HTTP-error while saving pdb structure")

    return response.text.splitlines(keepends=True)


def uniprot_to_pdb(seqs: "Seqs"):
    """Searches the pdb for any UniProt accession numbers for each descriptor in Seqs object.
    Updates the list of structures for a descriptor with missing pdb numbers.
    Overwrites any non-pdb structures stored in Seqs object if pdb structures were found.
    https://search.rcsb.org/index.html#search-example-8"""

    # Template for searching by uniprot
    query: dict = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "operator": "exact_match",
                    "value": "",
                    "attribute": "rcsb_polymer_entity_container_identifiers."
                                "reference_sequence_identifiers.database_accession"
                }
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "operator": "exact_match",
                    "value": "UniProt",
                    "attribute": "rcsb_polymer_entity_container_identifiers."
                                "reference_sequence_identifiers.database_name"
                }
            }
            ]
        },
        "return_type": "entry"
    }


    # Loop over all descriptors
    for dscrptr, accessions in seqs.accessions.items():
        # Get already known pdb entries
        pdbs = seqs.get_structures(dscrptr)
        # No previously known (pdb-type) structures
        if not pdbs or pdbs[0] != "pdb":
            pdbs = ["pdb"]

        # Loop over the list of accession numbers for the current decsriptor
        for accession in accessions:

            # Query the rscb pdb
            query["query"]["nodes"][0]["parameters"]["value"] = accession
            try:
                response = requests.get(SURL,
                                        params={
                                            "json": json.dumps(query)
                                        },
                                        timeout=10)
                response.raise_for_status()
            except (requests.HTTPError, requests.exceptions.Timeout):
                sys.exit("HTTP-error while querying RSCB PDB Search API.")

            # No PDB entries found
            if response.status_code == 204:
                break

            # Analyse the result
            try:
                results = response.json()
            except requests.exceptions.JSONDecodeError:
                sys.exit("Unexpected response format RSCB PDB Search API.")

            # PDB entries found
            try:
                for result in results["result_set"]:
                    # Add any unique pdbs
                    if result["identifier"].split("_", 1)[0] not in pdbs:
                        pdbs.append(result["identifier"].split("_", 1)[0])
            except KeyError:
                sys.exit("Unexpected JSON structure in RSCB PDB Search response.")

        seqs.add_structures(dscrptr, pdbs, verify=False)
