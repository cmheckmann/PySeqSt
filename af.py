"""Functions for AlphaFold"""
import sys
import requests

from seq import Seqs

# Define the EBI AlphaFold PSB API
EBI_URL = "https://alphafold.ebi.ac.uk/api/prediction/"


def download(url: str) -> list:
    """Download the AlphaFold PDBx/mmCIF file at provided url"""

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except (requests.HTTPError, requests.exceptions.Timeout):
        sys.exit("HTTP-error while saving EBI AlphaFold structure")

    return response.text.splitlines(keepends=True)

def from_uniprot(seqs: "Seqs"):
    """Find AlphaFold models in the EBI AlphaFold PSB for all sequences in a Seqs object
    which do not yet have a structure, using their UniProt Accession number.
    https://alphafold.ebi.ac.uk/api-docs"""

    for dscrptr, accessions in seqs.accessions.items():
        # Skip if there are already structures
        if seqs.get_structures(dscrptr):
            continue

        # Query the EBI AlphaFold PSB API with each accession number
        structures = ["EBI-AF"]
        for accession in accessions:
            try:
                response = requests.get(EBI_URL + accession, timeout=10)
                # Unfortunately the API raises 404 if it doesn't have a structure...
                if response.status_code == 404:
                    break
                response.raise_for_status()
            except (requests.HTTPError, requests.exceptions.Timeout):
                sys.exit("HTTP-error while querying EBI AlphaFold PSB API.")

            # Analyse the result
            try:
                structures.append(response.json()[0]["cifUrl"])
            except requests.exceptions.JSONDecodeError:
                sys.exit("Unexpected response format from EBI AlphaFold PSB API.")

        seqs.add_structures(dscrptr, structures, verify=False)
