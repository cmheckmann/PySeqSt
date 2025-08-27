"""Functions for file processing"""
import json
import os
import re
import sys

import af
import pd
from seq import Seqs

# Default output path
PATH = "./out"

PDB_URL = "https://files.rcsb.org/download/"

def new_output_folder():
    """Create a new folder for the output"""

    global PATH
    while True:
        try:
            os.mkdir(PATH)
            break
        except FileExistsError:
            p = PATH.rsplit("_", 1)
            if len(p) == 1 or not p[1].isdecimal():
                p.append("2")
            else:
                p[1] = str(int(p[1]) + 1)
            PATH = "_".join(p)


def open_blast(file: str) -> dict:
    """Open a user provide BLAST.json file and return the json object."""

    # Attempt to open file
    try:
        with open(file) as f:
            return json.load(f)
    except OSError:
        sys.exit(f"'{file}' could not be opened. ")
    except json.JSONDecodeError:
        sys.exit(f"Error reading '{file}', invalid JSON.")


def process_fasta(file: str) -> "Seqs":
    """Process user provided FASTA file of protein sequence(s).
    Return a new Seqs object."""

    # Attempt to open file
    try:
        with open(file) as f:
            sequences = f.readlines()
    except OSError:
        sys.exit(f"'{file}' could not be opened. ")

    # Grab all the descriptors and stitch together sequences spread over multiple lines
    seqs = Seqs()
    seq = ""
    dscrptr = ""
    invalid = 0
    for n, row in enumerate(sequences):
        # Stitch together seqeuence over multiple lines
        if not row.strip().startswith(">"):
            seq += row.strip()
            # Onto next row unless we reached the last row
            if n != len(sequences) - 1:
                continue

        # End of sequence; add to sequence object
        if seq != "":
            match seqs.add_seq(dscrptr, seq):
                case 0:
                    pass
                case 1:
                    print(f"Sequence\n'{dscrptr}':\n'{seq}'\nis not a valid protein sequence.\n"
                        "Verify it is contains only single-letter canonical (20) amino acids.\n")
                    invalid += 1
                case 2:
                    print(f"Sequence\n'{dscrptr}':\n'{seq}'\nis a duplicate sequence\n")
                case _:
                    sys.exit(f"Unkown error trying to add '{dscrptr}': '{seq}'")

            # Re-initialize seq
            seq = ""

        # Update descriptor
        if row.strip().startswith(">"):
            dscrptr = row.strip().removeprefix(">").strip()

    # Print number of invalid sequences and prompt user to continue
    if invalid > 0:
        print(f"{invalid} sequences were excluded from analysis due to invalid protein sequence.")
        _cont()

    return seqs


def save_blast(blast: dict):
    """Save a BLAST dictionary to BLAST.json"""

    out = f"{PATH}/BLAST.json"
    try:
        with open(out, "w") as f:
            json.dump(blast, f, indent=2)
    except OSError:
        print(f"Warning: BLAST output could not be saved to {out}.")
    else:
        print(f"BLAST output saved to {out}.")


def set_output_folder(path: str):
    """Set the PATH for any output files"""

    global PATH
    PATH = path
    # ensure relative path
    if not PATH.startswith("./"):
        PATH = "./" + PATH

    # Make the directory if it doesn't exist
    try:
        os.makedirs(PATH)
    except FileExistsError:
        print(f"WARNING: any existing files in {PATH} may be overwritten.")
        _cont()


def save_structures(seqs: "Seqs"):
    """Save all structures in Seqs object"""

    # Create folders to save the structures
    _make_outputdirs(seqs)
    # Save all structures
    for dscrptr, structures in seqs.structures.items():
        # PDB numbers
        if structures[0] == "pdb":
            for pdb in structures[1:]:
                print(f"Saving: {dscrptr}: {pdb}.cif")
                with open(f"{PATH}/{dscrptr}/{pdb}.cif", "w") as file:
                    file.writelines(pd.download(pdb))
        # EBI AlphaFold structures
        elif structures[0] == "EBI-AF":
            for url in structures[1:]:
                name = url.rsplit("/", 1)[1]
                print(f"Saving: {dscrptr}: {name}")
                with open(f"{PATH}/{dscrptr}/{name}", "w") as file:
                    file.writelines(af.download(url))


def _make_outputdirs(seqs: "Seqs"):
    """Create subfolders to store structures"""

    for dscrptr in seqs.dscrptrs:
        try:
            os.mkdir(PATH + "/" + dscrptr)
        except FileExistsError:
            pass

def _cont():
    """Ask the user if they want to continue"""

    if re.match(r"n",
            input("Press 'n'+enter to abort or 'any key'+enter to continue: "),
            flags=re.IGNORECASE):
        sys.exit("User request: exit programm")
