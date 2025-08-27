# Python Sequence to Structure (PySeqSt)
#### Video Demo:  https://youtu.be/_S9PuQmHux8
#### Description: Takes amino acid sequences and finds corresponding structures in the RCSB Protein Data Bank or EBI AlphaFold Protein Structure Database
## Requirements
requests https://requests.readthedocs.io/en/latest/

## How to run the program
|Args          |  Function                                       |
|--------------|-------------------------------------------------|
|-h  --help    |  Show help message                              |
|-f  --fasta   |  Provide input as a `.fasta` file               |
|-s, --seq     |  OR provide a single amino acid sequence        |
|-b, --blast   |  Pick up from a previous BLAST result           |
| --only_blast |  OR only perform the BLAST but do not process it|
|-o, --out     |  Specify output folder                          |

Without command-line arguments (or with only `--only_blast`), the user will be prompted for a SINGLE sequence (without descriptors, the program will set the descriptor to `seq_1`).\
Alternatively, a single sequence can be provided using `-s`.\
To run multiple sequences and/or allow for user defined descriptors, sequence(s) have to be provided as a `.fasta` file.\
Duplicate sequences will be rejected; duplicate descriptors will be appended with `_#` where `#` is incremented as needed to be unique. \
In all cases, invalid sequences (i.e. sequences that contain any symbol not matching a single letter, canonical amino acid or * as stop codon) will be rejected. The user will be given an option to exit the program to correct these sequences.

The default output folder is `.\out` or `.\out_#` where `#` is incremented as needed to be unique.\
When picking up from a previous BLAST result, the default output folder is whichever folder contains the BLAST result; files will be over-written, but the user will be prompted with an option to exit the program.\
A custom output folder can be provided with `-o`, if it already exists, files will be over-written, but the user will be prompted with an option to exit the program.

By default, after processing the user input, the sequences will be submitted to the NCBI BLAST API, which is a public resource and should be used responsibly.\
https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo for more info.

The BLAST result will then be analysed; to avoid rerunning identical BLAST queries, this analysis can also be started from a previous BLAST result.\
When starting from a prior BLAST, the user may still provide sequences to only process a sub-set of the BLAST result. Care should be taken to ensure descriptors AND sequences match.\
If no sequences are provided, they will be extracted from the BLAST result and the entire BLAST will be analysed.

The program will try to find all pdb entries matching the user provided sequence(s) (or mutants thereof), but results may be incomplete due to partial matches within any tag-regions.
To compensate, the program will also use any accession codes found to double check with the pdb, using the RCSB PDB Search API: https://www.rcsb.org/docs/programmatic-access/web-apis-overview#search-api \
The user can also modify parameters in ´bl.py´ to be more or less permissive with regard to the allowed cover, identity, and number of gaps.
To allow for tags, up to one consecutive run of gaps and 90% cover are allowed by default, as well as 100% identity outside of any gaps.

If no pdb structures are found, the program will instead try to retrieve an AlphaFold model from the EBI AlphaFold Protein Structure Database: https://alphafold.ebi.ac.uk/api-docs.

PDB files will be downloaded using the RCSB PDB download service: https://www.rcsb.org/docs/programmatic-access/file-download-services \
AlphaFold structures will be downloaded using a link obtained from the EBI AlphaFold Protein Structure Database API.

Structures will be saved in PDBx/mmCIF format, using a folder tree whereby each folder corresponds to the FASTA descriptor provided by the user.

## Files included
```
af.py
bl.py
fp.py
pd.py
PySeqSt.py
README.md
seq.py
```
## File descriptions
### af.py
#### AlphFold
Functions: \
`def download(url: str) -> list:` \
"""Download the AlphaFold PDBx/mmCIF file at provided url"""\
`def from_uniprot(seqs: "Seqs"):`\
"""Find AlphaFold models in the EBI AlphaFold PSB for all sequences in a Seqs object
    which do not yet have a structure, using their UniProt Accession number.
    https://alphafold.ebi.ac.uk/api-docs"""
### bl.py
#### BLAST
Change the following parameters to tune what PySeqSt counts as a match to the query sequence when looking through the BLAST output
```
IDENT = 1
COV = 0.9
GAPS = 1
```
Functions: \
`def extractseq(blast: dict) -> "Seqs":` \
"""Process user provided BLAST file to extract all query sequences.\
    Return a new Seqs object."""\
`def process(blast: dict, seqs: "Seqs") -> int:`\
"""Go through a BLAST json dictionary and find pdb entries for each query.\
    Ignore descriptors not in Seqs object.\
    Add a list of those pdb entries to the Seqs object.\
    The first entry of the list indicates the type of structure.\
    Return the number of sequences that do not match the Seqs object."""
> inner functions:\
> `def _checkhit() -> bool:`\
"""Determines if a hit sequence matches the query sequence"""\
> `def _acc_pdb_from_hit():`\
"""Appends pdb number and accession number to the respective lists for that query"""

`def run(seqs: "Seqs") -> dict:`\
"""Format and submit request to ncib BLAST. Check for usage limits at:\
https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo

Dumps BLAST results to a file, and returns the json dictionary.

Adapted from https://blast.ncbi.nlm.nih.gov/docs/web_blast.pl and\
https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html#urlapi"""

Helper functions:\
`def _convert_to_uniprot(a: list) -> list:`\
    """Convert accession numbers to UniProt\
    Adapted from https://www.uniprot.org/help/id_mapping_prog\
    and https://www.uniprot.org/api-documentation/idmapping"""\
`def _b_request(payload :dict) -> "requests.Response":`\
"""Send a request to the BLAST API"""\
`def _outcome(s: str) -> int:`\
"""Process the BLAST status."""\
`def _up_post(payload :dict) -> str:`\
"""Send a POST request to the UniProt ID Mapping API\
    https://www.uniprot.org/help/id_mapping"""\
`def _up_get(path: str) -> dict:`\
"""Send a GET request to the UniProt ID Mapping API\
    https://www.uniprot.org/help/id_mapping"""
### fp.py
#### File Processing
Default output path:
```
PATH = "./out"
```
Functions: \
`def new_output_folder():`\
"""Create a new folder for the output"""\
`def open_blast(file: str) -> dict:`\
"""Open a user provide BLAST.json file and return the json object."""\
`def process_fasta(file: str) -> "Seqs":`\
"""Process user provided FASTA file of protein sequence(s).\
    Return a new Seqs object."""\
`def save_blast(blast: dict):`\
"""Save a BLAST dictionary to BLAST.json"""\
`def set_output_folder(path: str):`\
"""Set the PATH for any output files"""\
`def save_structures(seqs: "Seqs"):`\
"""Save all structures in Seqs object"""

Helper functions:\
`def _make_outputdirs(seqs: "Seqs"):`\
"""Create subfolders to store structures"""\
`def _cont():`\
"""Ask the user if they want to continue"""
### pd.py
#### PDB
Functions: \
`def download(pdb: str) -> list:`\
"""Download the PDBx/mmCIF file corresponsing to a pdb entry"""\
`def uniprot_to_pdb(seqs: "Seqs"):`\
"""Searches the pdb for any UniProt accession numbers for each descriptor in Seqs object.\
    Updates the list of structures for a descriptor with missing pdb numbers.\
    Overwrites any non-pdb structures stored in Seqs object if pdb structures were found.\
    https://search.rcsb.org/index.html#search-example-8"""
### PySeqSt.py
#### Main Program
Functions:\
`def main():`\
"""Main program flow"""\
`def parse_cl_args() -> "argparse.Namespace":`\
"""Parses command line arguments"""\
`def process_input(seq: str) -> "Seqs":`\
"""Process user provided protein sequence.\
    Return a new Seqs object"""

Helper functions:\
`def _regex_type(pattern: str | re.Pattern):`\
"""Argument type for matching a regex pattern.\
    Modified from https://stackoverflow.com/a/78770201"""
### README.md
This file :)
### seq.py
#### Seqs class
"""An object of unique aa sequences and associated strcutures & accession numbers"""
Methods:\
`def __init__(self):`\
"""Initializes a Seqs object,\
        where self._seqs is a dictionary of descriptor: aa sequence\
        and self._structures is a dictionary of descriptor: structures\
        and self._accessions is a dictionary of descriptor: accession numbers"""\
`def add_accessions(self, dscrptr: str, accessions: list, seq="", *, verify=True) -> bool:`\
"""Adds a list of accession numbers for a descriptor in object, unless empty,\
        after verifying the descriptor and sequence, returns outcome of that verification\
        Verification can be skipped by setting verify=False."""\
`def add_structures(self, dscrptr: str, structures: list, seq="", *, verify=True) -> bool:`\
"""Adds a list of structures for a descriptor in object, unless empty,\
        after verifying the descriptor and sequence, returns outcome of that verification".\
        Verification can be skipped by setting verify=False.\
        The first entry of the list is expected to indicate the type of structure(s)"""\
`def add_seq(self, dscrptr: str, seq: str) -> int:`\
"""Adds a sequence with a descriptor to object. Rejects duplicate sequences (returns 2).\
        Returns 1 if sequence is not a valid amino acid sequence.\
        Empty descriptors will be set to seq_# and duplicate descriptors will be appended with _#\
        Returns 0 if sequence could be added i.e. is valid and unique."""

`@property`\
`def accessions(self) -> dict:`\
"""Return a dictionary for all descriptor: accession numbers in object"""\
`@property`\
`def dscrptrs(self) -> KeysView :`\
"""Return a list of all descriptors in object"""\
`@property`\
`def seqs(self) -> set:`\
"""Return a set of all sequences in object"""\
`@property`\
`def structures(self) -> dict:`\
"""Return a dictionary for all descriptor: structures in object"""

`def get_fasta(self) -> str:`\
"""Return all descriptors, sequences as a fasta formated string."""\
`def get_seq(self, dscrptr: str) -> str:`\
"""Return sequence matching the descriptor"""\
`def get_structures(self, dscrptr: str) -> list:`\
"""Return list of structures matching the descriptor"""

Helper methods:\
`@staticmethod`\
`def _validate_seq(seq: str) -> bool:`\
"""Verify that an amino acid sequence is valid.\
        For now not included: B: D/N, J: L/I, O: Pyrrolysine, U: Sec, Z: E/Q, X: any"""\
`def _verify_dscrptr(self, d: str, s: str) -> str:`\
"""Verify a descriptor and sequence match what is stored in object.\
        Also check possible (de)incremented versions."""


Folder `testing_results` contains the output folder structure that should be obtained following the commands.
