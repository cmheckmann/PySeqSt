"""Given a protein sequences, this program performs a BLAST search
and finds known structures, or AlphaFold models"""
import argparse
import re
import sys

import af
import bl
import fp
import pd
from seq import Seqs

#Define long messages
E1 = "Please provide a valid sequence, in the form of single-letter canonical (20) amino acids.\n" \
     "To include '>...' statements, provide fasta file with -f"
W1 = "WARNING: {} BLAST-queries didn't match any user-provided descriptor+sequence & were skipped."


def main():
    """Main program flow"""
    # Parse command line arguments
    args = parse_cl_args()

    # Switch functionality depending on command line args
    if not (args.seq or args.file):
        seqs = process_input(input("Protein sequence: ").strip())
        print("Sequence successfully loaded from user input.")
    elif args.seq:
        seqs = process_input(args.seq.strip())
        print("Sequence successfully loaded from provided commandline input.")
    elif args.file:
        seqs = fp.process_fasta(args.file.strip())
        print(f"Sequences successfully loaded from '{args.file.strip()}'")

    # Set user specified output folder
    if args.out:
        fp.set_output_folder(args.out.strip())


    # Load prior blast results; or run a Blast on the Seqs object
    if args.blast:
        blast = fp.open_blast(args.blast.strip())

        # Set the output folder to where the BLAST result is
        if not args.out:
            fp.set_output_folder(args.blast.strip().rsplit("/")[0])

        # If seqs is not yet set, set it from BLAST result
        try:
            seqs
        except NameError:
            seqs = bl.extractseq(blast)
            print(f"Successfully loaded sequences from {args.blast.strip()}")

    else:
        # Create a new output folder and run BLAST
        if not args.out:
            fp.new_output_folder()
        blast = bl.run(seqs)

    if args.only_blast:
        print("Finished")
        sys.exit()
    # Process blast and add any pdb entry to the Seqs object
    print("Processing BLAST...")
    c = bl.process(blast, seqs)
    if c > 0:
        print(W1.format(c))
    else:
        print("Sucessfully processed BLAST")

    # Check UniProt accession numbers in pdb database
    # and get any pdb numbers that were missed due to e.g. tags, point mutations
    print("Using UniProt accession number to find missed pdb entries...")
    pd.uniprot_to_pdb(seqs)

    # Query EBI Alphafold for UniProt accessions where no pdb
    print("Retrieving AlphaFold models from EBI AlphaFold Protein Structure Database...")
    af.from_uniprot(seqs)

    # TODO Create Alphafold model where there is no structure

    # Save all structures
    print("Saving structures...")
    fp.save_structures(seqs)
    print("Finished")
    sys.exit()


def parse_cl_args() -> "argparse.Namespace":
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Given a protein sequence, or a file of sequences, "
                    "this program performs a BLAST search, "
                    "finds known structures, or AlphaFold models. "
                    "Default mode queries user for input of a sequence, "
                    "optionally a single sequence or a file "
                    "of sequences can be provided as a command line argument. "
                    "Program can be restarted at intermediate steps. "
                    "If a BLAST.json file is provided without sequences, "
                    "all query sequences found in that file will be extracted."
        )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "--file", metavar="*.fasta",
                       help="Run on a file of FASTA protein sequences",
                       type=_regex_type(r"\.(?:fasta|fa)$"))
    group.add_argument("-s", "--seq", metavar="SEQ",
                       help="Run on a single protein sequence",
                       type=str)
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("-b", "--blast", metavar="*.json",
                        help="Use results from priot BLAST. Saves results in same folder as BLAST. "
                        "WARNING: files may be overwritten. -o to specify alternative directory",
                        type=_regex_type(r"\.json$"))
    group2.add_argument("--only_blast", action='store_true',
                        help="Stop after running BLAST")
    parser.add_argument("-o", "--out", metavar="./path",
                        help="Relative path to save output. WARNING: files may be overwritten",
                        type=str)
    return parser.parse_args()


def process_input(seq: str) -> "Seqs":
    """Process user provided protein sequence.
    Return a new Seqs object"""
    seqs = Seqs()
    if seqs.add_seq("", seq) != 0:
        sys.exit(E1)
    return seqs


def _regex_type(pattern: str | re.Pattern):
    """Argument type for matching a regex pattern.
    Modified from https://stackoverflow.com/a/78770201"""

    def closure_check_regex(arg_value):
        if not re.search(pattern, arg_value, flags=re.IGNORECASE):
            raise argparse.ArgumentTypeError("Invalid file extension")
        return arg_value

    return closure_check_regex


if __name__ == "__main__":
    main()
