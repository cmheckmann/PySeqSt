"""Defining the Seqs class"""
from collections.abc import KeysView
import re

class Seqs:
    """An object of unique aa sequences and associated strcutures ? accession numbers"""
    def __init__(self):
        """Initializes a Seqs object,
        where self._seqs is a dictionary of descriptor: aa sequence
        and self._structures is a dictionary of descriptor: structures
        and self._accessions is a dictionary of descriptor: accession numbers"""

        # Descriptor: sequence
        self._seqs = {}
        # Set of all unique sequences
        self._useqs = set()
        # Descriptor: strcutures (ids)
        self._structures = {}
        # Descriptor: accession numbers
        self._accessions = {}

    def add_accessions(self, dscrptr: str, accessions: list, seq="", *, verify=True) -> bool:
        """Adds a list of accession numbers for a descriptor in object, unless empty,
        after verifying the descriptor and sequence, returns outcome of that verification
        Verification can be skipped by setting verify=False."""

        # Verify descriptor and sequence match what is in object
        if verify:
            dscrptr = self._verify_dscrptr(dscrptr, seq)
            if not dscrptr:
                return False
        if len(accessions) > 0:
            self._accessions[dscrptr] = accessions
        return True

    def add_structures(self, dscrptr: str, structures: list, seq="", *, verify=True) -> bool:
        """Adds a list of structures for a descriptor in object, unless empty,
        after verifying the descriptor and sequence, returns outcome of that verification".
        Verification can be skipped by setting verify=False.
        The first entry of the list is expected to indicate the type of structure(s)"""

        # Valid structure tpyes:
        valid = ["pdb", "EBI-AF"]

        # Verify descriptor and sequence match what is in object
        if verify:
            dscrptr = self._verify_dscrptr(dscrptr, seq)
            if not dscrptr:
                return False

        # Verify structure type
        try:
            if structures[0] not in valid:
                raise IndexError
        except IndexError:
            raise ValueError(f"invalid structure type in first element, expect {valid}.")

        # Add structures if there are any
        if len(structures) > 1:
            self._structures[dscrptr] = structures
        return True

    def add_seq(self, dscrptr: str, seq: str) -> int:
        """Adds a sequence with a descriptor to object. Rejects duplicate sequences (returns 2).
        Returns 1 if sequence is not a valid amino acid sequence.
        Empty descriptors will be set to seq_# and duplicate descriptors will be appended with _#
        Returns 0 if sequence could be added i.e. is valid and unique."""

        seq = seq.strip().upper()
        dscrptr = dscrptr.strip()
        # Sequence not already in object
        if seq not in self.seqs:
            # Sequence is valid
            if self._validate_seq(seq):

                # Check for empty descriptor
                if dscrptr == "":
                    dscrptr = "seq_" + str(len(self._seqs) + 1)

                # Loop until desciptor is unique
                while True:
                    if dscrptr in self._seqs:
                        d = dscrptr.rsplit("_", 1)
                        if len(d) == 1 or not d[1].isdecimal():
                            d.append("2")
                        else:
                            d[1] = str(int(d[1]) + 1)
                        dscrptr = "_".join(d)
                    else:
                        break

                # add to self._seqs and self._useqs
                self._seqs[dscrptr] = seq
                self._useqs.add(seq)
                return 0

            # Invalid sequence
            return 1

        # Duplicate sequence
        return 2

    @property
    def accessions(self) -> dict:
        """Return a dictionary for all descriptor: accession numbers in object"""

        return self._accessions

    @property
    def dscrptrs(self) -> KeysView :
        """Return a list of all descriptors in object"""

        return self._seqs.keys()

    @property
    def seqs(self) -> set:
        """Return a set of all sequences in object"""
        return self._useqs

    @property
    def structures(self) -> dict:
        """Return a dictionary for all descriptor: structures in object"""

        return self._structures

    def get_fasta(self) -> str:
        """Return all descriptors, sequences as a fasta formated string."""

        out = ""
        for dscrptr, seq in self._seqs.items():
            out += f">{dscrptr}\n{seq}\n"
        return out.rstrip("\n")

    def get_seq(self, dscrptr: str) -> str:
        """Return sequence matching the descriptor"""

        return self._seqs.get(dscrptr)

    def get_structures(self, dscrptr: str) -> list:
        """Return list of structures matching the descriptor"""

        return self._structures.get(dscrptr)

    @staticmethod
    def _validate_seq(seq: str) -> bool:
        """Verify that an amino acid sequence is valid.
        For now not included: B: D/N, J: L/I, O: Pyrrolysine, U: Sec, Z: E/Q, X: any"""

        pattern = r"[ACDEFGHIKLMNPQRSTWVY-]+\*?"
        return re.fullmatch(pattern, seq) is not None

    def _verify_dscrptr(self, d: str, s: str) -> str:
        """Verify a descriptor and sequence match what is stored in object.
        Also check possible (de)incremented versions."""

        d = d.strip()
        s = s.strip().upper()

        # Remove any increments in supplied descriptor
        d_split = d.rsplit("_", 1)
        if len(d_split) == 1 or not d_split[1].isdecimal():
            pattern = rf"{d}(?:_\d+)?$"
        else:
            pattern = rf"{d_split[0]}(?:_\d+)?$"

        # Check if the deincremented descriptor (or incremented versions thereof)
        # match a descriptor in the object, with matching sequence
        for dscrptr in self.dscrptrs:
            if (re.match(pattern, dscrptr)
                and self.get_seq(dscrptr) == s):
                return dscrptr

        return ""
