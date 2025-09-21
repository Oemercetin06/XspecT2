"""
Mapping handler for the alignment-based misclassification detection.

Notes:
Developed by Oemer Cetin as part of a Bsc thesis at Goethe University Frankfurt am Main (2025).
(An Integration of Alignment-Free and Alignment-Based Approaches for Bacterial Taxon Assignment)
"""

import mappy, pysam, os, csv
from Bio import SeqIO
from xspect.definitions import fasta_endings

__author__ = "Cetin, Oemer"


class MappingHandler:
    """Handler class for all mapping related procedures."""

    def __init__(self, ref_genome_path: str, reads_path: str) -> None:
        """
        Initialise the mapping handler.

        This method sets up the paths to the reference genome and query sequences.
        Additionally, the paths to the output formats (SAM, BAM and TSV) are generated.

        Args:
            ref_genome_path (str): The path to the reference genome.
            reads_path (str): The path to the query sequences.
        """
        if not os.path.isfile(ref_genome_path):
            raise ValueError("The path to the reference genome does not exist.")

        if not os.path.isfile(reads_path):
            raise ValueError("The path to the reads does not exist.")

        if not ref_genome_path.endswith(tuple(fasta_endings)) and reads_path.endswith(
            tuple(fasta_endings)
        ):
            raise ValueError("The files must be FASTA-files!")

        stem = reads_path.rsplit(".", 1)[0] + "_mapped"
        self.ref_genome_path = ref_genome_path
        self.reads_path = reads_path
        self.sam = stem + ".sam"
        self.bam = stem + ".sorted.bam"
        self.tsv = stem + ".start_coordinates.tsv"

    def map_reads_onto_reference(self) -> None:
        """
        A Method that maps reads against the respective reference genome.

        This function creates a SAM file via Mappy and converts it into a BAM file.
        """
        # create header (entry = sequences of the reference genome)
        ref_seq = [
            {"SN": rec.id, "LN": len(rec.seq)}
            for rec in SeqIO.parse(self.ref_genome_path, "fasta")
        ]
        header = {"HD": {"VN": "1.0"}, "SQ": ref_seq}
        target_id = {sequence["SN"]: number for number, sequence in enumerate(ref_seq)}

        reads = list(SeqIO.parse(self.reads_path, "fasta"))
        if not reads:
            raise ValueError("Reads file is empty.")

        read_length = len(reads[0].seq)
        preset = "map-ont" if read_length > 150 else "sr"
        # create SAM-file
        aln = mappy.Aligner(self.ref_genome_path, preset=preset)
        with pysam.AlignmentFile(self.sam, "w", header=header) as out:
            for read in reads:
                read_seq = str(read.seq)
                for hit in aln.map(read_seq):
                    if hit.cigar_str is None:
                        continue
                    # add soft-clips so CIGAR length == len(read_seq) IMPORTANT!!
                    leftS = hit.q_st
                    rightS = len(read_seq) - hit.q_en
                    cigar = (
                        (f"{leftS}S" if leftS > 0 else "")
                        + hit.cigar_str
                        + (f"{rightS}S" if rightS > 0 else "")
                    )

                    mapped_region = pysam.AlignedSegment()
                    mapped_region.query_name = read.id
                    mapped_region.query_sequence = read_seq
                    mapped_region.flag = 16 if hit.strand == -1 else 0
                    mapped_region.reference_id = target_id[hit.ctg]
                    mapped_region.reference_start = hit.r_st
                    mapped_region.mapping_quality = (
                        hit.mapq or 255
                    )  # 0-60 (255 means unavailable)
                    mapped_region.cigarstring = cigar
                    out.write(mapped_region)
                    break  # keep only primary

        # create BAM-file
        pysam.sort("-o", self.bam, self.sam)
        pysam.index(self.bam)

    def get_total_genome_length(self) -> int:
        """
        Get the genome length from a BAM-file.

        This function opens a BAM-file and extracts the genome length information.

        Returns:
            int: The genome length.
        """
        with pysam.AlignmentFile(self.bam, "rb") as bam:
            return sum(bam.lengths)

    def extract_starting_coordinates(self) -> None:
        """
        Extract starting coordinates of mapped regions from a BAM-file.

        This function scans through a BAM-file and creates a TSV-file.
        The information that is extracted is the starting coordinate for each mapped read.
        """
        # create tsv-file with all start positions
        with open(self.tsv, "w") as tsv:
            tsv.write("reference_genome\tread\tmapped_starting_coordinate\n")
            try:
                with pysam.AlignmentFile(self.bam, "rb") as bam:
                    entry = {
                        i: seq["SN"] for i, seq in enumerate(bam.header.to_dict()["SQ"])
                    }
                    seen = set()
                    for ref_seq in bam.references:
                        for hit in bam.fetch(ref_seq):
                            if (
                                hit.is_unmapped
                                or hit.is_secondary
                                or hit.is_supplementary
                            ):
                                continue
                            key = (hit.reference_id, hit.reference_start)
                            if key in seen:
                                continue
                            seen.add(key)
                            tsv.write(
                                f"{entry[hit.reference_id]}\t{hit.query_name}\t{hit.reference_start}\n"
                            )
            except ValueError:
                tsv.write("dummy_reference\tdummy_read\t1000\n")

    def get_start_coordinates(self) -> list[int]:
        """
        Get the coordinates of a TSV-file.

        This function opens a TSV-file and saves all starting coordinates in a list.

        Returns:
            list[int]: The list containing all starting coordinates.

        Raises:
            ValueError: If no column with starting coordinates is found.
        """
        coordinates = []
        with open(self.tsv, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                val = row.get("mapped_starting_coordinate")
                if val is None:
                    raise ValueError("Column with starting coordinates not found.")
                coordinates.append(int(val))
        return coordinates
