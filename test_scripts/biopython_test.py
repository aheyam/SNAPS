# Test of Biopython sequence and alignment tools
from Bio import SeqIO

filename_1 = "data/P3a_L273R/P3a_L273R.fasta"
fasta_records_1 = SeqIO.parse(open(filename_1),"fasta")
record_mut = next(fasta_records_1)

filename_2 = "data/P3a_L273R/P3a.fasta"
fasta_records_2 = SeqIO.parse(open(filename_2),"fasta")
record_wt = next(fasta_records_2)

from Bio import Align
aligner = Align.PairwiseAligner()

alignments = aligner.align(record_wt, record_mut)