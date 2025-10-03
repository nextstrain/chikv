import argparse
from Bio import AlignIO, SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignment", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

alignment = AlignIO.read(args.alignment, "fasta") # TODO: replace this with full alignment later

gene_start = 9993
gene_end = 11310

slice = slice(gene_start, gene_end) # TODO: check that we get the right region -> done

seq_to_keep = []

for rec in alignment:
    segment = rec.seq[slice]
    print(segment)
    segment = str(segment)
    informative_bases_count = 0
    for base in segment:
        if base in "ATCG":
            informative_bases_count += 1
    coverage = informative_bases_count / len(segment)
    if coverage > 0.8:
        seq_to_keep.append(rec)

SeqIO.write(seq_to_keep, args.output, "fasta")