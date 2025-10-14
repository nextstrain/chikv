import argparse
from Bio import AlignIO, SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignment", required=True)
parser.add_argument("-o", "--output", required=True)
parser.add_argument("-s", "--segment-output", required=False)
args = parser.parse_args()

alignment = AlignIO.read(args.alignment, "fasta") # TODO: replace this with full alignment later

gene_start = 9993
gene_end = 11310

slice = slice(gene_start, gene_end) # TODO: check that we get the right region -> done

seq_to_keep = []

for rec in alignment:
    segment = rec.seq[slice]
    segment = str(segment)
    informative_bases_count = 0
    for base in segment:
        if base in "ATCG":
            informative_bases_count += 1
    coverage = informative_bases_count / len(segment)
    if coverage > 0.8:
        seq_to_keep.append(rec)


print(f"writing {len(seq_to_keep)} sequences to fasta file")
SeqIO.write(seq_to_keep, args.output, "fasta")

sliced_records = []
for rec in seq_to_keep:
    sliced_rec = rec[:]  # make a copy
    sliced_rec.seq = rec.seq[slice]
    sliced_rec.id = rec.id
    sliced_rec.description = f"{rec.description} [segment {gene_start}-{gene_end}]"
    sliced_records.append(sliced_rec)


print(f"Writing {len(sliced_records)} sliced sequences to {args.segment_output}")
SeqIO.write(sliced_records, args.segment_output, "fasta")