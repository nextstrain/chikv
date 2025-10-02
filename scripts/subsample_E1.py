from Bio import AlignIO, SeqIO

alignment = AlignIO.read("results/ingest/general/aligned.fasta", "fasta") # replace this with full alignment later

gene_start = 9993
gene_end = 11310

slice = slice(gene_start, gene_end) # check that we get the right region

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

SeqIO.write(seq_to_keep, "scripts/E1_sequences.fasta", "fasta")