from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, SimpleLocation

s = SeqIO.read("config/chikv_reference.gb", "genbank")

s.features = s.features[:1]

s.features.append(SeqFeature(SimpleLocation(76, 1681), type="CDS", qualifiers={"locus_tag": "nsp1"}))
s.features.append(SeqFeature(SimpleLocation(1681, 4075), type="CDS", qualifiers={"locus_tag": "nsp2"}))
s.features.append(SeqFeature(SimpleLocation(4075, 5665), type="CDS", qualifiers={"locus_tag": "nsp3"}))
s.features.append(SeqFeature(SimpleLocation(5665, 7498), type="CDS", qualifiers={"locus_tag": "nsp4"}))


s.features.append(SeqFeature(SimpleLocation(7566, 8349), type="CDS", qualifiers={"locus_tag": "C"}))
s.features.append(SeqFeature(SimpleLocation(8349, 8541), type="CDS", qualifiers={"locus_tag": "E3"}))
s.features.append(SeqFeature(SimpleLocation(8541, 9810), type="CDS", qualifiers={"locus_tag": "E2"}))
s.features.append(SeqFeature(SimpleLocation(9810, 9993), type="CDS", qualifiers={"locus_tag": "6K"}))
s.features.append(SeqFeature(SimpleLocation(9993, 11310), type="CDS", qualifiers={"locus_tag": "E1"}))


for rec in s.features:
    print(rec, "\n")

# s.features[4].qualifiers["translation"][0].find("WS")

SeqIO.write(s, "config/chikv_reference_adjusted.gb", "genbank")

# SeqIO.write(s, "config/chikv_reference_adjusted.gff", "gff")