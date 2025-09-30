from Bio import SeqIO
from Bio import SeqFeature, SimpleLocation

s = SeqIO.read("config/chikv_reference.gb", "genbank")

print(s.features)

s.features.append(SeqFeature)

s.features[4].qualifiers["translation"][0].find("WS")