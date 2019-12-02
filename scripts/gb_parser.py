#this transforms to fasta record with accession and species name +seq
from Bio import SeqIO
fast=open(snakemake.output[0],'w')

with open(snakemake.input[0], 'r') as tit:
    for rec in SeqIO.parse(tit, "genbank"):
        source = [f for f in rec.features if f.type == 'source'][0]
        sp=(" ".join(source.qualifiers['organism'][0].split()[0:2]))
        fast.write(">%s %s\n%s\n" %(rec.id,sp,rec.seq))
