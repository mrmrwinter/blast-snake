#Download genbank records:
from Bio import Entrez
idlist = sorted(list(set([line.strip() for line in open(snakemake.input[0])])))

f= open(snakemake.output[0],'w+')
Entrez.email = 'g.sellers@2011.hull.ac.uk'
handle = Entrez.efetch(db="nuccore", id=idlist, rettype="gb")
f.write(handle.read())
f.close()
