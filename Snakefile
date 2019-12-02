import os

# #configfile contains the sample wildcards and the paramater variables
configfile: "config.yaml"

#glob_wildcards
SAMPLES, = glob_wildcards("data/{sample}.fasta")

# run the following line of code first, straight into the bash terminal
# cd to the data directory, run the code, then cd bask into snakemake-blast
### rule samples_generation:
####	shell: "ls *.fasta > ../samplestest.txt"


# SAMPLE = list(set([line.strip() for line in open ("samples.txt")]))

### this is the target rule that snakemake will rely on if not given a target.
### also assists in buildings dags
rule all:
	input:
		expand("plots/{sample}.pdf", sample=SAMPLES)
		# expand("data/{sample}.fasta", sample=SAMPLES)


# ## this would be cool if you can ever get it working
# rule dag_creator:
# 	input:
# 		expand("plots/{sample}.pdf", sample=SAMPLES)
# 	script:
# 		"scripts/dagmaker.sh"


### this rule sends a query sequece to genbank and pulls the top 100 results
rule blast:
	input:
		"data/{sample}.fasta"
	output:
		"blast_out_accs/{sample}"
	shell:
		"blastn -db nt -query {input} -out {output} -evalue 0.001 -max_target_seqs 5 -outfmt '10 sacc' -num_threads 12"



rule acc_transform:
	input:
		"blast_out_accs/{sample}"
	output:
		"transformed_accs/{sample}"
	shell:
		"head -100 {input} | tr '\n' ',' | sed 's/,$//' > {output}"


#this next rule requires genbank login email
# this is currently grahams but needs to be made modular and eidtable from the configfile
rule genbank_pull:
	input:
		"blast_out_accs/{sample}"
	output:
		"pulled_gb/{sample}.gb"
	script:
		"scripts/gb_pull2.py"


# #this code will parse a genbank file and write a fasta containing GI, accession, species name, and seq.
rule gb_parsing:
 	input:
 		"pulled_gb/{sample}.gb"
	output:
		"blasted_fastas/{sample}"
	script:
		"scripts/gb_parser.py"


#this adds the query fastas to their respective msa to plot them in the tree with their blast hits
rule add_query_msa:
	input:
		"blasted_fastas/{sample}",
		"data/{sample}.fasta"
	output:
		"pre_mafft/{sample}"
	shell:
		"cat {input} > {output}"


rule mafft:
	input:
		"pre_mafft/{sample}"
	output:
		"mafft_out/{sample}"
	shell:
		"mafft --auto {input} > {output}"

# #auto - Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i and FFT-NS-2,
# #		according to data size. Default: off (always FFT-NS-2)

# #localpair - All pairwise alignments are computed with the Smith-Waterman algorithm.
# #		more accurate but slower than --6merpair. Suitable for a set of locally alignable
# #		sequences. Applicable to up to ~200 sequences. A combination with --maxiterate 1000 is
# #		recommended (L-INS-i). Default: off (6mer distance is used)


# #trimal cleans the ends of the alignment
rule trimal: # alignment curation
    input:
        "mafft_out/{sample}"
    output:
        "trimal/{sample}"
    shell:
        "trimal -in {input} -out {output} -gappyout -keepheader"


# #this rule processes an msa into a newick
rule fasttree:
	input:
		"trimal/{sample}"
	output:
		"newicks/{sample}.tree"
	shell:
		"FastTree -gtr -nt {input} > {output}"
# #GTR - GTR+CAT model. Generalised time-reversible
# # https://hpc-carpentry.github.io/hpc-python/15-snakemake-python/


# #this will write/draw a phylogeny from a newick file
rule ete3:
	input:
		"newicks/{sample}.tree"
	output:
		"plots/{sample}.pdf"
	script:
		"scripts/ete.py"
