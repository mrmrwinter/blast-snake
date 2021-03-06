                       ######### BLAST-SNAKE #########

################################################################################

configfile: "config.yaml"

#glob_wildcards
SAMPLES, = glob_wildcards("data/queries/{sample}.fasta")

### this is the target rule that snakemake will rely on if not given a target.
### also assists in buildings dags
rule all:
	input:
		expand("plots/{sample}.pdf", sample=SAMPLES)

################################################################################

#### BLASTING THE QUERY
# This rule sends a query sequence to genbank and pulls the top 'x' results
# in the form of accession numbers
rule blast:
	conda:
		"envs/blastsnake.yaml"
	input:
		"data/queries/{sample}.fasta"
	output:
		"data/blast_out_accs/{sample}"
	shell:
		"blastn -db nt -query {input} -out {output} -evalue 0.001 -max_target_seqs 5 -outfmt '10 sacc' -remote"

################################################################################

#### THIS STEP MIGHT NOT BE NECCESSARY

# #### DATA TRANSFORMATION
# ## This step converts newlines (\n) in the file to commas (,)
# ## This is the required format for the next rule
# rule acc_transform:
# 	input:
# 		"data/blast_out_accs/{sample}"
# 	output:
# 		"data/transformed_accs/{sample}"
# 	shell:
# 		"head -100 {input} | tr '\n' ',' | sed 's/,$//' > {output}"

################################################################################

#### PULLING FILES FROM GENBANK
# This rule requires genbank login email. This is currently grahams login but
# needs to be made modular and editable from the configfile.
# This rule pulls the genbank file that correlates with the accession numbers.
rule gb_pull:
	conda:
		"envs/blastsnake.yaml"
	input:
		"data/blast_out_accs/{sample}"
	output:
		"data/pulled_gb/{sample}.gb"
	params:
		"g.sellers@2011.hull.ac.uk"
	script:
		"scripts/gb_pull.py"

################################################################################

#### PARSING GENBANK FILES
# This code will parse a genbank file and write a fasta containing the sequence,
# and the GI, accession, and species name in the header.
rule gb_parsing:
	conda:
		"envs/blastsnake.yaml"
	input:
 		"data/pulled_gb/{sample}.gb"
	output:
		"data/parsed_gb/{sample}"
	script:
		"scripts/gb_parser.py"

################################################################################

#### TRANSFORM SPACES TO UNDERSCORES
# this will transform all spaces to underscores

rule sed_space2underscore:
	conda:
		"envs/blastsnake.yaml"
	input:
		"data/parsed_gb/{sample}"
	output:
		"data/trans_gb/{sample}"
	shell:
		"sed -e 's/ /_/g' {input} > {output}"

#### ADD QUERY TO MSA
# This adds the query sequence to their respective multi-fasta of hits
# This allows us to plot them in the tree with their respective BLASTN hits
rule add_query_msa:
	conda:
		"envs/blastsnake.yaml"
	input:
		"data/trans_gb/{sample}",
		"data/queries/{sample}.fasta"
	output:
		"data/pre_mafft/{sample}.msa"
	shell:
		"cat {input} > {output}"

################################################################################

#### ALIGNMENT WITH MAFFT
# MAFFT takes multiple FASTAs and aligns them based on homology
# This outputs an msa (multiple sequence alignment)
rule mafft:
	conda:
		"envs/blastsnake.yaml"
	input:
		"data/pre_mafft/{sample}.msa"
	output:
		"data/mafft_out/{sample}.aln"
	shell:
		"mafft --auto {input} > {output}"

################################################################################

#### TRIMAL TRIMMING
# Trimal cleans the ends of the alignment
# Uneven sequences in an alignment can disrupt the output
rule trimal:
	conda:
		"envs/blastsnake.yaml"
	input:
		"data/mafft_out/{sample}.aln"
	output:
		"data/trimal/{sample}.tml.aln"
	shell:
		"trimal -in {input} -out {output} -gappyout -keepheader"

################################################################################

#### FASTTREE
# This rule processes an msa into a newick file
# Newick files contain information about sequences relationships to each other
# This is necessary to plot a phylogeny
rule fasttree:
	conda:
		"envs/blastsnake.yaml"
	input:
		"data/trimal/{sample}.tml.aln"
	output:
		"data/newicks/{sample}.tree"
	shell:
		"FastTree -gtr -nt {input} > {output}"

## GTR - GTR+CAT model. Generalised time-reversible
## https://hpc-carpentry.github.io/hpc-python/15-snakemake-python/

################################################################################

#### ETE3 TREES
# This will write/draw a phylogeny from a newick file
# This is currently in a very barebones state
rule ete3:
	conda:
		"envs/blastsnake.yaml"
	input:
		"data/newicks/{sample}.tree"
	output:
		"plots/{sample}.pdf"
	script:
		"scripts/ete.py"

################################################################################

################################################################################







# BELOW HERE IS HALF FINISHED RULES THAT MIGHT OR MIGHT NOT BE USEFUL
# AS WELL AS OTHER WAYS OF PARSING BLAST

###############################################################################
# ## this would be cool if you can ever get it working
# rule dag_creator:
# 	input:
# 		expand("plots/{sample}.pdf", sample=SAMPLES)
# 	script:
# 		"scripts/dagmaker.sh"


# #parses a tabular file and outputs GI, accession, species name, and seq, writing them to fasta
# rule tab_to_fasta:
# 	input:
# 		"blast_out/{sample}"
# 	output:
# 		"blasted_fastas/{sample}"
# 	shell:
# 		"awk 'BEGIN { # }{ print '>'$1, $2 }' {input} > {output}"
# #
# # #-db nt = database nucleotide
# # #-remote = sends search to NCBI servers
# # #-outfmt = changes format of output. in this case xml is the favoured format
# #
# # #this will parse an xml file for seq and title
# #
# rule xml_to_fasta:
# 	run:
# 		blast_record = NCBIXML.read({input})
# 		write....
# # #
# #
# #
# # # #this rule converts sams to bams
# rule sam_to_bam:
# 	input:
# 		"blasted_sams/{sample}.sam"
# 	output:
# 		"blasted_bams/{sample}.bam"
# 	shell:
# 		"samtools view -Sb {input} > {output}"
# #
# #
# # #this rule converts bams to fasta
# rule bam_to_fasta:
# 	input:
# 		"blasted_bams/{sample}.bam"
# 	output:
# 		"blasted_fastas/{sample}.fasta"
# 	shell:
# 		"samtools fasta {input} -0 {output}"
# # #
# # #
# # # #this runs fastas through mafft to align them
#
# rule efetch_fasta_download_script:
# 	input:
# 		"transformed_accs/{sample}"
# 	output:
# 		"scripts/fasta_pull_{sample}.sh"
# 	shell:
# 		"sed 's/^/efetch -db nucleotide -format fasta -id /' {input} > {output}"
#
# rule fasta_from_acc:
# 	input:
# 		"scripts/fasta_pull_{sample}.sh"
# 	output:
# 		"pulled_fastas/{sample}.txt"
# 	shell:
# 		"bash {input} > {output}"
#
# rule gb_from_acc:
# 	input:
# 		"scripts/gb_pull_{sample}.sh"
# 	output:
# 		"pulled_gb/{sample}.gb"
# 	shell:
# 		"bash {input} > {output}"
#
# # #this is going to need the login thing
# rule efetch_genbank_download_script:
# 	input:
# 		"transformed_accs/{sample}"
# 	output:
# 		"scripts/gb_pull_{sample}.sh"
# 	shell:
# 		"sed 's/^/efetch -db nucleotide -id /' {input} > {output}"
