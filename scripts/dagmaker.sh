DAG_LINE='snakemake --dag (snakemake.input[0]) | dot -Tsvg > dag.svg' && echo $DAG_LINE > dag.sh
