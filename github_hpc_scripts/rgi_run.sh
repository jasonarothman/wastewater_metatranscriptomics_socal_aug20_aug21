#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=rgi_run      ## job name
#SBATCH -A katrine_lab     ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=60    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=slurm-%J.out ##output info file


module load anaconda/2020.07

eval "$(conda shell.bash hook)"

conda activate rgi

rgi load --card_json /dfs5/bio/rothmanj/databases/card.json
rgi load -i /dfs5/bio/rothmanj/databases/card.json --card_annotation /dfs5/bio/rothmanj/databases/card_database_v3.1.2.fasta
rgi load --wildcard_annotation /dfs5/bio/rothmanj/databases/wildcard_database_v3.0.2.fasta --wildcard_index /dfs5/bio/rothmanj/databases/wildcard/index-for-model-sequences.txt --card_annotation /dfs5/bio/rothmanj/databases/card_database_v3.1.2.fasta
rgi load --kmer_database /dfs5/bio/rothmanj/databases/wildcard/61_kmer_db.json --amr_kmers /dfs5/bio/rothmanj/databases/wildcard/all_amr_61mers.txt --kmer_size 61 --debug > kmer_load.61.log 2>&1

for f in $(ls *.fastq.1.gz | sed 's/.fastq.1.gz//' | sort -u)
do
rgi bwt -1 ${f}.fastq.1.gz -2 ${f}.fastq.2.gz -a bowtie2 -n 60 -o amr_analyses/${f} --clean --include_wildcard 
done