#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=humann3      ## job name
#SBATCH -A katrine_lab     ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=60    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=slurm-%J.out ##output info file


module load anaconda/2020.07

eval "$(conda shell.bash hook)"

conda activate biobakery3

for f in $(ls *.merged.fastq.gz | sed 's/.merged.fastq.gz//' | sort -u)
do
humann --input ${f}.merged.fastq.gz, -o ${f}.human.profile/ --threads 60
done