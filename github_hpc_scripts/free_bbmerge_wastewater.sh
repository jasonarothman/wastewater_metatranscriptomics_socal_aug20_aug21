#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=bbmerge_sewage      ## job name
#SBATCH -A katrine_lab     ## account to charge
#SBATCH -p free          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=30    ## number of cores the job needs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rothmanj@uci.edu
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=slurm-%J.out ##output info file

module load bbmap/38.87

for f in $(ls *.nohuman.fastq.1.gz | sed 's/.nohuman.fastq.1.gz//' | sort -u)
do
reformat.sh \
in1=${f}.nohuman.fastq.1.gz \
in2=${f}.nohuman.fastq.2.gz \
out=${f}.nohuman.merged.fastq.gz \
threads=30
done
