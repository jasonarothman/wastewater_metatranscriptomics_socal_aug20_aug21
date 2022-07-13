#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=braken_run      ## job name
#SBATCH -A katrine_lab     ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=36    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=slurm-%J.out ##output info file

module load anaconda/2020.07

eval "$(conda shell.bash hook)"

conda activate /data/homezvol0/rothmanj/.conda/envs/conda_software

for f in $(ls *.report | sed 's/.report//' | sort -u)
do
bracken -t 0 -d /dfs5/bio/rothmanj/databases/kraken_db_complete -i ${f}.report -o ${f}.species.bracken -w ${f}_report_estimates.species.report -r 100 -l S
done
