#! /bin/bash -login
#SBATCH -J snakemake-submission
#SBATCH -t 24:00:00
#SBATCH -p all
#SBATCH --mem 8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucas.penny@uhn.ca

. env.sh
module load samtools/1.9
module load picard/2.10.9
module load bwa/0.7.15
module load bamtools/2.3.0  
module load R
source activate frag-pipeline


echo 'Running on H4H cluster'
#snakemake -p --profile /cluster/home/t116306uhn/.config/snakemake/slurm
snakemake --rerun-incomplete --cluster "sbatch -p himem -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} -o slurm_log/{rule}_{wildcards}_out -e slurm_log/{rule}_{wildcards}_err" -p -j 20 --latency-wait 60