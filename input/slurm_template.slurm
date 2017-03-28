#!/bin/bash
     
#SBATCH --job-name=JOBNAME
#SBATCH --nodes=NODES
#SBATCH --ntasks=NTASKS
#SBATCH --cpus-per-task=NCPUS
#SBATCH --time=infinite
#SBATCH --nodelist nostromo-0-NODENUMBER

exit_handler() {
    date 
    echo "signal caught by shell script" 
    kill -SIGINT "${child}" 
}

trap exit_handler SIGINT INT TERM SIGTERM

 
# go to working directory
workdir=WORKDIR
cd $workdir

# Load necessary modules
module load foss

#set OMP Parameters
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=CLOSE 

#run
../../hybridsim_auto parameters.dat > log.out &
child=$!
wait "$child"
