#!/bin/bash 

#compile hybridsim 
cd ../build 
make 

#copy it to scratch
cd ..
#cp -f bin/hybridsim_auto /scratch-new/formanek/HYBRIDSIM

cd input 

# parameters for job submission
CPUsPerTask=1
CPUsPerNode=64
TasksPerNode=$[$CPUsPerNode/$CPUsPerTask]
echo "will submit $TasksPerNode per node"  
nodes=("5")

mkdir -p slurm
currentNodeIndex=0
currentNode=${nodes[$currentNodeIndex]}
runcount=0
while read rundir; do 
    echo $rundir
    jobname="$(echo $rundir | sed 's/\/scratch-new\/formanek\/HYBRIDSIM\/runs\///')"
    echo $jobname
    echo $currentNode
    sedjobname="$(echo $jobname | sed 's/\./\\\./g')"
    sed "s/JOBNAME/$sedjobname/" slurm_template.slurm > slurm/$jobname.slurm
    sed -i "s/NODES/1/" slurm/$jobname.slurm
    sed -i "s/NTASKS/1/" slurm/$jobname.slurm  
    sed -i "s/NCPUS/$CPUsPerTask/" slurm/$jobname.slurm  
    sed -i "s/NODENUMBER/$currentNode/" slurm/$jobname.slurm
    sedrundir="$(echo $rundir | sed 's/\./\\\./g')"
    sedrundir="$(echo $sedrundir | sed 's/\//\\\//g')"
    sed -i "s/WORKDIR/$sedrundir/" slurm/$jobname.slurm
    runcount=$[$runcount+1]
    if [ "$runcount" -eq "$TasksPerNode" ]; then
        runcount=0
        currentNodeIndex=$[$currentNodeIndex+1]
        currentNode=${nodes[$currentNodeIndex]}
    fi 
    sbatch slurm/$jobname.slurm 
    echo $rundir >>submitted_runs.dat
    sleep 30
done <runs_to_submit.dat
