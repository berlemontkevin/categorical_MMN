#!/bin/bash
# This line tells the shell how to execute this script, and is unrelated
# to SLURM.
   
# at the beginning of the script, lines beginning with "#SBATCH" are read by
# SLURM and used to set queueing options. You can comment out a SBATCH
# directive with a second leading #, eg:
##SBATCH --nodes=1
   
# we need 1 node, will launch a maximum of one task and use one cpu for the task: 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
   
# we expect the job to finish within 5 hours. If it takes longer than 5
# hours, SLURM can kill it:
#SBATCH --time=24:00:00
   
# we expect the job to use no more than 2GB of memory:
#SBATCH --mem=2GB
   
# we want the job to be named "myTest" rather than something generated
# from the script name. This will affect the name of the job as reported
# by squeue:
#SBATCH --job-name=oddball-task-integrator
 
# when the job ends, send me an email at this email address.
#SBATCH --mail-type=END
#SBATCH --mail-user=kevin.berlemont@nyu.edu
   
# both standard output and standard error are directed to the same file.
# It will be placed in the directory I submitted the job from and will
# have a name like slurm_12345.out
SBATCH --output=slurm_%j.out
 
# once the first non-comment, non-SBATCH-directive line is encountered, SLURM
# stops looking for SBATCH directives. The remainder of the script is  executed
# as a normal Unix shell script
  
# first we ensure a clean running environment:
module purge
# and load the module for the software we are using:
module load julia/1.6.1
  
# next we create a unique directory to run this job in. We will record its
# name in the shell variable "RUNDIR", for better readability.
# SLURM sets SLURM_JOB_ID to the job id, ${SLURM_JOB_ID/.*} expands to the job
# id up to the first '.' We make the run directory in our area under $SCRATCH, because at NYU HPC
# $SCRATCH is configured for the disk space and speed required by HPC jobs.
RUNDIR=$SCRATCH/1-project-categorical-MMN/run-${SLURM_JOB_ID/.*}
mkdir $RUNDIR
  
# we will be reading data in from somewhere, so define that too:
DATADIR=$SCRATCH/1-project-categorical-MMN/data/scripts
  
# the script will have started running in $HOME, so we need to move into the
# unique directory we just created
cd $RUNDIR
  
# now start the julia job:
julia $DATADIR/script2-oddball-task-integrator-time-constant.jl
