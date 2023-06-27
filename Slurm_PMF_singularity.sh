#!/bin/bash

#SBATCH --partition=normal
#SBATCH --job-name=wine_run

## NOTE: %u=userID, %x=jobName, %N=nodeID, %j=jobID, %A=arrayID, %a=arrayTaskID
#SBATCH --output=pmf_noGUI_try_%N_%j.out # output file
#SBATCH --error=pmf_noGUI_try_%N_%j.err # error file

## Email info for updates from Slurm
#SBATCH --mail-type=BEGIN,END,FAIL # ALL,NONE,BEGIN,END,FAIL,REQUEUE,..
#SBATCH --mail-user=tzhang23@gmu.edu

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1                # up to 128;
#SBATCH --mem-per-cpu=3500M                # memory per CORE

#SBATCH --export=NONE
#SBATCH --time=0-01:00:00                  # set to 1hr; please choose carefully

##SBATCH --array=1-12                       # Runs the program 12 times

module load singularity

SINGULARITY_BASE=/containers/hopper/Containers
CONTAINER=${SINGULARITY_BASE}/wine/wine.sif
SINGULARITY_RUN="singularity exec  -B ${PWD}:/host_pwd --pwd /host_pwd"
SCRIPT=PMF_bs_6f8xx_sealed_GUI_MOD.ini

${SINGULARITY_RUN} ${CONTAINER} wine ${PWD}/ME-2.exe ${SCRIPT} 