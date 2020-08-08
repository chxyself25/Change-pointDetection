#!/bin/bash

#BATCH --nodes=1 # request one node

#SBATCH --cpus-per-task=51  # ask for 8 cpus

#SBATCH --mem=128G #128GB of ram PER NODE, if you need more granular control lookup the --mem-per-cpu argument in the man page for sbatch

#SBATCH --partition=biocrunch # specify which partition your job should be submitted to

#SBATCH --time=2-00:00:00 # ask that the job be allowed to run for 24 days

# optional nice stuff below

#SBATCH --error=job.%J.err # tell it to store the error messages to a file

#SBATCH --output=job.%J.out # tell it to store the console text output to a file

#SBATCH --jobname="simulate Kd distribution" # a nice name for the job to have

# let's load some modules
module load r-doparallel/1.0.11-py2-r3.5-tlbjucn
module load r-expm/0.999-2-py2-r3.5-p6ucwpc
module load r-hmisc/4.1-1-py2-r3.5-2wm5k3n
module load r-numderiv/2016.8-1-py2-r3.5-5pb7s6o
module load gcc
module load r/3.6.0-py2-fupx2uq

# let's make sure we're where we expect to be in the filesystem tree

cd /work/LAS/zhuz-lab/xchang/Change-pointDetection/Kd_simulation

 #the commands we're running are below

Rscript ./Kd_sim.R
