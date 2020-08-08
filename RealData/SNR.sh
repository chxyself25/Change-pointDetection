#!/bin/bash

#SBATCH --nodes=1 # request one node

#SBATCH --cpus-per-task=51  # ask for 8 cpus

#SBATCH --mem-per-cpu=2G #128GB of ram PER NODE, if you need more granular control lookup the --mem-per-cpu argument in the man page for sbatch

#SBATCH --partition=biocrunch # specify which partition your job should be submitted to

#SBATCH --time=0-12:00:00 # ask that the job be allowed to run for 24 days

# optional nice stuff below

#SBATCH --error=job.%J.err # tell it to store the error messages to a file

#SBATCH --output=job.%J.out # tell it to store the console text output to a file

#SBATCH --job-name="calculate SNR" # a nice name for the job to have

# let's load some modules
module load r-doparallel/1.0.11-py2-r3.5-tlbjucn
module load r-devtools/1.12.0-py2-r3.5-3zfj3n2
module load r-expm/0.999-2-py2-r3.5-p6ucwpc
module load r-hmisc/4.1-1-py2-r3.5-2wm5k3n
module load r-numderiv/2016.8-1-py2-r3.5-5pb7s6o
module load r-mass/7.3-47-py2-r3.5-xtsjvcy
module load gcc/7.3.0-xegsmw4
module load r/3.6.0-py2-fupx2uq

# here we use the srun command in place of mpirun, since it has much better integration with the scheduler.
# srun --ntasks=544 gmx_mpi mdrun -s science.tpr -maxh 0.80

# let's make sure we're where we expect to be in the filesystem tree

cd /work/LAS/zhuz-lab/xchang/Change-pointDetection/RealData/

# the commands we're running are below

Rscript ./SNR.R
