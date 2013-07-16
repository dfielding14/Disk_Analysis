#!/bin/bash
# Job name:
#SBATCH --job-name=shell_stella
#
# Partition
#SBATCH --partition=henyey
#
# QoS:
#SBATCH --qos=small
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#
# Wall Clock Time:
#SBATCH --time=12:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=dfielding@berkeley.edu

## Run Command
mpirun -np 3 python shell_analyzer_stella.py --parallel

