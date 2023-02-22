#!/bin/bash

#SBATCH --job-name=cgrid_box
#SBATCH --output=output.txt
#SBATCH --time=10:00:00

####### 4 MPI ranks
#SBATCH --ntasks=8

mpirun ./cgrid_box /home/users/spetrov/SIMUG/SIMUG_v0/examples/H-CgridBox/config.txt
