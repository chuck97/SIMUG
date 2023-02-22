#!/bin/bash

#SBATCH --job-name=agrid_box
#SBATCH --output=output.txt
#SBATCH --time=5:00:00

####### 4 MPI ranks
#SBATCH --ntasks=8

mpirun ./agrid_box /home/users/spetrov/SIMUG/SIMUG_v0/examples/G-AgridBox/config.txt