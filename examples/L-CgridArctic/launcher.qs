#!/bin/bash

#SBATCH --job-name=agrid_cut_sphere
#SBATCH --output=output.txt
#SBATCH --time=12:00:00

####### MPI ranks
#SBATCH --ntasks=8

mpirun ./agrid_cut_sphere /home/users/spetrov/SIMUG/SIMUG_v0/examples/I-AgridCutSphere/config.txt