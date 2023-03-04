#!/bin/bash

#SBATCH --job-name=cgrid_cut_sphere
#SBATCH --output=output.txt
#SBATCH --time=24:00:00

####### MPI ranks
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12

mpirun ./cgrid_cut_sphere /home/users/spetrov/SIMUG/SIMUG_v0/examples/J-CgridCutSphere/config.txt