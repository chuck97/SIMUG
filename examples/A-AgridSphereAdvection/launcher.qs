#!/bin/bash

#SBATCH --job-name=agrid_sphere_adv
#SBATCH --output=output.txt
#SBATCH --time=24:00:00

####### MPI ranks
#SBATCH --ntasks=8

mpirun ./agrid_sphere_advection /home/users/spetrov/SIMUG/SIMUG_v0/examples/A-AgridSphereAdvection/config.txt