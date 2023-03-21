#!/bin/bash

#SBATCH --job-name=cgrid_sphere_adv
#SBATCH --output=output.txt
#SBATCH --time=24:00:00

####### MPI ranks
#SBATCH --ntasks=8

mpirun ./cgrid_sphere_advection /home/users/spetrov/SIMUG/SIMUG_v0/examples/B-CgridSphereAdvection/config.txt