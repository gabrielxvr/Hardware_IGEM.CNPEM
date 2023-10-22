#!/bin/bash

#SBATCH --job-name=magsim
#SBATCH --nodes=2
#SBATCH --ntasks=24
#SBATCH --mem=24G

python3 mag_sim.py
