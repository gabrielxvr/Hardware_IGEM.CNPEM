#!/bin/bash

#SBATCH --nodes=1                      #Numero de N  s
#SBATCH --ntasks-per-node=8            #Numero de tarefas por N
#SBATCH --ntasks=8                     #Numero total de tarefas MPI
#SBATCH -p head                        #Fila (partition) a ser utilizada
#SBATCH -J cap                         #Nome job

python3 mag_sim.py
