#!/bin/bash
#SBATCH --no-requeue
#SBATCH --job-name="k2y-example-silicon"
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --output=_scheduler-stdout.txt

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PARA_PREFIX=mpirun -np $SLURM_NTASKS

bash run.sh
