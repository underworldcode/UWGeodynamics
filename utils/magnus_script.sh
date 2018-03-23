#!/bin/bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=00:10:00
#SBATCH --image=docker:rbeucher/underworld2_geodynamics:dev
#SBATCH --volume="/scratch/q97/rb5533:/workspace/user_data"

echo "PRINTING ENVIRONMENT"
env

echo "PRINTING SLURM SCRIPT"
scontrol show job ${SLURM_JOBID} -ddd

module load shifter

# you can ignore the next two lines.. this is just pulling a script from *within* the container.
#shifter cp /opt/underworld2/docs/examples/1_06_Rayleigh_Taylor.ipynb .
#shifter /opt/underworld2/utils/ipynb_to_py.sh 1_06_Rayleigh_Taylor.ipynb

srun -n 48 shifter python 1_06_Rayleigh_Taylor.py
