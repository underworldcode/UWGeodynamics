Running on Pawsey / Magnus
==========================

The recommended way to run Underworld / UWGeodynamics model is to use
Shifter. Shifter is a wrapper around Docker that allows us to run docker
containers on Magnus.

You can have a look at the `Pawsey
documentation <https://support.pawsey.org.au/documentation/display/US/Shifter>`__
if you want to know more about Shifter:

Pre-requisites
--------------

::

   ssh username@magnus-1.pawsey.org.au

**A UWGeodynamics docker image is already available on Magnus**

::

   user@magnus-1:~>module load shifter
   user@magnus-1:~>shifter images
   magnus     docker     READY    17cc3c02ba   2018-05-09T08:47:59 rbeucher/underworld2_geodynamics:magnus

The following command will pull the latest version of the image:

::

   shifter pull docker:rbeucher/underworld2_geodynamics:magnus

Setting up a job
----------------

Here we assume that we have a copy of the UWGeodynamics Tutorial 1 model
saved as a python file (*Tutorial_1_ThermoMechanical_Model.py*), inside
a folder *UWGeo_Tutorial1* located in the
/scratch/your-project-account/your-username folder:

::

   rb5533@magnus-1:/scratch/q97/rb5533/UWGeo_Tutorial1> ls 
   Tutorial_1_ThermoMechanical_Model.py

SLURM file
~~~~~~~~~~

Following is an example of a SLURM file (*job.slurm*) showing how to run
Tutorial 1 on 1 node using 4 cores:

::

   #!/bin/bash

   #SBATCH --nodes=1
   #SBATCH --time=00:10:00
   #SBATCH --account=q97

   echo "PRINTING ENVIRONMENT"
   env

   echo "PRINTING SLURM SCRIPT"
   scontrol show job ${SLURM_JOBID} -ddd

   module load shifter

   srun -n4 shifter run --mpi rbeucher/underworld2_geodynamics:magnus python Tutorial_1_ThermoMechanical_Model.py 

Running a job
-------------

After the above we should have the following files in our
*UWGeo_Tutorial1* folder:

::

   rb5533@magnus-1:/scratch/q97/rb5533/UWGeo_Tutorial1> ls 
   Tutorial_1_ThermoMechanical_Model.py    job.slurm

The job can now be submitted to the queue system using:

::

   sbatch job.slurm

That’s it!!!
