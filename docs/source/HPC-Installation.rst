HPC Installation of Underworld / UWGeodynamics
==============================================

Installing Underworld2 and UWGeodynamics
----------------------------------------

Python 2:
~~~~~~~~~

1. Install a python-2.7.11 virtual environment in your HOME:

.. code:: bash

  module purge
  module load python/2.7.11 openmpi/3.1.2 hdf5/1.10.2p 
  cd ~
  pip install --user virtualenv
  ~/.local/bin/virtualenv python-2.7.11-venv
  source python-2.7.11-venv/bin/activate
  export HDF5_VERSION=1.10.2
  pip install --no-binary=mpi4py mpi4py
  CC="mpicc" HDF5_MPI="ON" HDF5_DIR=$HDF5_DIR pip install --no-binary=h5py h5py


2. Install Underworld

.. code:: bash
   
   UW_DIR=path-to-your-install
   git clone https://github.com/underworldcode/underworld2.git $UW_DIR
   
   module purge
   module load gcc/5.2.0 hdf5/1.10.2p petsc/3.8.4 swig/3.0.12 python/2.7.11 openmpi/3.1.2
   source path-to-your-python-venv
   export HDF5_VERSION=1.10.2
   
   cd $UW_DIR/libUnderworld
   ./configure.py --with-debugging=0
   ./compile.py -j4

    cd ~/python-2.7.11/lib/python2.7/site-packages
    echo "/short/q97/Underworld/underworld_latest" > underworld.pth
    echo "/short/q97/Underworld/underworld_latest/glucifer" > glucifer.pth

Now try to import underworld

.. code:: bash
 
    cd ~
    source python-2.7.11/bin/activate
    python -c "import underworld"

If it succeeds, you can install UWGeodynamics:

.. code:: bash

    git clone https://github.com/underworldcode/UWGeodynamics.git
    pip install UWGeodynamics/
    python -c "import UWGeodynamics
    rm -rf UWGeodynamics

**PBS script minimal example**

.. code:: bash

    #PBS -P q97
    #PBS -q express
    #PBS -l walltime=00:10:00
    #PBS -l mem=1GB
    #PBS -l jobfs=10MB
    #PBS -l ncpus=10
    #PBS -l software=underworld
    #PBS -l wd
    #PBS -N test
    
    module purge
    module load pbs dot gcc/5.2.0 hdf5/1.10.2p petsc/3.8.4 swig/3.0.12 python/2.7.11 openmpi/3.1.2
    source /short/q97/Underworld/python-2.7.11-venv/bin/activate
    
    MODELNAME="test"
    OUTPUTPATH=`pwd`
    SCRIPT="your-script.py"
    
    mpiexec --mca mpi_warn_on_fork 0 --mca opal_abort_print_stack 1 --mca mpi_param_check 1 \
     --mca mpi_add_procs_cutoff 256 python ./$SCRIPT 1> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.log 2> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.err


Running Underworld / UWGeodynamics on Pawsey MAGNUS using Shifter
=================================================================

The recommended way to run Underworld / UWGeodynamics model is to use
Shifter. Shifter is a wrapper around Docker that allows us to run docker
containers on Magnus.

You can have a look at the `Pawsey
documentation <https://support.pawsey.org.au/documentation/display/US/Shifter>`__
if you want to know more about Shifter:

Pre-requisites
--------------

.. code:: bash

   ssh username@magnus-1.pawsey.org.au

**A UWGeodynamics docker image is already available on Magnus**

.. code:: bash

   user@magnus-1:~>module load shifter
   user@magnus-1:~>shifter images
   magnus     docker     READY    17cc3c02ba   2018-05-09T08:47:59 underworldcode/uwgeodynamics:magnus

The following command will pull the latest version of the image:

.. code:: bash

   shifter pull docker:underworldcode/uwgeodynamics:magnus

Setting up a job
----------------

Here we assume that we have a copy of the UWGeodynamics Tutorial 1 model
saved as a python file (*Tutorial_1_ThermoMechanical_Model.py*), inside
a folder *UWGeo_Tutorial1* located in the
/scratch/your-project-account/your-username folder:

.. code:: bash

   rb5533@magnus-1:/scratch/q97/rb5533/UWGeo_Tutorial1> ls 
   Tutorial_1_ThermoMechanical_Model.py

SLURM file
~~~~~~~~~~

Following is an example of a SLURM file (*job.slurm*) showing how to run
Tutorial 1 on 1 node using 4 cores:

.. code:: bash

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

.. code:: bash

   rb5533@magnus-1:/scratch/q97/rb5533/UWGeo_Tutorial1> ls 
   Tutorial_1_ThermoMechanical_Model.py    job.slurm

The job can now be submitted to the queue system using:

.. code:: bash

   sbatch job.slurm

That’s it!!!
