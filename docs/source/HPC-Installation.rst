HPC Installation of Underworld / UWGeodynamics
==============================================

We provide several compiled versions of Underworld to users of the q97
group in the `/short/q97/Underworld` folder.

The name of the installations use the following pattern
`underworld_%branch_%date-of-compilation`

The folder also contains 2 symlinks called `underworld_latest` which
point to the latest release build available and
`underworld_development` which points to the latest development build
available.

Installing Underworld2
----------------------

If you are not a member of the q97 group of if you need to build your
own version of Underworld, we provide a script to perform the
installation automatically using the latest working libraries.


Python 2:
~~~~~~~~~

.. code:: bash
   
   #!/bin/bash
   # update_underworld.sh script

   # How to use this script
   # export UW_PREFIX="your-underworld-installation-directory"
   # source update_underworld.sh
   # update_underworld name-of-branch

   function update_underworld()
   {
   set -e
   today=$(date +"%m_%d_%Y")
   UW_DIR=/short/q97/underworld_master_$today
   git clone https://github.com/underworldcode/underworld2.git $UW_DIR
   cd $UW_DIR
   git checkout $1
   
   # setup modules
   module purge
   RUN_MODS='pbs dot python/2.7.11 python/2.7.11-matplotlib openmpi/1.10.2 mpi4py/2.0.0'
   module load gcc/5.2.0 hdf5/1.8.14p petsc/3.7.4 swig/3.0.12 $RUN_MODS # can use petsc/3.8.3 as well
   
   # setup environment
   H5PY_DIR=/apps/underworld/opt/h5py/2.7.1-python_2.7/lib/python2.7/site-packages/
   export PYTHONPATH=$H5PY_DIR:$PYTHONPATH
   export OPENGL_LIB=/apps/underworld/opt/mesa/17.1.5/build/linux-x86_64/gallium/targets/libgl-xlib/
   export OPENGL_INC=/apps/underworld/opt/mesa/17.1.5/include/GL
   
   # build and install code
   cd libUnderworld
   ./configure.py --python-dir=$PYTHON_ROOT --with-debugging=0  #--opengl-inc-dir=$OPENGL_INC --opengl-lib-dir=$OPENGL_LIB
   ./compile.py -j4
   }


Python 3:
~~~~~~~~~

The bleeding edge version of Underworld is only compatible with python 3
Here is a function to build a working python 3 version:


.. code:: bash

   #!/bin/bash
   # update_underworld.sh script

   # How to use this script
   # export UW_PREFIX="your-underworld-installation-directory"
   # source update_underworld.sh
   # update_underworld name-of-branch

   function update_underworld()
   {
   set -e
   today=$(date +"%m_%d_%Y")
   UW_DIR=$UW_PREFIX/underworld_development_$today
   git clone https://github.com/underworldcode/underworld2.git --depth 1 --branch development $UW_DIR
   cd $UW_DIR

   # setup modules
   module purge
   RUN_MODS='pbs dot python3/3.6.2 openmpi/3.1.2'
   PYTHON_ROOT=/apps/python3/3.6.2

   module load gcc/5.2.0 hdf5/1.10.2p petsc/3.8.4 swig/3.0.12 scons $RUN_MODS # can use petsc/3.8.3 as well

   export HDF5_VERSION=1.10.2

   # setup environment
   H5PY_DIR=/short/q97/Underworld/opt/lib/python3.6/site-packages
   MPI4PY_DIR=/short/q97/Underworld/opt/lib/python3.6/site-packages

   export PYTHONPATH=$MPI4PY_DIR:$H5PY_DIR:$PYTHONPATH
   export OPENGL_LIB=/apps/underworld/opt/mesa/17.1.5/build/linux-x86_64/gallium/targets/libgl-xlib/
   export OPENGL_INC=/apps/underworld/opt/mesa/17.1.5/include/GL

   # build and install code
   cd libUnderworld

   python3 configure.py --python-dir=/apps/python3/3.6.2 --with-debugging=0  #--opengl-inc-dir=$OPENGL_INC --opengl-lib-dir=$OPENGL_LIB
   python3 compile.py -j4

   }

To use the script you will have to first define the ``UW_PREFIX``
environment variable

.. code:: bash

   export UW_PREFIX="path-to-your-installation"

You can then source the script (Assuming it is called
update_underworld.sh)

.. code:: bash

   source update_underworld.sh

Finally run

.. code:: bash

   update_underworld branch-name

with `branch-name` being any of the Underworld GitHub branch or tag
available.

A PBS script example
~~~~~~~~~~~~~~~~~~~~

python 2
^^^^^^^^

.. code:: bash

  #PBS -P q97
  #PBS -q normal
  #PBS -l walltime=24:00:00
  #PBS -l mem=800GB
  #PBS -l jobfs=10GB
  #PBS -l ncpus=256
  #PBS -l software=underworld
  #PBS -l wd
  #PBS -N 3D_ExtQ_HR
  #PBS -M romain.beucher@unimelb.edu.au
  #PBS -m abe
  
  module purge
  module load pbs dot python/2.7.11 python/2.7.11-matplotlib openmpi/1.10.2 mpi4py/2.0.0 gcc/5.2.0
  
  export PYTHONPATH=/apps/underworld/opt/h5py/2.7.1-python_2.7/lib/python2.7/site-packages/:$HOME/programs/underworld2_development:$HOME/programs/underworld2_development/glucifer:/apps/underworld/opt/h5py/2.7.1-python_2.7/lib/python2.7/site-packages/:/apps/mpi4py/2.0.0/lib/python2.7/site-packages/
  export PYTHONPATH=$PYTHONPATH:$HOME/opt/UWGeodynamics
  
  #MODELNAME=$(git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/')
  MODELNAME="3D_EQHR"
  OUTPUTPATH=`pwd`
  SCRIPT="3D_Quick-restart.py"
  
  mpiexec --mca mpi_warn_on_fork 0 --mca opal_abort_print_stack 1 --mca mpi_param_check 1 --mca mpi_add_procs_cutoff 256 python ./$SCRIPT 1> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.log 2> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.err


python 3
^^^^^^^^

.. code:: bash

  #PBS -P q97
  #PBS -q normal
  #PBS -l walltime=24:00:00
  #PBS -l mem=500GB
  #PBS -l jobfs=10GB
  #PBS -l ncpus=128
  #PBS -l software=underworld
  #PBS -l wd
  #PBS -N 3D_ExtQ_HR
  #PBS -M romain.beucher@unimelb.edu.au
  #PBS -m abe
  
  module purge
  module load pbs dot python3/3.6.2 openmpi/3.1.2 gcc/5.2.0
  
  export PYTHONPATH=/short/q97/Underworld/opt/lib/python3.6/site-packages:$HOME/programs/underworld_development:$HOME/programs/underworld_development/glucifer
  export PYTHONPATH=$PYTHONPATH:$HOME/opt/UWGeodynamics
  
  #MODELNAME=$(git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/')
  MODELNAME="3D_Rot"
  OUTPUTPATH=`pwd`
  SCRIPT="3D_Rift-rotational.py"
  
  mpiexec --mca mpi_warn_on_fork 0 --mca opal_abort_print_stack 1 
  --mca mpi_param_check 1 --mca mpi_add_procs_cutoff 256 python3 
  ./$SCRIPT 1> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.log 2> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.err




Installing UWGeodynamics
------------------------

The installation on Raijin requires loading the python3 module

We highly recommend using the same python version used to build your
Underworld installation. The code has been tested for python >=3.5, we
recommend using python3/3.6.2

python 2
~~~~~~~~
.. code:: bash

  module load python python/2.7.11

Once the python module loaded, you can download the source of the
module:

.. code:: bash

  git clone https://github.com/underworldcode/UWGeodynamics.git 
  pip install -e UWGeodynamics --user


python 3
~~~~~~~~

.. code:: bash

  module load python python3/3.6.2

Once the python module loaded, you can download the source of the
module:

.. code:: bash

  git clone https://github.com/underworldcode/UWGeodynamics.git 
  pip3 install -e UWGeodynamics --user

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
