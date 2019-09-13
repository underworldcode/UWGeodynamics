Installation
============

Docker_ installation
--------------------

Docker containers provide and easy-way to set up and distribute
applications. They also provide a safe and consistent environment which
facilitate debugging and reproducibility of models. The image we provide
contains all the dependencies and configuration files required to run
Underworld models. Users can start developping model as soon as they
have downloaded the image, independently of the operating system running
on their machine.

We strongly encourage users to run UWGeodynamics using the docker images
we provide on `Docker Hub`_

Different version of the `underworldcode/uwgeodynamics` image can be
pulled using a tag:

1. The *latest* tag points to the github master branch and uses the latest
   *underworld* release.
2. The *dev* tag points to the github development and uses the development
   branch of *underworld*.
3. release tags such as *0.9.8* points to the specified version.

Command line
~~~~~~~~~~~~

Once you have installed docker on your system you can *pull* the
*UWGeodynamics* official image as follow:

.. code:: bash

  docker pull underworldcode/uwgeodynamics

You can list all the images available on your system as follow:

.. code:: bash

  docker images

An image can be deleted as follow:

.. code:: bash

  docker rmi underworldcode/uwgeodynamics

You can then start a docker container. (An instance of
an image).

.. code:: bash

  docker run -d \
     --name my_container \
     -p 8888:8888 \
     --mount source=myvol,target=/workspace/user-data \
     underworldcode/uwgeodynamics

You can access the container via your browser at the following
address: http://localhost:8888


It is also possible to ssh into the container as follow:

.. code:: bash

  docker exec -it my_container /bin/bash

You can list the containers currently existing on your machine by running:

.. code:: bash

  docker ps -a

The "a" means "all container". The :code:`docker ps` command only list
running containers.

Docker containers can be stop (so that they do not use CPU or RAM ressource):

.. code:: bash

  docker stop my_container

They can also be deleted:

.. code:: bash

  docker rm my_container

.. warning::

  It's a good idea to keep track of how many containers have been created as
  they can rapidly take a lot of space on your machine.

Kitematic_
~~~~~~~~~~

Kitematic_ is a program that provides a graphical user interface to
the *docker* daemon and to Docker Hub.
The software is available for Windows, MacOsx and Linux. Be aware that on
linux the installation may differ depending on the distribution you
are running.

1. Download and Install Kitematic_
2. Open Kitematic and search for the **uwgeodynamics** image.
3. Create a container by clicking on the create button.

You should now have a container appearing on the left side of your
kitematic window. The first thing to do now is to create a link between
a local directory (A directory on your physical hard drive) and a volume
directory inside the docker container. A volume is a special directory
that can be accessed from outside the container. It is the location you
will use to save your results.

Local Installation
------------------

This is not recommended and involves installing *Underworld* and all
its dependencies. Docker is highly recommended!!!

Requirements
~~~~~~~~~~~~

-  Python >= 2.7
-  A Working version of Underworld2 >=2.6.0 (Please refer to the
   Underworld documentation)
-  pint >= 0.8

**Note on Python 3 compatibility**:
The bleeding edge version of *Underworld* (development branch)
is now python 3 compatible only.
*UWGeodynamics* is python 3 ready and can thus be used with it.

Install
~~~~~~~

**from Pip**

The UWGeodynamics module can be installed directly from the Python
Package Index:

.. code:: bash

  pip install UWGeodynamics

**from sources**

The module source files are available through github_

.. code:: bash

  git clone https://github.com/underworldcode/UWGeodynamics.git

It can then be installed globally on your system using

.. code:: bash

  pip install UWGeodynamics/


HPC Installation
----------------

requirements
^^^^^^^^^^^^

- python >= 3.6
- openmpi >= 1.8
- hdf5 >= 1.8

NCI Raijin
~~~~~~~~~~

The following script install Underworld and UWGeodynamics on raijin:

.. code:: bash

  #!/bin/sh
  # This script installs underworld on raijin.nci.org.au
  # 
  # Usage:
  #  sh ./nci_raijin.sh <destfolder>
  #
  # exit when any command fails
  set -e
  
  DATE=`date +%d%b%Y` # could be used to date checkout eg,
  INSTALLPATH=`pwd`/$1
  UW_DIR=$INSTALLPATH/underworld-$DATE
  
  mkdir $INSTALLPATH
  
  cd $INSTALLPATH
  git clone https://github.com/underworldcode/underworld2.git $UW_DIR
  cd $UW_DIR
  git checkout master  # checkout the requested version
  
  # setup modules
  module purge
  RUN_MODS='pbs dot mpi4py/3.0.2-py36-ompi3'
  module load hdf5/1.10.2p petsc/3.9.4 gcc/5.2.0 mesa/11.2.2 swig/3.0.12 scons/3.0.1 $RUN_MODS
  echo "*** The module list is: ***"
  module list -t
  
  # setup python environment with preinstalled packages (h5py, lavavu, pint)
  export PYTHONPATH=/apps/underworld/opt/h5py/2.9.0-py36-ompi3/lib/python3.6/site-packages/h5py-2.9.0-py3.6-linux-x86_64.egg/:/apps/underworld/opt/lavavu/1.4.1_rc/:/apps/underworld/opt/pint/0.9_py36/lib/python3.6/site-packages/:$PYTHONPATH
  echo "*** New PYTHONPATH: $PYTHONPATH ***"
  
  # build and install code
  cd libUnderworld
  CONFIG="./configure.py  --python-dir=`python3-config --prefix` --with-debugging=0"
  echo "*** The config line is: ***"
  echo "$CONFIG"
  echo ""
  
  $CONFIG
  ./compile.py -j4
  cd .. ; source updatePyPath.sh 
  cd $INSTALLPATH
  
  # UWGeodynamics
  
  pip3 install git+https://github.com/underworldcode/uwgeodynamics --prefix=$INSTALLPATH
  
  cd $INSTALLPATH
  touch module_paths.sh
  
  echo "#!/bin/bash" >> module_paths.sh
  echo "source $UW_DIR/updatePyPath.sh" >> module_paths.sh
  echo "module purge" >> module_paths.sh
  echo "module load $RUN_MODS" >> module_paths.sh
  echo "" >> module_paths.sh
  echo "export PYTHONPATH=$UW_DIR:$UW_DIR/glucifer:$PYTHONPATH" >> module_paths.sh
  echo "" >> module_paths.sh
  echo "export PATH=$PATH" >> module_paths.sh
  echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> module_paths.sh
  
  echo "#####################################################################"
  echo "Underworld2 built successfully at:                                   "
  echo "  $UW_DIR                                                            "
  echo "#####################################################################"


We provide a minimal PBS script:

.. code:: bash

  #PBS -P project
  #PBS -q normal
  #PBS -l walltime=1:00:00
  #PBS -l mem=1GB
  #PBS -l jobfs=10MB
  #PBS -l ncpus=6
  #PBS -l software=underworld
  #PBS -l wd
  #PBS -N name
  
  source path-to-your-underworld-installation/module_paths.sh
  
  MODELNAME="Model"
  OUTPUTPATH=`pwd`
  SCRIPT="Model.py"
  
  mpiexec --mca mpi_warn_on_fork 0 --mca opal_abort_print_stack 1 --mca mpi_param_check 1 \
   --mca mpi_add_procs_cutoff 256 python ./$SCRIPT 1> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.log 2> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.err
  

Pawsey MAGNUS
-------------

The recommended way to run Underworld / UWGeodynamics model is to use
Shifter. Shifter is a wrapper around Docker that allows us to run docker
containers on Magnus.

You can have a look at the `Pawsey
documentation <https://support.pawsey.org.au/documentation/display/US/Shifter>`__
if you want to know more about Shifter:

Pre-requisites
~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~

Here we assume that we have a copy of the UWGeodynamics Tutorial 1 model
saved as a python file (*Tutorial_1_ThermoMechanical_Model.py*), inside
a folder *UWGeo_Tutorial1* located in the
/scratch/your-project-account/your-username folder:

.. code:: bash

   rb5533@magnus-1:/scratch/q97/rb5533/UWGeo_Tutorial1> ls 
   Tutorial_1_ThermoMechanical_Model.py

SLURM file
^^^^^^^^^^

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

   srun -n4 shifter run --mpi underworldcode/uwgeodynamics:magnus python Tutorial_1_ThermoMechanical_Model.py 

Running a job
~~~~~~~~~~~~~

After the above we should have the following files in our
*UWGeo_Tutorial1* folder:

.. code:: bash

   rb5533@magnus-1:/scratch/q97/rb5533/UWGeo_Tutorial1> ls 
   Tutorial_1_ThermoMechanical_Model.py    job.slurm

The job can now be submitted to the queue system using:

.. code:: bash

   sbatch job.slurm

That’s it!!!

.. _Jupyter: http://jupyter.org/
.. _Docker: https://www.docker.com
.. _Docker Hub: https://hub.docker.com/r/underworldcode/uwgeodynamics
.. _Kitematic: https://kitematic.com/
.. _github: https://github.com/underworldcode/UWGeodynamics.git
.. _Pint: https://pint.readthedocs.io/en/latest
