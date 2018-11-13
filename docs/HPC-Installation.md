
We provide several compiled versions of Underworld to users of the q97 group 
in the`/short/q97/Underworld`folder.

The name of the installations use the following pattern `underworld_%branch_%date-of-compilation`

The folder also contains 2 symlinks called `underworld_latest` which
point to the latest release build available and `underworld_development` which points to
the latest development build available.

#### Installing Underworld2

If you are not a member of the q97 group of if you need to build your own version of
Underworld, we provide a script to perform the installation automatically using the latest
working libraries.

```bash
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
UW_DIR=$UW_PREFIX/underworld_$1_$today
git clone https://github.com/underworldcode/underworld2.git --depth 1 --branch $1 $UW_DIR
cd $UW_DIR

# setup modules
module purge
RUN_MODS='pbs dot python3/3.6.2 openmpi/3.1.2'
PYTHON_ROOT=/apps/python3/3.6.2

module load gcc/5.2.0 hdf5/1.10.2p petsc/3.8.4 swig/3.0.12 scons $RUN_MODS # can use petsc/3.8.3 as well

export HDF5_VERSION=1.10.2

module list

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
```
To use the script you will have to first define the `UW_PREFIX` environment variable

```bash
export UW_PREFIX="path-to-your-installation"
```
You can then source the script (Assuming it is called update_underworld.sh)

```
source update_underworld.sh
```

Finally run

```
update_underworld branch-name
```

with `branch-name` being any of the Underworld GitHub branch or tag available.


#### Installing UWGeodynamics

The installation on Raijin requires loading the python3 module

We highly recommend using the same python version used to build
your Underworld installation.
The code has been tested for python >=3.5, we recommend using python3/3.6.2

```
      module load python python3/3.6.2
```

Once the python module loaded, you can download the source of the module:

```
    git clone https://github.com/rbeucher/UWGeodynamics.git 
    pip3 install -e UWGeodynamics --user
```

