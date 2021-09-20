#!/bin/bash

module purge
module load openmpi/4.0.2 hdf5/1.10.5p python3/3.8.5 scons/3.1.1

export GROUP=q97
export USER=
export INSTALL_NAME=UWGeodynamics_2.10.0

export CODES_PATH=/scratch/$GROUP/$USER/codes
export UW_OPT_DIR=$CODES_PATH/opt
export INSTALL_PATH=$CODES_PATH/$INSTALL_NAME

export SWIG_VERSION=3.0.12
export SWIG_PATH=$UW_OPT_DIR/swig-$SWIG_VERSION
export PATH=$SWIG_PATH/bin:$PATH

export OMPI_MCA_io=ompio

export CDIR=$PWD
export LD_PRELOAD=$OPENMPI_ROOT/lib/libmpi_usempif08_GNU.so.40:$OPENMPI_ROOT/lib/libmpi_usempi_ignore_tkr_GNU.so.40:$OPENMPI_ROOT/lib/libmpi_cxx.so.40

export UW_BRANCH=v2.10.0b
export UWGEO_BRANCH=development

install_swig() {
	tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
	cd $tmp_dir
	wget http://prdownloads.sourceforge.net/swig/swig-$SWIG_VERSION.tar.gz
	tar -xvzf swig-$SWIG_VERSION.tar.gz
	cd swig-$SWIG_VERSION
	./configure --prefix=$SWIG_PATH
	make
	make install
	rm $tmp_dir
	cd $CDIR
}


install_petsc(){
	source $INSTALL_PATH/bin/activate
        export PETSC_CONFIGURE_OPTIONS="--with-debugging=0 --prefix=/usr/local \
                --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3' \
                --with-zlib=1                   \
                --with-hdf5=1                   \
                --download-mumps=1              \
                --download-parmetis=1           \
                --download-metis=1              \
                --download-superlu=1            \
                --download-hypre=1              \
                --download-scalapack=1          \
                --download-superlu_dist=1       \
                --useThreads=0                  \
                --download-superlu=1            \
                --with-shared-libraries         \
                --with-cxx-dialect=C++11        \
		--prefix=/scratch/q97/codes/opt/petsc_3.12.3\
                --with-make-np=4"

       CC=mpicc CXX=mpicxx FC=mpif90 pip install petsc==3.12.3 -vvv
}

install_python_dependencies(){
	source $INSTALL_PATH/bin/activate
	pip3 install Cython
	pip3 install mpi4py
	unset HDF5_DIR
        export HDF5_VERSION=1.10.5
	export HDF5_LIBDIR=/apps/hdf5/1.10.5p/lib/ompi3
	export HDF5_INCLUDEDIR=/apps/hdf5/1.10.5p/include
        CC=mpicc HDF5_MPI="ON" pip3 install --no-cache-dir --no-binary=h5py h5py

}

install_underworld(){
	export PETSC_DIR=/scratch/q97/codes/opt/petsc_3.12.3
	source $INSTALL_PATH/bin/activate
	tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
	cd $tmp_dir
        git clone --branch $UW_BRANCH https://github.com/underworldcode/underworld2.git $tmp_dir
        pip3 install .
        rm -rf $tmp_dir	
	cd $CDIR
}

install_uwgeodynamics(){
	source $INSTALL_PATH/bin/activate
	tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
	cd $tmp_dir
    git clone --branch $UWGEO_BRANCH https://github.com/underworldcode/uwgeodynamics.git $tmp_dir
    pip3 install .
    rm -rf $tmp_dir	
	cd $CDIR
}

check_underworld_exists(){
	source $INSTALL_PATH/bin/activate
	return $(python3 -c "import underworld") 
}

check_uwgeodynamics_exists(){
	source $INSTALL_PATH/bin/activate
	return $(python3 -c "import UWGeodynamics") 
}

check_badlands_exists(){
	source $INSTALL_PATH/bin/activate
	return $(python3 -c "import badlands") 
}

install_badlands(){
       source $INSTALL_PATH/bin/activate
       pip3 install badlands
}

install_full_stack(){
    if ! command -v swig 2>/dev/null; then
           install_swig  
    else
           echo "Found swig"
    fi
    
    install_python_dependencies
    
    if ! check_underworld_exists; then
          install_underworld
    fi
    
    if ! check_uwgeodynamics_exists; then
          install_uwgeodynamics
    fi
    
    if ! check_badlands_exists; then
          install_badlands
    fi
}


if [ ! -d "$INSTALL_PATH" ]
then
    echo "Environment not found, creating a new one"
    mkdir $INSTALL_PATH
    python3 --version
    python3 -m venv $INSTALL_PATH
else
    echo "Found Environment"
    source $INSTALL_PATH/bin/activate
fi
