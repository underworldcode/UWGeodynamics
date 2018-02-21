#FROM underworldcode/underworld2_untested:v2.5
FROM myuw2_new
MAINTAINER romain.beucher@unimelb.edu

USER root

# required dependency for pyBadlands and friends
RUN apt-get update -qq && \
      DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
      apt-utils \
      libavcodec-dev \
      libavformat-dev \
      libavutil-dev \
      libswscale-dev \
      openmpi-bin \
      libhdf5-dev \
      liblapack-dev \
      llvm \
      libedit-dev \
      gfortran \
      libnetcdf-dev \ 
      libgeos-dev \
      libgeos++-dev \
      wget && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* && \
# dependency python packages
pip install \
        pillow \
        enum34 \
        pyvirtualdisplay \
        scipy \
        Cython==0.20 \
        markupsafe \
        zmq \
        singledispatch \
        backports_abc \
        certifi \
        jsonschema \
        path.py \
        git+https://github.com/badlands-model/triangle \
        pandas \
        plotly \ 
        numba==0.23.1 \
        ez_setup \
        netcdf4 \
        colorlover \
        cmocean \
        scikit-fuzzy \
        pyevtk \
	git+https://github.com/awickert/gFlex.git \ 
	shapely \
	descartes \
        jupyter_contrib_nbextensions \
        tqdm && \
# install jupyter extensions
jupyter contrib nbextension install --system && \
# basemap needs compilation :((
# note geos is install via apt
wget http://downloads.sourceforge.net/project/matplotlib/matplotlib-toolkits/basemap-1.0.7/basemap-1.0.7.tar.gz && \
  tar -zxvf basemap-1.0.7.tar.gz && \
  cd basemap-1.0.7 && \
  python setup.py build && \
  python setup.py install && \
  cd .. && \
  rm -rf basemap-1.0.7.tar.gz && \
  rm -rf basemap-1.0.7


# switch user
ENV NB_USER jovyan
USER $NB_USER
WORKDIR /opt
ENV BLAND_DIR /opt/pyBadlands
ENV BCOMP_DIR /opt/pyBadlands-Companion
ENV BWORK_DIR /opt/pyBadlands-workshop

RUN git clone --depth 1 https://github.com/rbeucher/UWGeodynamics.git
USER root
RUN pip install -e /opt/UWGeodynamics
USER $NB_USER

RUN mkdir /workspace/UWGeodynamics                          && \
rsync -av /opt/UWGeodynamics/examples/* /workspace/UWGeodynamics/examples/   && \
rsync -av /opt/UWGeodynamics/tutorials/* /workspace/UWGeodynamics/tutorials/ && \
rsync -av /opt/UWGeodynamics/manual/* /workspace/UWGeodynamics/manual/ 

# download pyBadland, companion and workshop
WORKDIR /opt
RUN git clone --depth 1 https://github.com/badlands-model/pyBadlands_serial.git $BLAND_DIR    && \
    git clone --depth 1 https://github.com/badlands-model/pyBadlands-Companion.git $BCOMP_DIR && \
    git clone --depth 1 https://github.com/badlands-model/pyBadlands-workshop.git $BWORK_DIR

# compile pyBadlands and companion
WORKDIR $BLAND_DIR/pyBadlands/libUtils
USER root
RUN make                         && \
    pip install -e $BLAND_DIR    && \
    cd /opt 
RUN pip install -e $BCOMP_DIR    && \
# update jupyter extensions while here
jupyter nbextension enable hide_input/main --system && \
jupyter nbextension enable init_cell/main --system
USER $NB_USER
RUN mkdir /workspace/volume      && \
    mkdir /workspace/companion   && \
    mkdir /workspace/workshop    && \
    mkdir /workspace/LavaVu

# Copy test files to workspace. Is this required??? takes 500Mb
RUN cp -av $BWORK_DIR/* /workspace/workshop/        && \
cp -av $BLAND_DIR/Examples/* /workspace/examples/   && \
cp -av $BCOMP_DIR/notebooks/* /workspace/companion/

# NOT SURE ABOUT /workspace/volume???
ENV LD_LIBRARY_PATH /workspace/volume/pyBadlands_serial/pyBadlands/libUtils:/$BLAND_DIR/pyBadlands/libUtils

# change user and update pythonpath
ENV PYTHONPATH $PYTHONPATH:$UW2_DIR

# move back to workspace directory
WORKDIR /workspace

# Trust underworld notebooks
RUN find -name \*.ipynb  -print0 | xargs -0 jupyter trust

# launch notebook
CMD ["jupyter", "notebook", "--ip='*'", "--no-browser"]
