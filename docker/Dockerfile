FROM underworldcode/underworld2_untested:v2.5
MAINTAINER romain.beucher@unimelb.edu

USER root

WORKDIR /opt
RUN git clone https://github.com/rbeucher/UWGeodynamics.git
RUN pip install -e /opt/UWGeodynamics
RUN mkdir /workspace/UWGeodynamics
RUN rsync -av /opt/UWGeodynamics/examples/* /workspace/UWGeodynamics/examples/
RUN rsync -av /opt/UWGeodynamics/tutorials/* /workspace/UWGeodynamics/tutorials/
RUN rsync -av /opt/UWGeodynamics/manual/* /workspace/UWGeodynamics/manual/

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
rm -rf /var/lib/apt/lists/*

# dependency python packages
RUN pip install \
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
        tqdm

# install jupyter extensions
RUN jupyter contrib nbextension install --system

# basemap needs compilation :((
# note geos is install via apt
RUN wget http://downloads.sourceforge.net/project/matplotlib/matplotlib-toolkits/basemap-1.0.7/basemap-1.0.7.tar.gz && \
  tar -zxvf basemap-1.0.7.tar.gz && \
  cd basemap-1.0.7 && \
  python setup.py build && \
  python setup.py install && \
  cd .. && \
  rm -rf basemap-1.0.7.tar.gz && \
  rm -rf basemap-1.0.7

# download pyBadland, companion and workshop
WORKDIR /opt
ENV BLAND_DIR /opt/pyBadlands
ENV BCOMP_DIR /opt/pyBadlands-Companion
ENV BWORK_DIR /opt/pyBadlands-workshop
RUN git clone https://github.com/badlands-model/pyBadlands_serial.git $BLAND_DIR && \
    git clone https://github.com/badlands-model/pyBadlands-Companion.git $BCOMP_DIR && \
    git clone https://github.com/badlands-model/pyBadlands-workshop.git $BWORK_DIR

# compile pyBadlands and companion
WORKDIR $BLAND_DIR/pyBadlands/libUtils
RUN make && \
    pip install -e $BLAND_DIR && \
    pip install -e $BCOMP_DIR && \
    mkdir /workspace/volume && \
    mkdir /workspace/companion && \
    mkdir /workspace/workshop && \
    mkdir /workspace/LavaVu

# Copy test files to workspace
RUN cp -av $BWORK_DIR/* /workspace/workshop/
RUN cp -av $BLAND_DIR/Examples/* /workspace/examples/
RUN cp -av $BCOMP_DIR/notebooks/* /workspace/companion/

# NOT SURE ABOUT /workspace/volume
ENV LD_LIBRARY_PATH=/workspace/volume/pyBadlands_serial/pyBadlands/libUtils:/$BLAND_DIR/pyBadlands/libUtils

# update all permissions for jovyan user
ENV NB_USER jovyan
RUN chown -R $NB_USER /workspace $UW2_DIR /home/$NB_USER

# change user and update pythonpath
USER $NB_USER
ENV PYTHONPATH $PYTHONPATH:$UW2_DIR

# update jupyter extensions 
RUN jupyter nbextension enable hide_input/main && \
    jupyter nbextension enable init_cell/main

# move back to workspace directory
WORKDIR /workspace

# Trust underworld notebooks
RUN find -name \*.ipynb  -print0 | xargs -0 jupyter trust

# launch notebook
CMD ["jupyter", "notebook", "--ip='*'", "--no-browser"]
