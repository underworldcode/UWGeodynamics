FROM underworldcode/underworld2_untested:dev
MAINTAINER romain.beucher@unimelb.edu

USER root

RUN apt-get update -qq
RUN DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommend cgdb sudo #python-mpltoolkits.basemap

# UWGeodynamics
WORKDIR /opt
RUN git clone -b development https://github.com/rbeucher/UWGeodynamics.git
RUN pip install -e /opt/UWGeodynamics
RUN mkdir /workspace/UWGeodynamics
RUN mkdir /workspace/MODELS_RESULTS
RUN rsync -av /opt/UWGeodynamics/examples/* /workspace/UWGeodynamics/examples/
RUN rsync -av /opt/UWGeodynamics/tutorials/* /workspace/UWGeodynamics/tutorials/
RUN rsync -av /opt/UWGeodynamics/manual/* /workspace/UWGeodynamics/manual/

# Badlands dependency
RUN pip install pandas
RUN pip install ez_setup
RUN pip install git+https://github.com/badlands-model/triangle.git
RUN pip install git+https://github.com/awickert/gFlex.git

# pyBadlands serial
WORKDIR /opt
RUN git clone https://github.com/rbeucher/pyBadlands_serial.git
RUN pip install -e pyBadlands_serial/
WORKDIR /opt/pyBadlands_serial/pyBadlands/libUtils
RUN make
RUN pip install -e /opt/pyBadlands_serial

# pyBadlands dependencies
RUN pip install cmocean
RUN pip install colorlover
# Force matplotlib 2.1.2 (Bug Badlands), Temporary
RUN pip install matplotlib==2.1.2

# pyBadlands companion
WORKDIR /opt
RUN git clone https://github.com/badlands-model/pyBadlands-Companion.git
RUN pip install -e /opt/pyBadlands-Companion
RUN mkdir /workspace/BADLANDS
RUN rsync -av /opt/pyBadlands-Companion/notebooks/* /workspace/BADLANDS/companion/

ENV PATH $PATH:/opt/pyBadlands_serial/pyBadlands/libUtils
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/pyBadlands_serial/pyBadlands/libUtils

# memory profiler and jupyterlab
RUN pip install memory_profiler
RUN pip install jupyterlab
RUN jupyter serverextension enable --py jupyterlab --sys-prefix

# update all permissions for jovyan user
ENV UW2_DIR /opt/underworld2
#RUN useradd -m -s /bin/bash -N jovyan
RUN echo "jovyan:docker" | chpasswd && adduser jovyan sudo
ENV NB_USER jovyan

# copy this file over so that no password is required
COPY jupyter_notebook_config.json /home/$NB_USER/.jupyter/jupyter_notebook_config.json

# update all permissions for user
RUN chown -R $NB_USER:users /workspace $UW2_DIR /home/$NB_USER /opt/pyBadlands_serial

# change user and update pythonpath
USER $NB_USER
ENV PYTHONPATH $PYTHONPATH:$UW2_DIR
ENV PYTHONPATH /workspace/user_data/UWGeodynamics:$PYTHONPATH

# move back to workspace directory
WORKDIR /workspace

# launch notebook
CMD ["jupyter", "lab", "--ip='*'", "--no-browser"]

