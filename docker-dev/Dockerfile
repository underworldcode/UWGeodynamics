FROM underworldcode/underworld2_untested:dev
MAINTAINER romain.beucher@unimelb.edu

USER root

RUN apt-get update -qq
RUN DEBIAN_FRONTEND=noninteractive apt-get install cgdb sudo

# UWGeodynamics
WORKDIR /opt
RUN git clone -b development https://github.com/rbeucher/UWGeodynamics.git && \
    pip install -e /opt/UWGeodynamics
RUN mkdir /workspace/UWGeodynamics && \
    mkdir /workspace/MODELS_RESULTS && \
    mv /opt/UWGeodynamics/examples /workspace/UWGeodynamics/. && \
    mv /opt/UWGeodynamics/tutorials /workspace/UWGeodynamics/. && \
    mv /opt/UWGeodynamics/manual /workspace/UWGeodynamics/.

# Badlands dependency
RUN pip install --no-cache-dir pandas \
                ez_setup \ 
                git+https://github.com/badlands-model/triangle.git \
                git+https://github.com/awickert/gFlex.git

# pyBadlands serial
RUN git clone https://github.com/rbeucher/pyBadlands_serial.git
RUN cd /opt/pyBadlands_serial/pyBadlands/libUtils && make
RUN pip install -e /opt/pyBadlands_serial

# pyBadlands dependencies
RUN git clone https://github.com/matplotlib/cmocean.git /tmp/cmocean && \
    pip install /tmp/cmocean/ && \ 
    rm -rf /tmp/cmocean
RUN pip install colorlover matplotlib==2.1.2

# pyBadlands companion
RUN git clone https://github.com/badlands-model/pyBadlands-Companion.git && \
    pip install -e /opt/pyBadlands-Companion && \
    mkdir /workspace/BADLANDS && \
    rsync -av /opt/pyBadlands-Companion/notebooks/* /workspace/BADLANDS/companion/

ENV PATH $PATH:/opt/pyBadlands_serial/pyBadlands/libUtils
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/pyBadlands_serial/pyBadlands/libUtils

# memory profiler
RUN pip install memory_profiler

# update all permissions for jovyan user
ENV UW2_DIR /opt/underworld2
RUN echo "jovyan:docker" | chpasswd && adduser jovyan sudo
ENV NB_USER jovyan

# copy this file over so that no password is required
COPY jupyter_notebook_config.json /home/$NB_USER/.jupyter/jupyter_notebook_config.json

# update all permissions for user
RUN chown -R $NB_USER:users /workspace $UW2_DIR /home/$NB_USER /opt/pyBadlands_serial

# change user and update pythonpath
USER $NB_USER
ENV PYTHONPATH $PYTHONPATH:$UW2_DIR
ENV PYTHONPATH /workspace/user_data/Underworld:/workspace/user_data/UWGeodynamics:$PYTHONPATH

# move back to workspace directory
WORKDIR /workspace

# launch notebook
CMD ["jupyter", "notebook", "--ip='*'", "--no-browser"]

