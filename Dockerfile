FROM underworldcode/underworld2:latest

USER root
ENV NB_USER jovyan
ENV NB_WORK /home/$NB_USER

WORKDIR /opt

RUN git clone -b master https://github.com/underworldcode/UWGeodynamics.git UWGeodynamics&& \
    pip3 install --no-cache-dir UWGeodynamics/ && \
    rm -rf $NB_WORK/UWGeodynamics && \
    mkdir $NB_WORK/UWGeodynamics && \
    mv ./UWGeodynamics/examples $NB_WORK/UWGeodynamics/. && \
    mv ./UWGeodynamics/tutorials $NB_WORK/UWGeodynamics/. && \
    mv ./UWGeodynamics/benchmarks $NB_WORK/UWGeodynamics/. && \
    mv ./UWGeodynamics/docs $NB_WORK/UWGeodynamics/ && \
    rm -rf UWGeodynamics && \
    chown -R $NB_USER:users $NB_WORK/UWGeodynamics

# Badlands dependency
RUN pip3 install --no-cache-dir pandas \
                ez_setup \ 
                tribad \
                git+https://github.com/awickert/gFlex.git \
                git+https://github.com/matplotlib/legacycontour.git \
                git+https://github.com/matplotlib/cmocean.git \
                colorlover matplotlib

# Compile and install pyBadlands
RUN git clone https://github.com/rbeucher/pyBadlands_serial.git pyBadlands &&\
    cd pyBadlands/pyBadlands/libUtils && \
    make && \
    cd /opt && \
    pip3 install --no-cache-dir -e pyBadlands && \
    rm -rf pyBadlands/Examples

USER jovyan
ENV PATH $PATH:/opt/pyBadlands/pyBadlands/libUtils
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/pyBadlands/pyBadlands/libUtils

WORKDIR $NB_WORK
CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--no-browser"]
