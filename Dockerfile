FROM underworldcode/underworld2:latest

USER root
ENV NB_USER jovyan
ENV NB_WORK /home/$NB_USER

WORKDIR /opt

# UWGeodynamics
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

#Badlands
RUN pip3 install badlands

USER jovyan

WORKDIR $NB_WORK
CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--no-browser"]
