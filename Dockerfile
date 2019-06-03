FROM underworldcode/underworld2:dev

USER root

# Badlands dependency
RUN pip3 install --no-cache-dir pandas \
                ez_setup \ 
                git+https://github.com/badlands-model/triangle.git \
                git+https://github.com/awickert/gFlex.git \
                git+https://github.com/matplotlib/legacycontour.git \
                git+https://github.com/matplotlib/cmocean.git \
                colorlover matplotlib

WORKDIR /opt

# Compile and install pyBadlands
RUN git clone https://github.com/badlands-model/pyBadlands_serial.git pyBadlands &&\
    cd pyBadlands/pyBadlands/libUtils && \
    make && \
    cd /opt && \
    pip3 install --no-cache-dir -e pyBadlands && \
    rm -rf pyBadlands/Examples

USER jovyan
ENV PATH $PATH:/opt/pyBadlands/pyBadlands/libUtils
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/pyBadlands/pyBadlands/libUtils

CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--no-browser"]
