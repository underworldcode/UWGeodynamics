# Stage 1: Inherit from underworldcode/underworld2:v29 and install dependency packages for Badlands
##########
FROM underworldcode/underworld2:v29 as base_runtime
MAINTAINER https://github.com/underworldcode/
# install runtime requirements
USER root
RUN apt-get update -qq \
&&  DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
        libxml2 \
        libpython3.7
RUN PYTHONPATH= /usr/bin/pip3 install --no-cache-dir setuptools scons 
# setup further virtualenv to avoid double copying back previous packages (h5py,mpi4py,etc)
RUN /usr/bin/python3 -m virtualenv --python=/usr/bin/python3 ${VIRTUAL_ENV}

# Stage 2: Build and install Badlands
##########
FROM base_runtime AS build_base
# install build requirements
RUN apt-get update -qq 
RUN DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
        build-essential \
        gfortran \
        python3-dev \
        swig \
        libxml2-dev
RUN PYTHONPATH= /usr/bin/pip3 install --no-cache-dir setuptools scons 
# setup further virtualenv to avoid double copying back previous packages (h5py,mpi4py,etc)
RUN /usr/bin/python3 -m virtualenv --python=/usr/bin/python3 ${VIRTUAL_ENV}
# Compile and install Badlands
RUN pip3 install badlands


# Stage 3: Resultant images
##########
FROM base_runtime
COPY --from=build_base ${VIRTUAL_ENV} ${VIRTUAL_ENV}
# Record Python packages, but only record system packages! 
# Not venv packages, which will be copied directly in.
RUN PYTHONPATH= /usr/bin/pip3 freeze >/opt/requirements.txt
# Record manually install apt packages.
RUN apt-mark showmanual >/opt/installed.txt
USER $NB_USER
WORKDIR $NB_WORK
CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--no-browser"]
