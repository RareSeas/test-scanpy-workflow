FROM imperialgenomicsfacility/scanpy-notebook-image
LABEL maintainer="RareSeas"
LABEL version="0.0.1"
LABEL description="Docker image for running Scanpy based single cell analysis"
ENV NB_USER vmuser
ENV NB_UID 1000
USER root
WORKDIR /
RUN apt-get -y update &&   \
    apt-get install --no-install-recommends -y \
      libfontconfig1 \
      libxrender1 \
      libreadline-dev \
      libreadline7 \
      libicu-dev \
      libc6-dev \
      icu-devtools \
      libjpeg-dev \
      libxext-dev \
      libcairo2 \
      libicu60 \
      libicu-dev \
      gcc \
      g++ \
      make \
      libgcc-5-dev \
      gfortran \
      git  && \
    apt-get purge -y --auto-remove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
USER $NB_USER
WORKDIR /home/$NB_USER
ENV TMPDIR=/home/$NB_USER/.tmp
ENV PATH=$PATH:/home/$NB_USER/miniconda3/bin/
ENV NUMBA_CACHE_DIR=/tmp
ENV CBDATA=/tmp
RUN rm -f /home/$NB_USER/environment.yml && \
    rm -f /home/$NB_USER/Dockerfile
COPY environment.yml /home/$NB_USER/environment.yml
COPY Dockerfile /home/$NB_USER/Dockerfile
USER root
RUN chown ${NB_UID} /home/$NB_USER/environment.yml && \
    chown ${NB_UID} /home/$NB_USER/Dockerfile && \
USER $NB_USER
WORKDIR /home/$NB_USER
RUN conda update -n base -c defaults conda && \
    conda env update -q -n notebook-env --file /home/$NB_USER/environment.yml && \
    conda clean -a -y && \
    rm -rf /home/$NB_USER/.cache && \
    rm -rf /tmp/* && \
    rm -rf ${TMPDIR} && \
    mkdir -p ${TMPDIR} && \
    mkdir -p /home/$NB_USER/.cache && \
    find miniconda3/ -type f -name *.pyc -exec rm -f {} \; 

# Add Data

ADD data /data

# Add experiment scripts

RUN git clone https://github.com/RareSeas/test-scanpy-workflow.git /scripts
RUN chmod a+rx /scripts/*

ENV PATH=/scripts:$PATH

WORKDIR /experiment

#EXPOSE 8888
#EXPOSE 8080

SHELL ["/bin/bash", "-c"]
ENTRYPOINT ["cat /scripts/testscript.py"]