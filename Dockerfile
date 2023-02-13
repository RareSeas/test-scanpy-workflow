FROM python:slim
RUN apt-get update && apt-get -y upgrade \
  && apt-get install -y --no-install-recommends \
    git \
    wget \
    g++ \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && echo "Running $(conda --version)" && \
    conda init bash && \
    . /root/.bashrc && \
    conda update conda && \
    conda create -n python-app && \
    conda activate python-app && \
    conda install python=3.8 pip && \
    echo 'print("Hello World!")' > python-app.py && \
    pip install anndata==0.7.6  && \
    pip install fa2==0.3.5  && \  
    pip install scanpy==1.9.1  && \
    pip install pandas==1.3.1  && \
    pip install scipy==1.6.3  && \
    pip install scrublet==0.2.3  && \
    pip install scikit-learn==0.24.2  && \
    pip install phate==1.0.7  && \
    pip install PhenoGraph==1.5.7  && \
    pip install leidenalg==0.8.4  && \
    pip install plotly==4.14.3  && \
    pip install iplotter==0.4.3  && \
    pip install cellbrowser==1.0.1  && \
    pip install scvelo==0.2.3  && \
    pip install pybind11==2.6.2  && \
    pip install velocyto==0.17.17  && \
    pip install umap-learn==0.4.6  && \
    pip install matplotlib==3.4.1  && \
    pip install matplotlib-inline==0.1.2  && \
    pip install ipykernel==5.5.5  && \
    pip install nbconvert==6.0.7  && \
    pip install scrublet==0.2.3  && \
    pip install pybiomart==0.2.0 && \
    conda install -c conda-forge h5py=3.2.1 && \
    conda install -c conda-forge llvmlite=0.34.0 && \
    conda install -c conda-forge numba=0.51.2 && \
    conda install -c conda-forge networkx=2.5.1 && \
    conda install -c conda-forge joblib=1.0.1 && \
    conda install -c conda-forge python-igraph=0.9.1 && \
    conda install -c conda-forge louvain=0.7.0 && \
    conda install -c conda-forge cython=0.29.23 && \
    conda install -c conda-forge click=7.1.2 && \
    conda install -c bioconda samtools=1.12


#RUN echo 'conda activate python-app \n\
#alias python-app="python python-app.py"' >> /root/.bashrc
#ENTRYPOINT [ "/bin/bash", "-l", "-c" ]
#CMD ["python python-app.py"]


ADD scanpydemo.py /root/miniconda3/bin
RUN chmod a+rx /root/miniconda3/bin/scanpydemo.py


RUN echo 'conda activate python-app \n\
alias scanpydemo="python scanpydemo.py"' >> /root/.bashrc
ENTRYPOINT [ "/bin/bash", "-l", "-c" ]
CMD ["python scanpydemo.py"]
