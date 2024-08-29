FROM continuumio/miniconda3:latest


RUN apt-get update --fix-missing && \
	yes|apt-get upgrade && \
	apt-get install -y \
		wget \
		bzip2 \
		unzip \
		make \
		g++ 

# Install conda environment
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict

RUN conda create -c egluckthaler -n starfish python=3.8 -y

ENV PATH /opt/conda/envs/starfish/bin:$PATH
ENV CONDA_DEFAULT_ENV starfish

RUN /bin/bash -c "source activate starfish" && \
	conda install mummer4=4.0.0 && \
	conda install bedtools=2.30 && \
	conda install sourmash=4.6.1 && \
	conda install blast=2.12.0 && \
	conda install mcl=14.137 && \
	conda install circos=0.69.8 && \
	conda install minimap2=2.24 && \
	conda install hmmer=3.3.2 && \
	conda install metaeuk=6.a5d39d9 && \
	conda install mafft && \
	conda install mmseqs2=14.7e284 && \
	conda install samtools=1.6 && \
	conda install eggnog-mapper=2.1.12 && \
	conda install scikit-image && \
	conda install gxx && \
	rm -rf /opt/conda/pkgs/* && \
	rm -rf /root/.cache/pip/* && \
	echo "source activate starfish" >> ~/.bashrc

RUN cd /opt/ && \
	git clone https://github.com/egluckthaler/starfish.git && \
    cd starfish/CNEFinder && \
    ./pre-install.sh && \
    make -f Makefile

WORKDIR /opt/starfish
COPY examples/run_starfish_example.sh .
ENV PATH /opt/starfish/CNEFinder:$PATH
ENV PATH /opt/starfish/bin:$PATH