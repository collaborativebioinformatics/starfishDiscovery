# Setup [starfish](https://github.com/egluckthaler/starfish/wiki/Installation) on DNAnexus ttyd using conda

## Get miniconda  

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`  

`bash Miniconda3-latest-Linux-x86_64.sh` 

`export PATH="/home/dnanexus/miniconda3/bin:$PATH"` # this is specific for DNAnexus ttyd

## Set config channels

`conda config --add channels bioconda`  

`conda config --add channels conda-forge`

`conda config --set channel_priority strict`

## Create environment

`conda create -c egluckthaler -n starfish python=3.8 starfish`

## Activate environment

`conda activate starfish`