FROM continuumio/miniconda3

# install GCC:
# for some reason this does won't install properly through conda
# it's required to build python packages from pip
RUN apt-get update --fix-missing && \
    apt-get install -y gcc vim && \
    apt-get clean;

# download and install conda environment
ADD https://raw.githubusercontent.com/koelling/amplimap/master/environment.yml environment.yml
RUN conda env create --file environment.yml

# set up
RUN echo "source activate amplimap" > ~/.bashrc
ENV PATH /opt/conda/envs/amplimap/bin:$PATH