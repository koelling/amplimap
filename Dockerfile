FROM continuumio/miniconda3

ADD https://raw.githubusercontent.com/koelling/amplimap/master/environment.yml environment.yml

# install GCC:
# for some reason this does won't install properly through conda
# it's required to build python packages from pip
# then download and install conda environment
RUN apt-get update --fix-missing && \
    apt-get install -y gcc && \
    apt-get clean && \
    conda env create --file environment.yml && \
    conda clean -i -l -t -y && \
    apt-get remove -y --purge gcc && \
    apt-get autoremove -y && \
    rm -rf /root/.cache/pip/*;

# set up
RUN echo "source activate amplimap" > ~/.bashrc
ENV PATH /opt/conda/envs/amplimap/bin:$PATH

# TO BUILD:
# cd docker; rm Dockerfile; cp ../Dockerfile .; docker build --no-cache -t amplimap .; docker run amplimap amplimap --version;
# docker tag amplimap koelling/amplimap:latest; docker push koelling/amplimap:latest