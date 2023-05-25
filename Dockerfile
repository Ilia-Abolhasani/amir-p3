FROM ubuntu:22.04
LABEL version="1.0.0"

WORKDIR /amir-p3

RUN apt-get update
RUN apt-get install -y autoconf
RUN apt-get install -y gcc-multilib
RUN apt-get install -y git
RUN apt-get install -y make
RUN apt-get install -y openjdk-8-jdk 
RUN apt-get install -y perl
RUN apt-get install -y python3.10
RUN apt-get install -y sqlite3
RUN apt-get install -y tar
RUN apt-get install -y unzip
RUN apt-get install -y wget 
RUN apt-get install -y ncbi-blast+
RUN apt-get install -y bedtools
RUN apt-get install -y texlive-font-utils
RUN apt-get install -y build-essential
RUN apt-get install -y python3-pip


COPY ./data /amir-p3/data
COPY ./config /amir-p3/config
COPY ./software /amir-p3/software
COPY ./requirements.txt /amir-p3/

WORKDIR /amir-p3/software
RUN find . -type f -wholename '*.tar.gz' -exec sh -c 'tar -xzvf {} && rm {}' \;

# unafold
RUN cd ./unafold-3.7 && ./configure && make && make install

# mfold
RUN cd ./mfold-3.6 && ./configure && make && make install

# contrafold
RUN cd ./contrafold/src && make clean && make

# mxfold
RUN wget https://github.com/keio-bioinformatics/mxfold2/releases/download/v0.1.1/mxfold2-0.1.1.tar.gz
RUN pip3 install mxfold2-0.1.1.tar.gz
RUN rm mxfold2-0.1.1.tar.gz


# vieena
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.1.tar.gz
RUN tar -xzvf ViennaRNA-2.5.1.tar.gz && rm ViennaRNA-2.5.1.tar.gz 
RUN cd ./ViennaRNA-2.5.1 && ./configure && make && make install


COPY ./src /amir-p3/src
RUN pip3 install ./requirements.txt

WORKDIR /amir-p3
