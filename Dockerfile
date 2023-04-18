FROM ubuntu:22.04
LABEL version="1.0.0"

WORKDIR /AmiRML_p/

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

#diamond
RUN mkdir diamond 
RUN cd diamond     
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz 
RUN tar -xzf diamond-linux64.tar.gz 
RUN rm diamond-linux64.tar.gz diamond_manual.pdf 
RUN cd ..


RUN apt-get install ncbi-blast+
RUN apt-get install bedtools

# mfold
RUN wget http://www.unafold.org/download/mfold-3.6.tar.gz
RUN tar -xvf mfold-3.6.tar.gz
#RUN rm mfold-3.6.tar.gz
RUN cd mfold-3.6
RUN ./configure
RUN make
RUN make install
RUN cd ..
RUN sudo apt install texlive-font-utils


# vieena
#!wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_20_04/viennarna_2.4.18-1_amd64.deb -O viennarna.deb
#!sudo dpkg -i ./viennarna.deb
#!sudo apt-get -f install
#!rm viennarna.deb

# mxfold
#!wget https://github.com/keio-bioinformatics/mxfold2/releases/download/v0.1.1/mxfold2-0.1.1.tar.gz
#!pip3 install mxfold2-0.1.1.tar.gz
#!rm mxfold2-0.1.1.tar.gz

# contrafold
#!wget http://contra.stanford.edu/contrafold/contrafold_v2_02.tar.gz
#!tar -xvzf contrafold_v2_02.tar.gz && rm contrafold_v2_02.tar.gz
#%cd contrafold/src
#!make clean
#!make 
# to file must changed to be complieable # utility.hpp and optimization.c++ files


WORKDIR /AmiRML_p/