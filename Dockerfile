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
RUN apt-get install -y ncbi-blast+
RUN apt-get install -y bedtools
RUN apt-get install -y texlive-font-utils


#diamond
RUN mkdir diamond 
RUN cd diamond     
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz 
RUN tar -xzf diamond-linux64.tar.gz 
RUN rm diamond-linux64.tar.gz diamond_manual.pdf 
RUN cd ..


# mfold
RUN wget http://www.unafold.org/download/mfold-3.6.tar.gz
RUN tar -xvf mfold-3.6.tar.gz
#RUN rm mfold-3.6.tar.gz
RUN cd mfold-3.6
RUN ./configure
RUN make
RUN make install
RUN cd ..

# unafold
RUN cd ./software/unafold-3.7
RUN ./configure
RUN make
RUN make install
RUN cd ..

# vieena
RUN wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_22_04/viennarna_2.5.1-1_amd64.deb -O viennarna.deb
RUN dpkg -i ./viennarna.deb
RUN apt-get -f install
RUN rm viennarna.deb

# mxfold
RUN wget https://github.com/keio-bioinformatics/mxfold2/releases/download/v0.1.1/mxfold2-0.1.1.tar.gz
RUN pip3 install mxfold2-0.1.1.tar.gz
RUN rm mxfold2-0.1.1.tar.gz

# contrafold
RUN wget http://contra.stanford.edu/contrafold/contrafold_v2_02.tar.gz
RUN tar -xvzf contrafold_v2_02.tar.gz && rm contrafold_v2_02.tar.gz
RUN cd software/contrafold/src
RUN make clean
RUN make 
# to file must changed to be complieable # utility.hpp and optimization.c++ files


WORKDIR /AmiRML_p/