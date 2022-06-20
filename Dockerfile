FROM ubuntu:21.04
MAINTAINER Ilia Abolhasani <Abolhasani.eliya@gmail.com>
LABEL version="1.0.0"

WORKDIR /CTAnalizer/

RUN apt-get update && apt-get install -y \    
    autoconf \        
    gcc-multilib \
    git \
    make \
    openjdk-8-jdk \
    perl \
    python \
    sqlite3 \
    tar \
    unzip \
    wget 

#diamond
RUN mkdir diamond && \
   cd diamond && \
    wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz && \
    tar -xzf diamond-linux64.tar.gz && \
    rm diamond-linux64.tar.gz diamond_manual.pdf && \
    cd /NGStools
    
WORKDIR /CTAnalizer/