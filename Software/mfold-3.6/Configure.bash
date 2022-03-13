#!/bin/bash

PREFIX=${1:-}
if [ "$PREFIX" ] ; then
    PREFIX="--prefix=$PREFIX"
fi

./configure $PREFIX CC=cc CXX=c++
