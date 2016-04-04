#!/bin/bash
cd $(dirname "$0")
CC=gcc
OS=`uname`
NODE=`uname -n`
FLAGS="-std=gnu99 -fopenmp -lmpi"
if [[ $NODE =~ h2ologin ]]; then
  CC=cc
  FLAGS=""
fi
if command -v clang-format >/dev/null 2>&1; then
  echo "Linting..."
  clang-format -i *.c *.h
fi
echo "Compiling..."
if [ $OS == "Darwin" ]; then
  CC=gcc-5
fi
echo "... using $CC."
rm *.out
$CC $FLAGS -c $(find . -name \*.c) -lm
if [ $OS == "Darwin" ]; then
  rm benchmark-native.o
fi
echo "Building..."

$CC $FLAGS $(find . -name \*.o -not -name k2.o) -o test.out -lm
$CC $FLAGS $(find . -name \*.o -not -name test.o) -o k2.out -lm
rm $(find . -name \*.o)
