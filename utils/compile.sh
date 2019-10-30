#/bin/bash

g++ -I/usr/include/eigen3 -O3 -mtune=native -fopenmp -o PE PE.cpp
