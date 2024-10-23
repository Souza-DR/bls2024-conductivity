#!/bin/bash

# gfortran -c -O4 drand.f90
# gfortran -c -O4 mergesort.f90
# gfortran -c -O4 geometry.f90
# gfortran -c -O4 geompack2.f90
# gfortran -c -O4 vorintpols.f90
# gfortran -c -Wall -O4 voro.f90
# gfortran -c -Wall -O4 voroma.f90

# gfortran mergesort.o drand.o geometry.o geompack2.o vorintpols.o voro.o voroma.o -o voroma

gfortran -fPIC -shared -O4 vorintpols.f90
gfortran -fPIC -shared -O4 -c voro.f90
gfortran -fPIC -shared -O4 -o libvoro.so mergesort.f90 geometry.f90 geompack2.f90 vorintpols.f90 voro.f90