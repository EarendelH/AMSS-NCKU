
## filein  = -I/usr/include -I/usr/lib/x86_64-linux-gnu/mpich/include -I/usr/lib/x86_64-linux-gnu/openmpi/lib/ -I/usr/lib/gcc/x86_64-linux-gnu/13/ -I/usr/include/c++/13/

## CPU-only (no CUDA), e.g. for AMD W7900D without GPU build
## filein: from mpic++ --showme:compile + system C++/gfortran
filein  = -I/usr/include/ -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/gcc/x86_64-linux-gnu/13/ -I/usr/include/c++/13/

## LDLIBS: from mpic++ --showme:link (-lmpi_cxx -lmpi) + gfortran
LDLIBS  = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib64 -L/usr/lib/gcc/x86_64-linux-gnu/13 -lmpi_cxx -lmpi -lgfortran

CXXAPPFLAGS  = -O3 -Wno-deprecated -Dfortran3 -Dnewc
#f90appflags = -O3 -fpp
f90appflags  = -O3 -x f95-cpp-input
f90          = gfortran 
f77          = gfortran 
CXX          = g++
CC           = gcc
CLINKER      = mpic++ 

## GPU (CUDA) not used for CPU-only build; uncomment and set paths if building ABEGPU
# Cu = nvcc
# CUDA_LIB_PATH = -L/usr/lib/cuda/lib64 -I/usr/include -I/usr/lib/cuda/include
# CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -Dfortran3 -Dnewc
