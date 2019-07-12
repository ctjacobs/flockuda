# Flockuda: A numerical model of predator-prey dynamics based on the Molecular Dynamics approach of Lee et al. (2006).

#  Copyright (C) 2019 Christian T. Jacobs

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# HDF5
HDF5_INCLUDE=/usr/include/hdf5/serial/
HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial/

# C compiler
CC=nvcc -arch=sm_35
CCFLAGS=-I. -L. -I$(HDF5_INCLUDE) -L$(HDF5_LIB) -I"/home/christian/Downloads/H5Part-1.6.6/build/include" -L"/home/christian/Downloads/H5Part-1.6.6/build/lib" -lcudart -lcurand -lhdf5 -lH5Part -DDEBUG=0

all: flockuda

flockuda: configuration.o forces.o prey.o predator.o main.o 
	$(CC) -o flockuda configuration.o forces.o prey.o predator.o main.o $(CCFLAGS)

main.o: main.cu
	$(CC) -dc main.cu $(CCFLAGS)

configuration.o: configuration.cu
	$(CC) -dc configuration.cu $(CCFLAGS)

forces.o: forces.cu
	$(CC) -dc forces.cu $(CCFLAGS)

predator.o: predator.cu
	$(CC) -dc predator.cu $(CCFLAGS)

prey.o: prey.cu
	$(CC) -dc prey.cu $(CCFLAGS)

clean:
	rm -rf *.o flockuda *.h5part
