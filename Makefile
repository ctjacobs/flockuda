# HDF5
HDF5_INCLUDE=/usr/include/hdf5/serial
HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial/

# C compiler
CC=nvcc -gencode arch=compute_35,code=sm_35
CCFLAGS=-I. -I$(HDF5_INCLUDE) -L$(HDF5_LIB)

all: flockuda

flockuda: main.o prey.o
	$(CC) -o flockuda main.o prey.o -lcudart -lcurand

main.o: main.cu
	$(CC) -dc main.cu $(CCFLAGS)

prey.o: prey.cu
	$(CC) -dc prey.cu $(CCFLAGS)

clean:
	rm -rf *.o flockuda
