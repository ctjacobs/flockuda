# HDF5
HDF5_INCLUDE=/usr/include/hdf5/serial/
HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial/

# C compiler
CC=nvcc -gencode arch=compute_35,code=sm_35
CCFLAGS=-I. -L. -I$(HDF5_INCLUDE) -L$(HDF5_LIB) -I"/home/christian/Downloads/H5Part-1.6.6/build/include" -L"/home/christian/Downloads/H5Part-1.6.6/build/lib" -lcudart -lcurand -lhdf5 -lH5Part

all: flockuda

flockuda: forces.o prey.o predator.o main.o 
	$(CC) -o flockuda forces.o prey.o predator.o main.o $(CCFLAGS)

main.o: main.cu
	$(CC) -dc main.cu $(CCFLAGS)

forces.o: forces.cu
	$(CC) -dc forces.cu $(CCFLAGS)

predator.o: predator.cu
	$(CC) -dc predator.cu $(CCFLAGS)

prey.o: prey.cu
	$(CC) -dc prey.cu $(CCFLAGS)

clean:
	rm -rf *.o flockuda *.h5part
