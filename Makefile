# -march=native -mtune=native -Ofast -funroll-loops 
FLAGS = 
DIR = heat-dat-files

# default serial
CC=gcc
TYPE=serial

# Ricorda: carica il modulo [devtoolset-10], nella versione gcc 4.* #pragma omp
# cancel non è supportata
ifeq ($(type), omp) 
	TYPE=omp
	FLAGS += -fopenmp -D_USE_OMP=1
endif
ifeq ($(type), mpi)
	TYPE=mpi
	CC=mpicc
endif

ifeq ($(out), y)
	FLAGS += -D_OUTPUT=1
endif

all:
	export OMP_CANCELLATION=true
	$(CC) -c -g utils.c -lm
	$(CC) -c -g heat2D-$(TYPE).c $(FLAGS) -lm
	$(CC) -g -o heat2D-$(TYPE) utils.o heat2D-$(TYPE).o $(FLAGS) -lm

clear:
	rm ./heat-dat-files/Heat*.dat

clear-dir:
	rm ./out/test*
	rm ./err/test*

# load-module:
#  	module load devtoolset-10
#  	module load cuda/10.2
#  	module load openmpi