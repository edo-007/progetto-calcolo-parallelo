
FLAGS = -lm 
DIR = heat-dat-files

ifeq ($(out), y)
	FLAGS += -D_OUTPUT=1
endif

ifeq ($(parflags), y)
	FLAGS += -fopenmp -D_USE_OMP=1
endif


serial:
	gcc -c -g utils.c
	gcc -c -g heat2D.c $(FLAGS)
	gcc -g -o heat2D utils.o heat2D.o $(FLAGS)

parallel:
	gcc -c -g utils.c $(FLAGS)
	gcc -c -g heat2D-parallel.c $(FLAGS)
	gcc -g -o heat2D-parallel utils.o heat2D-parallel.o $(FLAGS)

clear:
	rm ./heat-dat-files/Heat*.dat;\