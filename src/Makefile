# Makefile

CC = gcc
OUTPUT = -o $@
#OPTION = -llapacke -lcblas -llapack -lblas -lm 
OPTION = -llapacke -lcblas -llapack -lblas -lm -lgfortran
OBJ_CC = $(CC) $^ -c

all: simulation

simulation: simulation.o mur.o nt.o proposed.o utils.o 
	$(CC) $(OUTPUT) $^ $(OPTION)

simulation.o: simulation.c
	$(OBJ_CC)

mur.o: mur.c
	$(OBJ_CC)

nt.o: nt.c
	$(OBJ_CC)

proposed.o: proposed.c
	$(OBJ_CC)

utils.o: utils.c
	$(OBJ_CC) 

clean:
	rm *.o

distclean:
	rm *.o
	rm simulation

