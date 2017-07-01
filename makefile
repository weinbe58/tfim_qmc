
### Sets the particular extensions used in this makefile.
CC=g++
Opt=-O3 -std=c++11
CD=.

Inc= -I $(CD)


twoway.out: $(CD)/test.cpp $(CD)/qaqmc.h
	$(CC) $(Opt) $(Inc) -o twoway.out $(CD)/test.cpp 



