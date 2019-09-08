
### Sets the particular extensions used in this makefile.
CC=g++
Opt=-O3 -std=c++11 -Wall
CD=.

Inc= -I $(CD)


test_template_qaqmc.out: $(CD)/test_template_qaqmc.cxx $(CD)/template_qaqmc.h
	$(CC) $(Opt) $(Inc) -o test_template_qaqmc.out $(CD)/test_template_qaqmc.cxx 

test_template_proj.out: $(CD)/test_template_proj.cxx $(CD)/template_proj.h
	$(CC) $(Opt) $(Inc) -o test_template_proj.out $(CD)/test_template_proj.cxx 

all:
	make test_template_qaqmc.out
	make test_template_proj.out