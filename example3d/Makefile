CXX=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g
F77_LIBS=

libext = 
incs =    -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  

all: example

example: FT_libs example.cpp
	rm -f example
	$(CXX) -c -I../include $(incs) example.cpp
	$(CXX) example.o -I../include -L../lib/x86_64 $(libincs) -lFronTier -lm $(libs) -o example
clean:
	rm -rf *.o example 
tagsfile:
	ctags *.[chf] ../src/*/*.[chf] ../src/*/*/*.[chf]

-include ../devel-deps.inc
