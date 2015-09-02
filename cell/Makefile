CC=mpicxx -g -pedantic -Wno-long-long
CXX=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g
F77_LIBS=

libext = 
incs =    -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  

cell: FT_libs cell.cpp
	rm -f cell
	$(CXX) -c -I../include $(incs) cell.cpp
	$(CXX) cell.o -I../include -L../lib/x86_64  $(libincs) -lFronTier -lm $(libs) -o cell

-include ../devel-deps.inc

clean:
	rm -rf *.o cell 
tagsfile:
	ctags *.cpp ../src/*/*.[chf] ../src/*/*/*.[chf]
