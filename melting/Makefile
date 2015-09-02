CXX=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g
F77_LIBS=

libext = 
incs =   -I.. -I../util -g   -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  
build_cpu = x86_64
incsloc = -I. -I../iFluid -I../include

PETSC = yes
PETSC_ARCH=
PETSC_DIR=/usr/local/pkg/petsc-3.1-p7

Petsc_Include=-I/usr/local/pkg/petsc-3.1-p7/include
Petsc_Lib=-L/usr/local/pkg/petsc-3.1-p7/lib -lpetsc -llapack -lblas -ldl -lm -L/usr/X11R6/lib -lX11

melting : FT_libs

marray.o : melting.h marray.cpp
	${CXX} -c marray.cpp $(incsloc) $(incs) ${Petsc_Include}

solver.o : ../iFluid/solver.h ../iFluid/solver.cpp
	${CXX} -c ../iFluid/solver.cpp $(incsloc) $(incs) ${Petsc_Include}

mcartsn.o : melting.h mcartsn.cpp
	${CXX} -c mcartsn.cpp $(incsloc) $(incs) ${Petsc_Include}

melting: melting.h melting.cpp marray.o solver.o mcartsn.o
	${CXX} -c melting.cpp -I../include $(incs) ${Petsc_Include}
	$(CXX) -o melting melting.o -I../include -L../lib/${build_cpu}  solver.o mcartsn.o marray.o ${Petsc_Lib} -L../lib/${build_cpu} -lFronTier $(libincs)  $(libs) ${F77_LIBS}

-include ../devel-deps.inc

clean:
	rm -rf *.o melting 
output :
	rm -rf out/*

gif :
	hdf2gif out/intfc.hdf out/intfc.gif

TagsFile:
	ctags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

tagsfile:
	ctags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

etagsfile:
	etags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]
