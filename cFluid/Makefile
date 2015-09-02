CC=mpicxx -g -pedantic -Wno-long-long
CXX=mpicxx -g -pedantic -Wno-long-long -std=c++0x
F77=mpif77 -g
F77_LIBS=

libext = 
incs =   -I.. -I../util -g   -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  
build_cpu = x86_64
incsloc = -I. -I../include

PETSC = yes
PETSC_ARCH=
PETSC_DIR=/usr/local/pkg/petsc-3.1-p7

Petsc_Include=-I/usr/local/pkg/petsc-3.1-p7/include
Petsc_Lib=-L/usr/local/pkg/petsc-3.1-p7/lib -lpetsc -llapack -lblas -ldl -lm -L/usr/X11R6/lib -lX11
all: cFluid ../lib/libcF.a

cFluid :

cFtvd.o : cFluid.h cFtvd.cpp
	${CXX} -c cFtvd.cpp $(incsloc) $(incs) ${Petsc_Include}

cFweno.o : cFluid.h cFweno.cpp
	${CXX} -c cFweno.cpp $(incsloc) $(incs) ${Petsc_Include}

cFsub.o : cFluid.h cFsub.cpp
	${CXX} -c cFsub.cpp $(incsloc) $(incs) ${Petsc_Include}

cFcartsn.o : cFluid.h cFcartsn.cpp
	${CXX} -c cFcartsn.cpp $(incsloc) $(incs) ${Petsc_Include}

cFbase.o : cFluid.h cFbase.cpp
	${CXX} -c cFbase.cpp $(incsloc) $(incs) ${Petsc_Include}

cFphys.o : cFluid.h cFphys.cpp
	${CXX} -c cFphys.cpp $(incsloc) $(incs) ${Petsc_Include}

cFinit.o : cFluid.h cFinit.cpp
	${CXX} -c cFinit.cpp $(incsloc) $(incs) ${Petsc_Include}

cFeos.o : cFluid.h cFeos.cpp
	${CXX} -c cFeos.cpp $(incsloc) $(incs) ${Petsc_Include}

cFgriemann.o: cFluid.h cFgriemann.cpp
	${CXX} -c cFgriemann.cpp $(incsloc) $(incs) ${Petsc_Include}

cFluid: cFluid.h cFluid.cpp cFtvd.o cFweno.o cFsub.o cFcartsn.o cFbase.o cFphys.o cFinit.o cFeos.o cFgriemann.o 
	${CXX} -c cFluid.cpp -I../include $(incs) ${Petsc_Include}
	$(CXX) -o cFluid cFluid.o -g -I../include -L../lib/${build_cpu}  cFtvd.o cFweno.o cFsub.o cFcartsn.o cFbase.o cFphys.o cFinit.o cFeos.o cFgriemann.o ${Petsc_Lib} -L../lib/${build_cpu} -lFronTier $(libincs)  $(libs) ${F77_LIBS}

../lib/libcF.a: cFtvd.o cFweno.o cFsub.o cFcartsn.o cFbase.o cFphys.o cFinit.o cFeos.o
	$(AR) cru libcF.a $<
	ranlib libcF.a; mv libcF.a ../lib/$(build_cpu)/libcF.a

-include ../devel-deps.inc

clean:
	rm -rf *.o cFluid ../lib/$(build_cpu)/libcF.a 
output :
	rm -rf out/*

TagsFile:
	ctags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

tagsfile:
	ctags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

etagsfile:
	etags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]
