CXX=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g
F77_LIBS=
libext = 
incs =   -I.. -I../util -g   -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  
build_cpu = x86_64
incsloc = -I. -I../include -I../iFluid -I../crystal

PETSC = yes
PETSC_ARCH=
PETSC_DIR=/usr/local/pkg/petsc-3.1-p7
Petsc_Include=-I/usr/local/pkg/petsc-3.1-p7/include
Petsc_Lib=-L/usr/local/pkg/petsc-3.1-p7/lib -lpetsc -llapack -lblas -ldl -lm -L/usr/X11R6/lib -lX11

subsurf : FT_libs ../lib/libiF.a

iFbasic.o : ../iFluid/iFluid.h ../iFluid/iFbasic.cpp
	${CXX} -c ../iFluid/iFbasic.cpp $(incsloc) $(incs) ${Petsc_Include}

iFcartsn2d.o : ../iFluid/iFluid.h ../iFluid/iFcartsn2d.cpp
	${CXX} -c ../iFluid/iFcartsn2d.cpp $(incsloc) $(incs) ${Petsc_Include}

iFcartsn3d.o : ../iFluid/iFluid.h ../iFluid/iFcartsn3d.cpp
	${CXX} -c ../iFluid/iFcartsn3d.cpp $(incsloc) $(incs) ${Petsc_Include}

solver.o : ../iFluid/solver.h ../iFluid/solver.cpp
	${CXX} -c ../iFluid/solver.cpp $(incsloc) $(incs) ${Petsc_Include}
ellip.o : ../iFluid/solver.h ../iFluid/ellip.cpp
	${CXX} -c ../iFluid/ellip.cpp $(incsloc) $(incs) ${Petsc_Include}

ccartsn.o : ../crystal/crystal.h ../iFluid/solver.h ../crystal/ccartsn.cpp
	${CXX} -c ../crystal/ccartsn.cpp $(incsloc) $(incs) ${Petsc_Include}

subsub.o : subsurf.h ../iFluid/iFluid.h subsub.cpp
	${CXX} -c subsub.cpp $(incsloc) $(incs) ${Petsc_Include}

subsurf: ../iFluid/iFluid.h subsurf.h subsurf.cpp subsub.o ccartsn.o ../lib/libiF.a
	${CXX} -c subsurf.cpp $(incsloc) $(incs) ${Petsc_Include}
	$(CXX) -o subsurf subsurf.o subsub.o ccartsn.o -I../include -L../lib/${build_cpu} -liF ${Petsc_Lib} -L../lib/${build_cpu} -lFronTier $(libincs)  $(libs) ${F77_LIBS}

../lib/libiF.a: solver.o ellip.o iFbasic.o iFcartsn2d.o iFcartsn3d.o iFcylin3d.o iFdebug2d.o iFdebug3d.o
	$(AR) cru libiF.a $^
	ranlib libiF.a; mv libiF.a ../lib/$(build_cpu)/libiF.a

-include ../devel-deps.inc

clean:
	rm -rf *.o subsurf 
output :
	rm -rf out/*

gif :
	hdf2gif out/intfc.hdf out/intfc.gif

TagsFile:
	ctags *.h *.cpp ../iFluid/{iFluid.h,solver.h,solver.cpp,ellip.cpp,iFbasic.cpp,iFcartsn2d.cpp,iFcartsn3d.cpp,iFcylin3d.cpp} ../crystal/{crystal.h,ccartsn.cpp} ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

tagsfile:
	ctags *.h *.cpp ../iFluid/{iFluid.h,solver.h,solver.cpp,ellip.cpp,iFbasic.cpp,iFcartsn2d.cpp,iFcartsn3d.cpp,iFcylin3d.cpp} ../crystal/{crystal.h,ccartsn.cpp} ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

etagsfile:
	etags *.h *.cpp ../iFluid/{iFluid.h,solver.h,solver.cpp,ellip.cpp,iFbasic.cpp,iFcartsn2d.cpp,iFcartsn3d.cpp,iFcylin3d.cpp} ../crystal/{crystal.h,ccartsn.cpp} ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]
