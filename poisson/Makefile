CXX=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g
F77_LIBS=

libext = 
incs =    -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  
build_cpu = x86_64

PETSC = yes
PETSC_ARCH=
PETSC_DIR=/usr/local/pkg/petsc-3.1-p7

Petsc_Include=-I/usr/local/pkg/petsc-3.1-p7/include
Petsc_Lib=-L/usr/local/pkg/petsc-3.1-p7/lib -lpetsc -llapack -lblas -ldl -lm -L/usr/X11R6/lib -lX11

poisson : FT_libs
test_ebm2d : FT_libs
test_ebm3d : FT_libs

lpetsc.o : lpetsc.cpp poisson.h
	${CXX} -c lpetsc.cpp -I. -I../include $(incs) ${Petsc_Include}

poisson.o : poisson.cpp poisson.h
	${CXX} -c poisson.cpp -I. -I../include $(incs) ${Petsc_Include}

poisson: poisson.o lpetsc.o
	${CXX} -o poisson poisson.o -I../include -L../lib/${build_cpu}  lpetsc.o ${Petsc_Lib} -L../lib/${build_cpu} -lFronTier $(libincs)  $(libs) ${F77_LIBS}

# EBM
solver_petsc.o : solver_petsc.cpp solver_petsc.h solver.h
	${CXX} -c solver_petsc.cpp -I. -I../include $(incs) ${Petsc_Include}	
solver_lapack.o : solver_lapack.h solver_lapack.cpp solver.h
	${CXX} -c solver_lapack.cpp -I.

frontier_ppgrid.o : frontier_ppgrid.h frontier_ppgrid.cpp
	${CXX} -c frontier_ppgrid.cpp -I. -I../include $(incs) ${Petsc_Include}			
geometry.o : geometry.h geometry.cpp
	${CXX} -c geometry.cpp -I. -I../include $(incs) ${Petsc_Include}
ebm2d.o : ebm2d.cpp ebm2d.h solver_petsc.o frontier_ppgrid.o geometry.o
	${CXX} -c ebm2d.cpp -I. -I../include $(incs) ${Petsc_Include}			
ebm3d.o : ebm3d.cpp ebm3d.h solver_petsc.o solver_lapack.o frontier_ppgrid.o
	${CXX} -c ebm3d.cpp -I. -I../include $(incs) ${Petsc_Include}			
test_ebm2d.o : test_ebm2d.cpp ebm2d.o solver_petsc.o 
	${CXX} -c test_ebm2d.cpp -I. -I../include $(incs) ${Petsc_Include}				
test_ebm2d : test_ebm2d.o ebm2d.o solver_petsc.o frontier_ppgrid.o geometry.o 
	${CXX} -o test_ebm2d test_ebm2d.o ebm2d.o solver_petsc.o frontier_ppgrid.o geometry.o -I../include -L../lib/${build_cpu} ${Petsc_Lib} -L../lib/${build_cpu} -lFronTier $(libincs)  $(libs) ${F77_LIBS}

test_ebm3d.o : test_ebm3d.cpp ebm3d.o solver_petsc.o solver_lapack.o
	${CXX} -c test_ebm3d.cpp -I. -I../include $(incs) ${Petsc_Include}				
test_ebm3d : test_ebm3d.o ebm3d.o solver_petsc.o solver_lapack.o frontier_ppgrid.o
	${CXX} -o test_ebm3d test_ebm3d.o ebm3d.o solver_petsc.o solver_lapack.o frontier_ppgrid.o -I../include -L../lib/${build_cpu} ${Petsc_Lib} -L../lib/${build_cpu} -lFronTier $(libincs)  $(libs) ${F77_LIBS}

rungekutta_tga.o : rungekutta_tga.h rungekutta_tga.cpp
	${CXX} -c rungekutta_tga.cpp -I.

-include ../devel-deps.inc

clean:
	rm -rf *.o poisson test_ebm2d test_ebm3d 
clear:
	rm -rf *~
	rm -rf #*
plt :
	/bin/rm -rf *.plt *.vtk

backup :
	cp *.[ch] bak/
	cp makefile bak/
	cp *.txt bak/
save :
	cp *.[ch] ~/bak/poisson
	cp Makefile.in ~/bak/poisson
	cp *.txt ~/bak/poisson


tagsfile:
	ctags *.[ch] ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

etagsfile:
	etags *.[ch] ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

tag2d:
	etags test_ebm2d.[ch] ebm2d.[ch] frontier_ppgrid.[ch] ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]	
tag3d:
	etags test_ebm3d.[ch] ebm3d.[ch] frontier_ppgrid.[ch] ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]	

ebrowse:
	ebrowse *.[ch] ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

