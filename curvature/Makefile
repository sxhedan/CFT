CC=mpicxx -g -pedantic -Wno-long-long
CXX=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g
F77_LIBS=

libext = 
incs =    -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  

curvature: FT_libs curvature.cpp
	rm -f curvature
	$(CXX) -c -I../include $(incs) curvature.cpp
	$(CXX) curvature.o -I../include -L../lib/x86_64  $(libincs) -lFronTier -lm $(libs) -o curvature

-include ../devel-deps.inc

clean:
	rm -rf *.o curvature 
tagsfile:
	ctags *.cpp ../src/*/*.[chf] ../src/*/*/*.[chf]
