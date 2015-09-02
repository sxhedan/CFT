include ../iMesh-Defs.inc
#CXX=mpicxx -g -pedantic -Wno-long-long
#F77=mpif77 -g
#F77_LIBS=

#libext = 
incs =   -I.. -I../util -g   -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include  
#libincs =  -L/usr/local/pkg/openmpi/lib 
#libs =  -lmpi -lmpi_cxx  

iTaps: FT_libs iTaps.cpp
	rm -f iTaps
	$(IMESH_CXX) -c -I../include $(incs) $(IMESH_INCLUDES) iTaps.cpp
	$(IMESH_CXX) iTaps.o $(IMESH_LDFLAGS) $(IMESH_LIBS) -o iTaps

-include ../devel-deps.inc

clean:
	rm -rf *.o iTaps 
tagsfile:
	ctags *.cpp ../src/*/*.[chf] ../src/*/*/*.[chf]
