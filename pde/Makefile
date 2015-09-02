CXX=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g

libext = 
incs =   -I.. -I../util -g   -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  

wave: FT_libs wave.cpp
	rm -f wave
	$(CXX) -c -I../include $(incs) wave.cpp
	$(CXX) wave.o -I../include -L../lib/x86_64  $(libincs) -lFronTier -lm $(libs) -o wave

heat: FT_libs heat.cpp
	rm -f heat
	$(CXX) -c -I../include $(incs) heat.cpp
	$(CXX) heat.o -I../include -L../lib/x86_64  $(libincs) -lFronTier -lm $(libs) -o heat

bseq: FT_libs bseq.cpp
	rm -f bseq
	$(CXX) -c -I../include $(incs) bseq.cpp
	$(CXX) bseq.o -I../include -L../lib/x86_64  $(libincs) -lFronTier -lm $(libs) -o bseq

heat: FT_libs heat.cpp
-include ../devel-deps.inc

clean:
	rm -rf *.o wave
tagsfile:
	ctags *.cpp ../src/*/*.[chf] ../src/*/*/*.[chf]
