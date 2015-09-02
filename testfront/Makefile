CC=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g
libext=
incs=  -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs= -L/usr/local/pkg/openmpi/lib 
libs= -lmpi -lmpi_cxx  
build_cpu=x86_64

all:
	cd ../;make
	make testfront 



testfront: testfront.cpp
	rm -f testfront
	$(CC) -c -I../include $(incs) testfront.cpp
	$(CC) testfront.o -L../lib/$(build_cpu) $(libincs) -lFronTier -lm $(libs) -o testfront
clean:
	rm -rf *.o testfront 
tagsfile:
	ctags *.cpp ../src/*/*.[chf]
