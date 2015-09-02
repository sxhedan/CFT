CC=mpicc -g -pedantic -Wno-long-long
F77=mpif77 -g
name = util
libutil_a_OBJECTS =  \
cleanup.$(OBJEXT)  \
debug.$(OBJEXT)    \
error.$(OBJEXT)    \
fft.$(OBJEXT)\
fgetstrin.$(OBJEXT)\
machine.$(OBJEXT)  \
matrix.$(OBJEXT)   \
other.$(OBJEXT)    \
fsort.$(OBJEXT)    \
output.$(OBJEXT)   \
ppsub.$(OBJEXT)    \
quad.$(OBJEXT)     \
roots.$(OBJEXT)    \
runga.$(OBJEXT)    \
sphhar.$(OBJEXT)   \
screen.$(OBJEXT)   \
sigplot.$(OBJEXT)  \
simpleio.$(OBJEXT) \
times.$(OBJEXT)    \
uinit.$(OBJEXT)    \
vectormalloc.$(OBJEXT)\
vtk.$(OBJEXT)


libext = 
OBJEXT = o
ftlibincs =  -I.. -I../util -g   -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include  
build_cpu = x86_64
RANLIB = ranlib
VPATH = $(build_cpu)


all: FORCE
	$(MAKE) ../../lib/$(build_cpu)/lib$(name)$(libext).a

FORCE:
	


make.dep.u:../util/make.dep
	cp ../util/make.dep make.util.dep
	sed -e "1,$$ s,$(name)/,,g" make.util.dep > make.dep.temp
	grep hypre make.dep.temp -v  > make.dep.u
	rm make.dep.temp

../util/make.dep:
	cd ..;./util/mkfiles.pl
	cd ..;./util/mkdep.pl
	cd ..;mv make.dep util/make.dep


../../lib/$(build_cpu)/lib$(name)$(libext).a: make.dep.u $(build_cpu) $(libutil_a_OBJECTS)
	-mkdir ../../lib
	-mkdir ../../lib/$(build_cpu)
	-rm -rf libutil.a
	cd $(build_cpu);$(AR) cru libutil.a $(libutil_a_OBJECTS)
	cd $(build_cpu);$(RANLIB) libutil.a; mv libutil.a ../../../lib/$(build_cpu)/libutil$(libext).a

-include make.dep.u

diff:
	-diff -r . $(diffdir)/src/$(name) > DIFF


$(build_cpu):
	-mkdir $(build_cpu)
$(libutil_a_OBJECTS): %.o: %.c
	$(CC) -I.   -I../.. -I../  -g   -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include    -c $< -o $(build_cpu)/$@ 
 

clean:
	cd  $(build_cpu)/; rm -rf *.o
