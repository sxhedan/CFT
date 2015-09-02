# Sample Makefile for building libpcgnslib.a.
# 
# 11/07/2010. Duo Wang
build_cpu = x86_64
LIBPCGNS = libpcgnslib.a

# Option for HDF5. Note that HDF5 1.8 is required.
WITHHDF5 =
HDF5INC = 

# Compiler options
CC=mpicc -g -pedantic -Wno-long-long
CPPFLAGS = -I. 
CFLAGS = -O
incs =   -I.. -I../util -g   -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include  

AR = ar cr
ARFLAGS =
RANLIB = ranlib
vpath %.c adf adfh

OBJS = pcgnslib.o pcgns_util.o

# Rules
all: $(LIBPCGNS)

$(LIBPCGNS): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@
	cp $(LIBPCGNS) ../../../lib/$(build_cpu)/
$(OBJS) : %.o : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(incs) -c $< -o $@

clean:
	rm -f $(OBJS) $(LIBPCGNS)
