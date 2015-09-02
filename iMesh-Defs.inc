IMESH_CXX= mpicxx -g -pedantic -Wno-long-long
IMESH_CXXFLAGS = $(CXXFLAGS) -g
IMESH_CFLAGS = $(CFLAGS) -g
IMESH_LDFLAGS = -L$(FT_DIR)/lib/x86_64   -L/usr/local/pkg/openmpi/lib 
IMESH_LIBS = -lFronTier

IMESH_F77 = mpif77 -g
IMESH_F77_LIBS= 

IMESH_INCLUDES = -I$(FT_DIR)/include/FronTier  -I$(FT_DIR)/include/FronTier/intfc

IMESH_LIBS = -lFronTier -lm  -lmpi -lmpi_cxx  

IMESH_INCLUDEDIR = include/FronTier/intfc include/FronTier/front include/FronTier/util

FT_DIR=/nfs/user01/hyulim/tar_FTI/latest_tar/FTI
