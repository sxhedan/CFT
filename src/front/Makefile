CC=mpicc -g -pedantic -Wno-long-long
F77=mpif77 -g
name = front
FT_C_OBJECTS = fredist2d.$(OBJEXT) \
fadv.$(OBJEXT) \
fbdry1.$(OBJEXT) \
fbdry2.$(OBJEXT) \
fbdry3.$(OBJEXT) \
fbdry4.$(OBJEXT) \
fcheck3d.$(OBJEXT) \
fcorrspnd.$(OBJEXT) \
fcrosscur.$(OBJEXT) \
fcrossext.$(OBJEXT) \
fcrstatus.$(OBJEXT) \
fdiagnostic.$(OBJEXT) \
fgb2d.$(OBJEXT) \
fgb3dutil.$(OBJEXT) \
fgb3d.$(OBJEXT) \
fgb3comp.$(OBJEXT) \
fgrid.$(OBJEXT) \
finit.$(OBJEXT) \
fint.$(OBJEXT) \
fnode.$(OBJEXT) \
fnodesub.$(OBJEXT) \
fprint.$(OBJEXT) \
fprop2d.$(OBJEXT) \
fprop3d.$(OBJEXT) \
fredist.$(OBJEXT) \
fredist3d.$(OBJEXT) \
frp1.$(OBJEXT) \
frp2.$(OBJEXT) \
fscatter.$(OBJEXT) \
fredist1d.$(OBJEXT) \
fscat1d.$(OBJEXT) \
fscat2d.$(OBJEXT) \
fscat3d1.$(OBJEXT) \
fscat3d2.$(OBJEXT) \
fstate2d.$(OBJEXT) \
fstate.$(OBJEXT) \
ftop.$(OBJEXT) \
fsub.$(OBJEXT) \
funtan2d.$(OBJEXT) \
funtan3d.$(OBJEXT) \
fuserintfc.$(OBJEXT) \
fuserhooks.$(OBJEXT) \
foverture_patch.$(OBJEXT) \
foverture_adv.$(OBJEXT)   \
fvelo.$(OBJEXT)   \
fmap.$(OBJEXT)   \
fpatch3d.$(OBJEXT)   \
########################################################################################
# DO NOT MODIFY BELOW THIS LINE UNLESS YOUR POSITIVE YOU KNOW WHAT YOU WANT TO DO 
########################################################################################



libext = 
OBJEXT = o
ftlibincs =  -I.. -I../util -g   -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__MPI__ -I/usr/local/pkg/openmpi/include  
build_cpu = x86_64
RANLIB = ranlib
incs =   -D__MPI__ -I/usr/local/pkg/openmpi/include  
libincs =  -L/usr/local/pkg/openmpi/lib 
libs =  -lmpi -lmpi_cxx  

all: 
	cd ../intfc; $(MAKE)
	cd ../util; $(MAKE)
	-rm ../../lib/$(build_cpu)/libFronTier.a
	$(MAKE) ../../lib/$(build_cpu)/libFronTier.a


../../lib/$(build_cpu)/libFronTier.a: ../../lib/$(build_cpu)/lib$(name)$(libext).a
	cd ../intfc/; $(MAKE)
	cd ../../lib/$(build_cpu)/; rm -rf *.o; ar x libintfc.a; ar x libutil.a; ar x libfront.a; ar cru libFronTier.a *.o; ranlib libFronTier.a;  rm -rf *.o


######## FT BUILD INSTRUCTIONS ##############################
VPATH =  $(build_cpu):../
../../lib/$(build_cpu)/lib$(name)$(libext).a: make.dep $(build_cpu) $(FT_C_OBJECTS)
	-rm -rf lib$(name)$(libext).a
	-rm -rf ${build_cpu}/lib$(name)$(libext).a
	-mkdir ../../lib/
	-mkdir ../../lib/$(build_cpu)/
	cd $(build_cpu); $(AR) cru lib$(name).a $(FT_C_OBJECTS)
	cd $(build_cpu); $(RANLIB) lib$(name).a; mv lib$(name).a ../../../lib/$(build_cpu)/lib$(name)$(libext).a

$(build_cpu):
	-mkdir $(build_cpu)

make.dep:../util/make.dep
	cp ../util/make.dep make.util.dep
	sed -e "1,$$ s,$(name)/,,g" make.util.dep > make.dep

../util/make.dep:
	cd ..;./util/mkfiles.pl
	cd ..;./util/mkdep.pl
	cd ..;mv make.dep util/make.dep

include make.dep

diff:
	-diff -r . $(diffdir)/src/$(name) > DIFF

$(FT_C_OBJECTS):%.o: %.c
	$(CC) $(ftlibincs) -c $< -o $(build_cpu)/$@
######## END FT BUILD INSTRUCTIONS ##############################


clean:
	cd  $(build_cpu)/; rm -rf *.o
	cd ../intfc; $(MAKE) clean
tagsfile:
	ctags ../intfc/*.c ../util/*.c *.c ../intfc/*.h ../util/*.h *.h

$(build_cpu)/testfront.o:testfront.c
	$(CC) $(ftlibincs) $(incs) -c testfront.c -o $(build_cpu)/testfront.o
testfront: FORCE
	$(MAKE)
	$(MAKE) $(build_cpu)/testfront.o 
	$(CC) $(build_cpu)/testfront.o -L../../lib/$(build_cpu) $(libincs) -lfront$(libext) -lintfc$(libext) -lutil$(libext) $(libs) -lm -o testfront

TARFLIST = $$ME/README \
$$ME/build.sh \
$$ME/advbuild.sh \
$$ME/install-sh \
$$ME/config.status \
$$ME/config.sub \
$$ME/config.guess \
$$ME/devel-deps.inc.in \
$$ME/iMesh-Defs.inc.in \
$$ME/configure.in \
$$ME/Makefile.in \
$$ME/src/Makefile.in \
$$ME/README \
$$ME/install-sh \
$$ME/config.status \
$$ME/config.sub \
$$ME/config.guess \
$$ME/example/*.cpp \
$$ME/example/Makefile.in \
$$ME/example3d/*.cpp \
$$ME/example3d/Makefile.in \
$$ME/curvature/*.cpp \
$$ME/curvature/*.h \
$$ME/curvature/Makefile.in \
$$ME/subsurf/*.cpp \
$$ME/subsurf/*.h \
$$ME/subsurf/Makefile.in \
$$ME/airfoil/*.cpp \
$$ME/airfoil/*.h \
$$ME/airfoil/Makefile.in \
$$ME/crystal/*.cpp \
$$ME/crystal/*.h \
$$ME/crystal/Makefile.in \
$$ME/crystal/README \
$$ME/melting/*.cpp \
$$ME/melting/*.h \
$$ME/melting/Makefile.in \
$$ME/iFluid/*.cpp \
$$ME/iFluid/*.h \
$$ME/iFluid/Makefile.in \
$$ME/iFluid/in-* \
$$ME/cFluid/*.cpp \
$$ME/cFluid/*.h \
$$ME/cFluid/Makefile.in \
$$ME/cFluid/in-* \
$$ME/frgb/*.cpp \
$$ME/frgb/*.h \
$$ME/frgb/Makefile.in \
$$ME/frgb/in-* \
$$ME/finance/*.cpp \
$$ME/finance/*.h \
$$ME/finance/Makefile.in \
$$ME/cell/*.cpp \
$$ME/cell/Makefile.in \
$$ME/pde/*.cpp \
$$ME/pde/in-* \
$$ME/pde/Makefile.in \
$$ME/poisson/*.cpp \
$$ME/poisson/*.h \
$$ME/poisson/Makefile.in \
$$ME/iTaps/*.cpp \
$$ME/iTaps/Makefile.in \
$$ME/testfront/*.cpp \
$$ME/testfront/in/* \
$$ME/testfront/Makefile.in \
$$ME/src/util/mkdep.pl \
$$ME/src/util/mkfiles.pl \
$$ME/src/util/ftrules \
$$ME/src/util/gasrules \
$$ME/src/pak/Makefile.in \
$$ME/src/pak/odepack/Makefile.in \
$$ME/src/pak/odepack/*.f \
$$ME/src/pak/odepack/*.c \
$$ME/src/pak/dierckx/Makefile.in \
$$ME/src/pak/dierckx/*.f \
$$ME/src/pak/linpak/Makefile.in \
$$ME/src/pak/linpak/*.f \
$$ME/src/pak/blas/Makefile.in \
$$ME/src/pak/blas/*.f \
$$ME/src/intfc/*.c \
$$ME/src/intfc/*.h  \
$$ME/src/intfc/in \
$$ME/src/intfc/Makefile.in \
$$ME/src/intfc/in \
$$ME/src/util/*.h \
$$ME/src/util/*.c \
$$ME/src/util/Makefile.in \
$$ME/src/front/*.c \
$$ME/src/front/*.h \
$$ME/src/front/Makefile.in \
$$ME/src/front/in \
$$ME/devel-deps.inc.in \
$$ME/src/pak/cgnslib_3.0/Makefile* \
$$ME/src/pak/cgnslib_3.0/LICENSE* \
$$ME/src/pak/cgnslib_3.0/adf*/*.[ch] \
$$ME/src/pak/cgnslib_3.0/*.[ch]

tar:
	-mkdir in
	-cd ../intfc; mkdir in
	-cd testfront
	-rm *.sc 
	cd ../../ && ME=`basename $$PWD` && cd .. && tar -cf "`date +FronTier.Lite.%m_%d_%y.tar`" ${TARFLIST} && \
	gzip "`date +FronTier.Lite.%m_%d_%y.tar`" && mv `date +FronTier.Lite.%m_%d_%y.tar`.gz $$ME/.
