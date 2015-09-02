CC=mpicc -g -pedantic -Wno-long-long
F77=mpif77 -g
all:
	if [ 0 -eq 1 ]; then cd pak/cgnslib_3.0 && $(MAKE); fi
	# if [ 0 -eq 1 ]; then cd pak/pcgnslib_0.2.0 && $(MAKE); fi
	cd front  && $(MAKE)
	  
	if [ 1 -eq 1 ]; then cd pak/linpak && $(MAKE); fi
	if [ 1 -eq 1 ]; then cd pak/blas && $(MAKE); fi

clean:
	cd util && $(MAKE) clean
	cd front && $(MAKE) clean
	cd intfc && $(MAKE) clean
	 
	if [ 1 -eq 1 ]; then cd pak/linpak && $(MAKE) clean; fi
	if [ 1 -eq 1 ]; then cd pak/blas && $(MAKE) clean; fi
	if [ 0 -eq 1 ]; then cd pak/cgnslib_3.0 && $(MAKE) clean; fi
	# if [ 0 -eq 1 ]; then cd pak/pcgnslib_0.2.0 && $(MAKE) clean; fi

diff:
	cd util && $(MAKE) diff 
	cd front && $(MAKE) diff
	cd intfc && $(MAKE) diff
