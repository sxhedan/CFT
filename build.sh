#!/bin/bash

function parse_arguments {
 # Parse command-line arguments
    if [[ "$1" == "-h" ]]; then
        echo "$0    High-level build script for FronTier."
        echo "Usage: $0 [-d] [-n] [-g] [--with-cgns]  [--with-hdf]"
        echo
        echo "    -d          Enable debugging."
        echo "    -n          Just configure. Do not run make."
        echo "    -g          Enable compilation of Gas directory."
        echo "    --with-cgns Enable compilation with CGNS."
        echo "    --with-hdf5 Enable compilation with CGNS+HDF5."
        echo "    --with-hdf  Enable compilation with HDF4."
        echo "    --with-gd   Enable compilation with graphics drawing for 1-D post processing."
        echo "    --enable-itaps  Enable compilation with ITAPS."
        echo "    --enable-imesh  Enable compilation with IMESH."
	echo "    --with-***      Same parameter as old configure file."
        exit
    fi
    
    OPTS="-O3"
    for arg in $@ ; do
        if [[ "$arg" == "-d" ]]; then
	    OPTS="-g"
	    COPTS_GCC="-pedantic -Wno-long-long"
	    COPTS_ICC="-Wall"
        elif [[ "$arg" == "-n" ]]; then
	    NOMAKE=1
        elif [[ "$arg" == "--enable-itaps" ]]; then
	    CONF="$CONF --enable-itaps"
        elif [[ "$arg" == "--enable-imesh" ]]; then
	    CONF="$CONF --enable-imesh"
        elif [[ "$arg" == "-g" ]]; then
	    CONF="$CONF --with-gas"
        elif [[ "$arg" == "--with-cgns" ]]; then
            CONF="$CONF --with-cgns"
        elif [[ "$arg" == "--with-hdf5" ]]; then
            CONF="$CONF --with-cgns"
            WITHHDF5=1
        elif [[ "$arg" == "--with-hdf" ]]; then
            WITHHDF=1
        elif [[ "$arg" == "--with-gd" ]]; then
            WITHGD=1
	else
	    CONF="$CONF ${arg}"
        fi
    done
    if [[ "$OPTS" != "-g" ]]; then
        CONF="$CONF --with-no-debugging"
    fi
}

#################
function config_generic {
    # Generic platforms. 
    HASMPI=`which mpicxx`
    if [[ -n "$HASMPI" ]]; then
        # Assume OpenMPI.
        export CXX="mpicxx ${OPTS}"
        export F77="mpif77 ${OPTS}"
        export CC="mpicc ${OPTS}"
    else
        # Assume gcc and gfortran.
        export CXX="g++ ${OPTS}"
        export F77="gfortran ${OPTS}"
        export CC="gcc ${OPTS}"
    fi

    echo -n "Enter command for C compiler. (Hit RETURN to use default \"$CC\"): "
    read CC_
    if  [[ -n "$CC_" ]]; then export CC=$CC_; fi

    echo -n "Enter command for C++ compiler. (Hit RETURN to use default \"$CXX\"): "
    read CXX_
    if  [[ -n "$CXX2" ]]; then export CXX=$CXX_; fi

    echo -n "Enter command for F77 compiler. (Hit RETURN to use default \"$F77)\": "
    read F77_
    if  [[ -n "$F77_" ]]; then export F77=$F77_; fi

    if [[ ${CC:0:2} == "mp" ]] && [[ ${CXX:0:2} == "mp" ]] && [[ ${F77:0:2} == "mp" ]]; then
        CONF="$CONF --with-mpi"
    fi

    if [[ -n "$WITHHDF" ]]; then
        echo -n "Enter root directory for HDF4. (Hit RETURN to skip HDF4.): "
        read HDF4_DIR
        if [[ -n "$HDF4_DIR" ]]; then
            CONF="$CONF --with-hdf=$HDF4_DIR"
        fi
    fi
    if [[ -n "$WITHHDF5" ]]; then
        echo -n "Enter root directory for HDF5. (Hit RETURN to skip HDF5.): "
        read HDF5_DIR
        if [[ -n "$HDF5_DIR" ]]; then
            CONF="$CONF --with-hdf5=$HDF5_DIR"
        fi
    fi

    echo -n "Enter root directory for PETSC. PETSc is required to compile iFluid. Note that PETSc requires MPI, and the MPI wrappers you specified for CC and CXX must be the same as those used to compile PETSc. (Hit RETURN to skip PETSc.): "
    read PETSC_DIR
    export PETSC_DIR

    if [[ -n "$PETSC_DIR" ]]; then
        export CONF="--with-petsc=$PETSC_DIR --with-devel ${CONF}"

        echo -n "Enter PETSC_ARCH: "
        read PETSC_ARCH
        export PETSC_ARCH

        export PETSC_INCLUDE="-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include"

        export PETSC_LIB="-lblas -llapack -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -ldl -lm -lX11"
        echo -n "Enter PETSC_LIB. Note that you may also need to include libraries for LAPACK and BLAS also. (Hit return to use default $PETSC_LIB):"
        read PETSC_LIB_
        if  [[ -n "$PETSC_LIB_" ]]; then export PETSC_LIB=$PETSC_LIB_; fi
    fi

    PMAKE="-j2"
}

#################
function config_vogon {
    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    export PETSC_DIR=/usr/local/pkg/petsc-3.1-p7
    export PETSC_INCLUDE="-I${PETSC_DIR}/include"
    export PETSC_LIB="-L${PETSC_DIR}/lib -lpetsc -llapack -lblas -ldl -lm -L/usr/X11R6/lib -lX11"

    if [[ -n "$WITHHDF" ]]; then
        CONF="$CONF --with-hdf=/usr/local/pkg/HDF4"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/nfs/user01/duowang/hdf5-1.8.5"
        #CONF="$CONF --with-hdf5=/usr/local/pkg/HDF5"
    fi
    if [[ -n "$WITHGD" ]]; then
        CONF="$CONF --with-gd=/usr"
    fi

    export CONF="--with-openmpi=/usr/local/pkg/openmpi --with-petsc=${PETSC_DIR} -with-devel ${CONF}"
    PMAKE="-j2"
}

#################
function config_tulin {
    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    echo "export LD_LIBRARY_PATH=/usr/local/petsc/arch-linux2-c-debug/lib:\$LD_LIBRARY_PATH"

    export PETSC_DIR=/usr/local/petsc
    export PETSC_ARCH=arch-linux2-c-debug
    export PETSC_INCLUDE="-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include"
    export PETSC_LIB="-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -llapack -lblas -ldl -lm -L/usr/lib64 -lX11"

    if [[ -n "$WITHHDF" ]]; then
        echo "HDF not configured on TULIN"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        #CONF="$CONF --with-hdf5=/nfs/user01/duowang/hdf5-1.8.5"
        #CONF="$CONF --with-hdf5=/usr/local/pkg/HDF5"
	echo "HDF5 not configured on TULIN"
    fi
    if [[ -n "$WITHGD" ]]; then
        echo "GD not configured on TULIN"
    fi

    export CONF="--with-mpich=/usr/include/mpich2-x86_64 --with-petsc=${PETSC_DIR}"
    PMAKE="-j2"

}

#################
function config_hyun {
    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    #echo "export LD_LIBRARY_PATH=/usr/local/petsc/arch-linux2-c-debug/lib:\$LD_LIBRARY_PATH"

    export PETSC_DIR=
    export PETSC_ARCH=
    export PETSC_INCLUDE=
    export PETSC_LIB=

    if [[ -n "$WITHHDF" ]]; then
        echo "HDF not configured on TULIN"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        #CONF="$CONF --with-hdf5=/nfs/user01/duowang/hdf5-1.8.5"
        #CONF="$CONF --with-hdf5=/usr/local/pkg/HDF5"
	echo "HDF5 not configured on TULIN"
    fi
    if [[ -n "$WITHGD" ]]; then
        echo "GD not configured on TULIN"
    fi

    #export CONF="--with-mpich=/usr/include/openmpi-x86_64"
    export CONF="--with-mpi CC=mpicc CXX=mpicxx F77=mpif77"
    PMAKE="-j2"

}

#################
function config_ubuntu {
    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    export PETSC_DIR=/usr/lib/petscdir/3.0.0
    export PETSC_ARCH=linux-gnu-c-opt
    export PETSC_INCLUDE="-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include"
    export PETSC_LIB="-lblas -llapack -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -ldl -lm -L/usr/X11R6/lib -lX11"

    if [[ -n "$WITHHDF" ]]; then
        CONF="$CONF --with-hdf=/usr/lib"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/usr/lib"
    fi
    if [[ -n "$WITHGD" ]]; then
        CONF="$CONF --with-gd=/usr"
    fi

    export CONF="--with-mpi --with-petsc=${PETSC_DIR} -with-devel ${CONF}"
    PMAKE="-j2"
}

#################
function config_gentooORNL {
    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    export PETSC_DIR=/home/dealmeida/codes-tryout/petsc-3.1-p8_ompi-1.4.3_gcc-4.4.4_lpt1
    export PETSC_INCLUDE="-I${PETSC_DIR}/include"

    export PETSC_LIB="-L${PETSC_DIR}/lib -lpetsc -lflapack -lfblas -ldl -lm -L/usr/lib64 -lX11"

    if [[ -n "$WITHHDF" ]]; then
        CONF="$CONF --with-hdf=/usr/local/pkg/HDF4"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/nfs/user01/duowang/hdf5-1.8.5"
        #CONF="$CONF --with-hdf5=/usr/local/pkg/HDF5"
    fi
    if [[ -n "$WITHGD" ]]; then
        CONF="$CONF --with-gd=/usr"
    fi

    export CONF="--with-openmpi=/usr/local/ompi --with-petsc=${PETSC_DIR} --with-devel ${CONF}"
    PMAKE="-j2"
}

#################
function config_fission {
    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    export PETSC_DIR=/home/dealvf/libs/petsc-3.1-p8_ompi_gcc-1.4.3_fission
    export PETSC_INCLUDE="-I${PETSC_DIR}/include"

    export PETSC_LIB="-L${PETSC_DIR}/lib -lpetsc -lflapack -lfblas -ldl -lm -L/usr/X11R6/lib64 -lX11"

    if [[ -n "$WITHHDF" ]]; then
        CONF="$CONF --with-hdf=/usr/local/pkg/HDF4"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/nfs/user01/duowang/hdf5-1.8.5"
        #CONF="$CONF --with-hdf5=/usr/local/pkg/HDF5"
    fi
    if [[ -n "$WITHGD" ]]; then
        CONF="$CONF --with-gd=/usr"
    fi

    export CONF="--with-openmpi=/apps/local/openmpi/1.4.3/gcc-4.1.2/opt --with-petsc=${PETSC_DIR} --with-devel ${CONF}"
    PMAKE="-j2"
}

#################
function config_galaxy {
    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    export PETSC_DIR=/usr/local/pkg/petsc-3.1-p3
    export PETSC_ARCH=linux-gnu-c-debug
    export PETSC_INCLUDE="-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include"
    export PETSC_LIB="-lblas -llapack -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -ldl -lm -L/usr/X11R6/lib -lX11"

    if [[ -n "$WITHHDF" ]]; then
        CONF="$CONF --with-hdf=/usr/local/pkg/HDF4"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/nfs/user01/duowang/hdf5-1.8.5"
        #CONF="$CONF --with-hdf5=/usr/local/pkg/HDF5"
    fi
    if [[ -n "$WITHGD" ]]; then
        CONF="$CONF --with-gd=/usr"
    fi

    export CONF="--with-openmpi=/usr/local/pkg/openmpi --with-petsc=${PETSC_DIR} --with-devel ${CONF}"
    PMAKE="-j2"
}

#################
function config_seawulf {
    if [[ -z "$LD_LIBRARY_PATH" ]]; then
        echo "#### There seems to be a problem in your environment setting."
        echo "#### You need to set LD_LIBRARY_PATH to include /usr/local/pkg/openmpi/lib."
	exit
    fi
    echo "### Make sure you have the following environment variables set in .bash_profile or .cshrc:"
    echo "export PATH=$PATH:/usr/local/pkg/torque/bin:/usr/local/pkg/openmpi/bin:/usr/local/pkg/HDF4/bin"
    echo "export LD_LIBRARY_PATH=/usr/local/pkg/petsc-2.3.3-p11/lib:/usr/local/pkg/HDF4/lib:/usr/local/pkg/openmpi/lib"
    echo
    echo

    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    export PETSC_DIR=/usr/local/pkg/petsc-2.3.3-p11
    export PETSC_ARCH=linux-gnu-c-debug
    export PETSC_INCLUDE="-I${PETSC_DIR}/include -I${PETSC_DIR}/bmake/${PETSC_ARCH}"
    export PETSC_LIB="-lblas -llapack -L${PETSC_DIR}/lib/${PETSC_ARCH} -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -ldl -lm -lX11"

    if [[ -n "$WITHHDF" ]]; then
        CONF="$CONF --with-hdf=/usr/local/pkg/HDF4"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/nfs/user07/duowang/hdf5-1.8.5"
	#CONF="$CONF --with-hdf5=/usr/local/pkg/HDF5"
    fi
    export CONF="--with-openmpi --with-petsc=${PETSC_DIR} --with-devel ${CONF}"
    PMAKE="-j2"
}

#################
function config_nybluel {
    export PETSC_ARCH=bgl-ibm-opt
    export PETSC_DIR=/bgl/apps/petsc-2.3.3-p8/petsc-2.3.3-p8
    export PETSC_INCLUDE="-I${PETSC_DIR}/include -I${PETSC_DIR}/bmake/${PETSC_ARCH}"
    export HYPRE_DIR=/bgl/apps/petsc-2.3.3-p8/petsc-2.3.3-p8/externalpackages/hypre-2.0.0/bgl-ibm-opt/lib
    export SUPERLU_DIST_DIR=/bgl/apps/petsc-2.3.3-p8/petsc-2.3.3-p8/externalpackages/SuperLU_DIST_2.0-Jan_5_2006/bgl-ibm-opt
    export SUPERLU_DIR=/bgl/apps/petsc-2.3.3-p8/petsc-2.3.3-p8/externalpackages/SuperLU_3.0-Jan_5_2006/bgl-ibm-opt
    export NSS_FILES_DIR=/bgl/BlueLight/ppcfloor/blrts-gnu/powerpc-bgl-blrts-gnu/lib

    export PETSC_LIB="-L${PETSC_DIR}/lib/${PETSC_ARCH} -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lg2c -L/bgl/local/lib -L${HYPRE_DIR} -lHYPRE -L${SUPERLU_DIR} -lsuperlu_3.0 -L${SUPERLU_DIST_DIR} -lsuperlu_dist_2.0  -lc -L${NSS_FILES_DIR} -lnss_files -lnss_dns -lresolv -llapack.rts -lblas.rts -L/usr/lib -ldl -lm"

    export F77_LIBS="-L/opt/ibmcmp/xlsmp/bg/1.7/blrts_lib -L/opt/ibmcmp/xlmass/bg/4.4/blrts_lib -L/opt/ibmcmp/xlf/bg/11.1/blrts_lib -lxlf90 -lxlopt -lxlomp_ser -lxl -lxlfmath -lm -lc -lgcc"

    export CXX="mpixlcxx -DMPICH_IGNORE_CXX_SEEK ${OPTS}"
    export F77="mpixlf77 ${OPTS}"
    export CC="mpixlc ${OPTS}"

    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/bgl/apps/hdf5"
    fi

    export CONF="--with-mpi --with-extra-libs=-lmpich.rts --with-petsc=${PETSC_DIR} $CONF"
    PMAKE="-j8"
}

#################
function config_nybluep {
    export CXX="mpixlcxx_r -DMPICH_IGNORE_CXX_SEEK ${OPTS}"
    export F77="mpixlf77_r ${OPTS}"
    export CC="mpixlc_r ${OPTS}"

    export CONF="--with-mpi --with-extra-libs=-lmpich.cnk $CONF"
    PMAKE="-j8"
}

################
function config_pnnl {
    export CXX="mpiCC ${OPTS} ${COPTS_ICC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_ICC}"
    
    # Specify Petsc path
    export PETSC_ARCH=linux-hpmpi-intel
    export PETSC_DIR=/hptc_cluster/apps/libraries/petsc/2.3.3-p8/Intel
    export PETSC_INCLUDE="-I${PETSC_DIR}/include -I${PETSC_DIR}/bmake/${PETSC_ARCH}"
    export PETSC_LIB="-lblas -L${PETSC_DIR}/lib/${PETSC_ARCH} -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -ldl -lm -L/usr/X11R6/lib64"
    
    ./configure --with-mpi --with-petsc=$PETSC_DIR --with-devel ${CONF}
    PMAKE="-j8"
}

################
function config_abe {
    export CXX="mpicxx -DMPICH_IGNORE_CXX_SEEK -Wcheck ${OPTS} ${COPTS_ICC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_ICC}"

    # Specify Petsc path
    export PETSC_ARCH=abe-intel10-opt
    export PETSC_DIR=/usr/apps/math/petsc/petsc-2.3.3-p7
    export PETSC_INCLUDE="-I${PETSC_DIR}/include -I${PETSC_DIR}/bmake/${PETSC_ARCH}"
    export PETSC_LIB="-lblas -llapack -L${PETSC_DIR}/lib/${PETSC_ARCH} -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -L ${MKL_HOME}/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -lm -L/usr/X11R6/lib64"

    # add szip library
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/apps/hdf/szip/lib
    export scriptotherlibsinc="${scriptotherlibincs} -L/usr/apps/hdf/szip/lib"
    export scriptotherlibs="${scriptotherlibs} -lsz"

    if [[ -n "$WITHHDF" ]]; then
        CONF="$CONF --with-hdf=/usr/apps/hdf/hdf4/v423"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/usr/apps/hdf/phdf5/v185"
    fi

    export CONF="-with-mpich --with-petsc=$PETSC_DIR --with-devel ${CONF}"
    PMAKE="-j8"
}

#################
function config_icestorm {
    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"

    export PETSC_DIR=/home/dealvf/libs/petsc-3.1-p8_ompi_gcc-1.4.3_icestorm
    export PETSC_INCLUDE="-I${PETSC_DIR}/include"

    export PETSC_LIB="-L${PETSC_DIR}/lib -lpetsc -lflapack -lfblas -ldl -lm -L/usr/X11R6/lib64 -lX11"

    if [[ -n "$WITHHDF" ]]; then
        CONF="$CONF --with-hdf=/usr/local/pkg/HDF4"
    fi
    if [[ -n "$WITHHDF5" ]]; then
        CONF="$CONF --with-hdf5=/nfs/user01/duowang/hdf5-1.8.5"
        #CONF="$CONF --with-hdf5=/usr/local/pkg/HDF5"
    fi
    if [[ -n "$WITHGD" ]]; then
        CONF="$CONF --with-gd=/usr"
    fi

    export CONF="--with-openmpi=/usr/local/openmpi/openmpi-1.4.3/gcc-opt --with-petsc=${PETSC_DIR} --with-devel ${CONF}"
    PMAKE="-j2"
}


################# Main code starts here

parse_arguments $*

# Invoke machine specific config function
HOST=`uname -n`;
# Choose proper compilers based on machine
if [[ "${HOST}" == "fenp" ]]; then
    echo "Computer is recognized as NYBlue/P."
    config_nybluep
elif [[ "${HOST}" == "fen" ]]; then
    echo "Computer is recognized as NYBlue/L."
    config_nybluel
elif [[ "${HOST}" == "dontpanic" ]]; then
    echo "Computer is recognized as Vogon."
    config_vogon
elif [[ "${HOST}" == "sirius" ]]; then
    echo "Computer is recognized as Galaxy."
    config_galaxy
elif [[ "${HOST}" == "seawulf" ]]; then
    echo "Computer is recognized as Seawulf."
    config_seawulf
elif [[ "${HOST}" == "cu0login1" ]]; then
    echo "Computer is recognized as cu0login1."
    config_pnnl
elif [[ "${HOST}" == "jiao-pc3" ]]; then
    echo "Computer is recognized as Ubuntu."
    config_ubuntu
elif [[ "${HOST//[0-9]/}" == "honest.ncsa.uiuc.edu" ]]; then
    echo "Computer is recognized as NCSA Abe Linux cluster."
    config_abe
elif [[ "${HOST}" == "lpt1" ]]; then
    echo "Computer is recognized as Gentoo Linux for ORNL Laptop."
    config_gentooORNL
elif [[ "${HOST}" == "flogin1" ]]; then
    echo "Computer is recognized as INL Fission."
    config_fission
elif [[ "${HOST}" == "service0" ]]; then
    echo "Computer is recognized as INL Icestorm."
    config_icestorm
elif [[ "${HOST}" == "tulin.ams.sunysb.edu" ]]; then
    echo "Computer is recognized as Tulin."
    config_tulin
elif [[ "${HOST}" == "hyun.ams.sunysb.edu" ]]; then
    echo "Computer is recognized as Tulin."
    config_hyun
else
    echo "Computer was not recognized. Using generic configure options."
    config_generic
fi

echo "Configuring FronTier with the following commands:"
echo "export CC=\"$CC\""
echo "export CXX=\"$CXX\""
echo "export F77=\"$F77\""
if [[ -n "$PETSC_DIR" ]]; then
    echo "export PETSC_DIR=\"$PETSC_DIR\""
    echo "export PETSC_ARCH=\"$PETSC_ARCH\""
    echo "export PETSC_INCLUDE=\"$PETSC_INCLUDE\""
    echo "export PETSC_LIB=\"$PETSC_LIB\""
fi
echo "autoconf"
echo "./configure $CONF"

# Run autoconf to generate ./configure
autoconf

# Run ./configure.
./configure $CONF

# Finally, invoke make
if [[ -z "$NOMAKE" ]]; then
    make ${PMAKE}
fi
