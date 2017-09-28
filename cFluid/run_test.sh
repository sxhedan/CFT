#!/bin/bash
    #PBS -lnodes=1:ppn=12:towel
    #These two lines link the standard version of openmpi 
    LD_LIBRARY_PATH=/usr/local/pkg/openmpi/lib
    export LD_LIBRARY_PATH 
    #This chooses the mpirun script for our version of openmpi 
    MPIRUN=/usr/local/pkg/openmpi/bin/mpirun
    #These two lines just define the paths for input/output files for  
    #ease of typing 
    HOME=/nfs/scratch/dshe
    CFLUID=$HOME/CFT/FTI_REPO/cFluid
    #edit this line to change input/output file name 
    NAME=in-rt3d
    #Run the job!  
    $MPIRUN -np 8 $CFLUID/cFluid -d 3 -p 2 2 2 -i $CFLUID/$NAME -o $CFLUID/out-test-mp
