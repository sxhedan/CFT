#!/bin/bash
    #PBS -lnodes=1:ppn=12:towel
    #These two lines link the standard version of openmpi 
    LD_LIBRARY_PATH=/usr/local/pkg/openmpi/lib
    export LD_LIBRARY_PATH 
    #This chooses the mpirun script for our version of openmpi 
    MPIRUN=/usr/local/pkg/openmpi/bin/mpirun
    #These two lines just define the paths for input/output files for  
    #ease of typing 
    HOME=/nfs/user01/dshe
    CFLUID=$HOME/FTI_20150807/cFluid
    #edit this line to change input/output file name 
    NAME=in-vst-test
    OUTPUT=/nfs/scratch/dshe/FTI_20150807
    #Run the job!  
    $MPIRUN -np 1 $CFLUID/cFluid -d 3 -p 1 1 1 -i $CFLUID/$NAME -o $OUTPUT/out-vst-test-20150824-sp
