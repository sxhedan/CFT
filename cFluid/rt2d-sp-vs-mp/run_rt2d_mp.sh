#!/bin/bash
    #PBS -lnodes=2:ppn=12:towel
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
    NAME=in-rt2d
    #Run the job!  
    $MPIRUN -np 16 $CFLUID/cFluid -d 2 -p 2 8 -i $CFLUID/$NAME -o /nfs/scratch/dshe/FTI_20150807/out-rt2d-2x8-20150818
