#!/bin/bash
    #PBS -lnodes=6:ppn=12:towel
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
    NAME=in-rt2d-f-400
    #Run the job!  
    $MPIRUN -np 64 /nfs/user01/dshe/FTI_20150807/cFluid/cFluid -d 2 -p 4 16 -i $CFLUID/$NAME -o /nfs/scratch/dshe/FTI_20150807/out-rt2d-f-400-20150810
