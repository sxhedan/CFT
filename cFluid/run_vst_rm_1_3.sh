#!/bin/bash
    #PBS -lnodes=3:ppn=12:towel
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
    NAME=in-vst-rm-1-3
    OUTPUT=/nfs/scratch/dshe/FTI_20150807
    #Run the job!  
    $MPIRUN -np 32 $CFLUID/cFluid -d 3 -p 2 2 8 -i $CFLUID/$NAME -o $OUTPUT/out-vst-rm-1-3-20150821
