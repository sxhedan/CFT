#! /bin/bash

# Create archive/ if it doesn't exist
if [ ! -d "archive" ]
then
	mkdir archive
fi

rm -rf archive/*
mkdir archive/intfc

# Copy files to archive
list=`ls out`
echo "Copying .vtk files from out/ to archive/..."
for i in $list
do
      	mkdir archive/intfc/$i
	cp -r out/$i/intfc/vtk.* archive/intfc/$i/
done

exit 0
