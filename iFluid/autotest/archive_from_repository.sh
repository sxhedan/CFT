#! /bin/bash

archive2d="in-bubble in-drop in-emerge-bubble in-flow2d in-kh2d in-rt2d"
archive3d="in-rt3d in-rt3dBC in-stokes3d in-stokes3dBC"

# Download archive2d files  

cd autotest2d

if [ ! -d "archive" ]
then
    mkdir archive
fi
rm -rf archive/*
cd archive

echo "downloading 2d archived results..."
for i in ${archive2d}
do
    wget -q --user=sitsec --password=FronTier \
	http://localhost:8080/trac/raw-attachment/wiki/Autotest/${i}.tar.gz
done

# Extract archive2d files  

echo "extracting 2d archived results..."
for i in ${archive2d}
do
#tar xvfz ${i}.tar.gz > /dev/null
tar xvf ${i}.tar.gz > /dev/null         
done
cd ..
rm -f archive/*.tar.gz
cd ..

# Download archive3d files 

cd autotest3d

if [ ! -d "archive" ]
then
    mkdir archive
fi
rm -rf archive/*
cd archive

echo "downloading 3d archived results..."
for i in ${archive3d}
do
    wget -q --user=sitsec --password=FronTier \
	http://localhost:8080/trac/raw-attachment/wiki/Autotest/${i}.tar.gz
done

# Extract archive3d files 

echo "extracting 3d archived results..."
for i in ${archive3d}
do
#tar xvfz ${i}.tar.gz > /dev/null
tar xvf ${i}.tar.gz > /dev/null        
done
cd ..
rm -f archive/*.tar.gz
cd ..