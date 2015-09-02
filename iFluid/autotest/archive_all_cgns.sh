#! /bin/bash 

# Run all autotests' archive.sh scripts:
tests="autotest2d autotest3d"
for i in $tests
do
	cd $i
	echo "Running ./$i/archive_cgns.sh:"
	./archive_cgns.sh
	cd ..
done
