#! /bin/bash 

# Run all autotests' archive.sh scripts:
tests="autotest2d autotest3d"
for i in $tests
do
	cd $i
	echo "Running ./$i/archive.sh:"
	./archive.sh
	cd ..
done
