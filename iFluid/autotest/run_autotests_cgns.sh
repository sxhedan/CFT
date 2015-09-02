#! /bin/bash

# List of test directories. When a new test set is complete, add the name to
# this list.
tests="autotest2d autotest3d"

rm -f *out

HOST=`uname -n`; 
case "$HOST" in   
	p*) host_name=galaxy ;; 
        t*) host_name=vogon ;;
	w*) host_name=seawulf ;;
	a*) host_name=ncsa ;;
	*) echo "ERROR: Must run run_autotests_cgns.sh in interactive mode." 
	exit 0 ;;
esac 

# Check that $host_name is a recognized machine.
if [ "$host_name" == "galaxy" ] || [ "$host_name" == "seawulf" ] || \
	[ "$host_name" == "ncsa" ] || [ "$host_name" == "vogon" ]
then
	echo "Started run_autotests_cgns.sh on $host_name:"
	date
	echo ""
else
	echo "ERROR: The host machine is not recognized."
   	exit 0
fi

# Run tests
for i in $tests
do
	cd $i
	./${i}_cgns.sh $host_name > ../${i}_cgns.out 2>&1
	cd ..
	echo "Completed $i:"
	date
	echo ""
	echo "Note: Please refer to $i_cgns.out for test results."
	echo "Please refer to ./$i/out/test_name/test_name.out"
	echo "for individual test outputs."
	echo ""
done
echo "Finished run_autotests_cgns.sh on $host_name:"
date
echo ""
