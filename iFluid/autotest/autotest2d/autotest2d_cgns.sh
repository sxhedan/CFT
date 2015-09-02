#! /bin/bash
# autotest2d.sh - runs the tests provided by the iFluid executable.

HOST=$1

# Check if iFluid executable has been built.
if [ ! -e "../../iFluid" ]
then
    echo "$0 ERROR: iFluid executable cannot be found." >&2
    exit 1
fi

# Create out/ if it doesn't exist
if [ ! -d "out" ]
then
    mkdir out
fi

echo "========================================================================="
echo ".$HOST."
echo "autotest2d - performs tests provided in the \"in\" directory, using the "
echo "executable \"iFluid\"."
echo "========================================================================="
echo ""
echo "Execution tests - These tests assess only if a run has completed or not."
echo "The output is grepped and then the results are reported as below."
echo ""
echo "passed  = test ran to end without crashing."
echo "warning = passed but warnings found in output."
echo "FAILED  = test crashed."
echo ""
echo "Regression tests - These test, using diff, the output against an archived"
echo "result. Floating point rounding is a problem for these tests and may"
echo "cause them to fail superfluously. If the archived results were generated"
echo "on a different machine this may be a problem. The results are reported as"
echo "shown below."
echo ""
echo "passed  = output identical to archive"
echo "FAILED  = diff found differences"
echo "skipped = a regression test is not possible for this test"
echo ""
echo "-------------------------------------------------------------------------"
echo "Input name               Execution           Regression"
echo "-------------------------------------------------------------------------"

exepass=0
exewarn=0
exefail=0
regpass=0
regfail=0
regskipped=0

rm -rf in/*.sc
rm -rf in/*~
rm -rf out/*

# Finds list of input files
FILES=`ls in`
PATH_=`pwd`
if [[ -z "$MEXCGNS" ]]; then MEXCGNS=$HOME/mexcgns; fi

# Set max timestep value
tsmx="200"

# Launch all jobs
for i in $FILES
do
    mkdir out/$i
    cd out/$i
    
    if [ "$HOST" == "galaxy" ] || [ "$HOST" == "seawulf" ] || [ "$HOST" == "vogon" ]

    then
	mpirun -np 4 $PATH_/../../iFluid -d 2 -p 2 2 -i $PATH_/in/$i \
	    > $i.out 2>/dev/null
    elif [ "$HOST" == "ncsa" ]
    then
	mpirun $PATH_/../../iFluid -d 2 -i $PATH_/in/$i -np 4 \
	    > $i.out 2>/dev/null 
    else
	echo "The hostname provided is not recognized."
    	exit 0
    fi
    cd ../..
done

# Examine output
cores="0 1 2 3"
for i in $FILES
do
	# GOOD=""
	# WARN=""
	# FAIL=""
    printf "%-25s" "$i"
	# for j in $cores
        # do
    output=`cat out/$i/$i.out | grep CLEAN_UP`
    if [ "$HOST" == "ncsa" ]
    then
	goodresult="		CLEAN_UP, error = 0"
    else
	goodresult="		CLEAN_UP, error = 0
		CLEAN_UP, error = 0
		CLEAN_UP, error = 0
		CLEAN_UP, error = 0"
    fi
    if [ "$output" = "$goodresult" ]
    then
	warning=`cat out/$i/$i.out | grep WARNING`
	if [ "$warning" = "$NULL" ]
	then
				# GOOD="true"
	    printf "%-20s" "passed"
	    ((exepass++))
	else
			 	# WARN="true"   
	    printf "%-20s" "warning"
	    ((exewarn++))
	fi
    fi
    
    if [ "$output" != "$goodresult" ]
    then
			# FAIL="true"
	printf "%-20s" "FAILED"
	echo ""
	((exefail++))
	((regskipped++))
	continue
    fi
	# done
    
    
        # GOOD, WARN and FAIL might equal true, therefore put in this order.
        # if [ "$FAIL" = "true" ]
        # then
        #        printf "%-20s" "FAILED"
        #        echo ""
        #        ((exefail++))
        #        ((regskipped++))
        #        continue
        # elif [ "$WARN" = "true" ]
        # then
        #        printf "%-20s" "warning"
        #        ((exewarn++))
	# elif [ "$GOOD" = "true" ]
	# then
        #        printf "%-20s" "passed"
        #        ((exepass++))
	# fi
	
	# Compare results to archive
    if [ "$HOST" == "galaxy" ] || [ "$HOST" == "seawulf" ] || [ "$HOST" == "vogon" ]

    then
	cd
	cd $PATH_
	GOOD=""
	SKIP=""
	FAIL=""
	if [ ! -d "archive/${i}" ]
        then
            cd ..
	    echo "FAILED: no archive"
            cd autotest2d
	    ((regskipped++))
	    SKIP="true"
        else
	    for j in $cores
	    do	    
		regression1=`octave -q --eval "addpath('$MEXCGNS'); startup_mexcgns;\
		compare_cgns('$PATH_/out/${i}/intfc/cgns.ts0000${tsmx}-nd000${j}/2d-intfc.cgns',\
		'$PATH_/archive/${i}/intfc/cgns.ts0000${tsmx}-nd000${j}/2d-intfc.cgns');"`
		regression2=`octave -q --eval "addpath('$MEXCGNS'); startup_mexcgns; \
		    compare_cgns('$PATH_/out/${i}/intfc/cgns.ts0000${tsmx}-nd000${j}/velo.cgns', \
		    '$PATH_/archive/${i}/intfc/cgns.ts0000${tsmx}-nd000${j}/velo.cgns');"`
		if [ "$regression1" == "Files are similar." ]  \
		    && [ "$regression2" = "Files are similar." ]
		then
		    GOOD="true"
				# echo "passed"
				# ((regpass++))
		#	elif [ "$regression" = "skipped" ]
		#	then
			#	SKIP="true"
				# echo "skipped"
				# ((regskipped++))
		else
		    FAIL="true"
				# echo "FAILED"
				# ((regfail++))
		fi
	    done
	fi
	    
		# GOOD, WARN and FAIL might equal true. 
		# Therefore put in this order.
   	if [ "$FAIL" = "true" ]
	then
	    echo "FAILED"
	    ((regfail++))
	elif [ "$GOOD" = "true" ]
	then
	    echo "passed"
		((regpass++))
	fi
    else
	regression=`diff --normal \
	    out/$i/intfc/vtk.* \
	    archive/$i/vtk.*`
        if [ "$regression" = "$NULL" ]
        then
	    echo "passed"
	    ((regpass++))
        elif [ "$regression" = "skipped" ]
        then
	    echo "skipped"
	    ((regskipped++))
        else
	    echo "FAILED"
	    ((regfail++))
        fi
    fi
done

rm -rf in/*.sc
rm -rf intfc*
rm -rf FT-test*

echo ""
echo "Execution Tests                 Regression Tests"
echo "---------------                 ----------------"
printf "Number passed  = %-15s" "$exepass"             
echo "Number passed  = $regpass"
printf "Number warning = %-15s" "$exewarn"
echo "Number skipped = $regskipped"
printf "Number FAILED  = %-15s" "$exefail"               
echo "Number FAILED  = $regfail"

exit 0
