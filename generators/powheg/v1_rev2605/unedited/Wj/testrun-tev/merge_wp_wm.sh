#!/bin/bash

# 2 input files and total number of events
in1=$1
in2=$2
in3=$3

if [ $# -ne 3 ] ; then
    echo ' SYNTAX ERROR'
    echo '  usage: sh merge_wpwm.sh [arg1] [arg2] [arg3]'
    echo "  where [arg1] and [arg2] have to be the wp and the wm input files' prefixes"
    echo '  and [arg3] is the total number of wp + wm events'
    echo '  numevts entries in [arg1] and [arg2] files are NOT relevant'
    echo $1 $2 $3
    exit
fi

# compile needed executables
PWD=$(pwd)
#echo $PWD
echo ' * compile executable(s) needed by merging procedure:'
(cd ..; make merge_wp_wm >tmp.dat 2>&1; make pwhg_main >tmp.dat 2>&1; rm -f tmp.dat ; cd $PWD)
echo '   merging executable(s) compiled'
echo ' * start merging procedure'

#!: option to produce more events if grids already present?

if [ $in3 -gt 0 ] ; then
    echo ' * CREATION OF wp AND wm SAMPLES'

    file_wp=$1'-powheg.input'
    file_wm=$2'-powheg.input'

# check that input files are present
    if [ ! -f $file_wp ] ; then
	echo " Error: $file_wp NOT present"
	exit
    fi
    if [ ! -f $file_wm ] ; then
	echo " Error: $file_wm NOT present"
	exit
    
    fi

# it's responsibility of the user to check that the 2 input files
# are consistent

# temporary input files with numevts=1
    cp $file_wp  temp_wp-powheg.input
    cp $file_wm  temp_wm-powheg.input
    string_wp='numevts 1'
    cat temp_wp-powheg.input | sed -e '/numevts/ cZZZ' > tmp.dat
    cat tmp.dat | sed 's/ZZZ/'"$string_wp"'/g' > temp_wp-powheg.input
    string_wm='numevts 1'
    cat temp_wm-powheg.input | sed -e '/numevts/ cZZZ' > tmp.dat
    cat tmp.dat | sed 's/ZZZ/'"$string_wm"'/g' > temp_wm-powheg.input


# run pwhg_main
    is_wp_present='N'
    is_wm_present='N'
    if [ -e temp_wp-events.lhe ] ; then
	echo '   auxiliary file temp_wp-events.lhe already present. Use it? [Y/N]'
	read is_wp_present
	if [ $is_wp_present = 'Y' -o $is_wp_present = 'y' ] ; then
	    echo '   using temp_wp-events.lhe'
	elif [ $is_wp_present = 'N' -o $is_wp_present = 'n' ] ; then
	    echo '   Then please remove or rename the file temp_wp-events.lhe'
	    exit
	else
	    echo '   Invalid value, only Y or N accepted'
	    exit
	fi
    fi
    if [ -e temp_wm-events.lhe ] ; then
	echo '   auxiliary file temp_wm-events.lhe already present. Use it? [Y/N]'
	read is_wm_present
	if [ $is_wm_present = 'Y' -o $is_wm_present = 'y' ] ; then
	    echo '   using temp_wm-events.lhe'
	elif [ $is_wm_present = 'N' -o $is_wm_present = 'n' ] ; then
	    echo '   Then please remove or rename the file temp_wm-events.lhe'
	    exit
	else
	    echo '   Invalid value, only Y or N accepted'
	    exit
	fi
    fi

    PID_wp="-1"
#echo $PID_wp
    if [ $is_wp_present = 'N' -o $is_wp_present = 'n' ] ; then
	echo temp_wp | ../pwhg_main > temp_wp.out &
	PID_wp=$!
#echo $PID_wp
    fi

    PID_wm="-1"
#echo $PID_wm
    if [ $is_wm_present = 'N' -o $is_wm_present = 'n' ] ; then
	echo temp_wm | ../pwhg_main > temp_wm.out &
	PID_wm=$!
#echo $PID_wm
    fi

    if [ $PID_wp -ge 0 ] ; then
#    echo 'waiting t'
	wait $PID_wp
    fi
    if [ $PID_wm -ge 0 ] ; then
#    echo 'waiting tb'
	wait $PID_wm
    fi

    echo ' * wp and wm grids generation performed:'
    echo '   shell outputs stored in temp_wp.out and temp_wm.out'

# run program that compare wp and wm total cross section
    echo $3 | ../merge_wp_wm > totals.dat
    echo ' * totals.dat file created'

# parse totals.dat to extract the number of needed events
    grep '1nev_wp' totals.dat | awk '{print $1 > "nevwp"}'
    grep '2nev_wm' totals.dat | awk '{print $1 > "nevwm"}'
    nev_wp=$(cat nevwp)
    nev_wm=$(cat nevwm)
    rm -f nevwp nevwm

# update 2 input files with proper number of events
# and run pwhg_main to produce the wp and the wm samples
# If event files are already present, they are assumed to
# be already OK and ready for the merging procedure.
    if [ $is_wp_present = 'N' -o $is_wp_present = 'n' ] ; then
	rm -f temp_wp-events.lhe
	string_wp='numevts '$nev_wp
	cat temp_wp-powheg.input | sed -e '/numevts/ cZZZ' > tmp.dat
	cat tmp.dat | sed 's/ZZZ/'"$string_wp"'/g' > temp_wp-powheg.input
	echo temp_wp | ../pwhg_main > temp_wp_run.out &
    fi

    if [ $is_wm_present = 'N' -o $is_wm_present = 'n' ] ; then
	rm -f temp_wm-events.lhe
	string_wm='numevts '$nev_wm
	cat temp_wm-powheg.input | sed -e '/numevts/ cZZZ' > tmp.dat
	cat tmp.dat | sed 's/ZZZ/'"$string_wm"'/g' > temp_wm-powheg.input
	echo temp_wm | ../pwhg_main > temp_wm_run.out
    fi
    rm -f tmp.dat

    echo ' * wp and wm events generation performed:'
    echo '   shell outputs stored in temp_wp_run.out and temp_wm_run.out'
    echo ' * wp AND wm SAMPLES CREATED'
    echo ' '
    echo ' NOW CALL sh merge_wpwm.sh [arg1] [arg2] [-arg3]'
    echo ' with the SAME [arg1] [arg2] and the opposite'
    echo ' [arg3] IN ORDER TO PERFORM THE MERGING'
    echo ' '
else
    echo ' * MERGING wp AND wm SAMPLES'

    if [ -e temp_wp-events.lhe ] ; then
	echo '   temp_wp-events.lhe will be used'
    else
	echo ' Error: temp_wp-events.lhe missing'
	exit
    fi
    if [ -e temp_wm-events.lhe ] ; then
	echo '   temp_wm-events.lhe will be used'
    else
	echo ' Error: temp_wm-events.lhe missing'
	exit
    fi

# NB: special value for iprup (666)
    echo $3 | ../merge_wp_wm &
    PID_merge=$!

    wait $PID_merge
    echo 'numevts     1' > wp_wm_sample-powheg.input
    rm -f temp_wp*top #tmp.dat
    echo ' * OUTPUT FILE IS wp_wm_sample-events.lhe'
fi