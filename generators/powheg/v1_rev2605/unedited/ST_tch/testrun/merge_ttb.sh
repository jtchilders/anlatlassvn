#!/bin/bash

# 2 input files and total number of events
in1=$1
in2=$2
in3=$3

if [ $# -ne 3 ] ; then
    echo ' SYNTAX ERROR'
    echo '  usage: sh merge_ttb.sh [arg1] [arg2] [arg3]'
    echo "  where [arg1] and [arg2] have to be the t and the tb input files' prefixes"
    echo '  and [arg3] is the total number of t+tb events'
    echo '  numevts entries in [arg1] and [arg2] files are NOT relevant'
    echo $1 $2 $3
    exit
fi

# compile needed executables
PWD=$(pwd)
#echo $PWD
echo ' * compile executable(s) needed by merging procedure:'
(cd ..; make merge_t_tb >tmp.dat 2>&1; make pwhg_main >tmp.dat 2>&1; rm -f tmp.dat ; cd $PWD)
echo '   merging executable(s) compiled'
echo ' * start merging procedure'

#!: option to produce more events if grids already present?

if [ $in3 -gt 0 ] ; then
echo ' * CREATION OF t AND tbar SAMPLES'

file_t=$1'-powheg.input'
file_tb=$2'-powheg.input'

# check that input files are present
if [ ! -f $file_t ] ; then
    echo " Error: $file_t NOT present"
    exit
fi
if [ ! -f $file_tb ] ; then
    echo " Error: $file_tb NOT present"
    exit
fi

# it's responsibility of the user to check that the 2 input files
# are consistent

# temporary input files with numevts=1
cp $file_t  temp_t-powheg.input
cp $file_tb temp_tb-powheg.input
string_t='numevts 1'
cat temp_t-powheg.input | sed -e '/numevts/ cZZZ' > tmp.dat
cat tmp.dat | sed 's/ZZZ/'"$string_t"'/g' > temp_t-powheg.input
string_tb='numevts 1'
cat temp_tb-powheg.input | sed -e '/numevts/ cZZZ' > tmp.dat
cat tmp.dat | sed 's/ZZZ/'"$string_tb"'/g' > temp_tb-powheg.input


# run pwhg_main
is_t_present='N'
is_tb_present='N'
if [ -e temp_t-events.lhe ] ; then
    echo '   auxiliary file temp_t-events.lhe already present. Use it? [Y/N]'
    read is_t_present
    if [ $is_t_present = 'Y' ] ; then
	echo '   using temp_t-events.lhe'
    elif [ $is_t_present = 'N' ] ; then
	echo '   Then please remove or rename the file temp_t-events.lhe'
	exit
    else
	echo '   Invalid value, only Y or N accepted'
	exit
    fi
fi
if [ -e temp_tb-events.lhe ] ; then
    echo '   auxiliary file temp_tb-events.lhe already present. Use it? [Y/N]'
    read is_tb_present
    if [ $is_tb_present = 'Y' ] ; then
	echo '   using temp_tb-events.lhe'
    elif [ $is_tb_present = 'N' ] ; then
	echo '   Then please remove or rename the file temp_tb-events.lhe'
	exit
    else
	echo '   Invalid value, only Y or N accepted'
	exit
    fi
fi

PID_t="-1"
#echo $PID_t
if [ $is_t_present = 'N' ] ; then
echo temp_t | ../pwhg_main > temp_t.out &
PID_t=$!
#echo $PID_t
fi

PID_tb="-1"
#echo $PID_tb
if [ $is_tb_present = 'N' ] ; then
echo temp_tb | ../pwhg_main > temp_tb.out &
PID_tb=$!
#echo $PID_tb
fi

if [ $PID_t -ge 0 ] ; then
#    echo 'waiting t'
    wait $PID_t
fi
if [ $PID_tb -ge 0 ] ; then
#    echo 'waiting tb'
    wait $PID_tb
fi

echo ' * t and tbar grids generation performed:'
echo '   shell outputs stored in temp_t.out and temp_tb.out'

# run program that compare top and antitop total cross section
echo $3 | ../merge_t_tb > totals.dat
echo ' * totals.dat file created'

# parse totals.dat to extract the number of needed events
grep '1nev_t' totals.dat | awk '{print $1 > "nevt"}'
grep '2nev_tb' totals.dat | awk '{print $1 > "nevtb"}'
nev_t=$(cat nevt)
nev_tb=$(cat nevtb)
rm -f nevt nevtb

# update 2 input files with proper number of events
# and run pwhg_main to produce the t and the tbar samples
# If event files are already present, they are assumed to
# be already OK and ready for the merging procedure.
if [ $is_t_present = 'N' ] ; then
    rm -f temp_t-events.lhe
    string_t='numevts '$nev_t
    cat temp_t-powheg.input | sed -e '/numevts/ cZZZ' > tmp.dat
    cat tmp.dat | sed 's/ZZZ/'"$string_t"'/g' > temp_t-powheg.input
    echo temp_t | ../pwhg_main > temp_t_run.out &
fi

if [ $is_tb_present = 'N' ] ; then
    rm -f temp_tb-events.lhe
    string_tb='numevts '$nev_tb
    cat temp_tb-powheg.input | sed -e '/numevts/ cZZZ' > tmp.dat
    cat tmp.dat | sed 's/ZZZ/'"$string_tb"'/g' > temp_tb-powheg.input
    echo temp_tb | ../pwhg_main > temp_tb_run.out
fi
rm -f tmp.dat

echo ' * t and tbar events generation performed:'
echo '   shell outputs stored in temp_t_run.out and temp_tb_run.out'
echo ' * t AND tbar SAMPLES CREATED'

else
echo ' * MERGING t AND tbar SAMPLES'

if [ -e temp_t-events.lhe ] ; then
    echo '   temp_t-events.lhe will be used'
else
    echo ' Error: temp_t-events.lhe missing'
    exit
fi
if [ -e temp_tb-events.lhe ] ; then
    echo '   temp_tb-events.lhe will be used'
else
    echo ' Error: temp_tb-events.lhe missing'
    exit
fi

# NB: special value for iprup (666)
echo $3 | ../merge_t_tb &
PID_merge=$!

wait $PID_merge
echo 'numevts     1' > t_tb_sample-powheg.input
rm -f temp_t*top #tmp.dat
echo ' * OUTPUT FILE IS t_tb_sample-events.lhe'

fi