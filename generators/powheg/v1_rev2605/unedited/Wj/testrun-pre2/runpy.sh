#!/bin/bash


prg=../main-PYTHIA-lhef

if [ z$1 = z ]
then
    tune=320
else
    tune=$1
fi
label=-$tune

#----------------------------------------------------
for mode in nohad
do

    cat powheg.input-save > powheg.input
    
    echo pythiatune $tune >> powheg.input

    case $mode in
	nohad) echo nohad 1 >>  powheg.input ;;
	nompi) echo nompi 1 >>  powheg.input ;;
    esac

    for file in pwgevents-*.lhe
    do
	log=log`echo $file | sed 's/pwgevents// ; s/lhe//'`log
	echo $file | $prg > $log &
    done

    wait
    
    for W in W1 W2 W3 W4 W5 W6 W7
    do
	(echo 1 ; ls -c1 pwgPOWHEG+PYTHIA-output-00*-$W.top ; echo "") | mergedata
	
	mv fort.12 pwgPOWHEG+PYTHIA-output-$mode$label-ZJ-$W.top
    done
    (echo 4 ; ls -c1 pwgPOWHEG+PYTHIA-output-$mode$label-ZJ-W[1-7].top ; echo "") | mergedata
    mv fort.12 pwgPOWHEG+PYTHIA-output-$mode$label-ZJ-Whigh.top
    (echo 5 ; ls -c1 pwgPOWHEG+PYTHIA-output-$mode$label-ZJ-W[1-7].top ; echo "") | mergedata
    mv fort.12 pwgPOWHEG+PYTHIA-output-$mode$label-ZJ-Wlow.top
done
