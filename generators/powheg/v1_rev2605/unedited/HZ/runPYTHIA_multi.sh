#!/bin/bash


prg=../main-PYTHIA-lhef


#----------------------------------------------------
for mode in nohad
do

    cat powheg.input-save > powheg.input
    
#    echo pythiatune $tune >> powheg.input

#    case $mode in
#	nohad) echo nohad 1 >>  powheg.input ;;
#	nompi) echo nompi 1 >>  powheg.input ;;
#    esac

    for file in pwgevents-*.lhe
    do
	log=log`echo $file | sed 's/pwgevents// ; s/lhe//'`log
	echo $file | $prg > $log &
    done

    wait
    
    for W in W1 W2 W3 W4 W5 W6 W7
    do
	(echo 1 ; ls -c1 pwgPOWHEG+PYTHIA-output-00*-$W.top ; echo "") | mergedata
	
	mv fort.12 pwgPOWHEG+PYTHIA-output-$W.top
    done
    (echo 4 ; ls -c1 pwgPOWHEG+PYTHIA-output-W[1-7].top ; echo "") | mergedata
    mv fort.12 pwgPOWHEG+PYTHIA-output-max.top
    (echo 5 ; ls -c1 pwgPOWHEG+PYTHIA-output-W[1-7].top ; echo "") | mergedata
    mv fort.12 pwgPOWHEG+PYTHIA-output-min.top
done


for W in W1 W2 W3 W4 W5 W6 W7
do
    case $W in
	W1) ext=11;;
	W2) ext=HH;;
	W3) ext=H1;;
	W4) ext=1H;;
	W5) ext=21;;
	W6) ext=12;;
	W7) ext=22;;    
    esac    
    mv pwgPOWHEG+PYTHIA-output-$W.top pwgPOWHEG+PYTHIA-output-$ext.top
done

