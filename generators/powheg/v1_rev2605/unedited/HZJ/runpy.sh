#!/bin/bash


prg=../main-PYTHIA-lhef

cp powheg.input-save powheg.input

#echo 'nohad 1' >> powheg.input
#echo 'changescalup 1' >> powheg.input


for file in pwgevents-*.lhe
do
log=log`echo $file | sed 's/pwgevents// ; s/lhe//'`log
echo $file | $prg > $log &
done

wait

(echo 1 ; ls -c1 pwg-00*POWHEG+PYTHIA-output.top ; echo "") | mergedata

mv fort.12 pwgPOWHEG+PYTHIA-output.top


