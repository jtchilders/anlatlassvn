#!/bin/bash


prg=../lhef_analysis

cp powheg.input-save powheg.input

for file in pwgevents-*.lhe
do
log=log`echo $file | sed 's/pwgevents// ; s/lhe//'`log
echo $file | $prg > $log &
done

wait

(echo 1 ; ls -c1 pwgLHEF_analysis-*.top ; echo "") | mergedata

mv fort.12 pwgLHEF_analysis.top


