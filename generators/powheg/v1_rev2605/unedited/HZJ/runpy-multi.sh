#!/bin/bash


prg=../main-PYTHIA-lhef

for label in 11 1H H1 HH 12 21 22
do
case $label in
11) rfac=1 ; ffac=1 ;;
1H) rfac=1 ; ffac=0.5 ;;
H1) rfac=0.5 ; ffac=1 ;;
HH) rfac=0.5 ; ffac=0.5 ;;
12) rfac=1 ; ffac=2 ;;
21) rfac=2 ; ffac=1 ;;
22) rfac=2 ; ffac=2 ;;
esac

cat powheg.input-save | sed "s/renscfact.*/renscfact $rfac/ ; s/facscfact.*/facscfact $ffac/ " > powheg.input

#echo nohad 1 >>  powheg.input

for file in pwgevents-*.lhe
do
log=log`echo $file | sed 's/pwgevents// ; s/lhe//'`log
echo $file | $prg > $log &
done

wait

(echo 1 ; ls -c1 pwg-00*POWHEG+PYTHIA-output.top ; echo "") | mergedata

mv fort.12 pwgPOWHEG+PYTHIA-output-$label.top

done
