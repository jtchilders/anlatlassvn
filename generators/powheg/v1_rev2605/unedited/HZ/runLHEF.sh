#!/bin/bash


prg=../lhef_analysis

cp powheg.input-save powheg.input

for file in pwgevents-*.lhe
do
    log=log`echo $file | sed 's/pwgevents// ; s/lhe//'`log
    echo $file | $prg > $log &
done

wait

for W in W1 W2 W3 W4 W5 W6 W7
do
    (echo 1 ; ls -c1 pwgLHEF_analysis-00*-$W.top ; echo "") | mergedata   
    mv fort.12 pwgLHEF_analysis-$W.top
done

(echo 4 ; ls -c1 pwgLHEF_analysis-W[1-7].top ; echo "") | mergedata
mv fort.12 pwgLHEF_analysis-max.top
(echo 5 ; ls -c1 pwgLHEF_analysis-W[1-7].top ; echo "") | mergedata
mv fort.12 pwgLHEF_analysis-min.top

	# 1) facscfact=1   ; renscfact=1   ;;
	# 2) facscfact=0.5 ; renscfact=0.5 ;;
	# 3) facscfact=0.5 ; renscfact=1   ;;
	# 4) facscfact=1   ; renscfact=0.5 ;;
	# 5) facscfact=2   ; renscfact=1   ;;
	# 6) facscfact=1   ; renscfact=2   ;;
	# 7) facscfact=2   ; renscfact=2   ;;

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
    mv pwgLHEF_analysis-$W.top pwgLHEF_analysis-$ext.top
done

