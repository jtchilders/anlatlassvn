#!/bin/bash

> Timings.txt


# First compile the pwhg_main executable in the ../ directory
#

# two stages of importance sampling grid calculation
for igrid in {1..2}
do

(echo -n st1 xg$igrid ' ' ; date ) >> Timings.txt

cat powheg.input-save | sed "s/xgriditeration.*/xgriditeration $igrid/ ; s/parallelstage.*/parallelstage 1/" > powheg.input

for i in {1..48}
do
echo $i | ../pwhg_main > run-st1-xg$igrid-$i.log 2>&1 &
done
wait

done



# compute NLO and upper bounding envelope for underlying born comfigurations
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 2/ ' > powheg.input
(echo -n st2 ' ' ; date ) >> Timings.txt
for i in {1..48}
do
echo $i | ../pwhg_main > run-st2-$i.log 2>&1 &
done
wait


# compute upper bounding coefficients for radiation
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 3/' > powheg.input
(echo -n st3 ' ' ; date ) >> Timings.txt
for i in {1..48}
do
echo $i | ../pwhg_main > run-st3-$i.log 2>&1 &
done
wait



# generate events 
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 4/' > powheg.input
(echo -n st4 ' ' ; date ) >> Timings.txt
for i in {1..48}
do
echo $i | ../pwhg_main > run-st4-$i.log 2>&1 &
done
wait

(echo -n end ' ' ; date ) >> Timings.txt

# Now all events are available. This is an example of using the reweighting
# feature. We compute reweighting information for the 7 points scale variation

function char {
	case $1 in
	    [1-9])  echo 000$1 ;;
	    [1-9][0-9])  echo 00$1 ;;
	esac
}


for iscales in {1..6}
do
    case $iscales in
	1) facscfact=0.5 ; renscfact=0.5 ;;
	2) facscfact=0.5 ; renscfact=1   ;;
	3) facscfact=1   ; renscfact=0.5 ;;
	4) facscfact=2   ; renscfact=1   ;;
	5) facscfact=1   ; renscfact=2   ;;
	6) facscfact=2   ; renscfact=2   ;;
    esac
# must be at the parallel stage 4, and the files generated at the previous stages
# must be present.
    cat  powheg.input-save | sed "s/parallelstage.*/parallelstage 4/ ; s/storeinfo_rwgt/compute_rwgt/ ; s/facscfact.*/facscfact $facscfact/ ; s/renscfact.*/renscfact $renscfact/ " > powheg.input
    for i in {1..48}
    do
	ch=`char $i`

	../pwhg_main <<EOF > run-$i.log 2>&1 &
$i
pwgevents-$ch.lhe
EOF
	    
    done

    wait
	
    for i in {1..48}
    do
	ch=`char $i`
	\mv -f pwgevents-rwgt-$ch.lhe pwgevents-$ch.lhe 
	
    done
    
done




