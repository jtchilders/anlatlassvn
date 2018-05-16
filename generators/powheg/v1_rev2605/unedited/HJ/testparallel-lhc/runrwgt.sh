#!/bin/bash

function char {
	case $1 in
	    [1-9]) echo 000$1 ;;
	    [1-9][0-9]) echo 00$1 ;;
	    [1-9][0-9][0-9]) echo 0$1 ;;
	    *) echo invalid ; exit -1 ;;
	esac
}


for iscales in {1..7}
do
    case $iscales in
	1) facscfact=1   ; renscfact=1   ;;
	2) facscfact=0.5 ; renscfact=0.5 ;;
	3) facscfact=0.5 ; renscfact=1   ;;
	4) facscfact=1   ; renscfact=0.5 ;;
	5) facscfact=2   ; renscfact=1   ;;
	6) facscfact=1   ; renscfact=2   ;;
	7) facscfact=2   ; renscfact=2   ;;
    esac

    cat  powheg.input-save | sed "s/storeinfo_rwgt/compute_rwgt/ ; s/facscfact.*/facscfact $facscfact/ ; s/renscfact.*/renscfact $renscfact/ ; s/parallelstage.*/parallelstage 4/ ; s/xgriditeration/#xgriditeration/; s/fastbtlbound/#fastbtlbound/" > powheg.input
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




