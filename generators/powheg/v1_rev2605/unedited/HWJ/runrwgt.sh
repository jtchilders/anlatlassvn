#!/bin/bash

function char {
	case $1 in
	    1)  echo 0001 ;;
	    2)  echo 0002 ;;
	    3)  echo 0003 ;;
	    4)  echo 0004 ;;
	    5)  echo 0005 ;;
	    6)  echo 0006 ;;
	    7)  echo 0007 ;;
	    8)  echo 0008 ;;
	    9)  echo 0009 ;;
	    10) echo 0010 ;;
	esac
}


for iscales in {1..7}
do
    case $iscales in
	1) facscfact=0.5 ; renscfact=0.5 ;;
	2) facscfact=0.5 ; renscfact=1   ;;
	3) facscfact=1   ; renscfact=0.5 ;;
	4) facscfact=2   ; renscfact=1   ;;
	5) facscfact=1   ; renscfact=2   ;;
	6) facscfact=2   ; renscfact=2   ;;
	7) facscfact=1   ; renscfact=1   ;;
    esac

    cat  powheg.input-save | sed "s/storeinfo_rwgt/compute_rwgt/ ; s/facscfact.*/facscfact $facscfact/ ; s/renscfact.*/renscfact $renscfact/ ; s/parallelstage/#parallelstage/ ; s/xgriditeration/#xgriditeration/; s/fastbtlbound/#fastbtlbound/" > powheg.input
    for i in {1..8}
    do
	ch=`char $i`

	../pwhg_main <<EOF > run-$i.log 2>&1 &
$i
pwgevents-$ch.lhe
EOF
	    
    done

    wait
	
    for i in {1..8}
    do
	ch=`char $i`
	\mv -f pwgevents-rwgt-$ch.lhe pwgevents-$ch.lhe 
	
    done
    
done




