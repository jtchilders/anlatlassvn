#!/bin/bash

prg=../pwhg_main1


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
    
    cat  powheg.input-save | sed "s/storeinfo_rwgt/compute_rwgt/ ; s/facscfact.*/facscfact $facscfact/ ; s/renscfact.*/renscfact $renscfact/ ; s/parallelstage/#parallelstage/ ; s/xgriditeration/#xgriditeration/; s/fastbtlbound/#fastbtlbound/" > powheg.input
    
    
    for i in {1..48}
    do
	ch=`char $i`
	
	$prg <<EOF > run-$i.log 2>&1 &
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
	11) echo 0011 ;;
	12) echo 0012 ;;
	13) echo 0013 ;;
	14) echo 0014 ;;
	15) echo 0015 ;;
	16) echo 0016 ;;
	17) echo 0017 ;;
	18) echo 0018 ;;
	19) echo 0019 ;;
	20) echo 0020 ;;
	21) echo 0021 ;;
	22) echo 0022 ;;
	23) echo 0023 ;;
	24) echo 0024 ;;
	25) echo 0025 ;;
	26) echo 0026 ;;
	27) echo 0027 ;;
	28) echo 0028 ;;
	29) echo 0029 ;;
	30) echo 0030 ;;
	31) echo 0031 ;;
	32) echo 0032 ;;
	33) echo 0033 ;;
	34) echo 0034 ;;
	35) echo 0035 ;;
	36) echo 0036 ;;
	37) echo 0037 ;;
	38) echo 0038 ;;
	39) echo 0039 ;;
	40) echo 0040 ;;
	41) echo 0041 ;;
	42) echo 0042 ;;
	43) echo 0043 ;;
	44) echo 0044 ;;
	45) echo 0045 ;;
	46) echo 0046 ;;
	47) echo 0047 ;;
	48) echo 0048 ;;
	49) echo 0049 ;;
    esac
}

