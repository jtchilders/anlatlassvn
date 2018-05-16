#!/bin/bash

# This runs powheg storing reweighting information

cat powheg.input-save | sed 's/#storeinfo_rwgt.*/storeinfo_rwgt 1/' > powheg.input

../pwhg_main

# This runs powheg storing reweighting information for 6 more scale combinations

for i in {1..6}
do

case $i in
1) renfac=0.5; facfac=0.5  ;;
2) renfac=0.5; facfac=1    ;;
3) renfac=1  ; facfac=0.5  ;;
4) renfac=1  ; facfac=2    ;;
5) renfac=2  ; facfac=1    ;;
6) renfac=2  ; facfac=2    ;;
esac

cat powheg.input-save | sed "s/#compute_rwgt.*/compute_rwgt 1/ ; s/#renscfact.*/renscfact $renfac/ ; s/#facscfact.*/facscfact $facfac/ ; " > powheg.input

../pwhg_main

mv -f pwgevents-rwgt.lhe pwgevents.lhe

done
