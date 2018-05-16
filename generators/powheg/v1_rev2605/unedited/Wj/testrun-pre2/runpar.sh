#!/bin/bash

dir=`dirname $0`

function char {
	case $1 in
	    [1-9])  echo 000$1 ;;
	    [1-9][0-9])  echo 00$1 ;;
	    [1-9][0-9][0-9])  echo 0$1 ;;
	esac
}


cp powheg.input-save powheg.input
echo 'parallelstage 1' >> powheg.input
echo 'xgriditeration 1' >> powheg.input

for i in {1..48}
do
    echo $i | ../pwhg_main > run-1-1-`char $i`.log &
done
wait

cp powheg.input-save powheg.input
echo 'parallelstage 1' >> powheg.input
echo 'xgriditeration 2' >> powheg.input

for i in {1..48}
do
    echo $i | ../pwhg_main > run-1-2-`char $i`.log &
done
wait

for stage in 2 3 4
do

cp powheg.input-save powheg.input
echo "parallelstage $stage" >> powheg.input
for i in {1..48}
do
    echo $i | ../pwhg_main > run-$stage-`char $i`.log &
done
wait

done

# Uncomment the following to compute new weights for all scale variations
# $dir/reweight.sh
