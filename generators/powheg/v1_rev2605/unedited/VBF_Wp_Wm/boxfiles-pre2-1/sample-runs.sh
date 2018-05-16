#!/bin/bash


# compute grid 1st iter 

cat powheg.input-save | sed 's/xgriditeration.*/xgriditeration 1/ ; s/parallelstage.*/parallelstage 1/' > powheg.input
for i in {1..4}
do
echo $i | ../pwhg_main > run-1a-$i.log 2>&1 &
done
wait


# compute grid 2nd iter 
cat powheg.input-save | sed 's/xgriditeration.*/xgriditeration 2/; s/parallelstage.*/parallelstage 1/' > powheg.input
for i in {1..4}
do
echo $i | ../pwhg_main > run-1b-$i.log 2>&1 &
done
wait


# compute NLO
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 2/ ' > powheg.input
for i in {1..4}
do
echo $i | ../pwhg_main > run-2-$i.log 2>&1 &
done
wait


# compute in parallel upper bounding grid
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 3/' > powheg.input
for i in {1..4}
do
echo $i | ../pwhg_main > run-3-$i.log 2>&1 &
done
wait



# generate events 
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 4/' > powheg.input
for i in {1..4}
do
echo $i | ../pwhg_main > run-4-$i.log 2>&1 &
done

