#!/bin/bash

# this function does a single step of grid improvement
function doxgrid {

cat powheg.input-save | sed "s/#multigrid.*/multigrid $1/ ; s/ncall2.*/ncall2 0/ " > powheg.input

for i in {1..20}
do
echo $i | ../pwhg_main > run-$i.log 2>&1 &
done

wait

\rm pwgxgrid.dat

cat powheg.input-save | sed 's/#multigrid.*/multigrid 0/ ; s/ncall2.*/ncall2 0/ ' > powheg.input
\rm pwgxgrid.dat

../pwhg_main <<EOF
1
EOF

\rm pwggridinfo*
}


# first grid improvment step

doxgrid 1
mv pwgbtlgrid.top pwgbtlgrid-1.top 

# second grid improvment step

doxgrid 2
mv pwgbtlgrid.top pwgbtlgrid-2.top 

# third grid improvment step

doxgrid 3
mv pwgbtlgrid.top pwgbtlgrid-3.top 

cat powheg.input-save | sed 's/nubound.*/nubound 0/' > powheg.input

# compute in parallel upper bounding grid
for i in {1..20}
do
echo $i | ../pwhg_main > run-$i.log 2>&1 &
done

wait

cat powheg.input-save | sed 's/numevts.*/numevts 0/' > powheg.input

# compute in parallel upper bounds for event generation
for i in {1..20}
do
echo $i | ../pwhg_main > run-$i.log 2>&1 &
done

wait

cat powheg.input-save > powheg.input

# Generate events in parallel
for i in {1..20}
do
echo $i | ../pwhg_main > run-$i.log 2>&1 &
done



