#!/bin/bash

svn export ../MadGraphStuff MadTMP

if [ -e proc_card.dat ]
then
\cp proc_card.dat MadTMP/Cards/
else
echo proc_card.dat missing!
exit -1
fi

cd MadTMP


./NewProcess.sh $*


\cp -a MadGraph_POWHEG/my_proc/Source/DHELAS .
\cp -a MadGraph_POWHEG/my_proc/Source/MODEL .
\cp -a MadGraph_POWHEG/my_proc/SubProcesses Madlib

\rm Madlib/coupl.inc

cp MODEL/coupl.inc Madlib/

echo editing DHELAS/Makefile
sed -i 's/FC[\t ]*=.*// ; s/^DEST[\t ]*=.*/DEST = ..\// ;  s/^LIBRARY[\t ]*=.*/LIBRARY = ..\/libdhelas3.a/'   DHELAS/Makefile

echo editing MODEL/makefile
sed -i 's/^F77[\t ]*=.*// ; s/-ffixed-line-length-132// ; s/^LIBDIR[\t ]*=.*/LIBDIR = ..\// ' MODEL/makefile
echo editing Madlib/makefile
sed -i 's/^F77[\t ]*=.*// ; s/^LIBDIR[\t ]*=.*/LIBDIR = ..\// ' Madlib/makefile


# -n: do not overwrite existing files

echo > ../MGfiles.list

for i in *
do
# -a : exists, file or directory
if ! [ -a ../$i ]
then
    mv $i ../
    echo $i >> ../MGfiles.list
fi
done
