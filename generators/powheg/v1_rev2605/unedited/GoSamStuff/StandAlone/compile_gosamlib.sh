#!/bin/sh

cd $1/../GoSamlib/

FILES=$(echo *.f *.f90 ' ' | sed 's/qlonshellcutoff.f// ; s/qlconstants.f// ; s/.f /.o /g ; s/.f90 /.o /g ; ')

cd $1

ar cru gosamlib.a $FILES
