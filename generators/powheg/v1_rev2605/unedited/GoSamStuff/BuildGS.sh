#!/bin/sh

# ####################################################################
# To successfully employ this script, make sure that your bash is in #
# the bin - directory or adjust accordingly.                         #
######################################################################

RUNDIR=${PWD}
GOSAMSTUFFDIR=${PWD}/../GoSamStuff
GOSAMDIR=${PWD}/GoSam_POWHEG
GOSAMLIBDIR=${PWD}/GoSamlib
VIRDIR=$GOSAMDIR/Virtual
CONTRIBDIR=$GOSAMSTUFFDIR/GoSam-contrib

# If input is 'HELP':
if [ "$1" = "help" ] || [ "$1" = "" ]
then
    echo "******************************************************************************************"
    echo "-- HELP menu for the script buildvirt --                                                  "
    echo "Possible running options:                                                                 "
    echo "./BuildGS virtual       : generates the full virtual code,                                "
    echo "./BuildGS interface      : creates the new files needed to run the virtual amplitude,      "
    echo "./BuildGS standalone    : makes standalone virtual code,                                  "
    echo "./BuildGS allvirt       : makes virtual, interface and standalone all in one go,           "
    echo "./BuildGS help          : shows this menu.                                                "
    echo "******************************************************************************************"
    exit
fi

# If input is 'CLEANVIRT':
if [ "$1" = "cleanvirt" ]
then
    echo "This will delete the files for the virtual amplitude."
    echo "Are you sure you want to proceed? (Yes/No)"
    read ANS
    if [ "$ANS" = "Yes" ]
    then
	if [ -d $RUNDIR/GoSamlib ]
	then
	    echo "A standalone directory /GoSamlib also exists."
	    echo "Do you want to remove this as well? (Yes/No)"
	    read ANSLIB
	    if [ "$ANSLIB" = "Yes" ]
	    then
		echo "---> Removing/GoSamlib directory ..."
		rm -fr $RUNDIR/GoSamlib
	    else
		echo "Directory /GoSamlib will not be removed!"
	    fi
	fi
	echo "---> Making backup copy of generation files and cards ..."
	cp -f $GOSAMDIR/gosam.rc $GOSAMDIR/filter.py $GOSAMDIR/orderfile.lh $RUNDIR/
	echo "---> Removing GoSam_POWHEG directory ..."
	rm -f $GOSAMDIR
	exit
    else
	echo "Aborted!"
    fi
fi

# GENVIRT
if [ "$1" = "virtual" ] || [ "$1" = "allvirt" ]
then
    if [ ! -d $GOSAMDIR ]
    then
	mkdir $GOSAMDIR
    fi 
    
    if [ -f $RUNDIR/MadGraph_POWHEG/my_proc/SubProcesses/orderfile.lh ]
    then
	cp $RUNDIR/MadGraph_POWHEG/my_proc/SubProcesses/orderfile.lh $GOSAMDIR/orderfile.lh
    else
	echo "No orderfile found!!"
	exit
    fi
    
    if [ -f $RUNDIR/gosam.rc ]
    then
	:
    else
	echo "------------------------------------------"
	echo "File 'gosam.rc' not found."
	read -p "Take the one in '../GoSamStuff/gosam.rc'? [Y/n]: " REPLY
	if [ $REPLY = "Y" ]
	then
	    cp $GOSAMSTUFFDIR/Templates/gosam.rc $RUNDIR/gosam.rc
	else
	    echo "Please provide a 'gosam.rc' file in "${PWD}
	    exit
	fi
    fi
    echo "---> GoSam is writing the code for virtual part ..."
    cd $RUNDIR
    gosam.py --olp --mc=powhegbox --config=gosam.rc --ignore-unknown --force --destination=$VIRDIR $GOSAMDIR/orderfile.lh
    cp $RUNDIR/gosam.rc $GOSAMDIR
    cp -f $RUNDIR/filter.py $GOSAMDIR
    
    cd $VIRDIR
    echo "---> Generating code for virtual amplitudes ..."
    make source
    if [ "$1" = "virtual" ]
    then
	exit
    fi
fi


# INTERFACE
if [ "$1" = "interface" ] || [ "$1" = "allvirt" ]
then
    if [ ! -f write_pwhg_files.f ]
    then
	cp $GOSAMSTUFFDIR/Templates/write_pwhg_files.f $GOSAMDIR/write_pwhg_files.f
    fi

    echo "---> Generating new files to run the virtual amplitude ..."
    cd $GOSAMDIR
    gfortran -o write_pwhg_files write_pwhg_files.f
    ./write_pwhg_files

    echo "---> Copying new files to process directory ..."
    echo $RUNDIR
    cd $RUNDIR
    mv $RUNDIR/virtual.f $RUNDIR/virtual.f.dummy
    mv $RUNDIR/init_couplings.f $RUNDIR/init_couplings.f.old
    cp $GOSAMDIR/virtual_new.f $RUNDIR/virtual.f
    cp $GOSAMDIR/init_couplings_new.f $RUNDIR/init_couplings.f
    cp $GOSAMSTUFFDIR/Templates/PhysPars.h $RUNDIR/PhysPars.h
    cp $GOSAMSTUFFDIR/Templates/pwhg_gosam.f $RUNDIR/pwhg_gosam.f

    echo "---> Removing old executables ..."
    rm -f $GOSAMDIR/write_pwhg_files
    rm -f $RUNDIR/Madlib/write_proc_labels
    rm -f $RUNDIR/Madlib/write_proc_labels_real

    echo "*******************************************************"
    echo "The new files:                                         "
    echo "- virtual.f                                            "
    echo "- init_couplings.f                                     "
    echo "were created. The old files were renamed to:           "
    echo "- virtual.f.dummy                                      "
    echo "- init_couplings.f.old                                 "
    echo "                                                       "
    echo "The new file                                           "
    echo "-PhysPars.h                                            "
    echo "-pwhg_gosam.f                                          "
    echo "was copied from the 'GoSamStuff/Templates' folder      "
    echo "                                                       "
    echo "Please before compiling generate a standalone version  "
    echo "of the virtual code.                                   "
    echo "*******************************************************"
    if [ "$1" = "interface" ]
    then
	exit
    fi
fi

# STANDALONE
if [ "$1" = "standalone" ] || [ "$1" = "allvirt" ]
then
    MAKEDEPF90=$(command -v makedepf90 | gawk '{ if($1=="") printf "0"; else printf "present"}'); 
    if [ $MAKEDEPF90 = "present" ]
    then
	if [ ! -d $RUNDIR/GoSamlib ]
	then
	    mkdir $RUNDIR/GoSamlib
	fi
	echo "---> Generating standalone copy of virtual code ..."
	for file in `find $GOSAMDIR -path "*.f90"`
	do    
	    nmodules=`grep -i '^ *end *module ' $file | wc -l`  
	    if [ "$nmodules" = "1" ]
		echo -n "."
	    then
		modname=`grep -i '^ *end *module ' $file | sed 's/^ *end *module //' | tr -d ' '`.f90
		\cp $file $GOSAMLIBDIR/$modname
	    elif  [ "$nmodules" = "0" ]
	    then
		:
	    else
		echo $file
		echo 'More than 1 module per file! Exiting ...'
		exit -1
	    fi
	done

        # copy all needed fortran and include files from the contrib directory
        # to the GoSamlib subdirectory
	for file in  `find $CONTRIBDIR -path "*.f.in"`  `find $CONTRIBDIR -path "*.f90.in"`
	do
	    filename="${file##*/}"
	    
	    sed -e 's/\@DATADIR\@/GoSamlib/g' \
		-e 's/\@PACKAGE\@/gosam-contrib/g' \
		-e 's/\@VERSION\@/2.1.1/g' \
		-e 's/\@case_with_lt\@/!AC!/g' \
		-e 's/\@case_with_golem\@/!AC!/g' \
		-e 's/\@case_with_avh\@/    /g' \
		-e 's/\@case_with_ql\@/    /g' \
		-e 's/\@case_wout_lt\@/    /g' \
		-e 's/\@case_wout_golem\@/    /g' \
		-e 's/\@case_wout_avh\@/!AC!/g' \
		-e 's/\@case_wout_ql\@/!AC!/g' \
		-e 's/\@lt_real_kind\@/kind(1.0d0)/g' \
		-e 's/\@fortran_real_kind\@/kind(1.0d0)/g' $file > $GOSAMLIBDIR/"${filename%.*}"
	    
	done

	for file in  `find $CONTRIBDIR -path "*.[hf]"`  `find $CONTRIBDIR -path "*.[hf]90"`
	do
	    \cp $file $GOSAMLIBDIR/"${file##*/}"
	done

        # create the dependencies file
	cd $GOSAMLIBDIR
	
	makedepf90 *.f90 > deps.txt
	
	cat $GOSAMSTUFFDIR/StandAlone/Makefile.nodeps deps.txt > $GOSAMLIBDIR/Makefile.virt.dep
	
	cp -f $GOSAMSTUFFDIR/StandAlone/compile_gosamlib.sh $GOSAMLIBDIR/compile_gosamlib.sh
	
	echo "done"
	echo "---> Copying generation files to:"
	echo $RUNDIR"/GoSamlib/ ..."
	cp -f $GOSAMDIR/gosam.rc $GOSAMLIBDIR/gosam.rc
	cp -f $GOSAMDIR/orderfile.lh $GOSAMLIBDIR/orderfile.lh
	cp -f $GOSAMDIR/orderfile.olc $GOSAMLIBDIR/orderfile.olc
	if [ -f $GOSAMDIR/filter.py ]
	then
	    cp -f $GOSAMDIR/filter.py $GOSAMLIBDIR/filter.py
	fi
    else
	echo "*******************************************************"	
	echo "The program 'makedepf90' was not found on your machine."
	echo "You need to install it to produce a working standalone "
	echo "version for the virtual amplitude code. A tarball with "
	echo "the code can be found in the folder:                   "
	echo $GOSAMSTUFFDIR/StandAlone
	echo "*******************************************************"
    fi
    if [ "$1" = "standalone" ]
    then
	exit
    fi
fi

exit

