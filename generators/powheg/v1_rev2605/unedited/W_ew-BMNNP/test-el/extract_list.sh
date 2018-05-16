#!/bin/sh
if [ "$1" = "--help" ] || [ "$1" = "" ]; then
	echo "Usage: ./extract_list --list"
	echo "        list distributions"
	echo 
	echo "       ./extract_list <name of distribution> <file>"
	echo "       <file> = pwgNLO.top"
	echo "              = pwgpwhgalone-output.top"
	echo "              = POWHEG+<shower>.top"
else
	if [ "$1" = "--list" ]; then
		grep "#" pwgNLO.top
	else
		grep -n "#" $2 > powheg_tmp 
		tmp1=`grep -wn $1 powheg_tmp`
		line1=`echo $tmp1 | awk -F":" '{ print $2 }'`
		tmp1=`echo $tmp1 | awk -F":" '{ print $1 }'`
		tmp2=`echo $tmp1 + 1 | bc`
		line2=`awk "NR==$tmp2" powheg_tmp | awk -F":" '{ print $1 }'`
		unset tmp1
		unset tmp2
		rm powheg_tmp
		line1=`echo $line1 + 1 | bc`
		line2=`echo $line2 - 1 | bc`
		awk "NR==$line1,NR==$line2" $2 > $1_$2_gnuplot
	fi
fi
