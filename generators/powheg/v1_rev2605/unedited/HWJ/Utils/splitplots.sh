#!/bin/bash 

function usage {
echo 'usage: splitplots file.gp'
echo 'split the plots separated by "reset"'
exit -1
}

case a$1 in
*.gp) if ! [ -e $1 ] ;  then usage ; fi ;;
*) usage ;;
esac

i=0
(cat $1 ; echo EOF) | while :
do
read line
if [ "a$line" = aEOF ]
then
exit 0
fi
if [ "a$line" = areset ]
then
    if [ -e singleplot.gp ]
	then
	i=$[$i+1]
	gnuplot singleplot.gp
#	if [ $i = 2 ]
#	then exit
#	fi
    fi
    echo $line > singleplot.gp
else
    echo $line >> singleplot.gp
fi
done


#\rm singleplot.gp


