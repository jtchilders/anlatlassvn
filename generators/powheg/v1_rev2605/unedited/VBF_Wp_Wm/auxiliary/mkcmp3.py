#!/usr/bin/python

# regular expressions
import math
import re
import sys

# debugger
#import pdb
#pdb.set_trace()

ifil=0
fin1=open(sys.argv[1])
fin2=open(sys.argv[2])
fin3=open(sys.argv[3])
val=''
if len(sys.argv) == 5:
	val=sys.argv[4]
cmp='cmp3'
rat='rat'
fou=open(cmp+'.top','w')


def findplot1():
	while True:
		line=fin1.readline()
		if line=='':
			return 1
		if (re.search('TITLE TOP',line)!=None or re.search('TITLE BOTTOM',line)!=None) and  re.search('""',line)==None :
			title=line
			fou.write(line)
			while True:
				line=fin1.readline()
				fou.write(line)
				if re.search('SET ORDER',line)!=None:
					line='set order x y '+val+' dy\n'
					break
			findplot2(re.sub('  *',' ',title))
			findplot3(re.sub('  *',' ',title))
		fou.write(line)


def findplot2(title):
	fin2.seek(1)
	while True:
		line=fin2.readline()
		if line=='':
			return 1
		if (re.sub('  *',' ',line)==title):
			while True:
				line=fin2.readline()
				if re.search('SET ORDER',line)!=None:
					fou.write(' set color red\n')
					break
			while True:
				line=fin2.readline()
				if re.search('HIST',line)!=None:
					break
				fou.write(line)
			fou.write(line)
			fou.write(' PLOT\n')
			fou.write('set color white\n')
			return 0


def findplot3(title):
	fin3.seek(1)
	while True:
		line=fin3.readline()
		if line=='':
			return 1
		if (re.sub('  *',' ',line)==title):
			while True:
				line=fin3.readline()
				if re.search('SET ORDER',line)!=None:
					fou.write(' set color green\n')
					break
			while True:
				line=fin3.readline()
				if re.search('HIST',line)!=None:
					break
				fou.write(line)
			fou.write(line)
			fou.write(' PLOT\n')
			fou.write('set color white\n')
			return 0



fou.write('set bar y size=0.02 permanent')

findplot1()

fin1.close()
fin2.close()
fin3.close()
fou.close()
print 'bye'
