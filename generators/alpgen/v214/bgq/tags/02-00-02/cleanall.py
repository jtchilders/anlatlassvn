#!/usr/bin/env python

import os,glob

file_list = glob.glob('*/Makefile')

for file in file_list:
   dir = file.split('/')[0]
   os.system('make -C ' + dir + ' cleanall')

exe_list = glob.glob('*/*gen')
for exe in exe_list:
   os.remove(exe)
