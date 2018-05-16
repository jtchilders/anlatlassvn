#!/usr/bin/env python

import os,glob

file_list = glob.glob('*/Makefile')

for file in file_list:
   dir = file.split('/')[0]
   os.system('make -C ' + dir + ' cleanall')

rm_list = glob.glob('*/*gen')
rm_list += glob.glob('*/*.mod')
rm_list += glob.glob('*work/*gen')
rm_list += glob.glob('*work/*gen_*')
rm_list += glob.glob('*work/*gen90')
rm_list += glob.glob('*work/*gen90_*')
rm_list += glob.glob('*lib/XXX.f90')
for item in rm_list:
   os.remove(item)
