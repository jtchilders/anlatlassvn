#!/usr/bin/env python
import json

files_per_directory = {}
libraries_per_directory = {}
current_directory = ''
for line in open('make.log'):

   if 'Entering directory' in line:
      start_index = line.find('`') + 1
      end_index = line.find("'",start_index)
      current_directory = line[start_index:end_index]
      parts = current_directory.split('/')
      print 'Entering ',current_directory
      if current_directory not in files_per_directory:
         files_per_directory[current_directory] = []

      if current_directory not in libraries_per_directory:
         libraries_per_directory[current_directory] = {}

   if line.startswith('libtool: compile:') and ' -c ' in line:
      
      end_index = line.find('.C')
      print line
      if end_index == -1:
         #print 'ERROR: no .C found'
         end_index = line.find('.c')
         if end_index == -1:
            #print 'ERROR: no .c found'
            end_index = line.find('.f')
            if end_index == -1:
               print 'ERROR on line: ',line[:-1]
               continue

      start_index = line.rfind(' ',0,end_index)

      compile_file = line[start_index:end_index+2]
      #print line
      #print current_directory
      print 'Found file: ',compile_file
      files_per_directory[current_directory].append(compile_file)
      
      #print files_per_directory

   if line.startswith('libtool: link:') and '-o' in line:
      start_index = line.find('-o') + 3
      end_index = line.find(' ',start_index)

      compile_library = line[start_index:end_index].replace('.libs/','')
      print line
      print 'found library: ',compile_library



      libraries_per_directory[current_directory][compile_library] = []

      parts = line.split()
      for part in parts:
         if '.o' in part or part.startswith('-l') or part.endswith('.so'):
            libraries_per_directory[current_directory][compile_library].append(part)

      print 'found files for library',compile_library,libraries_per_directory[current_directory]
      

      

print files_per_directory
print ' '
print libraries_per_directory


json.dump(files_per_directory,open('files_per_directory.json','w'))
json.dump(libraries_per_directory,open('libraries_per_directory.json','w'))

