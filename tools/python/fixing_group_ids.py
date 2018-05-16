#!/usr/bin/env python
import os,sys,subprocess,shlex,glob
from mysql.mysql_wrapper import MySQL



def run(cmd,ignore_error_code=False,shell=False):
   command = shlex.split(cmd)
   if shell:
      command = cmd
   print command
   os.system(command)
   p = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=False)
   lines = ''
   for line in p.stdout:
      print line[:-1]
      lines += line
   p.wait()
   if p.returncode != 0 and not ignore_error_code:
      raise Exception(' cmd: ' + cmd + ' failed with exit code: ' + str(p.returncode) + '\n' + lines)
   return lines

def main():
   DEFAULT_TABLE_NAME='argo_core_argodbentry'
   DEFAULT_SERVER_NAME='localhost'
   DEFAULT_USERNAME='root'
   DEFAULT_PASSWORD=''
   DEFAULT_DB_NAME='argo_cluster_production'

   mysql = MySQL(DEFAULT_SERVER_NAME,DEFAULT_USERNAME,DEFAULT_PASSWORD,DEFAULT_DB_NAME)

   where_conditions = "group_identifier='None'"

   entries,labels = mysql.select(where_conditions=where_conditions,db_table_name=DEFAULT_TABLE_NAME)

   #group_id_postfix = '_HPC.TXT.mc15_v1_i11'

   zdecmode_names = {1:'ee', 2:'mumu',3:'tautau',4:'ll' }
   wdecmode_names = {1:'enu',2:'munu',3:'taunu', 4:'lnu'}
   ihvy_names     = {4:'c',5:'b',6:'t'}
   
   for entry in entries:
      folder = entry[labels['output_url']]
      folder = folder.split('//')[1].replace('atlasgridftp02.hep.anl.gov','')

      njets = -1
      wjets = False
      wdecmode = -1
      zjets = False
      zdecmode = -1
      heavy_jets = False
      heavy_qq = False
      heavy_q = False
      ihvy = -1

      print 'id = ' + str(entry[labels['id']]) + '  folder = ' + folder

      if not os.path.exists(folder):
         print folder, 'does not exit'
         group_identifier = 'job folder missing'
         set_command = 'group_identifier="'+group_identifier+'"'
         where_conditions = 'id="'+str(entry[labels['id']])+'"'
         print set_command,where_conditions
         #mysql.update(set_command=set_command,where_conditions=where_conditions)
      elif os.path.exists(os.path.join(folder,'alpout.input.1')):
         lines = run('grep njets ' + os.path.join(folder,'alpout.input.1'))
         if len(lines) > 0:
            njets = int(float(lines.split()[1]))

         try:
            lines = run('grep iwdecmod ' + os.path.join(folder,'alpout.input.1'))
            if len(lines) > 0:
               wdecmode = int(float(lines.split()[1]))
               wjets = True
         except Exception,e:
            try:
               lines = run('grep izdecmod ' + os.path.join(folder,'alpout.input.1'))
               if len(lines) > 0:
                  zdecmode = int(float(lines.split()[1]))
                  zjets = True
            except:
               raise Exception(' cannot find izdecmod or iwdecmod in input file.')
         
         try:
            lines = run('grep ihvy ' + os.path.join(folder,'alpout.input.1'))
            if len(lines) > 0:
               ihvy = int(float(lines.split()[1]))
               heavy_jets = True
            
               if wjets:
                  try:
                     lines = run('grep wqqgen90 ' + os.path.join(folder,glob.glob(os.path.join(folder,'*.cobaltlog'))[0]))
                     if len(lines) > 0: heavy_qq = True
                  except:
                     try:
                        lines = run('grep wqgen90 ' + os.path.join(folder,glob.glob(os.path.join(folder,'*.cobaltlog'))[0]))
                        if len(lines) > 0: heavy_q = True
                     except: pass
               elif zjets:
                  try:
                     lines = run('grep zqqgen90 ' + os.path.join(folder,glob.glob(os.path.join(folder,'*.cobaltlog'))[0]))
                     if len(lines) > 0: heavy_qq = True
                  except: pass
         except: pass


         group_identifier = 'alpgen'
         if zjets: 
            group_identifier += '.Z' + zdecmode_names[zdecmode] 
         elif wjets:
            group_identifier += '.W' + wdecmode_names[wdecmode]

         if heavy_q:
            group_identifier += ihvy_names[ihvy]
         elif heavy_qq:
            group_identifier += ihvy_names[ihvy] + ihvy_names[ihvy]

         group_identifier += ('%ijets' % njets)

         print folder,'   ',group_identifier
         set_command = 'group_identifier="'+group_identifier+'"'
         where_conditions = 'id="'+str(entry[labels['id']])+'"'
         print set_command,where_conditions
         #mysql.update(set_command=set_command,where_conditions=where_conditions)
      else:
         try:
            files = glob.glob(folder + '/*sherpa*')
            if len(files) > 0:
               set_command = 'group_identifier="sherpa"'
               where_conditions = 'id="'+str(entry[labels['id']])+'"'
               print folder
               print set_command,where_conditions
               #mysql.update(set_command=set_command,where_conditions=where_conditions)
         except:
            print folder
         



def old_stuff_for_sherpa():
   search_folder = '/grid/atlas/hpc/data/sherpa/inputcards/Wenu/joboptions'

   run_file = os.path.join(folder,'Run.py')
   jo_file = ''
   if os.path.exists(run_file):
      stdout = run('/users/hpcusers/svn/tools/bash_scripts/find_identical_files.sh ' + search_folder + ' ' + run_file,True)
      words = stdout.split()
      jo_file = words[3]
      #file_start_index = stdout.find(search_folder)
      #file_end_index = stdout.find(' ',file_start_index+1)
      #print file_start_index,file_end_index
      #if file_start_index >= 0:
      #   jo_file = stdout[file_start_index:file_end_index]
   else: pass
   print 'found file: ' + jo_file
   jo_file = jo_file.split('/')[-1].replace('.py','').replace('MC15.','').replace('_BFilter','') # get rid of path

   if len(jo_file) > 10:

      gid = 'group.phys-gener.sherpa020101.' + jo_file + '_HPC.TXT.mc15_v1_i12'

      print gid

      set_command = "group_identifier='" + gid + "'"
      where_conditions = "id='" + str(entry[labels['id']]) + "'"
      #mysql.update(set_command=set_command,where_conditions=where_conditions,db_table_name=DEFAULT_TABLE_NAME)





if __name__ == '__main__':
   main()
