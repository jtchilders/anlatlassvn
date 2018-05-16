#!/usr/bin/env python
import os,subprocess
from mysql import mysql_wrapper as mdb

mysql = mdb.MySQL(db_server_password='1w@ANLHEP')

input_file_template = {
                       'wjet':'alpgen_inputs/vjets_%s_syst/w%s%ijets_%s.txt',
                       'zjet':'alpgen_inputs/vjets_%s_syst/z%s%ijets_%s.txt',
                       'zqq': 'alpgen_inputs/vjets_%s_syst/z%s%ijets_%s.txt',
                       'wqq': 'alpgen_inputs/vjets_%s_syst/w%s%ijets_%s.txt',
                      }
group_id_template = 'alpgen.%s.%s.%s.%ijets.%s'

systematics = ['ktfac2','ktfac0p5','qfac2','qfac0p5']
njets = [0,1,2,3,4,5]
flavors = {'wjet':['enu','munu'],
           'zjet':['ee','mumu'],
           'zqq':['eebb','eecc','mumubb','mumucc'],
           'wqq':['enubb','enucc','munubb','munucc'],
          }
processes = ['wjet','zjet','wcjet','wqq','zqq']
energies = ['13TeV','8TeV']

event_settings = { 'wjet':{ 
                           0: {'e':32e6,'i':30,'w':300e6,'n': 80e3},
                           1: {'e':16e6,'i':30,'w':300e6,'n':100e3},
                           2: {'e':8e6, 'i':30,'w':160e6,'n':150e3},
                           3: {'e':4e6, 'i':30,'w':90e6, 'n':300e3},
                           4: {'e':2e6, 'i':30,'w':50e6, 'n':350e3},
                           5: {'e':1e6, 'i':30,'w':50e6, 'n':  1e6},
                          },
                   'zjet':{
                           0: {'e':96e6,'i':10,'w':32e7, 'n':  5e4},
                           1: {'e':48e6,'i':10,'w':16e7, 'n':  5e4},
                           2: {'e':24e6,'i':10,'w':8e7,  'n':  4e5},
                           3: {'e':12e6,'i':10,'w':4e7,  'n':  2e6},
                           4: {'e':6e6, 'i':10,'w':2e7,  'n':  1e6},
                           5: {'e':3e6, 'i':10,'w':1e7,  'n':  5e5},
                          },
                   'zqq':{
                           0: {'e':24e6,'i':10,'w':9e7,  'n':  4e5},
                           1: {'e':12e6,'i':10,'w':5e7,  'n':  2e6},
                           2: {'e':6e6, 'i':10,'w':3e7,  'n':  1e6},
                           3: {'e':3e6, 'i':10,'w':2e7,  'n':  5e5},
                          },
                   'wqq':{
                           0: {'e':64e6, 'i':10,'w':16e7,  'n':  4e5},
                           1: {'e':32e6, 'i':10,'w':8e7,   'n':  2e6},
                           2: {'e':16e6, 'i':10,'w':4e7,   'n':  1e6},
                           3: {'e':8e6,  'i':10,'w':2e7,   'n':  5e5},
                          },

                 }
job_shape_settings = { 
                        0: { 'o':512,  'c':4, 't':30 },
                        1: { 'o':512,  'c':32,'t':30 },
                        2: { 'o':512,  'c':64,'t':30 },
                        3: { 'o':2048, 'c':64,'t':60 },
                        4: { 'o':4096, 'c':64,'t':60 },
                        5: { 'o':16384,'c':64,'t':60 },
                     }

njet = 0
process = 'wqq'
energy = '13TeV'
machine = 'mira'

for flavor in flavors[process]:
   print 'flavor: ',flavor
   for syst in systematics:
      print '  syst: ',syst
      group_id = group_id_template % (process,flavor,energy,njet,syst)
      print '   group id: ',group_id
      input_file = input_file_template[process] % (energy,flavor,njet,syst)
      print '  input_file: ',input_file
      
      # check for entry in DB
      entries,labels = mysql.select(where_conditions='group_identifier=\'%s\''%group_id)
      
      evtset = event_settings[process][njet]
      js = job_shape_settings[njet]
      
      if len(entries) > 0:
         entry = entries[-1]
         regenerate_path = entry[labels['output_url']].replace('gsiftp://atlasgridftp02.hep.anl.gov','')
         print 'found previous job, regenerate_path:',regenerate_path
         
         cmd = '/users/hpcusers/svn/jobpayloads/submit_alpgen.py -r %s -n %i -p %s -o %i -c %i -a %s -t %i -s %s -g %s' % (
            regenerate_path,evtset['n'],
            process,js['o'],js['c'],input_file,js['t'],machine,group_id
            )
      else:
         cmd = '/users/hpcusers/svn/jobpayloads/submit_alpgen.py -e %i -i %i -w %i -n %i -p %s -o %i -c %i -a %s -t %i -s %s -g %s' % (
            evtset['e'],evtset['i'],evtset['w'],evtset['n'],
            process,js['o'],js['c'],input_file,js['t'],machine,group_id
            )
      
      print '   cmd: ',cmd
      
      if True:
         p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,env=os.environ)
         stdout,stderr = p.communicate()
         
         print p.returncode
         if p.returncode == 0:
            print stdout
            print stderr
         else:
            print stdout
            print stderr
            raise Exception('error occurred, returncode = ' + str(p.returncode))
         

      #raw_input('press enter to continue...')





