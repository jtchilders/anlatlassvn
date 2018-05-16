#!/usr/bin/env python
import os,sys,subprocess,glob


def main():
   print 'usage: <script> "bash glob string"'
   print '       you must use the quotations or the glob string will not parse properly'
   glob_string = sys.argv[1]
   filenames = glob.glob(glob_string)
   


   script = """#!/usr/bin/env bash
source /users/jchilders/scripts/setupAthena.sh 17.2.11.12
checkFile.py -f  \\\n"""
   for filename in filenames:
      script += '    ' + filename + '  \\\n'
   
   print script
   script_file = open('tmp_script.sh','w')
   script_file.write(script)
   script_file.close()
   script_file = None

   p = subprocess.Popen('bash tmp_script.sh',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
   #line = p.stdout.readline()
   #while line:
   #   print line[0:-1]
   #   line = p.stdout.readline()
   
   stdout,stderr = p.communicate()
   #print stderr
   print 'loop lines'
   total_events = 0
   stdout_lines = stdout.splitlines()
   for line in stdout_lines:
      if line.find('Nbr Events:') >= 0:
         print line
         total_events += int(line.split()[2])

   print 'total number of events = ' + str(total_events)


if __name__ == "__main__":
   main()



