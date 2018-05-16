#!/usr/bin/env python
import os,sys,optparse,logging,multiprocessing
import alpgen_reorganize_unw
import MySQLdb
logger = logging.getLogger(__name__)

class OneTask:
   def __init__(self,
                db_group_pattern,
                collection_name,
                pythia_inclusive,
                max_output=1e12,
               ):
      self.db_group_pattern      = db_group_pattern
      self.collection_name       = collection_name
      self.pythia_inclusive      = pythia_inclusive
      self.max_output            = max_output

tasks = [ 
   
   #OneTask('alpgen.wenuc0jets' , 'group.phys-gener.alpgen214.361860.AlpgenPythia_P2012_WenucNp0_HPC.TXT.mc15_v3',False,12e6),
   #OneTask('alpgen.wmunuc0jets' , 'group.phys-gener.alpgen214.361865.AlpgenPythia_P2012_WmunucNp0_HPC.TXT.mc15_v3',False,12e6),
   #OneTask('alpgen.wtaunuc0jets' , 'group.phys-gener.alpgen214.361870.AlpgenPythia_P2012_WtaunucNp0_HPC.TXT.mc15_v3',False,12e6),

   #OneTask('alpgen.wenuc1jets' , 'group.phys-gener.alpgen214.361861.AlpgenPythia_P2012_WenucNp1_HPC.TXT.mc15_v3',False,4.6e6),
   #OneTask('alpgen.wmunuc1jets' , 'group.phys-gener.alpgen214.361866.AlpgenPythia_P2012_WmunucNp1_HPC.TXT.mc15_v3',False,4.6e6),
   #OneTask('alpgen.wtaunuc1jets' , 'group.phys-gener.alpgen214.361871.AlpgenPythia_P2012_WtaunucNp1_HPC.TXT.mc15_v3',False,4.6e6),

   OneTask('alpgen.wenuc2jets' , 'group.phys-gener.alpgen214.361862.AlpgenPythia_P2012_WenucNp2_HPC.TXT.mc15_v3',False,2.4e6),
   OneTask('alpgen.wmunuc2jets' , 'group.phys-gener.alpgen214.361867.AlpgenPythia_P2012_WmunucNp2_HPC.TXT.mc15_v3',False,2.4e6),
   OneTask('alpgen.wtaunuc2jets' , 'group.phys-gener.alpgen214.361872.AlpgenPythia_P2012_WtaunucNp2_HPC.TXT.mc15_v5',False,2.4e6),

   #OneTask('alpgen.wenuc3jets' , 'group.phys-gener.alpgen214.361863.AlpgenPythia_P2012_WenucNp3_HPC.TXT.mc15_v3',False,0.6e6),
   #OneTask('alpgen.wmunuc3jets' , 'group.phys-gener.alpgen214.361868.AlpgenPythia_P2012_WmunucNp3_HPC.TXT.mc15_v3',False,0.6e6),
   #OneTask('alpgen.wtaunuc3jets' , 'group.phys-gener.alpgen214.361873.AlpgenPythia_P2012_WtaunucNp3_HPC.TXT.mc15_v3',False,0.6e6),

   #OneTask('alpgen.wenuc4jets' , 'group.phys-gener.alpgen214.361864.AlpgenPythia_P2012_WenucNp4_HPC.TXT.mc15_v3',True,0.2e6),
   #OneTask('alpgen.wmunuc4jets' , 'group.phys-gener.alpgen214.361869.AlpgenPythia_P2012_WmunucNp4_HPC.TXT.mc15_v3',True,0.2e6),
   #OneTask('alpgen.wtaunuc4jets' , 'group.phys-gener.alpgen214.361874.AlpgenPythia_P2012_WtaunucNp4_HPC.TXT.mc15_v3',True,0.2e6),
]





def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(process)d %(levelname)s:%(name)s:%(message)s')

   '''
   parser = optparse.OptionParser(description='')
   parser.add_option('-i','--input',dest='input',help='input')
   options,args = parser.parse_args()

   
   manditory_args = [
                     
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   '''
   
   
   pool = multiprocessing.Pool(5)
   
   pool.map(run_reorg,tasks)


   logger.info('Done')
   

def run_reorg(task):
   
   logger.info(' Running task with: ')
   logger.info('    db_group_pattern: ' + str(task.db_group_pattern))
   logger.info('    collection_name:  ' + str(task.collection_name))
   logger.info('    pythia_inclusive: ' + str(task.pythia_inclusive))
   logger.info('    max_output:       ' + str(task.max_output))

   # access database and get input_basenames
   input_basenames,ids = get_input_basenames(task.db_group_pattern)
   try:
      alpgen_reorganize_unw.alpgen_reorganize_unw(
           input_basenames,
           task.collection_name,
           task.max_output,
           pythia_inclusive = task.pythia_inclusive,
           files_per_subset = 99999,
           pythia_tmp_output = '/tmp/pythia_output.hepmc.p' + str(os.getpid())
          )
   except Exception,e:
      if 'No Input Files provided' in str(e):
         logger.info(' failed: ' + str(task.collection_name))
         return
      else:
         raise
    
   # change group_identifiers
   change_group_identifiers(task.collection_name,ids)

   logger.info(' done with: ' + str(task.collection_name))


   
def get_input_basenames(db_group_pattern):
   input_basenames = []
   ids = []
   try:
      connection = MySQLdb.connect('localhost','root','1w@ANLHEP','argo_cluster_production')
   except:
      logger.exception(' error connecting to database ')
      raise
   
   try:
      cursor = connection.cursor()
   except:
      logger.exception(' error getting cursor to database connection ')
      connection.close()
      raise

   try:
      #logger.info(' group_identifier: ' + db_group_pattern)
      cursor.execute("SELECT * FROM argo_core_argodbentry WHERE group_identifier LIKE '%" + db_group_pattern + "%' AND state_current='HISTORY';")
      
      entries = cursor.fetchall()
      #logger.info(' entries: ' + str(len(entries)))
      
      description = cursor.description
      label = {}
      for i in range(len(description)):
         label[description[i][0]] = i

      for entry in entries:
         #logger.info(' Entry id = ' + str(entry[label['id']]) + '  group_identifier = ' + entry[label['group_identifier']] + '  state_current = ' + entry[label['state_current']] )

         work_dir = entry[label['output_url']]
         id = entry[label['id']]
         work_dir = work_dir.replace('gsiftp://atlasgridftp02.hep.anl.gov','')

         #logger.info('work_dir: ' + work_dir)
         
         
         # check if input folder exists
         if not os.path.exists(work_dir):
            work_dir = os.path.join('/atlasfs/atlas/hpc/grid/argo',os.path.basename(work_dir))
            if not os.path.exists(work_dir):
               logger.error(' Could not find work directory for job id = ' + str(id) + ' in ' + work_dir)
               raise Exception(' Could not find work directory for job id = ' + str(id) + ' in ' + work_dir)
            
         input_basename = os.path.join(work_dir,'alpout')

         input_basenames.append(input_basename)

         ids.append(id)
      
      connection.close()
   except:
      logger.exception(' error access database ')
      connection.close()
      raise

   #logger.info('input_basenames: ' + str(input_basenames))

   return input_basenames,ids

def change_group_identifiers(new_group_identifier,ids=[]):
   try:
      connection = MySQLdb.connect('localhost','root','1w@ANLHEP','argo_cluster_production')
   except:
      logger.exception(' error connecting to database ')
      raise
   
   try:
      cursor = connection.cursor()
   except:
      logger.exception(' error getting cursor to database connection ')
      connection.close()
      raise

   try:
      for id in ids:
         cursor.execute("UPDATE argo_core_argodbentry SET group_identifier='" + new_group_identifier + "' WHERE id='" + str(id) + "';")

   except:
      logger.exception(' error updating group identifier ')
      connection.close()
      raise



if __name__ == "__main__":
   main()
