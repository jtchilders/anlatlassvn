#!/usr/bin/env python
import os,sys,optparse,logging,copy,json,time
import MySQLdb
from AlpgenUnwParFile import AlpgenUnwParFile
from AlpgenInputFile import AlpgenInputFile
from submit_alpgen import submit_alpgen
logger = logging.getLogger(__name__)

DEFAULT_DB_SERVER       = 'localhost'
DEFAULT_DB_USERNAME     = 'root'
DEFAULT_DB_PASSWORD     = '1w@ANLHEP'
DEFAULT_DB_NAME         = 'argo_cluster_production'
DEFAULT_TABLE_NAME      = 'argo_core_argodbentry'
DEFAULT_WALLTIME        = -1
DEFAULT_NODES           = -1

QUEUE_LIMIT             = 15 # Mira Queue Limit (-1 because Balsam seems to veto last job, perhaps it is counting the running job against the queue limit)

MAX_EVENTS_REACHED      = 1 # signal that the total number of events requested has been reached

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('','--db-server',dest='db_server',help='The server name on which the database is hosted.',default=DEFAULT_DB_SERVER)
   parser.add_option('','--db-username',dest='db_username',help='The username for the MySQL server',default=DEFAULT_DB_USERNAME)
   parser.add_option('','--db-password',dest='db_password',help='The password for the MySQL server',default=DEFAULT_DB_PASSWORD)
   parser.add_option('','--db-name',dest='db_name',help='The name of the MySQL DB.',default=DEFAULT_DB_NAME)
   parser.add_option('-t','--tablename',dest='tablename',help='The MySQL DB table name.',default=DEFAULT_TABLE_NAME)
   parser.add_option('-g','--group-identifier',dest='group_identifier',help='The group identifier with which to screen db entries.')
   parser.add_option('-n','--num-evts',dest='num_evts',help='Continue submitting events until this number of events is reached.')
   parser.add_option('-i','--itr-evts',dest='itr_evts',help='The number of events per run.')
   parser.add_option('-x','--no-submit',dest='submit',help='do not submit the message to ARGO. For testing purposes.',action='store_false',default=True)
   parser.add_option('','--walltime',dest='walltime',help='Override the walltime of the previous job.',default=DEFAULT_WALLTIME,type='int')
   parser.add_option('','--nodes',dest='nodes',help='Override the number of nodes for each job.',type='int',default=DEFAULT_NODES)
   options,args = parser.parse_args()

   
   manditory_args = [
                     'db_server',
                     'db_username',
                     'db_name',
                     'tablename',
                     'group_identifier',
                     'num_evts',
                     'itr_evts',
                     'submit',
                     'walltime',
                     'nodes',
                  ]

   for man in manditory_args:
      if not options.__dict__[man]:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)

   logger.info(' group_identifier = %30s' % options.group_identifier)
   logger.info(' num_evts         = %30s' % options.num_evts)
   logger.info(' itr_evts         = %30s' % options.itr_evts)
   logger.info(' db_server        = %30s' % options.db_server)
   logger.info(' db_username      = %30s' % options.db_username)
   logger.info(' db_name          = %30s' % options.db_name)
   logger.info(' tablename        = %30s' % options.tablename)
   logger.info(' submit           = %30s' % str(options.submit))
   logger.info(' walltime         = %30s' % str(options.walltime))
   logger.info(' nodes            = %30s' % str(options.nodes))

   while True:
      retval = batch_submit_alpgen(
                           group_identifier  = options.group_identifier,
                           num_evts          = int(options.num_evts),
                           itr_evts          = int(options.itr_evts),
                           db_server         = options.db_server,
                           db_username       = options.db_username,
                           db_password       = options.db_password,
                           db_name           = options.db_name,
                           tablename         = options.tablename,
                           submit            = options.submit,
                           walltime          = options.walltime,
                           nodes             = options.nodes,
                         )
      if retval == MAX_EVENTS_REACHED:
         logger.info('reached max events.')
         break
      time.sleep(10)
   


def batch_submit_alpgen(
                        group_identifier,
                        num_evts,
                        itr_evts,
                        db_server         = DEFAULT_DB_SERVER,
                        db_username       = DEFAULT_DB_USERNAME,
                        db_password       = DEFAULT_DB_PASSWORD,
                        db_name           = DEFAULT_DB_NAME,
                        tablename         = DEFAULT_TABLE_NAME,
                        submit            = True,
                        walltime          = DEFAULT_WALLTIME,
                        nodes             = DEFAULT_NODES,
                       ):
   
   total_evts = 0
   n_total_evts = 0
   latest_row = None
   last_history_row  = None
   total_unfinished = 0
   queued_jobs = 0
   con = MySQLdb.connect(db_server,db_username,db_password,db_name)
   with con:
      cur = con.cursor(MySQLdb.cursors.DictCursor)
      # first find the number of queued jobs
      cur.execute('SELECT * FROM ' + tablename)
      rows = cur.fetchall()
      
      # count the number of queued jobs
      for row in rows:
         if 'SUBJOB_QUEUED' in row['state_current'] or 'SUBJOB_SUBMITTED' in row['state_current']:
            queued_jobs += 1

      
      # get all the jobs for this group id
      cur.execute('SELECT * FROM ' + tablename + ' WHERE group_identifier="' + group_identifier + '"')

      rows = cur.fetchall()
      
      for row in rows:

         if 'HISTORY' in row['state_current']:
            # extract events generated
            job_path = row['output_url'].replace('gsiftp://atlasgridftp02.hep.anl.gov','').replace('//','/')
            unwPar = AlpgenUnwParFile.read_file(os.path.join(job_path,'alpout_unw.par'))
            if unwPar is not None:
               total_evts += unwPar.event_count
               n_total_evts += 1
            if last_history_row is None:
               last_history_row = copy.deepcopy(row)
         elif 'FAILED' in row['state_current']:
            continue
         else:
            total_unfinished += 1

         # check for latest job
         if latest_row is None:
            latest_row = copy.deepcopy(row)
         elif latest_row['id'] < row['id']:
            latest_row = copy.deepcopy(row)

      logger.info('total_evts       = ' + str(total_evts))
      logger.info('total_jobs       = ' + str(n_total_evts))
      logger.info('total_unfinished = ' + str(total_unfinished))
      #con.close()
   
   events_queued_and_done = (total_evts + (1.*total_evts/n_total_evts)*total_unfinished)
   logger.info('Total events in queue and completed: ' + str(events_queued_and_done))
   if events_queued_and_done > num_evts:
      return MAX_EVENTS_REACHED

   logger.info(' Jobs Queued: ' + str(queued_jobs))
   if queued_jobs >= QUEUE_LIMIT:
      return
   
   # regenerate latest row
   if n_total_evts > 0:
      # only regenerate if there is room on the queue
      if ( total_unfinished < QUEUE_LIMIT ):

         # extract regen path from info from latest row
         job_path = latest_row['output_url'].replace('gsiftp://atlasgridftp02.hep.anl.gov','').replace('//','/')
         json_text = latest_row['job_list_text'].replace('\n','')
         balsam_jobs = json.loads(json_text)
         if len(balsam_jobs) > 0:
            evgen_job = balsam_jobs[-1]

            process = None
            if 'zjetgen' in evgen_job['executable_args']:
               process = 'zjet'
            elif 'wjetgen' in evgen_job['executable_args']:
               process = 'wjet'
            else:
               logger.error('Process could not be identified = ' + str(evgen_job['executable_args']))
               return
            logger.info('process = ' + str(process))
            # extract iseed
            alpIn1 = AlpgenInputFile()
            alpIn1.read(os.path.join(job_path,'alpout.input.1'))
            iseed1 = int(float(alpIn1.options['iseed1'].value)) + 1
            logger.info('iseed1 = ' + str(iseed1))
            # set walltime
            if walltime == DEFAULT_WALLTIME:
               walltime = evgen_job['wall_minutes']
            # set number of nodes
            numnodes = evgen_job['nodes']
            if nodes != DEFAULT_NODES:
               numnodes = nodes
            # submit regeneration job
            submit_alpgen(
                           evts_per_iter        = 0,
                           numiters             = 0,
                           num_warmup           = 0,
                           num_weighted         = itr_evts,
                           process              = process,
                           numnodes             = numnodes,
                           cpus_per_node        = evgen_job['processes_per_node'],
                           alpgen_input_file    = '',
                           pdf_filename         = 'cteq6l1.tbl',
                           walltime             = walltime,
                           site                 = evgen_job['target_site'],
                           regen_path           = job_path,
                           iseed1               = iseed1,
                           group_identifier     = group_identifier,
                           submit               = submit,   
                         )
   # if we have reached the end send a signal so as to stop the loop
   
   
   return 0
         


if __name__ == "__main__":
   main()
