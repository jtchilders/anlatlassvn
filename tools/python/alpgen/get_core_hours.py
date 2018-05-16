#!/usr/bin/env python
import os,sys,optparse,logging,glob,time
from mysql.mysql_wrapper import MySQL

from argo_core.user_job_analysis.alpgen.AlpgenInputFile import AlpgenInputFile
from argo_core.user_job_analysis.alpgen.AlpgenUnwParFile import AlpgenUnwParFile
from argo_core.user_job_analysis.alpgen.AlpgenComboRamdiskOutputLog import AlpgenComboRamdiskOutputLog
from argo_core.user_job_analysis.logfiles.CobaltLog import CobaltLog
from argo_core.user_job_analysis.logfiles.EdisonLog import EdisonLog
from argo_core.user_content.alpgen.tools import sizeof_fmt,get_condor_runtime,get_cobalt_queue_run_time,get_timezone_offset,get_num_unw_evt,get_unw_file_size
import common_core.Serializer as Serializer

logger = logging.getLogger(__name__)

DEFAULT_TABLE_NAME='argo_core_argodbentry'
DEFAULT_SERVER_NAME='localhost'
DEFAULT_USERNAME='root'
DEFAULT_PASSWORD='1w@ANLHEP'
DEFAULT_DB_NAME='argo_cluster_production'

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('--db-table',dest='db_table',help='DB Table name. [default='+DEFAULT_TABLE_NAME+']',default=DEFAULT_TABLE_NAME)
   parser.add_option('--db-server',dest='db_server',help='DB Server name. [default='+DEFAULT_SERVER_NAME+']',default=DEFAULT_SERVER_NAME)
   parser.add_option('--db-user',dest='db_user',help='DB user name. [default='+DEFAULT_USERNAME+']',default=DEFAULT_USERNAME)
   parser.add_option('--db-password',dest='db_password',help='DB password. [default='+DEFAULT_PASSWORD+']',default=DEFAULT_PASSWORD)
   parser.add_option('--db-name',dest='db_name',help='DB name. [default='+DEFAULT_DB_NAME+']',default=DEFAULT_DB_NAME)
   options,args = parser.parse_args()

   
   manditory_args = [
                     'db_table',
                     'db_server',
                     'db_user',
                     'db_password',
                     'db_name',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   
   mysql = MySQL(options.db_server,options.db_user,options.db_password,options.db_name)

   where_conditions = "id > 820 AND id < 2633 AND state_current='HISTORY'"
   entries,label = mysql.select(where_conditions=where_conditions)

   logger.info('got ' + str(len(entries)) + ' entries.')
   total_cputime = 0
   total_events = 0
   for entry in entries:
      output_url = entry[label['output_url']]
      job_path = output_url.replace('gsiftp://atlasgridftp02.hep.anl.gov','')
      if not os.path.exists(job_path):
         job_path = job_path.replace('/grid/','/atlasfs/')

      if not os.path.exists(job_path):
         logger.info('no job path found')
         continue

      logger.info('job_path: ' + job_path)

      balsam_jobs = Serializer.deserialize(entry[label['job_list_text']])
      evtgen_job = balsam_jobs[-1]
      numnodes = evtgen_job['nodes']
      
      if os.path.exists(os.path.join(job_path,'alpout_unw.par')):
         alpgenUnwParFile = AlpgenUnwParFile.read_file(os.path.join(job_path,'alpout_unw.par'))
         total_events += alpgenUnwParFile.event_count

      alpgen_cobalt_output_log_filenames   = sorted(glob.glob(os.path.join(job_path,'*.output')))
      alpgen_torque_output_log_filenames   = sorted(glob.glob(os.path.join(job_path,'alpgenCombo.sh.o*')))
      alpgen_slurm_output_log_filenames    = sorted(glob.glob(os.path.join(job_path,'slurm-*.out')))
      alpgen_condor_output_log_filenames   = sorted(glob.glob(os.path.join(job_path,'condor_stdout.txt.*')))
      alpOutputLog = None
      alpgen_output_log_filename = ''

      if ( len(alpgen_cobalt_output_log_filenames) > 0 and 
            any(evtgen_job['target_site'] == x for x in ['mira','cetus','vesta'])):
         alpgen_output_log_filename = alpgen_cobalt_output_log_filenames[-1]
      elif ( len(alpgen_torque_output_log_filenames) > 0 and 
             evtgen_job['target_site'] == 'edison'):
         alpgen_output_log_filename = alpgen_torque_output_log_filenames[-1]
      elif (len(alpgen_slurm_output_log_filenames) > 0 and
             evtgen_job['target_site'] == 'edison'):
         alpgen_output_log_filename = alpgen_slurm_output_log_filenames[-1]
      elif (len(alpgen_condor_output_log_filenames) > 0 and
             evtgen_job['target_site'] == 'argo_cluster'):
         alpgen_output_log_filename = alpgen_condor_output_log_filenames[-1]
      else:
         logger.error('no output log file found')
         alpgen_output_log_filename = ''
      
      if os.path.exists(alpgen_output_log_filename):
         alpOutputLog = AlpgenComboRamdiskOutputLog(alpgen_output_log_filename)
         if alpOutputLog.evtgen_failed:
            evtgen_run_time = 'Failed'
         else:
            evtgen_run_time  = time.strftime('%H:%M:%S', time.gmtime(alpOutputLog.evtgen_total_secs))

            if alpOutputLog.evtgen_gen_secs > 0:
               evtgen_gen_time = time.strftime('%H:%M:%S', time.gmtime( alpOutputLog.evtgen_gen_secs ))
            else:
               evtgen_gen_time = 'N/A'

         if alpOutputLog.unw_failed:
            unw_run_time = 'Failed'
         else:
            unw_run_time     = time.strftime('%H:%M:%S', time.gmtime(alpOutputLog.unw_total_secs))
            if alpOutputLog.unw_gen_secs > 0:
               unw_gen_time = time.strftime('%H:%M:%S', time.gmtime( alpOutputLog.unw_gen_secs ))
            else:
               unw_gen_time = 'N/A'
         if alpOutputLog.agg_failed:
            agg_run_time = 'Failed'
         else:
            agg_run_time     = time.strftime('%H:%M:%S', time.gmtime(alpOutputLog.agg_total_secs))


         if alpOutputLog.total_secs > 0:
            total_cputime += float(numnodes)*(float(alpOutputLog.total_secs)/60./60.)*16.
            total_time = time.strftime('%H:%M:%S', time.gmtime( alpOutputLog.total_secs ))
         else:
            total_time = str(alpOutputLog.total_secs)
      
      logger.info('total_cputime: ' + str(total_cputime))
      logger.info('total_events: ' + str(total_events))


if __name__ == "__main__":
   main()
