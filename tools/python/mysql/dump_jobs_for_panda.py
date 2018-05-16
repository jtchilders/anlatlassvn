#!/usr/bin/env python
import os,sys,optparse,logging
from mysql.mysql_wrapper import MySQL
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
   
   where_conditions = "group_identifier='group.phys-gener.alpgen214.361725.AlpgenPythia_P2012_ZtautauNp5_HPC.TXT.mc15_v1' AND id > 1450 AND state_current='HISTORY'"
   entries,label = mysql.select(where_conditions=where_conditions)

   logger.info('got ' + str(len(entries)) + ' entries.')

   for entry in entries:
      set_command = "group_identifier='alpgen.ztautau5jets'"
      where_conditions = "id='" + str(entry[label['id']]) + "'"

      mysql.update(set_command=set_command,where_conditions=where_conditions)

   
      
   


if __name__ == "__main__":
   main()
