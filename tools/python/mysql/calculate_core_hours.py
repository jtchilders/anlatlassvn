#!/usr/bin/env python
import os,sys,optparse,logging,imp
from mysql.mysql_wrapper import MySQL
logger = logging.getLogger(__name__)


DEFAULT_TABLE_NAME='argo_core_argodbentry'
DEFAULT_SERVER_NAME='localhost'
DEFAULT_USERNAME='root'
DEFAULT_PASSWORD=''
DEFAULT_DB_NAME='argo_cluster_production'

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-g','--grp-filter',dest='grp_filter',help='The filter string to apply to the group_identifier')
   parser.add_option('-p','--pk-filter',dest='pk_filter',help='Only use jobs above pk number')

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

   change_db_attribute(
            grp_filter     = options.grp_filter,
            pk_filter      = options.pk_filter,
            db_table       = options.db_table,
            db_server      = options.db_server,
            db_user        = options.db_user,
            db_password    = options.db_password,
            db_name        = options.db_name,
         )


def change_db_attribute(
         grp_filter,
         pk_filter,
         db_table,
         db_server,
         db_user,
         db_password,
         db_name,
      ):
   
   logger.info('%-12s %-s' % ('grp_filter',str(grp_filter)))
   logger.info('%-12s %-s' % ('pk_filter',str(pk_filter)))
   logger.info('%-12s %-s' % ('db_table',db_table))
   logger.info('%-12s %-s' % ('db_server',db_server))
   logger.info('%-12s %-s' % ('db_user',db_user))
   logger.info('%-12s %-s' % ('db_password',db_password))
   logger.info('%-12s %-s' % ('db_name',db_name))

   mysql = MySQL(db_server,db_user,db_password,db_name)

   where_conditions = ''
   if grp_filter is not None:
      where_conditions = "group_identifier='" + grp_filter + "' "
   if pk_filter is not None:
      where_conditions += " pk='" + pk
   entries,label = mysql.select(where_conditions=where_conditions)

   logger.info('retrieved ' + str(len(entries)) + ' entries.')

   for entry in entries:
      attr_value = entry[label[attribute]]
      new_attr_value = replace
      if transform is not None:
         new_attr_value = transform(attr_value)
      logger.info('old: ' + attr_value + ' new: ' + new_attr_value)

      set_command = attribute + "='" + new_attr_value +"'"
      where_conditions = "id='" + str(entry[label['id']]) + "'"

      mysql.update(set_command=set_command,where_conditions=where_conditions)
      


if __name__ == "__main__":
   main()
