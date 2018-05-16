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
   parser.add_option('-a','--attribute',dest='attribute',help='Search MySQL DB using this attribute name')
   parser.add_option('-m','--attr-match',dest='attr_match',help='Search MySQL DB for attributes that match this string.')
   parser.add_option('-r','--replace',dest='replace',help='Replace found attribute with this string. replace and transform cannot both be used.')
   parser.add_option('-l','--like',dest='use_like',action='store_true',help='Instead of requiring the attr-match string to match exactly, allow wildcards "%" ',default=False)
   parser.add_option('-t','--transform',dest='transform',help='A file that contains a transform(attr_in) function which transforms each attr-match into a new attribute which will be set in the DB. The function should return the transformed string. replace and transform cannot both be used.')
   parser.add_option('--db-table',dest='db_table',help='DB Table name. [default='+DEFAULT_TABLE_NAME+']',default=DEFAULT_TABLE_NAME)
   parser.add_option('--db-server',dest='db_server',help='DB Server name. [default='+DEFAULT_SERVER_NAME+']',default=DEFAULT_SERVER_NAME)
   parser.add_option('--db-user',dest='db_user',help='DB user name. [default='+DEFAULT_USERNAME+']',default=DEFAULT_USERNAME)
   parser.add_option('--db-password',dest='db_password',help='DB password. [default='+DEFAULT_PASSWORD+']',default=DEFAULT_PASSWORD)
   parser.add_option('--db-name',dest='db_name',help='DB name. [default='+DEFAULT_DB_NAME+']',default=DEFAULT_DB_NAME)
   options,args = parser.parse_args()

   
   manditory_args = [
                     'attribute',
                     'attr_match',
                     'use_like',
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
   
   if options.replace is not None and options.transform is not None:
      logger.error('replace and transform cannot be set at the same time. Only one at a time.')
      parser.print_help()
      sys.exit(-1)
   if options.replace is None and options.transform is None:
      logger.error('replace or transform must be specified, but not both.')
      parser.print_help()
      sys.exit(-1)

   transform = None
   if options.transform is not None:
      directory, module_name = os.path.split(options.transform)
      module_name = os.path.splitext(module_name)[0]

      path = list(sys.path)
      sys.path.insert(0, directory)
      try:
         module = __import__(module_name)
      finally:
         sys.path[:] = path # restore

      transform = module.transform

   change_db_attribute(
            attribute      = options.attribute,
            attr_match     = options.attr_match,
            use_like       = options.use_like,
            replace        = options.replace,
            transform      = transform,
            db_table       = options.db_table,
            db_server      = options.db_server,
            db_user        = options.db_user,
            db_password    = options.db_password,
            db_name        = options.db_name,
         )


def change_db_attribute(
         attribute,
         attr_match,
         use_like,
         replace,
         transform,
         db_table,
         db_server,
         db_user,
         db_password,
         db_name,
      ):
   
   logger.info('%-12s %-s' % ('attribute',attribute))
   logger.info('%-12s %-s' % ('attr_match',attr_match))
   logger.info('%-12s %-s' % ('use_like',str(use_like)))
   logger.info('%-12s %-s' % ('replace',replace))
   logger.info('%-12s %-s' % ('transform',str(transform)))
   logger.info('%-12s %-s' % ('db_table',db_table))
   logger.info('%-12s %-s' % ('db_server',db_server))
   logger.info('%-12s %-s' % ('db_user',db_user))
   logger.info('%-12s %-s' % ('db_password',db_password))
   logger.info('%-12s %-s' % ('db_name',db_name))

   mysql = MySQL(db_server,db_user,db_password,db_name)

   where_conditions = ''
   if use_like:
      where_conditions = attribute + " LIKE '" + attr_match + "'"
   else:
      where_conditions = attribute + "='" + attr_match + "'"
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
