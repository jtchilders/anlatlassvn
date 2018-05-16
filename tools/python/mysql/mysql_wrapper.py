import MySQLdb,logging
logger = logging.getLogger(__name__)

DEFAULT_TABLE_NAME='argo_core_argodbentry'
DEFAULT_SERVER_NAME='localhost'
DEFAULT_USERNAME='root'
DEFAULT_PASSWORD=''
DEFAULT_DB_NAME='argo_cluster_production'

class MySQL:
   ''' wrapper class for MySQLdb '''
   def __init__(self,
                db_server_name=DEFAULT_SERVER_NAME,
                db_server_login=DEFAULT_USERNAME,
                db_server_password=DEFAULT_PASSWORD,
                db_name=DEFAULT_DB_NAME,
                ):
      # open a connection
      try:
         self.connection = MySQLdb.connect(db_server_name,db_server_login,db_server_password,db_name)
      except Exception,e:
         logger.error('Received Exception when opening connection to MySQLdb, server_name='+db_server_name+', db_name='+db_name+': ' + str(e))
         raise
      # get a cursor from the connection.
      try:
         self.cursor = self.connection.cursor()
      except Exception,e:
         logger.error('Received Exception when getting cursor: ' + str(e))
         raise

   def select(self,
            field_list = '*',  # can be a list of the field you want output.
            where_conditions = None, # can be a conditional list used to select "X='5' AND Y='Time To Go'"
            db_table_name = DEFAULT_TABLE_NAME,
             ):
      try:
         cmd = 'SELECT ' + field_list + ' FROM ' + db_table_name
         if where_conditions is not None:
            cmd += ' WHERE ' + where_conditions
         cmd += ';'
         logger.debug('execute: ' + cmd)
         self.cursor.execute(cmd)
      except Exception,e:
         logger.error('SELECT Received Exception running execute with cmd= \n' + cmd + '\n Exception: ' + str(e))
         raise

      entries = self.cursor.fetchall()
      labels = self.get_column_labels(self.cursor)

      return entries,labels

   def update(self,
            set_command, # should be: "field='value'"
            db_table_name = DEFAULT_TABLE_NAME,
            where_conditions = None, # can be a conditional list used to select "X='5' AND Y='Time To Go'"
            ):
      try:
         cmd = 'UPDATE ' + db_table_name + ' SET ' + set_command
         if where_conditions is not None:
            cmd += ' WHERE ' + where_conditions
         cmd += ';'
         self.cursor.execute(cmd)
      except Exception,e:
         logger.error('UPDATE Received Exception running execute with cmd= \n' + cmd + '\n Exception: ' + str(e))
         raise

   @staticmethod
   def get_column_labels(cursor):
      description = cursor.description
      label = {}
      for i in range(len(description)):
         label[description[i][0]] = i
      return label