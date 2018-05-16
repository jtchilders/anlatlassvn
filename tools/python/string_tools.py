#!/usr/bin/env python
import os,sys,logging
logger = logging.getLogger(__name__)

class Set:
   def __init__(self,variable_name = None,value = None):
      self.variable_name = variable_name
      self.value = value

def get_variable_and_value(
                           string,
                           variable_name,
                           variable_value_separator = '=',
                           set_separator = ' ',
                          ):
   ''' Finds a variable and value pattern inside a string, i.e.  get_variable_and_value('x=5 y=6 z=7','x') returns '5'.
       get_variable_and_value('x=5,y=6,z=7','y','=',',') returns '6'.
   '''
   variable_start_index = string.find(variable_name)
   value_start_index = variable_start_index+len( variable_name + variable_value_separator)
   value_end_index  = string.find(set_separator,value_start_index) - 1
   
   if value_start_index == value_end_index:
      return string[value_start_index]
   return string[value_start_index:value_end_index]

def get_variables_and_values(
                             string,
                             variable_value_separator = '=',
                             set_separator = ' ',
                            ):
   ''' 
      Finds all variables and values inside a string, i.e. get_variables_and_values('x=6;y=7;z=8','=',';') returns
      {'x':'6','y':'7','z':'8'}.
   '''

   stripped_string = string.strip()
   sets = stripped_string.split(set_separator)
   output = {}
   for set in sets:
      var,val = split_variable_and_value(set,variable_value_separator)
      output[var] = val

   return output



def split_variable_and_value(
                             string,
                             variable_value_separator = '=',
                            ):
   '''
      returns the values on either side of the separator
   '''
   parts = string.split(variable_value_separator)
   
   return parts[0].strip(),parts[1].strip()



if (__name__ == '__main__'):
   
   logging.basicConfig()

   string = 'x=5 y=6 z=7'
   print get_variable_and_value(string,'y')


