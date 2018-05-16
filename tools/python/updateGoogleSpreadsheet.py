#!/usr/bin/env python
import os,sys,optparse,logging
import gspread.

def main():
   logging.basicConfig(level=logging.INFO)

   parser = optparse.OptionParser(description='Make updates to a google spreadsheet')
   parser.add_option('-k','--sheet-key',dest='spreadsheet_key',help='The key of the Google spreadsheet to edit.')
   parser.add_option('-r','--row-string',dest='row_string',help='This string is used to find the row to edit. THe first row which has this string will be edited, not consecutive matching rows.')
   options,args = parser.parse_args()

   if options.spreadsheet_key is None:
      parser.error('Must specify -k')
   if options.row_string is None:
      parser.error('Must specify -r')
   

   return 0


if __name__ == "__main__":
   sys.exit(main())
