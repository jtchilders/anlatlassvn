import pycurl,os,sys,logging,optparse,urllib
import xml.etree.ElementTree as ET
try:
    # Python 3
    from io import BytesIO
except ImportError:
    # Python 2
    from StringIO import StringIO as BytesIO
logger = logging.getLogger(__name__)

curl_buffer = None

class gsheets:
   ''' class for accessing Google Spreadsheets '''

   URL_LOGIN='https://www.google.com/accounts/ClientLogin'
   AUTH_KEY="Auth="
   LIST_FEED_REF='http://schemas.google.com/spreadsheets/2006#listfeed'
   CELLS_FEED_REF='http://schemas.google.com/spreadsheets/2006#cellsfeed'

   class XMLLinkElement:
      def __init__(self,element):
         self.element = element
         try:
            self.href = self.element.attrib['href']
         except: pass
         try:
            self.rel  = self.element.attrib['rel']
         except: pass
         try:
            self.type = self.element.attrib['type']
         except: pass
      def __str__(self):
         string = ''
         string += ' link: href = %-s\n' % self.href
         string += ' link: rel  = %-s\n' % self.rel
         string += ' link: type = %-s\n' % self.type
         return string

   class XMLEntryElement:
      def __init__(self,element):
         self.element = element
         self.read_element()
      def __str__(self):
         string = ''
         if 'id' in dir(self):
            string += ' id             = %-s\n' % self.id
         if 'updated' in dir(self):
            string += ' updated        = %-s\n' % self.updated
         if 'category_scheme' in dir(self):
            string += ' category: scheme = %-s\n' % self.category_scheme
         if 'category_term' in dir(self):
            string += ' category: term = %-s\n' % self.category_term
         if 'title' in dir(self):
            string += ' title          = %-s\n' % self.title
         if 'title_type' in dir(self):
            string += ' title_type     = %-s\n' % self.title_type
         if 'content' in dir(self):
            string += ' content        = %-s\n' % self.content
         if 'content_type' in dir(self):
            string += ' content_type   = %-s\n' % self.content_type
         if 'links' in dir(self):
            for link in self.links:
               string += str(link)
         return string
      def read_element(self):
         for child in self.element:
            if 'id' in child.tag:
               self.id = child.text
            elif 'updated' in child.tag:
               self.updated = child.text
            elif 'category' in child.tag:
               self.category_scheme = child.attrib['scheme']
               self.category_term = child.attrib['term']
            elif 'title' in child.tag:
               self.title = child.text
               self.title_type = child.attrib['type']
            elif 'content' in child.tag:
               self.content = child.text
               self.content_type = child.attrib['type']
            elif 'link' in child.tag:
               if 'links' in dir(self):
                  self.links.append(gsheets.XMLLinkElement(child))
               else:
                  self.links = [ gsheets.XMLLinkElement(child) ]
      def get_link_by_rel(self,string):
         for link in self.links:
            if string in link.rel:
               return link
         return None
      def get_link_by_href(self,string):
         for link in self.links:
            if string in link.href:
               return link
         return None

   class XMLSheetElement(XMLEntryElement):
      def __init__(self,element):
         self.element = element
         self.read_element()
         for child in self.element:
            if 'rowCount' in child.tag:
               self.row_count = int(child.text)
            elif 'colCount' in child.tag:
               self.col_count = int(child.text)
         self.link_listfeed = self.get_link_by_rel('listfeed')
         self.link_cellsfeed = self.get_link_by_rel('cellsfeed')
         self.link_visualizationApi = self.get_link_by_rel('virtualizationApi')
      def __str__(self):
         string = gsheets.XMLEntryElement.__str__(self)
         string += ' row_count      = %-10i    col_count = %-i\n' % (self.row_count,self.col_count)
         return string
      def get_rows(self,header,tag_name_prefix):
         global curl_buffer
         curl_buffer = ''
         curl = pycurl.Curl()
         if self.link_listfeed:
            curl.setopt(curl.URL,self.link_listfeed.href)
         curl.setopt(curl.HTTPHEADER,header)
         curl.setopt(curl.WRITEFUNCTION, write_function)
         curl.perform()
         root = ET.fromstring(curl_buffer)


         indent(root)
         ET.dump(root)


         self.rows = []
         for entry in root.findall(tag_name_prefix + 'entry'):
            row = gsheets.XMLRowElement(entry)
            self.rows.append(row)
         logger.info(' length rows: ' + str(len(self.rows)))


   class XMLRowElement(XMLEntryElement):
      def __init__(self,element):
         self.element = element
         self.read_element()
         self.parse_content()
      def read_element(self):
         gsheets.XMLEntryElement.read_element(self)
      def parse_content(self):
         self.columns = []
         underscore_index = self.content.find('_')
         while underscore_index >= 0:
            colon_index = self.content.find(':',underscore_index)
            #logger.info('content: ' + self.content)
            #logger.info('indices: ' + str(underscore_index) + ' ' + str(colon_index))
            column_string = self.content[underscore_index:colon_index]
            #logger.info('column_string: ' + column_string)
            column_elements = self.element.findall('{http://schemas.google.com/spreadsheets/2006/extended}' + column_string)
            #logger.info('elements:' + str(column_elements))
            if len(column_elements) > 0:
               self.columns.append(column_elements[0].text)
            #logger.info( column_string + ' ' + column_elements[0].text)
            underscore_index = self.content.find('_',colon_index)
      def __str__(self):
         string = gsheets.XMLEntryElement.__str__(self)
         string += ' column contents for this row: \n'
         for column in self.columns:
            string += ' ' + column
         string += '\n'
         return string




   def __init__(self,
                login_info_file,
                spreadsheet_key,
               ):
      self.login_info_file =login_info_file
      self.spreadsheet_key = spreadsheet_key
      self.auth_key = None
      self.worksheet_feed_url = 'https://spreadsheets.google.com/feeds/worksheets/' + str(self.spreadsheet_key) + '/private/full'
      self.header = None
      self.tag_name_prefix = None
      self.sheets = {}


   def login(self):
      # read in login information
      login_file = open(self.login_info_file)
      lines = login_file.readlines()
      login_file.close()

      # make file for output
      auth_filename = 'auth.out'
      out_file = open(auth_filename,'w')

      # setup command
      url = gsheets.URL_LOGIN
      parameters = {'Email':lines[0],'Passwd':lines[1]}
      post_data = {'accountType':'GOOGLE',
                   'source':'Google-cURL-Example',
                   'service':'wise'
                  }

      curl = pycurl.Curl()
      curl.setopt(curl.URL, url + '?' + urllib.urlencode(parameters))
      curl.setopt(curl.POSTFIELDS,urllib.urlencode(post_data))
      curl.setopt(curl.WRITEDATA,out_file)
      # run curl
      curl.perform()

      curl.close()
      out_file.close()

      # grab authentication key
      auth_file = open(auth_filename)
      for line in auth_file:
         if gsheets.AUTH_KEY in line:
            self.auth_key = line[len(gsheets.AUTH_KEY):-1]
            break
      auth_file.close()

      self.header = ['Authorization: GoogleLogin auth=' + str(self.auth_key) ]

      logger.debug('authentication key: ' + str(self.header))


   def get_sheets(self):
      global curl_buffer
      curl_buffer = ''
      curl = pycurl.Curl()
      curl.setopt(curl.URL,str(self.worksheet_feed_url))
      curl.setopt(curl.HTTPHEADER,self.header)
      curl.setopt(curl.WRITEFUNCTION, write_function)
      curl.perform()
      print '>>>> ',curl_buffer
      self.worksheet_feed_xml_root = ET.fromstring(curl_buffer)

      # Google XML tags are prefixed with the Atom version like this
      # {http://www.w3.org/2005/Atom}feed
      # or
      # {http://www.w3.org/2005/Atom}entry
      # I will extract this URL instead of assuming it is always the same
      root_tag = self.worksheet_feed_xml_root.tag

      index = root_tag.find('feed')
      self.tag_name_prefix = root_tag[0:index]

      # extract the URLs for the spreadsheet inside the file
      for spreadsheet in self.worksheet_feed_xml_root.findall(self.tag_name_prefix + 'entry'):

         sheet = gsheets.XMLSheetElement(spreadsheet)
         self.elementTree = spreadsheet
         self.sheets[sheet.title] = sheet

   def get_rows(self):
      for sheet in self.sheets.values():
         sheet.get_rows(self.header,self.tag_name_prefix)





def write_function(buffer):
   global curl_buffer
   if buffer:
      curl_buffer += buffer

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def main():
   logging.basicConfig(level=logging.DEBUG)

   parser = optparse.OptionParser(description='')
   parser.add_option('-l','--login-info',dest='login_info_file',help='Input file where first line is google login, second line is password.')
   parser.add_option('-k','--sheet-key',dest='sheet_key',help='Google Spreadsheet key from URL: http://docs.google.com/spreadsheet/ccc?key=XXXXX')
   options,args = parser.parse_args()

   if options.login_info_file is None or options.sheet_key is None:
      parser.error('must specify options.')

   get_sheets(options.login_info_file,options.sheet_key)


def get_sheets(login_info_file,sheet_key):
   gs = gsheets(login_info_file,sheet_key)
   gs.login()
   gs.get_sheets()
   gs.get_rows()

   for sheet in gs.sheets.values():
      logger.info(str(sheet))

      i = 0
      for row in sheet.rows:
         logger.info( str(i) + '\n' + str(row))
         i += 1




   return gs





if __name__ == '__main__':
   main()
