#!/usr/bin/env python
import os,sys,optparse,logging
logger = logging.getLogger(__name__)

def main():
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   parser.add_option('-i','--input',dest='input',help='input')
   parser.add_option('-l','--sherpa-libs',dest='sherpa_libs',help='location of Sherpa libraries')
   parser.add_option('-s','--sherpa-install',dest='sherpa_install',help='location of Sherpa')
   parser.add_option('-m','--hepmc-install',dest='hepmc_install',help='location of HepMC')
   parser.add_option('-q','--sqlite-install',dest='sqlite_install',help='location of sqlite')
   parser.add_option('-z','--zlib',dest='zlib',help='static zlib library',default='')
   options,args = parser.parse_args()

   
   manditory_args = [
                     'input',
                     'sherpa_libs',
                     'sherpa_install',
                     'hepmc_install',
                     'sqlite_install',
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   

   LIBRARY_SEARCH = ['/bin/sh','libtool','--tag=CXX','--mode=link']
   MAKE_DIRECTORY_SEARCH = ['make[','Entering directory']
   SKIP_LIBRARIES = ['HepMC','sqlite3','gfortran','m','quadmath',
                     'dl','z','mpichf90-gcc','mpich-gcc','opa-gcc',
                     'mpl-gcc','pami-gcc','SPI','SPI_cnk','rt',
                     'pthread','stdc++',
                    ]
   SKIP_LIBRARY_PATHS = [
                         #'/soft/compilers/gcc/4.8.4/powerpc64-bgq-linux/lib',
                         '/soft/libraries/alcf/current/gcc/ZLIB/lib',

                        ]
   libraries= {}
   current_directory = None
   for line in open(options.input):
      
      # keep track of which directory we are in with make
      if all(x in line for x in MAKE_DIRECTORY_SEARCH):
         start_index = line.find('`') + 1
         end_index = line.find("'")
         current_directory = line[start_index:end_index]

      elif all(x in line for x in LIBRARY_SEARCH):
         # found a linking command

         # extract all the '.lo' files
         # first split the command into words
         parts = line.split()
         # loop over the parts to check for things I want
         library_name = None
         object_files = []
         dependent_libraries = []
         dependent_lib_search_paths = []
         for i,part in enumerate(parts):
            # found the output
            if part == '-o':
               library_name = parts[i+1]
            # found a libtool object file
            elif '.lo' in part:
               object_files.append(part)
            # found a depedent library
            elif '-l' == part[0:2]:
               if not any(x == part[2:] for x in SKIP_LIBRARIES):
                  dependent_libraries.append(part[2:])
            # found a dependent library path
            elif '-L' == part[0:2]:
               dependent_lib_search_paths.append(part[2:])
         libraries[library_name] = {
               'name':library_name,
               'objects':object_files,
               'dependent_libraries':dependent_libraries,
               'dependent_lib_search_paths':dependent_lib_search_paths,
               'current_directory':current_directory,
               }

   #print libraries
   #print libraries['libSherpaMain.la']

   # now build the one shared library
   sherpaMain = libraries['libSherpaMain.la']

   # loop over library dependencies
   includes = get_dependencies(sherpaMain['dependent_libraries'],libraries)
   
   static_cmd = ''
   
   # creat cmd
   libtool_cmd = '/bin/sh ./libtool  --tag=CXX   --mode=link mpic++ -Wl,--no-as-needed -g -O2  -rdynamic -Wl,--no-as-needed  -o libSherpaLibs.la -rpath ' + options.sherpa_libs + ' -Wl,--no-undefined'
   
   for file in sherpaMain['objects']:
      if file not in libtool_cmd:
         libtool_cmd += ' ' + sherpaMain['current_directory'] + '/' + file

   for file in includes:
      if file not in libtool_cmd:
         libtool_cmd += ' ' + file

   for path in SKIP_LIBRARY_PATHS:
      libtool_cmd += ' -L' + path
   #for lib in SKIP_LIBRARIES:
   #   if lib != 'quadmath' and lib != 'HepMC' and lib != 'sqlite3':
   #      libtool_cmd += ' -l' + lib
   #libtool_cmd += ' -L/users/hpcusers/svn/generators/sherpa/v220/HepMC-2.06.08/install/lib'

   # add the object files from hepmc and sqlite
   libtool_cmd += ' ' + options.hepmc_install + '/src/*.lo'
   libtool_cmd += ' ' + options.sqlite_install + '/sqlite3.lo'

   # save this command for the static build
   static_cmd = libtool_cmd

   # add gzip libs and dl libs
   libtool_cmd += ' -lz'
   libtool_cmd += ' -ldl'

   if not os.path.exists('libSherpaLibs.la'):
      print libtool_cmd
      os.system(libtool_cmd)

   # now build dynamically loaded Sherpa 
   if not os.path.exists('Sherpa') and os.path.exists('libSherpaLibs.la'):
      libtool_cmd = '/bin/sh ./libtool  --tag=CXX   --mode=link mpic++ -Wl,--no-as-needed -dynamic -g -O2  -rdynamic -Wl,--no-as-needed  -o Sherpa ' + options.sherpa_install + '/SHERPA/Run/Main.o -rpath ' + options.sherpa_libs + ' -L./ -lSherpaLibs'
      print libtool_cmd
      os.system(libtool_cmd)

   # now build statically linked Sherpa
   static_cmd = static_cmd.replace('mpic++','mpic++ ' + options.sherpa_install + '/SHERPA/Run/Main.o')
   static_cmd += ' -static'
   static_cmd = static_cmd.replace('-o libSherpaLibs.la','-o SherpaStatic')
   if options.zlib != '':
      static_cmd += ' ' + options.zlib
   else:
      static_cmd += ' -lz'
   static_cmd += ' -ldl'
   print static_cmd
   os.system(static_cmd)
   


def get_dependencies(dependent_libraries,libraries):
   includes = []
   # loop over dependencies
   for library in dependent_libraries:
      # get the build information for this dependency
      try:
         dep_lib = libraries['lib'+library+'.la']
      except KeyError,e:
         logger.info('skipping library: ' + library)
         
      if len(dep_lib['dependent_libraries']) > 0:
         dep_includes = get_dependencies(dep_lib['dependent_libraries'],libraries)
         includes += dep_includes

      for obj in dep_lib['objects']:
         includes.append(dep_lib['current_directory'] + '/' + obj)

   return includes


if __name__ == "__main__":
   main()
