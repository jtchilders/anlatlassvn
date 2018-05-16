#!/usr/bin/env python
import os,sys,optparse
import submit_alpgen


def main():
   parser = optparse.OptionParser(description='submit many alpgen jobs to ARGO')
   parser.add_option('-e','--evts-per-iter',dest='evts_per_iter',help='number of events per warmup iteration',type='int')
   parser.add_option('-i','--num-iter',dest='numiters',help='number of iterations for the warmup',type='int')
   parser.add_option('-w','--warmup-weighted',dest='num_warmup',help='number of event to in the warmup, after the iterations complete',type='int')
   parser.add_option('-n','--num-weighted',dest='num_weighted',help='number of weighted events to generate.',type='int')
   parser.add_option('-p','--process',dest='process',help='define the process to generate, 2Q,4Q,hjet,top,wjet,zjet,Njet,etc.')
   parser.add_option('-o','--num-nodes',dest='numnodes',help='number of nodes to use on destination machine',type='int')
   parser.add_option('-c','--cpus-per-node',dest='cpus_per_node',help='number of CPUs per node to use on destination machine',type='int')
   parser.add_option('-a','--alpgen-input',dest='alpgen_input_file',help='The AlpGen input file which carries all the options for this generation job')
   parser.add_option('-f','--pdf-filename',dest='pdf_filename',help='The PDF Filename for Alpgen')
   parser.add_option('-t','--wall-time',dest='walltime',help='The wall time to submit to the queue in minutes.',type='int')
   parser.add_option('-s','--site',dest='site',help='Balsam site name on which to run the event generation')
   parser.add_option('-d','--dev',dest='dev',help='Run in development mode, means warmup is sent to argo_cluster_dev',action='store_true',default=False)
   parser.add_option('-r','--regen',dest='regen_path',help='Regenerate the events from a previous warmup, the argument is the path from which to copy the grid1/2 files. The input file should still be provided using the -a option.')
   parser.add_option('-x','--no-submit',dest='submit',help='do not submit the message to ARGO. For testing purposes.',action='store_false',default=True)
   parser.add_option('','--iseed1',dest='iseed1',help='override random number, iseed1, in alpgen input file',type='int')
   parser.add_option('','--iseed2',dest='iseed2',help='override random number, iseed2, in alpgen input file',type='int')
   parser.add_option('','--iseed3',dest='iseed3',help='override random number, iseed3, in alpgen input file',type='int')
   parser.add_option('','--iseed4',dest='iseed4',help='override random number, iseed4, in alpgen input file',type='int')
   
   multi_submit_alpgen(
                       options.evts_per_iter,
                       options.numiters,
                       options.num_warmup,
                       options.num_weighted,
                       options.process,
                       options.numnodes,
                       options.cpus_per_node,
                       options.alpgen_input_file,
                       options.pdf_filename,
                       options.walltime,
                       options.site,
                       options.dev,
                       options.regen_path,
                       options.submit,
                       options.iseed1,
                       options.iseed2,
                       options.iseed3,
                       options.iseed4,
                      )

def multi_submit_alpgen(
                       evts_per_iter,
                       numiters,
                       num_warmup,
                       num_weighted,
                       process,
                       numnodes,
                       cpus_per_node,
                       alpgen_input_file,
                       pdf_filename,
                       walltime,
                       site,
                       dev                = False,
                       regen_path         = None,
                       submit             = True,
                       iseed1             = None,
                       iseed2             = None,
                       iseed3             = None,
                       iseed4             = None,
                       ):
   



if __name__ == '__main__':
   sys.exit(main())
