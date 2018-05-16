import logging
logger = logging.getLogger(__name__)


class Computer:
# Collects information on HPC parameters

   def __init__(self,
                name,
                architecture, 
                min_nodes,
                max_nodes,
                cores_per_node,
               ):
   
      if architecture not in ('x86','bgp','bgq'):
         logger.error('Unknown computer architecture')
         return

      self.name = name
      self.architecture = architecture
      self.min_nodes = min_nodes
      self.max_nodes = max_nodes
      self.cores_per_node = cores_per_node
      self.cpu_hours = None

      self.shapes = []
      nodes = self.min_nodes
      while nodes <= self.max_nodes:
         self.shapes.append(
                            Shape(
                                  self,
                                  nodes,
                                  self.cores_per_node,
                                  0,
                                  1
                                 )
                           )
         nodes = nodes * 2



class Shape:
# A place to store information on the various
#  job shapes.  It's convenient to do some
#  XML building here.

   def __init__(self,
                computer,
                nodes,
                cores_used,
                events,
                sequence
               ):
      self.computer        = computer
      self.nodes           = nodes
      self.cores_used      = cores_used
      self.events          = events
      self.sequence        = sequence
      self.flag            = self.get_flag()
      self.prescript       = None
      self.prescript_args  = None
      self.postscript      = None
      self.postscript_args = None

   def get_flag(self):
      if self.computer.architecture == 'bgp':
         if self.cores_used == 4:
            return '--mode=vn'
         elif self.cores_used == 2:
            return '--mode=dual'
         elif self.cores_used == 1:
            return '--mode=smp'
      elif self.computer.architecture == 'bgq':
         if self.cores_used > 8:
            return '--mode c16'
         elif self.cores_used > 4:
            return '--mode c8'
         elif self.cores_used > 2:
            return '--mode c4'
         elif self.cores_used == 2:
            return '--mode c2'
         elif self.cores_used == 1:
            return '--mode c1'
         else:
            logger('flag not found for cores_used = ' + self.cores_used + ' and architecture = ' + self.computer.architecture)
            sys.exit(-1)
      
      return None

computers = {}
computers['frontend']      = Computer('frontend','x86',1,1,1)
computers['intrepid']      = Computer('intrepid','bgp',512,4096,4)
computers['challenger']    = Computer('challenger','bgp',16,1024,4)
computers['surveyor']      = Computer('surveyor','bgp',64,1024,4)
computers['mira']          = Computer('mira','bgq',512,4096,64)
computers['cetus']         = Computer('cetus','bgq',128,1024,64)
computers['vesta']         = Computer('vesta','bgq',1,1024,64)
computers['vesta-dev']     = Computer('vesta-dev','bgq',1,1024,64)
computers['argo_cluster']  = Computer('argo_cluster','x86',1,5,16)
computers['testmachine']   = Computer('testmachine','x86',1,5,16)

def found_computer(name):
   if name in computers.keys():
      return True
   return False

def check_computer_config(name,nodes,processes_per_node):
   computer = computers[name]
   if computer.min_nodes < nodes and nodes > computer.max_nodes:
      logger.error('Job is requesting more nodes ('+str(nodes)+') than allowed ('+str(computer.min_nodes)+'-'+str(computer.max_nodes)+') for this machine ('+name+')')
      return False
   if processes_per_node > computer.cores_per_node:
      logger.error('Job is requestiong more cores per node('+str(processes_per_node)+') than allowed ('+str(computer.cores_per_node)+') for this machine ('+name+')')
      return False
   return True

