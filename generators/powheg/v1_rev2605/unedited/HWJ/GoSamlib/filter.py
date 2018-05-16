# vim: ts=3:sw=3:expandtab

def zero_loop(d):
   """
   Loops which cancel due to Furry's theorem.
   """
   return d.chord(QUARKS) == d.loopsize() and \
          d.loopsize() == 3 and \
          d.loopvertices([A], QUARKS, QUARKS) == 1 and \
          d.loopvertices([g], QUARKS, QUARKS) == 2
          
def v_floop(d):
   """
   Closed quark loops with gauge boson attached to the loop
   """
   return d.chord(QUARKS) == d.loopsize() and \
          d.loopvertices([Z,A], QUARKS, QUARKS) >= 1 \
          and not zero_loop(d)

def no_hff(d):
   """
   No Higgs attached to massless fermions.
   """
   return d.vertices([H], [U,D,S,C,B], [Ubar,Dbar,Sbar,Cbar,Bbar]) == 0

def top_loop(d):
   return d.chord([T,Tbar]) == d.loopsize()

def b_loop(d):
   return d.chord([B,Bbar]) == d.loopsize()
