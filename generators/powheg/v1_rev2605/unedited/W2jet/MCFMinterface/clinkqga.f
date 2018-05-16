c Subroutine to link colours for
c quark - gluon - aquark
c in planar order
      subroutine clinkqga(ic1,ic2,ic3)
      integer iclabel
      common/ciclabel/iclabel
      integer ic1(2),ic2(2),ic3(2)
c ic1 is a quark: has colour, zero anticolor
      ic1(1)=iclabel+2
      ic1(2)=0
c ic2 is a gluon: link to quark
      ic2(1)=iclabel+3
      ic2(2)=iclabel+2
c ic3 is an antiquark
      ic3(1)=0
      ic3(2)=iclabel+3
      iclabel=iclabel+10
      end





