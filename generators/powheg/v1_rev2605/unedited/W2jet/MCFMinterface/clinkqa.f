c Subroutine to link colours for
c quark-antiquark
c in planar order
      subroutine clinkqa(ic1,ic2)
      integer iclabel
      common/ciclabel/iclabel
      integer ic1(2),ic2(2)
c ic1 is a quark: has colour, zero anticolor
      ic1(1)=iclabel+2
      ic1(2)=0
c ic2 is an antiquark: has anticolour, zero color
      ic2(1)=0
      ic2(2)=iclabel+2
      iclabel=iclabel+10
      end


