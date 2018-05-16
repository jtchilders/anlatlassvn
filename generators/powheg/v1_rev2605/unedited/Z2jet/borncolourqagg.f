      subroutine borncolourqagg(icol1,icol2,icol3,icol4)
c                                q     qbar  g     g
      implicit none
      include 'qq_cs.f'
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      double precision r,random
      
!      write(6,*) 'borncolourqagg: qq_cs',qq_cs

      r=random()*(qq_cs(1)+qq_cs(2))      
      
      if(r.lt.qq_cs(1)) then
         call colourjoin2g(icol1,icol2,icol3,icol4)
      else
         call colourjoin2g(icol1,icol2,icol4,icol3)
      endif
      end

      subroutine colourjoin2g(icol1,icol2,icol3,icol4)
c                             q     qbar  g     g
c perform a planar colour connection on the planar sequence
c q qbar g g
c  q[1](icol1),g[2](icol4)
c  g[1](icol4),g[2](icol3)
c  g[1](icol3),qbar[2](icol2)
      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      integer newcolor
      icol1(2)=0
      icol2(1)=0
      call getnewcolor(newcolor)
      icol1(1)=newcolor
      icol4(2)=newcolor
      call getnewcolor(newcolor)
      icol4(1)=newcolor
      icol3(2)=newcolor
      call getnewcolor(newcolor)
      icol3(1)=newcolor
      icol2(2)=newcolor
      end

      
