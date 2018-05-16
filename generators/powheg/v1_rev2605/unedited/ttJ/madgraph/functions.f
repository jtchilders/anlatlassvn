      subroutine switchmom(p1,p,ic,jc,n)
c**************************************************************************
c     Changes stuff for crossings
c**************************************************************************
      implicit none
      integer n
      integer jc(n),ic(n)
      real*8 p1(0:3,n),p(0:3,n)
      integer i,j
c-----
c Begin Code
c-----
      do i=1,n
         do j=0,3
            p(j,ic(i))=p1(j,i)
         enddo
      enddo
      do i=1,n
         jc(i)=1
      enddo
      jc(ic(1))=-1
      jc(ic(2))=-1
      end
     
      function flavequiv_perm(n,aflav,bflav,perm)
      implicit none
      logical flavequiv_perm
      integer n, aflav(n),bflav(n),perm(n)
      include '../nlegborn.h'
      include '../../include/pwhg_flst.h'
      integer j,k,kb,itmp,ib(nlegreal-2+6)
      call intassign(n,bflav,ib)
      do j=1,n
         perm(j)=j
      enddo
      do j=1,n
         if(aflav(j).ne.ib(j)) then
            if(j.le.2) then
               flavequiv_perm=.false.
               return
            endif
            do k=j+1,n
               if(aflav(j).eq.ib(k)) then
                  itmp=ib(j)
                  ib(j)=ib(k)
                  ib(k)=itmp
                  itmp=perm(j)
                  perm(j)=perm(k)
                  perm(k)=itmp
                  goto 10
               endif
            enddo
            flavequiv_perm=.false.
            return
         endif
 10      continue
      enddo
      flavequiv_perm=.true.
      end

