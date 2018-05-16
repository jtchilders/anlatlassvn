
      function random()
      real * 8 random
      real * 8 saverandom
      logical fixed
      COMMON/tmpfixed/fixed
      data fixed/.false./
      save saverandom
      if(fixed) then
         random=saverandom
         return
      endif
      call rm48(random,1)
      saverandom=random
      end


      subroutine resetrandom
      call RM48IN(54217137,0,0)
      end

      subroutine randomsave
      implicit none
      integer j,ipar(3,10)
      data j/0/
      save j,ipar
      j=j+1
      if(j.gt.10) then
         write(*,*) ' Too many recursive calls to randomsave'
         stop
      endif
      call rm48ut(ipar(1,j),ipar(2,j),ipar(3,j))
      return
      entry randomrestore
      if(j.le.0) then
         write(*,*) ' Too many calls to randomrestore'
         stop
      endif
      call rm48in(ipar(1,j),ipar(2,j),ipar(3,j))
      j=j-1
      return
      end

      subroutine readcurrentrandom(i1,n1,n2)
      implicit none
      integer i1,n1,n2
      call rm48ut(i1,n1,n2)
      end

      subroutine setrandom(i1,n1,n2)
      implicit none
      integer i1,n1,n2
      integer i1cur,n1cur,n2cur,j
      real * 8 tmp
      if(i1.eq.0) then
c This is used for complete initialization
         call resetrandom
         return
      endif
c Reinitializing the random number may be expensive;
c If we need a sequence number greater than the current status
c just call the generator enough times to get there. 
      call rm48ut(i1cur,n1cur,n2cur)
      if(i1.eq.i1cur.and.n2.eq.n2cur.and.n1.ge.n1cur) then
         do j=n1cur,n1-1
            call rm48(tmp,1)
         enddo
         call rm48ut(i1cur,n1cur,n2cur)
         if(i1.eq.i1cur.and.n2.eq.n2cur.and.n1.eq.n1cur) then
c Succeded
            return
         else
c Failed!
            write(*,*) ' setrandom: debug ...'
            call exit(-1)
         endif
      endif
c reinitialize from scratch
      if (I1.gt.0.and.n1.ge.0.and.n2.ge.0) then
c     restart a previous run
            call rm48in(I1,N1,N2)
      else
         write(*,*) 'ERROR: setrandom called with',i1,n1,n2
         call exit(-1)
      endif
      end


      subroutine savecurrentrandom
      implicit none
      integer ipar(3)
      common/crandom/ipar
      call rm48ut(ipar(1),ipar(2),ipar(3))
      end


      subroutine getcurrentrandom(i1,n1,n2)
      implicit none
      integer i1,n1,n2
      integer ipar(3)
      common/crandom/ipar
      i1 = ipar(1)
      n1 = ipar(2)
      n2 = ipar(3)
      end

      subroutine printcurrentrandom
      implicit none
      integer ipar(3)
      common/crandom/ipar
      write(*,*) 'Random number seeds: ',ipar(1),ipar(2), ipar(3)
      end
