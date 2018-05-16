
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


      subroutine setrandom(i1,n1,n2)
      implicit none
      integer i1,n1,n2
      
      if (I1.gt.0) then
         if (((N1.gt.0).and.(N2.ge.0)).or.(N1.ge.0).and.(N2.gt.0)) then
c     restart a previous run or start a new run with this initialization
            call rm48in(I1,N1,N2)
         else
c     just change the random seed
            call rm48in(I1,0,0)
         endif
      else
         call resetrandom
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
