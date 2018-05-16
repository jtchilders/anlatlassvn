      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_flst.h'
      include '../include/pwhg_kn.h'
      include '../include/pwhg_flg.h'
      integer i1,i2,i3,i4,i5,i6,i7,k,ii(7)
      equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),
     #  (i4,ii(4)),(i5,ii(5)),(i6,ii(6)),(i7,ii(7))
      integer j,tmp
      logical flavequiv,condition
      external flavequiv
      logical debug
      parameter (debug=.false.)

ccccccccccccccccccccccccccccc
c     !:
      integer nflav
      parameter (nflav=5)
cccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccc
c     from born flavour structure to jborn (for simplicity, here
c     use only MASSLESS colored particles, sorted as:
c     id_plus,id_minus,id_final1,id_final2)
c     Particle 3 and 4 are the 2 leptons
      integer bornflst2pwhgcode(-5:5,-5:5,-5:5,-5:5)
      common/cbornflst2pwhgcode/bornflst2pwhgcode
ccccccccccccccccccccccccccccccccccccc


c     check nlegborn. This is only a sanity check while we are TESTING 
c     the code and we change often from one process to the other
      if (nlegborn.ne.6) then
         write(*,*) ' ERROR: set nlegborn to the appropriate value'
         write(*,*) ' for this process in nlegborn.h'
         call exit(1)
      endif

      if(nflav.ne.5) then
         write(*,*) '***********************************'
         write(*,*) '*** RUNNING WITH nflav = ',nflav
         write(*,*) '***********************************'
      endif      

*********************************************************************
c     index of the first LIGHT coloured parton in the final state
      flst_lightpart=5
*********************************************************************

      i3=11
      i4=-11

*********************************************************************
***********            BORN SUBPROCESSES              ***************
*********************************************************************
      flst_nborn=0
      do i1=-nflav,nflav
         do i2=-nflav,nflav
            if(i1*i2.gt.0) then
c both quarks or antiquarks
               i5=i1
               i6=i2
               goto 10
            endif
c a quark and a gluon
            if(i1*i2.eq.0) then
               if(i1.ne.0) then
                  i5=i1
                  i6=i2
                  goto 10
               elseif(i2.ne.0) then
                  i5=i2
                  i6=i1
                  goto 10
               endif
            endif
c two gluons or quark antiquark from here on
            if(i1+i2.ne.0) then
c different quark antiquark
               i5=i1
               i6=i2
               goto 10
            endif
c two gluons or quark antiquark of the same flavour
            do i5=0,nflav
               i6=-i5
c     exclude gg -> Z gg
               if (i1.eq.0.and.i5.eq.0) goto 1234                  
               flst_nborn=flst_nborn+1
               if(flst_nborn.gt.maxprocborn) goto 999
               do k=1,6
                  flst_born(k,flst_nborn)=ii(k)
               enddo
 1234          continue
            enddo
            goto 11
 10         continue
            flst_nborn=flst_nborn+1
            if(flst_nborn.gt.maxprocborn) goto 999
            do k=1,6
               flst_born(k,flst_nborn)=ii(k)
            enddo
 11         continue
         enddo
      enddo
      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif

cccccccccccccccccccccccccccccccccccccccccc
c     Needed to use an OLP program
cccccccccccccccccccccccccccccccccccccccccc
      do j=1,flst_nborn
         i1=flst_born(1,j)
         i2=flst_born(2,j)
         i5=flst_born(5,j)
         i6=flst_born(6,j)
         bornflst2pwhgcode(i1,i2,i5,i6)=j
      enddo

*********************************************************************
***********            REAL SUBPROCESSES              ***************
*********************************************************************
      flst_nreal=0
c     4q1g - qq type
      i7=0
      do i1=-nflav,nflav
         do i2=-nflav,nflav
            do i5=-nflav,nflav
               do i6=-nflav,nflav
                  condition=
     $                 (i1.ne.0).and.
     $                 (abs(i1).eq.abs(i2)).and.
     $                 (abs(i1).eq.abs(i5)).and.
     $                 (abs(i1).eq.abs(i6))
                  condition=condition.and.
     $                 (((i1.eq.i5).and.(i2.eq.i6)).or.
     $                 ((i1.eq.i6).and.(i2.eq.i5)).or.
     $                 ((i1.eq.i2).and.(i1.eq.i5).and.(i1.eq.i6)))
                  if(condition) then
                     do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                        if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                       goto 20
                     enddo
                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     do k=1,nlegreal
                        flst_real(k,flst_nreal)=ii(k)
                     enddo
 20                  continue
                  endif
               enddo
            enddo
         enddo
      enddo
      if(debug) write(*,*) '# of 4q1g - qq type: ',flst_nreal
      tmp=flst_nreal
  
c     4q1g - gq type
      i1=0
      do i2=-nflav,nflav
         do i5=-nflav,nflav
            do i6=-nflav,nflav
               do i7=-nflav,nflav
                  condition=
     $                 (i2.ne.0).and.
     $                 (i2.eq.i5).and.
     $                 (abs(i6).eq.abs(i2)).and.
     $                 (i6+i7.eq.0)
                  if(condition) then
                     do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                        if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                       goto 21
                     enddo
                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     do k=1,nlegreal
                        flst_real(k,flst_nreal)=ii(k)
                     enddo
 21                  continue
                  endif
               enddo
            enddo
         enddo
      enddo
      if(debug) write(*,*) '# of 4q1g - gq type: ',flst_nreal-tmp
      tmp=flst_nreal

c     4q1g - qg type
      i2=0
      do i1=-nflav,nflav
         do i5=-nflav,nflav
            do i6=-nflav,nflav
               do i7=-nflav,nflav
                  condition=
     $                 (i1.ne.0).and.
     $                 (i1.eq.i5).and.
     $                 (abs(i6).eq.abs(i1)).and.
     $                 (i6+i7.eq.0)
                  if(condition) then
                     do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                        if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                       goto 22
                     enddo
                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     do k=1,nlegreal
                        flst_real(k,flst_nreal)=ii(k)
                     enddo
 22                  continue
                  endif
               enddo
            enddo
         enddo
      enddo
      if(debug) write(*,*) '# of 4q1g - qg type: ',flst_nreal-tmp
      tmp=flst_nreal

c     2q2qprime1g - qq type
      i7=0
      do i1=-nflav,nflav
         do i2=-nflav,nflav
            do i5=-nflav,nflav
               do i6=-nflav,nflav
                  condition=
     $                 (i1*i2*i5*i6.ne.0).and.
     $                 (((i1.eq.i5).and.(i2.eq.i6).and.(i1.ne.i2)).or.
     $                 ((i1.eq.i6).and.(i2.eq.i5).and.(i1.ne.i2)).or.
     $        ((i1+i2.eq.0).and.(i5+i6.eq.0).and.(abs(i1).ne.abs(i5))))
                  if(condition) then
                     do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                        if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                       goto 23
                     enddo
                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     do k=1,nlegreal
                        flst_real(k,flst_nreal)=ii(k)
                     enddo
 23                  continue
                  endif
               enddo
            enddo
         enddo
      enddo
      if(debug) write(*,*) '# of 2q2qprime1g - qq type: ',flst_nreal-tmp
      tmp=flst_nreal

c     2q2qprime1g - gq type
      i1=0
      do i2=-nflav,nflav
         do i5=-nflav,nflav
            do i6=-nflav,nflav
               do i7=-nflav,nflav
                  condition=
     $                 (i2.ne.0).and.
     $                 (i2.eq.i5).and.
     $                 (abs(i6).ne.abs(i2)).and.
     $                 (i6*i7.ne.0).and.
     $                 (i6+i7.eq.0)
                  if(condition) then
                     do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                        if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                       goto 24
                     enddo
                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     do k=1,nlegreal
                        flst_real(k,flst_nreal)=ii(k)
                     enddo
 24                  continue
                  endif
               enddo
            enddo
         enddo
      enddo
      if(debug) write(*,*) '# of 2q2qprime1g - gq type: ',flst_nreal-tmp
      tmp=flst_nreal

c     2q2qprime1g - qg type
      i2=0
      do i1=-nflav,nflav
         do i5=-nflav,nflav
            do i6=-nflav,nflav
               do i7=-nflav,nflav
                  condition=
     $                 (i1.ne.0).and.
     $                 (i1.eq.i5).and.
     $                 (abs(i6).ne.abs(i1)).and.
     $                 (i6*i7.ne.0).and.
     $                 (i6+i7.eq.0)
                  if(condition) then
                     do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                        if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                       goto 25
                     enddo
                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     do k=1,nlegreal
                        flst_real(k,flst_nreal)=ii(k)
                     enddo
 25                  continue
                  endif
               enddo
            enddo
         enddo
      enddo
      if(debug) write(*,*) '# of 2q2qprime1g - qg type: ',flst_nreal-tmp
      tmp=flst_nreal

c     2q3g - qq type
      i5=0
      i6=0
      i7=0
      do i1=-nflav,nflav
         do i2=-nflav,nflav
            condition=
     $           (i1*i2.ne.0).and.
     $           (i1+i2.eq.0)
            if(condition) then
               do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                  if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                 goto 26
               enddo
               flst_nreal=flst_nreal+1
               if(flst_nreal.gt.maxprocreal) goto 998
               do k=1,nlegreal
                  flst_real(k,flst_nreal)=ii(k)
               enddo
 26            continue
            endif
         enddo
      enddo
      if(debug) write(*,*) '# of 2q3g - qq type: ',flst_nreal-tmp
      tmp=flst_nreal

c     2q3g - qg/gq type
      i6=0
      i7=0
      do i1=-nflav,nflav
         do i2=-nflav,nflav
            do i5=-nflav,nflav
               condition=
     $              (i1*i2.eq.0).and.
     $              (i5.ne.0).and.
     $              (((i1.eq.0).and.(i2.eq.i5)).or.
     $              ((i2.eq.0).and.(i1.eq.i5)))
               if(condition) then
                  do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                     if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                    goto 27
                  enddo
                  flst_nreal=flst_nreal+1
                  if(flst_nreal.gt.maxprocreal) goto 998
                  do k=1,nlegreal
                     flst_real(k,flst_nreal)=ii(k)
                  enddo
 27               continue
               endif
            enddo
         enddo
      enddo
      if(debug) write(*,*) '# of 2q3g - qg/gq type: ',flst_nreal-tmp
      tmp=flst_nreal

c     2q3g - gg type
      i1=0
      i2=0
      i7=0
      do i5=-nflav,nflav
         do i6=-nflav,nflav
            condition=
     $           (i5*i6.ne.0).and.
     $           (i5+i6.eq.0)
            if(condition) then
               do j=1,flst_nreal
c     Check that an inequivalent configuration is generated
                  if(flavequiv(nlegreal,flst_real(1,j),ii(1)))
     $                 goto 28
               enddo
               flst_nreal=flst_nreal+1
               if(flst_nreal.gt.maxprocreal) goto 998
               do k=1,nlegreal
                  flst_real(k,flst_nreal)=ii(k)
               enddo
 28            continue
            endif
         enddo
      enddo
      if(debug) write(*,*) '# of 2q3g - gg type: ',flst_nreal-tmp
      tmp=flst_nreal

      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
         enddo
      endif
      return
 998  write(*,*) 'init_processes: increase maxprocreal'
      stop
 999  write(*,*) 'init_processes: increase maxprocborn'
      end
 




