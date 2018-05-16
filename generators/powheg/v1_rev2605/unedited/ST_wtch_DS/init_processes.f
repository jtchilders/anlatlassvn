      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'LesHouches.h'
      include 'PhysPars.h'

      integer ipl,imn,i1,i2,i3
      integer nflav
      parameter (nflav=5)

      integer three_ch(-6:6)
      data three_ch /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      character *3 qcd_name(-6:6)
      data qcd_name /'tx ','bx ','cx ','sx ','ux ','dx ',
     #         'g  ','d  ','u  ','s  ','c  ','b  ','t  '/

      logical debug
      parameter (debug=.false.)

      logical condition_b_wt
      logical condition_blike_wt,condition_1_wt,condition_2_wt
      logical condition_wt
      integer l(3)

      integer j,k
      real *8 powheginput
      external powheginput

      lprup(1)=10000


ccccccccccccccccccccccccccccccc
c     decide t or tbar process
      ttype=powheginput('ttype')
      if(abs(ttype).ne.1) then
         write(*,*) 'Unrecognised ttype in input file'
         write(*,*) 'admitted values: 1 for t, -1 for tbar'
         call exit(1)
      endif
ccccccccccccccccccccccccccccccc    


cccccccccccccccccccccccccccccc
c     wt-channel subprocesses
cccccccccccccccccccccccccccccc


cccccccccccccccccccc
c     index of the first MASSLESS coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=5
cccccccccccccccccccc


*********************************************************************
***********            BORN SUBPROCESSES              ***************
*********************************************************************
      write(*,*) 'POWHEG: looking for Born amplitudes'
      flst_nborn=0
      i1=-24
      i2=6
      do ipl=-nflav,nflav
         do imn=-nflav,nflav
            condition_b_wt=.false.
c     flavour conditions to yield a valid amplitude
            condition_b_wt=(
     $((three_ch(ipl) + three_ch(imn)).eq.
     $((-3) + three_ch(i2))).and.
     $(ipl*imn.eq.0))
            if(condition_b_wt) then
               flst_nborn=flst_nborn+1
               if(flst_nborn.gt.maxprocborn) goto 999
               flst_born(1,flst_nborn)=ipl
               flst_born(2,flst_nborn)=imn
               flst_born(3,flst_nborn)=i1
               flst_born(4,flst_nborn)=i2
            endif
         enddo
      enddo
      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif

*********************************************************************
***********            REAL SUBPROCESSES              ***************
*********************************************************************
      write(*,*) 'POWHEG: looking for real amplitudes'
      flst_nreal=0
      do ipl=-nflav,nflav
         do imn=-nflav,nflav
            do i3=-nflav,nflav
c     charge conservation
               if(three_ch(ipl)+three_ch(imn).eq.
     $(-3) + three_ch(i2) + three_ch(i3)) then
c     flavour conditions to yield a valid amplitude
c     born-like + outgoing gluon (leg 5, i.e. i3=0)
                  condition_blike_wt=(
     $                 (i3.eq.0).and.
     $                 (((three_ch(ipl) + three_ch(imn)).eq.
     $                 ((-3) + three_ch(i2))).and.
     $                 (ipl*imn.eq.0)))
c     (ipl+imn=0)
                  condition_1_wt=(
     $                 (ipl + imn.eq.0).and.
     $                 (three_ch(i2)+three_ch(i3)+ (-3).eq.0))
c     (ipl or imn is the same of i3, 
c     ipl*imn and ipl+imn is not zero)
                  condition_2_wt=(
     $                 (ipl*imn.ne.0).and.
     $                 (ipl+imn.ne.0).and.
     $                 ((ipl.eq.i3).or.
     $                 (imn.eq.i3)))

                  condition_wt=condition_blike_wt.or.condition_1_wt.or.
     $                 condition_2_wt
                  
                  if(condition_wt) then
c     check that the conditions of t-channel are reciprocally exclusive
                     do j=1,3
                        l(j)=0
                     enddo
                     if(condition_blike_wt) l(1)=1
                     if(condition_1_wt) l(2)=1
                     if(condition_2_wt) l(3)=1
                     if(l(1)+l(2)+l(3).ne.1) then
                        write(*,*) 
     #'problem in looking for real subprocesses (wt-ch)'
c$$$                        write(*,*) l
c$$$                        write(*,*) ipl,imn,i1,i2,i3
                        call exit(1)
                     endif

                     flst_nreal=flst_nreal+1
                     if(flst_nreal.gt.maxprocreal) goto 998
                     
                     flst_real(1,flst_nreal)=ipl
                     flst_real(2,flst_nreal)=imn
                     flst_real(3,flst_nreal)=i1
                     flst_real(4,flst_nreal)=i2
                     flst_real(5,flst_nreal)=i3
                  endif
               endif
            enddo
         enddo
      enddo

      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
         enddo
      endif

ccccccccccccccccccccccccccccccccccccccc
c     charge conjugation
      if(ttype.eq.-1) then
         do j=1,flst_nborn
            flst_born(1,j)=-flst_born(1,j)
            flst_born(2,j)=-flst_born(2,j)
            flst_born(3,j)=-flst_born(3,j)
            flst_born(4,j)=-flst_born(4,j)
         enddo
         do j=1,flst_nreal
            flst_real(1,j)=-flst_real(1,j)
            flst_real(2,j)=-flst_real(2,j)
            flst_real(3,j)=-flst_real(3,j)
            flst_real(4,j)=-flst_real(4,j)
            flst_real(5,j)=-flst_real(5,j)
         enddo
      endif
cccccccccccccccccccccccccccccccccccccc

      return
 998  write(*,*) 'init_processes: increase maxprocreal'
      call exit(1)
 999  write(*,*) 'init_processes: increase maxprocborn'
      call exit(1)
      end
 

