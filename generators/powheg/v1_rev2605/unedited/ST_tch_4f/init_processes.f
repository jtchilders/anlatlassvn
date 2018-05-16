      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'LesHouches.h'
      include 'PhysPars.h'
      
      integer ipl,imn,i1,i2,iborn,ileg,iflav
      integer nflav
      parameter (nflav=4)

      integer three_ch(-6:6)
      data three_ch /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      character *3 qcd_name(-6:6)
      data qcd_name /'tx ','bx ','cx ','sx ','ux ','dx ',
     #         'g  ','d  ','u  ','s  ','c  ','b  ','t  '/

      logical debug
      parameter (debug=.false.)

      logical condition_b_t,
     #condition_blike_t,condition_g1_t,condition_g2_t,condition_t 

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


c$$$      print*, '********************************************'
c$$$      print*, '********************************************'
c$$$      print*, '***********      WARNING       *************'
c$$$      print*, 'POSSIBLE SUBTLETIES WITH NFLAV AND st_nlight'
c$$$      print*, '               SEE WBB CODE '
c$$$      print*, '***********      WARNING       *************'
c$$$      print*, '********************************************'
c$$$      print*, '********************************************'


cccccccccccccccccccccccccccccc
c     t-channel subprocesses (4f)
c     3rd is the TOP, 4th the BOTTOM
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
      do ipl=-nflav,nflav
         do imn=-nflav,nflav
            do i1=-nflav,nflav
               condition_b_t=.false.
c     flavour conditions to yield a valid amplitude
c     gluon in position 1:
               condition_b_t=(
     $              (three_ch(ipl).eq.0.and.
     $              (three_ch(i1)-three_ch(imn)).eq.-3).or.
     $              (three_ch(imn).eq.0.and.
     $              (three_ch(i1)-three_ch(ipl)).eq.-3))
               if(condition_b_t) then
                  flst_nborn=flst_nborn+1
                  if(flst_nborn.gt.maxprocborn) goto 999
                  flst_born(1,flst_nborn)=ipl
                  flst_born(2,flst_nborn)=imn
                  flst_born(3,flst_nborn)=6
                  flst_born(4,flst_nborn)=-5
                  flst_born(5,flst_nborn)=i1
               endif
            enddo
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


      do iborn=1,flst_nborn
c     emission of an extra gluon
         flst_nreal=flst_nreal+1
         do ileg=1,nlegborn
            flst_real(ileg,flst_nreal)=flst_born(ileg,iborn)
         enddo
         flst_real(nlegreal,flst_nreal)=0
c     case Q G
         if(flst_born(2,iborn).eq.0) then
c     on the massless current, introduce a splitting
c     g->qq, that matches the q entering the Born diagram
            flst_nreal=flst_nreal+1
            flst_real(1,flst_nreal)=0
            flst_real(2,flst_nreal)=0
            do ileg=3,5
               flst_real(ileg,flst_nreal)=flst_born(ileg,iborn)
            enddo
            flst_real(6,flst_nreal)=-flst_born(1,iborn)
            if(flst_real(6,flst_nreal).lt.flst_real(5,flst_nreal)) then
               flst_nreal=flst_nreal-1
            endif
c     make the incoming gluon in the Born to be radiated
c     by a quark line
            do iflav=-nflav,nflav
               if(iflav.ne.0) then
                  flst_nreal=flst_nreal+1
                  flst_real(1,flst_nreal)=flst_born(1,iborn)
                  flst_real(2,flst_nreal)=iflav
                  do ileg=3,5
                     flst_real(ileg,flst_nreal)=flst_born(ileg,iborn)
                  enddo
                  flst_real(nlegreal,flst_nreal)=iflav
               endif
            enddo
         elseif(flst_born(1,iborn).eq.0) then
c     case G Q


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccc NOT NEEDED, ALREADY TAKEN INTO ACCOUNT BEFORE ccccccccc
c$$$c     on the massless current, introduce a splitting
c$$$c     g->qq, that matches the q entering the Born diagram
c$$$            flst_nreal=flst_nreal+1
c$$$            flst_real(1,flst_nreal)=0
c$$$            flst_real(2,flst_nreal)=0
c$$$            do ileg=3,5
c$$$               flst_real(ileg,flst_nreal)=flst_born(ileg,iborn)
c$$$            enddo
c$$$            flst_real(6,flst_nreal)=-flst_born(2,iborn)
c$$$            if(flst_real(6,flst_nreal).lt.flst_real(5,flst_nreal)) then
c$$$               flst_nreal=flst_nreal-1
c$$$            endif
ccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     make the incoming gluon in the Born to be radiated
c     by a quark line
c     If initial state is qq or qbarqbar, it was already included
c     in the Q G case !
            do iflav=-nflav,nflav
               if(iflav.ne.0) then
                  if(iflav.ne.flst_born(2,iborn)) then
                     flst_nreal=flst_nreal+1
                     flst_real(2,flst_nreal)=flst_born(2,iborn)
                     flst_real(1,flst_nreal)=iflav
                     do ileg=3,5
                        flst_real(ileg,flst_nreal)=flst_born(ileg,iborn)
                     enddo
                     flst_real(nlegreal,flst_nreal)=iflav
                  endif
               endif
            enddo


         endif
      enddo


      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
c$$$            if(flst_real(4,j).gt.flst_real(5,j)) then
c$$$               write(*,*) '                 |__ reverse order'
c$$$            endif
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
            flst_born(5,j)=-flst_born(5,j)
         enddo
         do j=1,flst_nreal
            flst_real(1,j)=-flst_real(1,j)
            flst_real(2,j)=-flst_real(2,j)
            flst_real(3,j)=-flst_real(3,j)
            flst_real(4,j)=-flst_real(4,j)
            flst_real(5,j)=-flst_real(5,j)
            flst_real(6,j)=-flst_real(6,j)
         enddo
      endif
cccccccccccccccccccccccccccccccccccccc

      return
 998  write(*,*) 'init_processes: increase maxprocreal'
      call exit(1)
 999  write(*,*) 'init_processes: increase maxprocborn'
      call exit(1)
      end
 

