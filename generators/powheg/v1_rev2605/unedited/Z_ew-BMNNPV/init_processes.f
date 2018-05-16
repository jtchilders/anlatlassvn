      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_par.h'
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      integer i1,i2,i3,i4,i5,k,ii(nlegreal)
      equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),
     #  (i4,ii(4)),(i5,ii(5))
      logical debug
      parameter (debug=.true.)
      integer j
      logical condition
      real * 8 powheginput
      external powheginput
c     vector boson id and decay
      integer vdecaymode,tmp
      common/cvecbos/vdecaymode
c     lepton masses
      real *8 ktmaxqed
      common/showerqed/ktmaxqed
      logical ifphotoninduced
      integer phind

c******************************************************
c     Choose the process to be implemented
c******************************************************
c   decay products of the vector boson
      vdecaymode=powheginput('vdecaymode')
c   consider photon induced processes
      phind=powheginput('#photoninduced')
      if ( phind.lt.0 ) phind = -100000

      if ( phind.eq.1 ) then 
          ifphotoninduced=.true.
          ifdis=.true.
c          ifdis=.false.
          write(*,*) 'photon induced processes on, only MRSTQED2004'
      else
          ifphotoninduced=.false.
          ifdis=.false.
      endif
 
      par_isrtinycsi = 1d-11
      par_isrtinyy = 1d-11

      par_fsrtinycsi = 1d-15
      par_fsrtinyy = 1d-15

      ktmaxqed = powheginput("#ktmaxqed")
      if (ktmaxqed.le.0d0) ktmaxqed  = 0.01d0     
c      flg_jacsing = .true.

      if ((vdecaymode.lt.11).or.(vdecaymode.gt.16)) then
         write(*,*) 'ERROR: The decay mode you selected'
     #  //' is not allowed (Up to now only leptonic decays)'
         stop
      endif
            
      write(*,*) 
      write(*,*) ' POWHEG: Single Z production and decay'
      if (vdecaymode.eq.11) then 
          write(*,*) '         to e- e+ '
      elseif (vdecaymode.eq.13) then 
          write(*,*) '         to mu- mu+ '
      elseif (vdecaymode.eq.15) then
          write(*,*) '         to tau- tau+ '
      else
c     not yet implemented
         write(*,*) 'non leptonic Z decays '//
     #        'not yet implemented'
         stop
      endif
      write(*,*) 
      
c     change the LHUPI id of the process according to vector boson id
c     and decay
      lprup(1)=10000+vdecaymode ! 10000+idup of charged decay product of the W
      
c*********************************************************     
c
c     index of the first light particle in the final state
c     that can give rise to collinear singularities;
c     The charged leptons are considered MASSIVE particles in this code,
c     and thus do not give rise to collinear singularities
      flst_lightpart=5
c     Z decay products
      i3=vdecaymode
      i4=-i3

*********************************************************************
***********            BORN SUBPROCESSES              ***************
*********************************************************************
      flst_nborn=0
      do i1=-5,5
         do i2=-5,5
            if(i1.ne.0.and.i1+i2.eq.0) then
c     q qbar
               flst_nborn=flst_nborn+1
               if(flst_nborn.gt.maxprocborn) goto 999
               do k=1,nlegborn
                  flst_born(k,flst_nborn)=ii(k)
               enddo
            endif
         enddo
      enddo
      if (ifphotoninduced) then
          i1=22
          i2=22
          flst_nborn=flst_nborn+1
          if(flst_nborn.gt.maxprocborn) goto 999
          do k=1,nlegborn
              flst_born(k,flst_nborn)=ii(k)
          enddo
      endif
      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif

      

*********************************************************************
***********            REAL SUBPROCESSES              ***************
*********************************************************************
      flst_nreal=0
      do i1=-5,5
         do i2=-5,5
            do i5=-5,5
               condition=.false.
               if(.not.(i1.eq.0.and.i2.eq.0)) then
c     exclude gg
                  if((i1.ne.0).and.(i1+i2.eq.0).and.(i5.eq.0)) then
c     q qbar -> g
                     condition=.true.
                  elseif((i1.eq.0).and.(i2.eq.i5)) then
c     g q
                     condition=.true.
                  elseif((i2.eq.0).and.(i1.eq.i5)) then
c     q g
                     condition=.true.
                  endif
               endif
               if (condition) then
                  flst_nreal=flst_nreal+1
                  if(flst_nreal.gt.maxprocreal) goto 998
                  do k=1,nlegreal
                     flst_real(k,flst_nreal)=ii(k)
                  enddo
               endif
            enddo
            condition=(i1+i2.eq.0.and.i1.ne.0)
            if (condition) then
               flst_nreal=flst_nreal+1
               if(flst_nreal.gt.maxprocreal) goto 998
               do k=1,nlegborn
                  flst_real(k,flst_nreal)=ii(k)
               enddo
c Photon in final state
               flst_real(nlegreal,flst_nreal)=22
            endif
 11         continue
         enddo
      enddo

c Photon in initial state
      if (ifphotoninduced) then
          i1=22
          do i2=-5,5
             do i5=-5,5
                if (i2.eq.i5.and.i2.ne.0) then
                   flst_nreal=flst_nreal+1
                   if(flst_nreal.gt.maxprocreal) goto 998
                   do k=1,nlegreal
                       flst_real(k,flst_nreal)=ii(k)
                   enddo
                endif
             enddo
          enddo

          i2=22
          do i1=-5,5
             do i5=-5,5
                if (i1.eq.i5.and.i1.ne.0) then
                   flst_nreal=flst_nreal+1
                   if(flst_nreal.gt.maxprocreal) goto 998
                   do k=1,nlegreal
                       flst_real(k,flst_nreal)=ii(k)
                    enddo
                endif
             enddo
          enddo
      endif

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
