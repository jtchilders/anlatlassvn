      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'constants.f' 
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_par.h'
      include 'LesHouches.h'
      integer i1,i2,i3,i4,i5,i6,i7,k,ii(7)
      equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),
     #  (i4,ii(4)),(i5,ii(5)),(i6,ii(6)),(i7,ii(7))
      logical debug
      parameter (debug=.true.)
      integer j
      integer charge3(-6:6)
      data charge3 /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      real * 8 powheginput
      external powheginput
c     vector boson id and decay
c     lepton masses
      include 'cvecbos.h'
      include 'vvsettings.f'
      real *8 lepmass(3),decmass,decmass1,decmass2
      common/clepmass/lepmass,decmass,decmass1,decmass2
      logical condition

c      par_isrtinycsi = 0
c      par_isrtinyy = 0
c      par_fsrtinycsi = 0
c      par_fsrtinyy = 0
c******************************************************
c     Choose the process to be implemented
c******************************************************

c     ID of vector boson produced
c     decay products of the vector boson
      vdecaymodeWp=powheginput('vdecaymodeWp')
      vdecaymodeWm=powheginput('vdecaymodeWm')

      if (powheginput("#zerowidth").eq.1) then 
         zerowidth = .true. 
         write(*,*) 'Zerowidth approximation' 
      else
         zerowidth = .false. 
         write(*,*) 'Generating off-shell W-bosons'
      endif         
      
      if (powheginput("#dronly").eq.1) then 
         dronly = .true. 
         write(*,*)  'Double resonant diagrams only'
c  cant have srdiags if zerowidth is true
      elseif (zerowidth) then
         write(*,*) 'Single resonant diagrams require off-shell W'
         write(*,*) 'Double resonant diagrams only'
         dronly = .true.
      else
         dronly = .false. 
         write(*,*)  'Including single resonant diagrams'
      endif

      if(powheginput("#withdamp").ne.0) then
         flg_withdamp = .true.
         flg_bornzerodamp = .true.
      else
         flg_withdamp = .false.
         flg_bornzerodamp = .false.
      endif

      if ((abs(vdecaymodeWp).ne.11).and.(abs(vdecaymodeWp).ne.13)
     .     .and.(abs(vdecaymodeWp).ne.15)
     .     .and.(abs(vdecaymodeWp).ne.1)   ! down quark
     .     .and.(abs(vdecaymodeWp).ne.3)   ! strange quark
     .     .and.(abs(vdecaymodeWp).ne.7)   ! all hadrons
     .                     ) then
         stop 'W decay mode should be put in terms of the charged 
     .   leptons: +/- 11, +/-13, +/-15 only, '//
     .   'or hsdrons : +/-1, +/-3, +/-7 for either +-1 or +-3'
      endif

      if ((abs(vdecaymodeWm).ne.11).and.(abs(vdecaymodeWm).ne.13)
     .     .and.(abs(vdecaymodeWm).ne.15)
     .     .and.(abs(vdecaymodeWm).ne.1)   ! down quark
     .     .and.(abs(vdecaymodeWm).ne.3)   ! strange quark
     .     .and.(abs(vdecaymodeWm).ne.7)   ! all hadrons
     .                     ) then
         stop 'W decay mode should be put in terms of the charged 
     .   leptons: +/- 11, +/-13, +/-15 only, '//
     .   'or hadrons : +/-1, +/-3, +/-7 for either +-1 or +-3'
      endif

c     change the LHUPI id of the process according to vector boson id
c     and decay
c     10000+idup of first decay product of W1 + decay product of W2
      lprup(1)=10000-100*vdecaymodeWp+vdecaymodeWm 

c we pretend that quarks are not quarks, otherwise POWHEG
c makes them radiate
      if(abs(vdecaymodeWm).le.7)
     . vdecaymodeWm = sign(1,vdecaymodeWm)*(abs(vdecaymodeWm) + 100)
      if(abs(vdecaymodeWp).le.7)
     . vdecaymodeWp = sign(1,vdecaymodeWp)*(abs(vdecaymodeWp) + 100)

      if (vdecaymodeWp.gt.0) then
         stop 'W+ decays inconsistent' 
      endif
      if (vdecaymodeWm.lt.0) then
         stop 'W- decays inconsistent' 
      endif

c -- this seems redundant, since we know the sign of vdecaymode
c -- are these even used?
      if (vdecaymodeWp.lt.0) then
         idvecbos1 = 24
         idvecbos2 = -24
      else
         idvecbos1 = -24
         idvecbos2 = 24
      endif

      write(*,*) 
      write(*,*) ' POWHEG: WW production and decay'
      if (vdecaymodeWp.eq.-11) write(*,*) '         to e+ ve  '
      if (vdecaymodeWp.eq.-13) write(*,*) '         to mu+ vm '
      if (vdecaymodeWp.eq.-15) write(*,*) '         to tau+ vt'
      if (vdecaymodeWp.eq.-101) write(*,*) '         to u dbar'
      if (vdecaymodeWp.eq.-103) write(*,*) '         to c sbar'
      if (vdecaymodeWm.eq.-107) write(*,*) '        to u dbar or c sbar'
      write(*,*)'            and'
      if (vdecaymodeWm.eq.11) write(*,*) '         to e- ve~  '
      if (vdecaymodeWm.eq.13) write(*,*) '         to mu- vm~ '
      if (vdecaymodeWm.eq.15) write(*,*) '         to tau- vt~'
      if (vdecaymodeWm.eq.101) write(*,*) '         to ubar d'
      if (vdecaymodeWm.eq.103) write(*,*) '         to cbar s'
      if (vdecaymodeWm.eq.107) write(*,*) '         to ubar d or cbar s'

      

c     index of the first coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=7
c     Z decay products
      i4=vdecaymodeWp
c     this is a convoluted way of doing things, but it ends up pairing 
c     lepton-antineutrino, and antilepton-neutrino
      i3=-sign(1,i4)*(abs(i4)+1)
      i5=vdecaymodeWm
      i6=-sign(1,i5)*(abs(i5)+1)

*********************************************************************
***********            BORN SUBPROCESSES              ***************
****  At the moment, CKM matrix is implicitly the unit matrix
****  This should be changed. 
*********************************************************************
      flst_nborn=0

      do i1=-nf,nf
         do i2=-nf,nf
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
      do i1=-nf,nf
         do i2=-nf,nf
            do i7=-nf,nf
               condition=.false.
               if(.not.(i1.eq.0.and.i2.eq.0)) then
c     exclude gg
                  if((i1.ne.0).and.(i1+i2.eq.0).and.(i7.eq.0)) then
c     q qbar -> g
                     condition=.true.
                  elseif((i1.eq.0).and.(i2.eq.i7)) then
c     g q
                     condition=.true.
                  elseif((i2.eq.0).and.(i1.eq.i7)) then
c     q g
                     condition=.true.
                  endif
               endif
               if(condition) then
                  flst_nreal=flst_nreal+1
                  if(flst_nreal.gt.maxprocreal) goto 998
                  do k=1,nlegreal
                     flst_real(k,flst_nreal)=ii(k)
                  enddo
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
      return
 998  write(*,*) 'init_processes: increase maxprocreal'
      call exit(-1)
 999  write(*,*) 'init_processes: increase maxprocborn'
      call exit(-1)
      end
 

      block data lepmass_data
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      data lepmass /0.51099891d-3,0.1056583668d0,1.77684d0/
      end
