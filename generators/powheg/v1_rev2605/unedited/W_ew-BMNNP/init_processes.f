      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_par.h'
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      include 'pwhg_physpar.h'
      include 'pwhg_st.h'
      integer i1,i2,i3,i4,i5,k,ii(nlegreal)
      equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),
     #  (i4,ii(4)),(i5,ii(5))
      logical debug
      parameter (debug=.true.)
      integer j
      integer charge3(-6:6)
      data charge3 /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      logical condition
      real * 8 powheginput
      external powheginput
c     vector boson id and decay
      integer idvecbos,vdecaymode,tmp
      common/cvecbos/idvecbos,vdecaymode
c     lepton masses
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      real *8 kt2minqed
      common/showerqed/kt2minqed
      real * 8 cmass, bmass
c******************************************************
c     Choose the process to be implemented
c******************************************************
c    ID of vector boson produced
      idvecbos=powheginput('idvecbos')
c   decay products of the vector boson
      vdecaymode=powheginput('vdecaymode')
 
      par_isrtinycsi = 1d-8
      par_isrtinyy = 1d-8

      par_fsrtinycsi = 1d-11
      par_fsrtinyy = 1d-11

      kt2minqed = powheginput("#kt2minqed")
      if (kt2minqed.le.0d0) kt2minqed  = 0.001d0**2
c      flg_jacsing = .true.

      if (lepmass(1).ne.0.51099892d-3) then
         write(*,*) 'block data lepmass not loaded. stop running' 
         stop
      endif
      
      if(idvecbos.eq.24) then
         write(*,*) 
         write(*,*) ' POWHEG: Single W+ production and decay ' 
         if (vdecaymode.eq.-11) then
            write(*,*) '         to e+ ve '
         elseif (vdecaymode.eq.-13) then
            write(*,*) '         to mu+ vmu'
         elseif (vdecaymode.eq.-15) then
            write(*,*) '         to tau+ vtau'
         else
            write(*,*) 'ERROR: The decay mode you selected' /
     $           /' is not allowed '
            call exit(-1)
         endif
      elseif(idvecbos.eq.-24) then
         write(*,*) 
         write(*,*) ' POWHEG: Single W- production and decay ' 
         if (vdecaymode.eq.11) then
            write(*,*) '         to e- ve~ '
         elseif (vdecaymode.eq.13) then
            write(*,*) '         to mu- vmu~'
         elseif (vdecaymode.eq.15) then
            write(*,*) '         to tau- vtau~'
         else
            write(*,*) 'ERROR: The decay mode you selected' /
     $           /' is not allowed '
            call exit(-1)
         endif
      else
         write(*,*) 'ERROR: The ID of vector boson you selected' 
     $        //' is not allowed (24: W+ -24: W-)'
         stop
      endif

c     change the LHUPI id of the process according to vector boson id
c     and decay
      lprup(1)=10000+vdecaymode ! 10000+idup of charged decay product of the W
      
      if(lprup(1).eq.10011) then
         decmass=lepmass(1)
         
      elseif(lprup(1).eq.(10000-11)) then
         decmass=lepmass(1)
        
      elseif(lprup(1).eq.10013) then
         decmass=lepmass(2)
         
      elseif(lprup(1).eq.(10000-13)) then
         decmass=lepmass(2)

      elseif(lprup(1).eq.10015) then
         decmass=lepmass(3)
         
      elseif(lprup(1).eq.(10000-15)) then
         decmass=lepmass(3) 
  
      else
c     not yet implemented
         write(*,*) 'non leptonic W decays '//
     #        'not yet implemented'
         stop
      endif   

c     Set here lepton and quark masses for momentum reshuffle in the LHE event file
      do j=1,st_nlight         
         physpar_mq(j)=0d0
      enddo
      do j=1,3
         physpar_ml(j)=lepmass(j)
      enddo
c     read eventual c and b masses from the input file
      cmass=powheginput("#cmass_lhe")
      if (cmass.gt.0d0) physpar_mq(4)=cmass
      bmass=powheginput("#bmass_lhe")
      if (bmass.gt.0d0) physpar_mq(5)=bmass



c*********************************************************     
c
c     index of the first light particle in the final state
c     that can give rise to collinear singularities;
c     The charged leptons are considered MASSIVE particles in this code,
c     and thus do not give rise to collinear singularities
      if(powheginput('#easlight').eq.1d0) then
         flst_lightpart=3
      else
         flst_lightpart=5
      endif
      i3=vdecaymode
      if ((idvecbos.eq.24).and.(vdecaymode.lt.0)) then
         i4=-vdecaymode+1
      elseif ((idvecbos.eq.-24).and.(vdecaymode.gt.0)) then
         i4=-(vdecaymode+1)
      endif
c     Born graphs
      flst_nborn=0
      condition=.false.
      do i1=-5,5
         do i2=-5,5
            condition=(charge3(i1)+charge3(i2)).eq.(sign(3,idvecbos))
            if(condition) then
c     q qbar'
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
     
c     Real graphs    
      flst_nreal=0
      condition=.false.
      do i1=-5,5
         do i2=-5,5
            if (abs(i1).eq.abs(i2)) goto 11
            do i5=-5,5
               condition=.false.
               if ((i1.eq.0).and.(i2.ne.0)) then
                  condition=(charge3(i2)-charge3(i5))
     $                 .eq.(sign(3,idvecbos))    
               endif
               if ((i2.eq.0).and.(i1.ne.0)) then
                  condition=(charge3(i1)-charge3(i5))
     $                 .eq.(sign(3,idvecbos))
               endif
               if (i5.eq.0) then
                condition=(charge3(i1)+charge3(i2))
     $                 .eq.(sign(3,idvecbos))
               endif   
               if(condition) then
                  flst_nreal=flst_nreal+1
                  if(flst_nreal.gt.maxprocreal) goto 998
                  do k=1,nlegreal
                     flst_real(k,flst_nreal)=ii(k)
                  enddo
               endif
            enddo
            condition=(charge3(i1)+charge3(i2))
     $           .eq.(sign(3,idvecbos))
            if(condition) then
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
      
      block data lepmass_data 
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      data lepmass /0.51099892d-3,0.105658369d0,1.77699d0/
      end
