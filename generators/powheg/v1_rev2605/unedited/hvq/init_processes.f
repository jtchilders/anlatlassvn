      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_par.h'
      include 'LesHouches.h'
      include 'pwhg_physpar.h'
      logical debug
      parameter (debug=.false.)
      integer i,ihvq,j
      character *3 flav(-5:5)
      data (flav(i),i=-5,5) 
     #     /'b~','c~','s~','u~','d~','g','d','u','s','c','b'/
      real * 8 powheginput,qmass
      real * 8 cmass, bmass
c     lepton masses
      real * 8 lepmass(3)
      data lepmass /0.51099891d-3,0.1056583668d0,1.77684d0/

      par_isrtinycsi = 1d-8
      par_isrtinyy = 1d-8
      par_fsrtinycsi = 1d-8
      par_fsrtinyy = 1d-8
c flag to do importance sampling in x variable in collinear remnants
      flg_collremnsamp=.true.      
c     number of light flavors
      qmass=powheginput("qmass")
      if(qmass.lt.3) then         
         st_nlight = 3
      elseif(qmass.lt.10) then
         st_nlight = 4
      else
         st_nlight = 5
      endif
      kn_masses(1)=0
      kn_masses(2)=0
      kn_masses(3)=qmass
      kn_masses(4)=qmass
      kn_masses(5)=0
      kn_minmass=2*qmass

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

      ihvq=st_nlight+1
      lprup(1)=1000+ihvq
c     index of the first light coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=5

      flst_nborn=1
      flst_nreal=1

      flst_born(1,flst_nborn)=0
      flst_born(2,flst_nborn)=0
      flst_born(3,flst_nborn)=ihvq
      flst_born(4,flst_nborn)=-ihvq

      flst_real(1,flst_nreal)=0
      flst_real(2,flst_nreal)=0
      flst_real(3,flst_nreal)=ihvq
      flst_real(4,flst_nreal)=-ihvq
      flst_real(5,flst_nreal)=0

      do i=-st_nlight,st_nlight
         if(i.ne.0) then
            flst_nborn=flst_nborn+1
            flst_born(1,flst_nborn)=i
            flst_born(2,flst_nborn)=-i
            flst_born(3,flst_nborn)=ihvq
            flst_born(4,flst_nborn)=-ihvq
            
            flst_nreal=flst_nreal+1
            flst_real(1,flst_nreal)=i
            flst_real(2,flst_nreal)=-i
            flst_real(3,flst_nreal)=ihvq
            flst_real(4,flst_nreal)=-ihvq
            flst_real(5,flst_nreal)=0
            
            flst_nreal=flst_nreal+1
            flst_real(1,flst_nreal)=i
            flst_real(2,flst_nreal)=0
            flst_real(3,flst_nreal)=ihvq
            flst_real(4,flst_nreal)=-ihvq
            flst_real(5,flst_nreal)=i
            
            flst_nreal=flst_nreal+1
            flst_real(1,flst_nreal)=0
            flst_real(2,flst_nreal)=i
            flst_real(3,flst_nreal)=ihvq
            flst_real(4,flst_nreal)=-ihvq
            flst_real(5,flst_nreal)=i
         endif
      enddo
      call init_top_dec(ihvq)
      end


      subroutine init_top_dec(nhvq)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'      
      include 'LesHouches.h'      
      integer nhvq
      integer itdec
      integer iwp1,iwp2,iwm1,iwm2
      real * 8 mdecwp1,mdecwp2,mdecwm1,mdecwm2,totbr
      real * 8 powheginput
      external powheginput
      if(nhvq.eq.6) then
         itdec=powheginput('#topdecaymode')
         if(itdec.eq.-1000000) itdec=0
         if(itdec.ne.0) then
            lprup(1)=300000+itdec
c first call to pickwdecay, to initialize and get back the branching fraction
            call  pickwdecays(iwp1,mdecwp1,iwp2,mdecwp2,
     #                 iwm1,mdecwm1,iwm2,mdecwm2,totbr)
            rad_branching=totbr
         endif
      endif
      end




