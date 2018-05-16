      subroutine init_phys
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_st.h'
      include 'pwhg_rad.h'
      include 'pwhg_dbg.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_physpar.h'
      character * 5 scheme
      character * 3 whichpdfpk
      real * 8 powheginput
      integer iorder,iret,iun,j
      external whichpdfpk,powheginput
c Initialization of default values for common block
c variables. These may be overridden by the user program
c init_processes.

      par_diexp=powheginput("#par_diexp")
      par_dijexp=powheginput("#par_dijexp")
      par_2gsupp=powheginput("#par_2gsupp")
      if(par_diexp.lt.0) par_diexp=1
      if(par_dijexp.lt.0) par_dijexp=1
      if(par_2gsupp.lt.0) par_2gsupp=1
      if(par_diexp.ne.par_dijexp) then
         write(*,*) 'par_dijexp not equal to par_diexp;'
         write(*,*) 'not possible at present!!! Program exits'
         call exit(-1)
      endif
c
      par_isrtinycsi = 1d-6
      par_isrtinyy = 1d-6
      par_fsrtinycsi = 1d-5
      par_fsrtinyy = 1d-6
c
      rad_branching=1
c this is set to true in processes where the FSR jacobian
c can become singular (massless recoil particle)
      flg_jacsing=.false.
c flag to use importance sampling in the x variable in
c collinear remnant generation. Needed for charm at LHC
      flg_collremnsamp=.false.
c End initialization of common block defaults.
      pdf_ih1=powheginput('ih1')
      pdf_ih2=powheginput('ih2')
      if(whichpdfpk().eq.'lha') then
         pdf_ndns1=powheginput('lhans1')
         pdf_ndns2=powheginput('lhans2')
      elseif(whichpdfpk().eq.'mlm') then
         pdf_ndns1=powheginput('ndns1')
         pdf_ndns2=powheginput('ndns2')
      else
         write(*,*) ' unimplemented pdf package',whichpdfpk()
         stop
      endif
      if(pdf_ndns1.ne.pdf_ndns2) then
         st_lambda5MSB=powheginput('QCDLambda5')
      else
         call genericpdfpar(pdf_ndns1,pdf_ih1,st_lambda5MSB,
     1                      scheme,iorder,iret)
         if(iret.ne.0) then
            write(*,*) ' faulty pdf number ',pdf_ndns1
            stop
         endif
      endif
      kn_beams(0,1)=powheginput('ebeam1')
      kn_beams(0,2)=powheginput('ebeam2')
      kn_beams(1,1)=0
      kn_beams(1,2)=0
      kn_beams(2,1)=0
      kn_beams(2,2)=0
      kn_beams(3,1)=kn_beams(0,1)
      kn_beams(3,2)=-kn_beams(0,2)
      kn_sbeams=4*kn_beams(0,1)*kn_beams(0,2)

c generation cut: see Gen_born_phsp.f
      kn_ktmin=powheginput("#bornktmin")
      if(kn_ktmin.lt.0) kn_ktmin=0

c masses for light fermions, used in momentum reshuffling
      do j=1,6
         physpar_mq(j)=0
      enddo
      do j=1,3
         physpar_ml(j)=0
      enddo

c thresholds 
      rad_ptsqmin=powheginput('#ptsqmin')
      if(rad_ptsqmin.lt.0) rad_ptsqmin=0.8d0
      rad_charmthr2=powheginput('#charmthr')
      if(rad_charmthr2.lt.0) rad_charmthr2=1.5d0
      rad_charmthr2=rad_charmthr2**2
      rad_bottomthr2=powheginput('#bottomthr')
      if(rad_bottomthr2.lt.0) rad_bottomthr2=5d0
      rad_bottomthr2=rad_bottomthr2**2
c scale factors
      st_renfact=powheginput('#renscfact')
      st_facfact=powheginput('#facscfact')
      if(st_facfact.lt.0) st_facfact=1
      if(st_renfact.lt.0) st_renfact=1

c     if true, perform the check that there are no coloured light
c     partons before flst_lightpart
      flg_lightpart_check=.true.
c
c initialize Lambda values for radiation
      call init_rad_lambda
c
      call init_processes

      call init_couplings

c initialize number of singular regions
      rad_nkinreg=1+(nlegborn-flst_lightpart+1)
      call genflavreglist


      dbg_softtest=.true.
      dbg_colltest=.true.
      if(powheginput("#softtest").eq.0) dbg_softtest=.false.
      if(powheginput("#colltest").eq.0) dbg_colltest=.false.
      if(flg_withdamp) then
         write(*,*) ' POWHEG: no soft tests if withdamp is set'
         dbg_softtest=.false.
      endif
      if(flg_bornonly) then
         write(*,*)
     $        ' POWHEG: no soft and coll. tests if bornonly is set'
         dbg_softtest=.false.
         dbg_colltest=.false.
      endif
      if (dbg_softtest.or.dbg_colltest) then         
         call newunit(iun)
         open(unit=iun,file='pwhg_checklimits')
         call checklims(iun)
         call flush(iun)
         write(*,*) ' POWHEG:  '
         write(*,*) ' Check of soft/collinear limits performed'
         write(*,*) ' Results in file pwhg_checklimits'
      endif
      end

