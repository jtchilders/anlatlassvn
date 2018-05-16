c  The next subroutines opens some histograms and prepares them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist_slp_fat
      implicit none
      include  'LesHouches.h'
      include 'pwhg_book.h'
      integer diag
      integer max_diag
      parameter (max_diag=300)
      real * 8 binsize(max_diag)
      common/pwhghistcommon/binsize
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG

c ========================
            
      call pwhginihist
 
C -----------------------------------------------------
C HISTOGRAMS WITH VBF CUTS in FAT-JET SETUP:

      diag = 1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(all VBF cuts)','LIN',
     &                binsize(diag),0d0,1d0)
      
c     cross section sanity check
      diag = 2
      binsize(diag) = 5d0
      call pwhgbookup(diag,'3 jet incl.- VBF CUTS','LIN',
     .     binsize(diag),2.5d0,72.5d0)

      diag = 3
      binsize(diag) = 5d0
      call pwhgbookup(diag,'2 jet incl.- VBF CUTS','LIN',
     .     binsize(diag),2.5d0,72.5d0)

      diag = 4
      binsize(diag) = 5d0
      call pwhgbookup(diag,'2 jet excl.- VBF CUTS','LIN',
     .     binsize(diag),2.5d0,72.5d0)

      diag = 5
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT jet 1- VBF CUTS','LOG',
     .     binsize(diag),0d0,880d0)

      diag = 6
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT jet 2- VBF CUTS','LOG',
     .     binsize(diag),0d0,520d0)

      diag = 7
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y jet 1- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 8
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y jet 2- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 9
      binsize(diag) = 20d0
      call pwhgbookup(diag,'M_reconstructed(2l2j)- VBF CUTS','LOG',
     .     binsize(diag),0d0,2000d0)

      diag = 10
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y j1j2(tag)- VBF CUTS','LIN',
     .     binsize(diag),-6d0,6d0)

      diag = 11
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT fat jet - VBF CUTS','LOG',
     .     binsize(diag),0d0,520d0)

      diag = 12
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT LEPT- VBF CUTS','LOG',
     .     binsize(diag),0d0,2000d0)

      diag = 13
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT MISS- VBF CUTS','LOG',
     .     binsize(diag),0d0,1000d0)

      diag = 14
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Eta_lept- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 15
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Eta_fat_jet- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 16
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Pt J3- VBF CUTS','LOG',
     .     binsize(diag),0d0,200d0)

      diag = 17
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y J3- VBF CUTS, kt3>pt_tag_min','LOG',
     .     binsize(diag),-5d0,5d0)

      diag = 18
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y J3- VBF CUTS, kt3>pt_cut','LOG',
     .     binsize(diag),-5d0,5d0)

      diag = 19
      binsize(diag) = 50d0
      call pwhgbookup(diag,'M_j1j2(tag) - VBF CUTS','LIN',
     .     binsize(diag),0d0,3000d0)

      diag = 20
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Phi(j1j2)- VBF CUTS','LIN',
     .     binsize(diag),0d0,3.2d0)

      diag = 21
      binsize(diag) = 2d0
      call pwhgbookup(diag,'M(fat_jet)- VBF CUTS','Lin',
     .     binsize(diag),40d0,120d0)
  
      diag = 22
      binsize(diag) = 2d0
      call pwhgbookup(diag,'m(34)- VBF CUTS','LIN',
     .     binsize(diag),40d0,120d0)

      diag = 23
      binsize(diag) = 2d0
      call pwhgbookup(diag,'m(56)- VBF CUTS','LIN',
     .     binsize(diag),60d0,100d0)
 
      diag = 24
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y*- VBF CUTS, kt3>pt_jet_min','LIN',
     .     binsize(diag),-4d0,4d0)
      
      diag = 25
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y*- VBF CUTS, kt3>pt_cut','LIN',
     .     binsize(diag),-4d0,4d0)
      
      end 
      
c ========================

      subroutine analysis_slp_fat(dsig0)
      implicit none
      real * 8 dsig0,dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'nlegborn.h'
      include 'LesHouches.h'
      include 'cvecbos.h'
      include 'PhysPars.h'
      include 'PhysPars_Higgs.h'

      logical ini
      data ini/.true./
      save ini
      integer diag
      integer max_diag
      parameter (max_diag=300)
      real * 8 binsize(max_diag)
      common/pwhghistcommon/binsize
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
c      data WHCPRG/'NLO   '/
      integer nleptons,nnu
      integer ileptons(200),inu(200)
      integer   maxjet
      parameter (maxjet=2048)
      real * 8  kt(maxjet),eta(maxjet),rap(maxjet),
     1    phi(maxjet),pj(4,maxjet),ptrel(maxjet)
      real * 8 ptel1,ptel2,etael1,etael2,ptmiss
      real * 8 pj12(4),y12, pel12(4),etael12
      real*8 ptl(100),etal(100)
      real * 8 mll,phill
      real * 8 invmass 
      real * 8 mjj,phijj
      real * 8 rap_max,rap_min
      real * 8 r,fphi
      integer ihep,j,mjets
      real * 8 etafromp,ptfromp 
      logical passcuts_vbf
     
c cut parameters:
      real*8 ptmin,pt_cut
      real*8 ptmiss_min,etal_max,ptl_min
      real*8 pt_tag_min,pt_dec_min,pt_dec_cut,yj_max
      real*8 yjj_min,Rjj_min,mjj_tag_min
      real*8 Rjl_min,Rll_min
      logical rap_gap,rap_sign

      common /jetcuts/ptmin,yj_max,Rjj_min

      real*8 Rjl_tmp0,Rjl_tmp1
      integer ij,il
      real*8 Rll_tmp
      real*8 Rjj_tmp,Rjj_tmp0

      integer itag1,itag2,itag3
      real*8 pt_1,pt_2,pt_3
      integer ltag1,ltag2
      real*8 ptl_1,ptl_2

      real*8 pmiss_sum1,pmiss_sum2,pl_sum(1:4)
      real*8 m34,m56,p34(0:3),p56(0:3)
      real*8 ystar

      integer iveto1,iveto2,nveto
      real*8 pt_v1,pt_v2
      real*8 pt_veto
      logical rap_veto

      real*8 ptl1_min,ptl2_min,mll_min,mll_max,phill_max

      real * 8 powheginput
      external powheginput

      real*8 mj
      real*8 mj_tmp(maxjet),mj_diff(maxjet)
      real*8 pjj(4,maxjet,maxjet)
      integer i,jdec

      real*8   kt_dec(2),eta_dec(2),rap_dec(2),
     1         phi_dec(2),pj_dec(4,2),ptrel_dec(2)
      real*8   kt_fat,eta_fat,rap_fat,
     1         phi_fat,pj_fat(4),ptrel_fat

      real*8   kt_j(maxjet),eta_j(maxjet),rap_j(maxjet),
     1         phi_j(maxjet),pj_j(4,maxjet),ptrel_j(maxjet)
      integer ii,mjets_j
      
      real*8 eta_dec_max,Rjd_min,Rjd_tmp0,Rjd_tmp1
      real*8 mdec_jj
      real*8 rdec(1:4),rmdec
      complex*16 sol1,sol2
      real*8 arg
      real*8 plep(1:4),pneu(1:4),rneu(1:4)
      real*8 rap_veto_eta
      real*8 pw1(1:4),pw2(1:4)      
      real*8 pt_boost_min

      real*8 pwjet(1:4)
      logical passjetcut
      integer wtag

c =======================================================

c initialize all cuts:
C
      etal_max   = 1d10 
      eta_dec_max   = 1d10
      ptl_min    = 0d0  
      ptmiss_min = 0d0

      pt_tag_min = 0d0
      pt_dec_min = 0d0
      ptmin = 0d0
      pt_cut = 0d0
      yj_max  = 1d10

      pt_boost_min = 0d0
      
      yjj_min = 0d0
      rap_sign = .false.
      rap_gap  = .false.

      Rjj_min = 0d0
      mjj_tag_min = 0d0

      Rjl_min = 0d0
      Rjd_min = 0d0

      pt_veto = 1d10
      rap_veto = .false.
      rap_veto_eta = 0d0

      passjetcut = .false.
      wtag = 0
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c VBF cuts (can be modified by user):
C
        Rjj_min = 1d0      !parameter for jet algorithm
        ptmin = 1d0        !parameter for jet algorithm

        pt_boost_min = 300d0

        ptl_min    = 320d0 
        ptmiss_min = 30d0

        pt_tag_min = 25d0  !tag jets
        pt_cut = 10d0      !parameter for y3 distributions
        yj_max  = 4.5d0

        yjj_min  = 3.0d0
        rap_sign = .true.
        rap_gap  = .true.

        mjj_tag_min = 600d0

        Rjl_min = 0.3d0
        Rjd_min = 0.3d0

        passcuts_vbf = .true.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      
c from pico to femto      
      dsig=dsig0*1000
      diag = 0 
      if (ini) then

         ph_Hmass   = powheginput('hmass')
         ph_Hwidth   = powheginput('hwidth')

         ph_Wmass   = powheginput('wmass')
         ph_Wwidth   = powheginput('wwidth')

         write (*,*)
         write (*,*) '********************************************'
         if(WHCPRG.eq.'NLO   ') then
            write (*,*) 'Fixed order ANALYSIS CALLED '
         
            if(((abs(vdecaymodeWp).eq.101.and.abs(vdecaymodeWm).eq.11)
     $           .or. 
     $        (abs(vdecaymodeWp).eq.103.and.abs(vdecaymodeWm).eq.11).or.
     $        (abs(vdecaymodeWp).eq.107.and.abs(vdecaymodeWm).eq.11).or.  
     $        (abs(vdecaymodeWp).eq.101.and.abs(vdecaymodeWm).eq.13).or.
     $        (abs(vdecaymodeWp).eq.103.and.abs(vdecaymodeWm).eq.13).or.
     $        (abs(vdecaymodeWp).eq.107.and.abs(vdecaymodeWm).eq.13)
     $        ).and.vdecaymodeWp*vdecaymodeWm.lt.0) then

              continue
            else
             write(*,*) '**************************************'
             write(*,*) ' template analysis works only for W-  '
             write(*,*) ' bosons decaying to electrons or muons'
             write(*,*) ' and W+ decaying hadronically         '
             write(*,*) '                 STOP                 '
             write(*,*) '**************************************'
             call exit(1)
            endif
      
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS CALLED     '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS CALLED     '
         endif  

c
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) '                ANALYSIS CUTS               '
         write(*,*) '*********************************************'
         write(*,*) '*********************************************'
            write(*,*) ''
            write(*,*) 'Higgs parameters:'
            write(*,*) 'Higgs mass [GeV]= ',ph_Hmass
            write(*,*) 'Higgs width [GeV]= ',ph_Hwidth
            write(*,*) 'W mass [GeV]= ',ph_Wmass
            write(*,*) 'W width [GeV]= ',ph_Wwidth
            write(*,*) ''
            write(*,*) 'jet cuts:'
            write(*,*) 'jet ptmin = ',ptmin
            write(*,*) 'jet3 pt_cut = ',pt_cut
            write(*,*) 'tag jets pt_tag_min = ',pt_tag_min
            write(*,*) 'yj_max  = ',yj_max
            write(*,*) 'yjj_min = ',yjj_min
            write(*,*) 'Rjj_min = ',Rjj_min
            write(*,*) 'mjj_tag_min = ',mjj_tag_min
            write(*,*) ''
            write(*,*) 'boosted jet: pt >',pt_boost_min
            write(*,*) ''
            write(*,*) 'lepton cuts:'
            write(*,*) 'ptl_min    = ',ptl_min 
            write(*,*) 'etal_max   = ',etal_max 
            write(*,*) ''
            write(*,*) 'ptmiss_min    = ',ptmiss_min
            write(*,*) ''
            write(*,*) 'extra cuts:'
            write(*,*) 'Rjl_min  = ',Rjl_min
            write(*,*) 'rap_gap  = ',rap_gap
            write(*,*) 'rap_sign  = ',rap_sign
            write(*,*) ''
            write(*,*) 'fat jet reconstructed around MW'
         write(*,*) '********************************************'
         write(*,*) '********************************************'        
         ini=.false.
      endif
         
c find electrons and muons and neutrinos
      nleptons=0
      nnu=0
      do ihep=1,nhep
         if(isthep(ihep).eq.1) then
            if((abs(idhep(ihep)).eq.11).or.
     %         (abs(idhep(ihep)).eq.13)) then
               nleptons=nleptons+1
               ileptons(nleptons)=ihep
            elseif((abs(idhep(ihep)).eq.12).or.
     %             (abs(idhep(ihep)).eq.14)) then
               nnu=nnu+1
               inu(nnu)=ihep
            endif
         endif
      enddo
            
      if(nleptons.gt.200.or.nnu.gt.200) then
         write(*,*) ' crazy event, too many leptons'
         return
      endif
      if(nleptons.lt.1.or.nnu.lt.1) then
         write(*,*) ' crazy event, missing leptons'
         return
      endif
c sort by pt
      call sortbypt(nleptons,ileptons)
      call sortbypt(nnu,inu)

c initialize:
	do j = 1,maxjet
	   kt(j) = 0d0
	   kt_j(j) = 0d0
	enddo   
        kt_fat = 0d0
        kt_dec(:) = 0d0

      mjets_j=10
      call buildjets_slp(mjets_j,kt_j,eta_j,rap_j,phi_j,
     &     pj_j,ptrel_j,wtag,pwjet,passjetcut)

      if (.not.passjetcut) return

      ii = 0
      do i = 1,mjets_j
c select fat jet from array of all jets:
         if (i.eq.wtag) then
            kt_fat    = kt_j(wtag)
            eta_fat   = eta_j(wtag)
            rap_fat   = rap_j(wtag)
            phi_fat   = phi_j(wtag)
            ptrel_fat = ptrel_j(wtag)

            pj_fat(1:4) = pj_j(1:4,wtag)
c new array with remaining jets:
         else   
            ii = ii+1
            kt(ii)    = kt_j(i)
            eta(ii)   = eta_j(i)
            rap(ii)   = rap_j(i)
            phi(ii)   = phi_j(i)
            ptrel(ii) = ptrel_j(i)

            pj(:,ii) = pj_j(:,i)
         endif
      enddo   
      mjets = ii

c invariant mass of the fat jet:
      mj = invmass(pj_fat)
	
cccccccccccccccccccccccccccccccccccccccc
c
c now start analysis:

c require one very hard jet:
c      if (kt_wjet.lt.pt_boost_min) return 
      if (kt_fat.lt.pt_boost_min) return 

c set to zero for debugging:
      pj_dec(:,:)=0d0

      if (kt(1).lt.pt_tag_min) return
      if (kt(2).lt.pt_tag_min) return

c************************************

	do ihep = 1,nhep
	if (idhep(ihep).eq.24) then
	 pw1(1:4) = phep(1:4,ihep)
	endif 
	if (idhep(ihep).eq.-24) then
	 pw2(1:4) = phep(1:4,ihep)	 
	endif 
	enddo   

      pmiss_sum1 = 0d0
      pmiss_sum2 = 0d0
      do il = 1,nnu ! neutrinos
         pmiss_sum1 = pmiss_sum1 + phep(1,inu(il))
         pmiss_sum2 = pmiss_sum2 + phep(2,inu(il))
      enddo   
      ptmiss = sqrt((pmiss_sum1)**2+
     1              (pmiss_sum2)**2)
	
c reconstructed neutrino momentum:

      plep(1:4) =   phep(1:4,ileptons(1))
      pneu(1) =   pmiss_sum1
      pneu(2) =   pmiss_sum2

      arg = (plep(4)**2*(ph_wmass**4 + 4d0*ph_wmass**2*(
     &   plep(1)*pneu(1) + plep(2)*pneu(2)) + 
     &   4d0*(plep(1)**2*pneu(1)**2 + plep(3)**2*pneu(1)**2 + 
     &   2d0*plep(1)*plep(2)*pneu(1)*pneu(2) + plep(2)**2*pneu(2)**2 + 
     &   plep(3)**2*pneu(2)**2-plep(4)**2*(pneu(1)**2+pneu(2)**2))))

      if (arg.lt.0) return

      sol1 = (plep(3)*ph_wmass**2 + 2d0*plep(1)*plep(3)*pneu(1) + 
     &     2d0*plep(2)*plep(3)*pneu(2) - 
     &     dsqrt(arg))/(2.d0*(plep(4)**2-plep(3)**2))


      sol2 =  (plep(3)*ph_wmass**2 + 2d0*plep(1)*plep(3)*pneu(1) + 
     &     2d0*plep(2)*plep(3)*pneu(2) + 
     &     dsqrt(arg))/(2.d0*(plep(4)**2 - plep(3)**2))

      if (abs(sol1).lt.abs(sol2)) then
         rneu(4) = dsqrt(ptmiss**2+dreal(sol1)**2)
         rneu(3) = dreal(sol1)
      else
         rneu(4) = dsqrt(ptmiss**2+dreal(sol2)**2)
         rneu(3) = dreal(sol2)
      endif
      rneu(1:2) = pneu(1:2)  

c compare to reconstructed quantities:
      rdec(1:4) = plep(1:4)+rneu(1:4)+ 
     &            pj_fat(1:4)
      rmdec = invmass(rdec)
       
      etael1=etafromp(phep(1,ileptons(1)))
      ptel1=ptfromp(phep(1,ileptons(1)))
c
cccccccccccccccccccccccccccccccccc
c
c check, if event passes VBF cuts:
c
c identify 2 hardest jets ("tag" jets) and 3rd hardest jet, 
c and check if they pass pt and rap cuts:
      pt_1 = 0d0
      pt_2 = 0d0
      pt_3 = 0d0
      itag1 = 0
      itag2 = 0
      itag3 = 0

      do j = 1, mjets
           if ( kt(j).gt.pt_1 .and. kt(j).gt.pt_tag_min 
     &          .and. abs(rap(j)).lt.yj_max) then
              itag1 = j   
              pt_1  = kt(j)
           end if
      end do
      do j = 1, mjets
           if ( kt(j).gt.pt_2 .and. kt(j).gt.pt_tag_min 
     &          .and. abs(rap(j)).lt.yj_max
     &          .and. j.ne.itag1 
     &         ) then
              itag2 = j    
              pt_2  = kt(j)
           end if         
      end do
c
      do j = 1, mjets
           if ( kt(j).gt.pt_3 .and. kt(j).gt.ptmin 
     &          .and. abs(rap(j)).lt.yj_max
     &          .and. j.ne.itag1 
     &          .and. j.ne.itag2
     &         ) then
              itag3= j    
              pt_3 = kt(j)
           end if         
      end do      
      if (itag1.eq.0.or.itag2.eq.0) passcuts_vbf = .false. 
      if (itag1.eq.0.or.itag2.eq.0) return

c identify hardest lepton: 
      ptl_1 = 0d0
      ltag1 = 0
      do j = 1,nleptons
         ptl(j) = ptfromp(phep(1,ileptons(j)))
         etal(j)= etafromp(phep(1,ileptons(j)))
      enddo
      do j = 1,nleptons 
           if ( ptl(j).gt.ptl_1 .and. ptl(j).gt.ptl_min 
     &          .and. abs(etal(j)).lt.etal_max) then
              ltag1 = j   
              ptl_1 = ptl(j)
           end if
      end do
      if (ltag1.eq.0) passcuts_vbf = .false.
      if (ltag1.eq.0) return 

      pmiss_sum1 = 0d0
      pmiss_sum2 = 0d0
      do il = 1,nnu !neutrinos
         pmiss_sum1 = pmiss_sum1 + phep(1,inu(il))
         pmiss_sum2 = pmiss_sum2 + phep(2,inu(il))
      enddo  
      ptmiss = sqrt((pmiss_sum1)**2+
     1              (pmiss_sum2)**2)

c
cccccccccccccccccccccccccccccccccccccc
c
C     now come histograms with VBF cuts 
c
      if(mjets.lt.2) passcuts_vbf = .false.      
      if (itag1.eq.0.or.itag2.eq.0) passcuts_vbf = .false. 
      if (ltag1.eq.0) passcuts_vbf = .false.

      if (kt(itag1).lt.pt_tag_min) passcuts_vbf = .false.   
      if (kt(itag2).lt.pt_tag_min) passcuts_vbf = .false. 

      if(abs(rap(itag1)).gt.yj_max) passcuts_vbf = .false. 
      if(abs(rap(itag2)).gt.yj_max) passcuts_vbf = .false.
      if(abs(rap(itag1)-rap(itag2)).lt.yjj_min) passcuts_vbf = .false.
      if (rap_sign.and.(rap(itag1)*rap(itag2).ge.0d0)) 
     &     passcuts_vbf = .false.

      if(abs(rap_fat).gt.yj_max) passcuts_vbf = .false. 
       
      Rjj_tmp = r(pj(1:4,itag1),pj(1:4,itag2))
      if (Rjj_tmp.lt.Rjj_min) passcuts_vbf = .false. 
     
      if(ptmiss.lt.ptmiss_min) passcuts_vbf = .false.  
     
      etael1=etal(ltag1)
      if(abs(etael1).gt.etal_max) passcuts_vbf = .false.  

      ptel1=ptl(ltag1)
      if(ptel1.lt.ptl1_min)  passcuts_vbf = .false.  

c separation of tagging jets from leading charged lepton:
      Rjl_tmp1 = 1d10
      do ij = 1,mjets
         if(ij.eq.itag1.or.ij.eq.itag2) then
         do il = 1,nleptons
            if(il.eq.ltag1) then
               Rjl_tmp0 = r(pj(1:4,ij),phep(1:4,ileptons(il)))
               Rjl_tmp1 = min(Rjl_tmp0,Rjl_tmp1)
            endif   
         enddo !il
         endif
      enddo !ij
      if (Rjl_tmp1.lt.Rjl_min) passcuts_vbf = .false. 

c separation of tagging jets from fat jet:
      Rjd_tmp1 = 1d10
      do ij = 1,mjets
         if(ij.eq.itag1.or.ij.eq.itag2) then
               Rjd_tmp0 = r(pj(1:4,ij),pj_fat(1:4))
               Rjd_tmp1 = min(Rjd_tmp0,Rjd_tmp1)
         endif
      enddo !ij
      if (Rjd_tmp1.lt.Rjd_min) passcuts_vbf = .false. 

      pj12 = pj(:4,itag1)+pj(:4,itag2)
      mjj = invmass(pj12)
      if(mjj.lt.mjj_tag_min)  passcuts_vbf = .false. 

      if (rap_gap.and.(min(rap(itag1),rap(itag2))
     &             .ge.min(etael1,eta_fat)))
     & passcuts_vbf = .false. 
      if (rap_gap.and.(max(rap(itag1),rap(itag2))
     &             .le.max(etael1,eta_fat)))
     & passcuts_vbf = .false. 
c
c*************************

      if (passcuts_vbf) then

      diag = 1   
      call pwhgfill(diag,0.5d0,dsig)

c     three jet inclusive.
      diag=2
      do j=1,3
         if(mjets.lt.j) kt(j)=0
      enddo
      do j=5,70,5
         if (itag3.ne.0) then 
         if(kt(itag1).gt.j.and.kt(itag2).gt.j.and.kt(itag3).gt.j) then
            call pwhgfill(diag,dble(j),dsig)
         endif
         endif 
      enddo

c     two jet inclusive and two jet exclusive.
c     It is exclusive if the third jet is below the cut
c     inclusive otherwise
      diag=3
      do j=1,3
         if(mjets.lt.j) kt(j)=0
      enddo
      do j=5,70,5
         if(kt(itag1).gt.dble(j).and.kt(itag2).gt.dble(j)) then
c inclusive:
            call pwhgfill(diag,dble(j),dsig)
c exclusive:
            if(itag3.eq.0) then 
               call pwhgfill(diag+1,dble(j),dsig)
            elseif (kt(itag3).lt.j) then 
               call pwhgfill(diag+1,dble(j),dsig)
            endif
         endif
      enddo
      diag=4

c     PT of jet 1
      diag=5
      call pwhgfill(diag,kt(itag1),dsig/binsize(diag))
c     PT of jet 2
      diag=6
      call pwhgfill(diag,kt(itag2),dsig/binsize(diag))
c     Y of jet 1
      diag=7
      call pwhgfill(diag,rap(itag1),dsig/binsize(diag))
c     Y of jet 2
      diag=8
      call pwhgfill(diag,rap(itag2),dsig/binsize(diag))

C     M_dec(reconstruced)
      diag=9
      call pwhgfill(diag,rmdec,dsig/binsize(diag))

c     Y jet1 - Y jet 2
      diag=10
      pj12 = pj(:4,itag1)+pj(:4,itag2)
      y12 = rap(itag1)-rap(itag2)
      call pwhgfill(diag,y12,dsig/binsize(diag))

C     Pt fat jet    
      diag=11
      call pwhgfill(diag,kt_fat,dsig/binsize(diag))

C     Pt lept     
      diag=12
      call pwhgfill(diag,ptel1,dsig/binsize(diag))

C     Pt miss
      diag=13
      call pwhgfill(diag,ptmiss,dsig/binsize(diag))

C     ETA lept
      diag=14
      call pwhgfill(diag,etael1,dsig/binsize(diag))

C     eta fat jet    
      diag=15
      call pwhgfill(diag,eta_fat,dsig/binsize(diag))

c     pt of third jet
      diag=16
      if(itag3.ne.0) call pwhgfill(diag,kt(itag3),dsig/binsize(diag))

c     y of third jet, kt3>pt_tag_min
      diag=17
      if(itag3.ne.0 .and. kt(itag3) > pt_tag_min ) 
     .     call pwhgfill(diag,rap(itag3),dsig/binsize(diag))
c     y of third jet, kt3>pt_cut
      diag=18
      if(itag3.ne.0 .and. kt(itag3) > pt_cut ) 
     .     call pwhgfill(diag,rap(itag3),dsig/binsize(diag))

C     M_j1j2
      diag=19
      mjj = invmass(pj12)
      if(mjets.ge.2) call pwhgfill(diag,mjj,dsig/binsize(diag))

c     Phi(j1j2) 
      diag=20 
      phijj=fphi(pj(1:4,itag1),pj(1:4,itag2))
      if(mjets.ge.2) call pwhgfill(diag,phijj,dsig/binsize(diag))
     
C     Minv(fat_jet) 
      diag=21
      call pwhgfill(diag,mj,dsig/binsize(diag))

C     M(34) 
      diag=22
      m34 = invmass(pw1)
      call pwhgfill(diag,m34,dsig/binsize(diag))

C     M(56) 
      diag=23
      m56 = invmass(pw2)
      call pwhgfill(diag,m56,dsig/binsize(diag))     

      if(itag3.ne.0 .and. kt(itag3) > pt_tag_min ) then 
          ystar     = rap(itag3)-(rap(itag1)+rap(itag2))/2d0
c         y* 
          diag=24
          call pwhgfill(diag,ystar,dsig/binsize(diag))
      endif

      if(itag3.ne.0 .and. kt(itag3) > pt_cut ) then 
          ystar     = rap(itag3)-(rap(itag1)+rap(itag2))/2d0
c         y* 
          diag=25
          call pwhgfill(diag,ystar,dsig/binsize(diag))
      endif

      endif !passcuts_vbf

      end

ccccccccccccccccccccccccccccccccccccc

      subroutine buildjets_slp(mjets,kt,eta,rap,phi,pjet,ptrel,
     &     wtag,wjet,passcuts)
c     arrays to reconstruct jets
      implicit none
      integer mjets
      real * 8  kt(mjets),eta(mjets),rap(mjets),phi(mjets),
     &          pjet(4,mjets),ptrel(mjets)
      include   'hepevt.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu
      real * 8 palg,pp,tmp
c
      real*8 ptmin,yjmax,R
      common /jetcuts/ptmin,yjmax,R
      logical passcuts 
      real *8 wjet(4) 
      integer wtag 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,mjets
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
      enddo
C - Extract final state particles to feed to jet finder
c      print*,'==============='
      do j=1,nhep
         if (isthep(j).eq.1.and.
     1        (abs(idhep(j)).lt.11.or.abs(idhep(j)).gt.16)) then
            if(ntracks.eq.maxtrack) then
               write(*,*) 'analyze: need to increase maxtrack!'
               write(*,*) 'ntracks: ',ntracks
               stop
            endif
            ntracks=ntracks+1
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,j)
            enddo
            itrackhep(ntracks)=j
         endif
      enddo
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------------------------------------------------------------- C
C - Inclusive jet pT and Y spectra are to be compared to CDF data:    - C    
C --------------------------------------------------------------------- C
c
c note: ptmin and R are in common block "jetcuts"
c
c C/A (c.f. ZZ code):
      palg=0
      call fastjetfatjet(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec,wtag,wjet,passcuts)
      
      mjets=min(mjets,njets)
      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4 .and. passcuts) then
         write(*,*) ' bug!'
      endif
c end check consistency
c
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         kt(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(kt(j)**2+pjet(3,j)**2)
         eta(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         rap(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      call computeptrel(ptrack,ntracks,rap,kt,phi,mjets,jetvec,ptrel)
      end
