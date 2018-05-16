c  The next subroutines opens some histograms and prepares them 
c      to receive data 
c  You can substitute these with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist_lep
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

c ========================
 
      diag = 1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(incl cuts)','LIN',binsize(diag),0d0,1d0)

C -----------------------------------------------------
C     -- HISTOGRAMS WITH VBF CUTS 

      diag = 2
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(all VBF cuts)','LIN',
     &                binsize(diag),0d0,1d0)
      
c  cross section sanity check
      diag = 3
      binsize(diag) = 5d0
      call pwhgbookup(diag,'3 jet incl.- VBF CUTS','LIN',
     .     binsize(diag),2.5d0,72.5d0)

      diag = 4
      binsize(diag) = 5d0
      call pwhgbookup(diag,'2 jet incl.- VBF CUTS','LIN',
     .     binsize(diag),2.5d0,72.5d0)

      diag = 5
      binsize(diag) = 5d0
      call pwhgbookup(diag,'2 jet excl.- VBF CUTS','LIN',
     .     binsize(diag),2.5d0,72.5d0)

      diag = 6
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT jet 1- VBF CUTS','LOG',
     .     binsize(diag),0d0,880d0)

      diag = 7
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT jet 2- VBF CUTS','LOG',
     .     binsize(diag),0d0,520d0)

      diag = 8
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y jet 1- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 9
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y jet 2- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 10
      binsize(diag) = 5d0
      call pwhgbookup(diag,'M(2l2v)- VBF CUTS','LOG',
     .     binsize(diag),0d0,500d0)

      diag = 11
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y j1j2(tag)- VBF CUTS','LIN',
     .     binsize(diag),-6d0,6d0)

      diag = 12
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT LEPT1- VBF CUTS','LOG',
     .     binsize(diag),0d0,440d0)

      diag = 13
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT LEPT2- VBF CUTS','LOG',
     .     binsize(diag),0d0,440d0)

      diag = 14
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT MISS- VBF CUTS','LOG',
     .     binsize(diag),0d0,600d0)

      diag = 15
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Eta_lept1- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 16
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Eta_lept2- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 17
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Eta_l1l2- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 18
      binsize(diag) = 20d0
      call pwhgbookup(diag,'M(l1l2)- VBF CUTS','LOG',
     .     binsize(diag),0d0,960d0)

      diag = 19
      binsize(diag) = 10d0
      call pwhgbookup(diag,'MT_WW- VBF CUTS','LOG',
     .     binsize(diag),0d0,700d0)

      diag = 20
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Phi(l1l2)- VBF CUTS','LIN',
     .     binsize(diag),0d0,3.2d0)

      diag = 21
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Pt J3- VBF CUTS','LOG',
     .     binsize(diag),0d0,200d0)

      diag = 22
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y J3- VBF CUTS, kt3>pt_jet_min','LOG',
     .     binsize(diag),-5d0,5d0)

      diag = 23
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y J3- VBF CUTS, kt3>pt_cut','LOG',
     .     binsize(diag),-5d0,5d0)

      diag = 24
      binsize(diag) = 50d0
      call pwhgbookup(diag,'M_j1j2(tag) - VBF CUTS','LIN',
     .     binsize(diag),0d0,3000d0)

      diag = 25
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Phi(j1j2)- VBF CUTS','LIN',
     .     binsize(diag),0d0,3.2d0)
  
      diag = 26
      binsize(diag) = 2d0
      call pwhgbookup(diag,'m(34)- VBF CUTS','LIN',
     .     binsize(diag),20d0,120d0)

      diag = 27
      binsize(diag) = 2d0
      call pwhgbookup(diag,'m(56)- VBF CUTS','LIN',
     .     binsize(diag),20d0,120d0)
 
      diag = 28
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y*- VBF CUTS, kt3>pt_jet_min','LIN',
     .     binsize(diag),-4d0,4d0)
   
      diag = 29
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y*- VBF CUTS, kt3>pt_cut','LIN',
     .     binsize(diag),-4d0,4d0)
   
  
      diag = 30
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'(y1+y2)/2- VBF CUTS','LIN',
     .     binsize(diag),0d0,4d0)

c
C XSEC WITH VBF CUTS and CJV CUTS
c
      diag = 31
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(CJV cuts)','LIN',
     &                binsize(diag),0d0,1d0)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      end 
      
      subroutine analysis_lep(dsig0)
      implicit none
      real * 8 dsig0,dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'nlegborn.h'
      include  'LesHouches.h'
      include 'cvecbos.h'
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
      data WHCPRG/'NLO   '/
      integer nleptons,nnu
      integer ileptons(200),inu(200)
c
c for mW-distributions:
      integer nel,nenu
      integer iel(200),ienu(200)
c
      integer nmu,nmnu
      integer imu(200),imnu(200)
c
      integer   maxjet
      parameter (maxjet=2048)
      real * 8  kt(maxjet),eta(maxjet),rap(maxjet),
     1    phi(maxjet),pj(4,maxjet),ptrel(maxjet)
      real * 8 ptel1,ptel2,etael1,etael2,ptmiss
      real*8 pen(4),pmn(4)
      real * 8 pj12(4),y12, pel12(4),etael12
      real*8 ptl(100),etal(100)
      real * 8  mll, mtww, rj1l, rj2l, phill
      real * 8 Etmisstilde, Etll,invmass 
      real * 8 mjj,phijj
      real * 8 r,fphi
      integer ihep,j,mjets
      real * 8 etafromp,ptfromp 
      logical passcuts_vbf
     
c cut parameters:
      real*8 ptmin,pt_cut
      real*8 ptmiss_min,etal_max,ptl_min
      real*8 pt_tag_min,pt_jet_min,yj_max
      real*8 yjj_min,Rjj_min,mjj_min
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
      real*8 ystar, rap_av,Rj12

      integer iveto1,iveto2,nveto
      real*8 pt_v1,pt_v2
      real*8 pt_veto
      logical rap_veto
      logical passcuts_cjv

      real*8 m4l
      real*8 ptl1_min,ptl2_min,mll_min,mll_max,phill_max

      real * 8 powheginput
      external powheginput
c
c =======================================================

c initialize all cuts:
C
      etal_max   = 1d10 
      ptl_min    = 0d0  
      ptmiss_min = 0d0

      ptl1_min = 0d0
      ptl2_min = 0d0
      mll_min  = 0d0
   
      pt_tag_min = 0d0
      pt_jet_min = 0d0
      pt_cut = 0d0
      yj_max  = 1d10
      
      yjj_min = 0d0
      rap_sign = .false.
      rap_gap  = .false.

      Rjj_min = 0d0
      mjj_min = 0d0

      Rjl_min = 0d0
      Rll_min = 0d0
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c VBF cuts (can be modified by user):
C
        etal_max   = 2.5d0
        ptl_min    = 15d0
        ptl2_min   = 15d0
        ptl1_min   = 25d0
        mll_min = 0d0   
        mll_max = 1d10
        phill_max = 1.8d0

        ptmiss_min = 25d0

        pt_tag_min = 25d0  !tag jets
        pt_jet_min = 25d0  !non-tag extra jets
        ptmin = 1d0        !parameter for jet algorithm
        pt_cut = 10d0       !parameter for y3 distributions
        yj_max  = 4.5d0

        yjj_min  = 3.8d0
        rap_sign = .true.
        rap_gap  = .true.

        Rjj_min = 0.4d0
        mjj_min = 500d0

        Rjl_min = 0.3d0
        Rll_min = 0.3d0
C
        passcuts_vbf = .true. 
        passcuts_cjv = .true.

C CJV cuts:
        pt_veto = 20d0
        rap_veto = .true.

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      
c from pico to femto      
      dsig=dsig0*1000
      diag = 0 
      if (ini) then

         ph_Hmass  = powheginput('hmass')
         ph_Hwidth = powheginput('hwidth')

         write (*,*)
         write (*,*) '********************************************'
         if(WHCPRG.eq.'NLO   ') then
            write (*,*) 'Fixed order ANALYSIS CALLED '
         
          if(((abs(vdecaymodeWp).eq.11.and.abs(vdecaymodeWm).eq.11).or. 
     $        (abs(vdecaymodeWp).eq.13.and.abs(vdecaymodeWm).eq.13).or.
     $        (abs(vdecaymodeWp).eq.11.and.abs(vdecaymodeWm).eq.13).or.
     $        (abs(vdecaymodeWp).eq.13.and.abs(vdecaymodeWm).eq.11)   
     $        ).and.vdecaymodeWp*vdecaymodeWm.lt.0) then

              continue
           else
             write(*,*) '**************************************'
             write(*,*) ' template analysis works only for W  '
             write(*,*) ' bosons decaying to electrons or muons'
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
         write (*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) '                ANALYSIS CUTS                     '
         write(*,*) '*********************************************'
         write(*,*) '*********************************************'
            write(*,*) ''
            write(*,*) 'Higgs parameters:'
            write(*,*) 'Higgs mass [GeV]= ',ph_Hmass
            write(*,*) 'Higgs width [GeV]= ',ph_Hwidth
            write(*,*) ''
            write(*,*) 'jet cuts:'
            write(*,*) 'jet ptmin = ',ptmin
            write(*,*) 'jet3 pt_cut = ',pt_cut
            write(*,*) 'tag jets pt_tag_min = ',pt_tag_min
            write(*,*) 'non-tag jets pt_jet_min = ',pt_jet_min
            write(*,*) 'yj_max  = ',yj_max
            write(*,*) 'yjj_min = ',yjj_min
            write(*,*) 'Rjj_min = ',Rjj_min
            write(*,*) 'mjj_min = ',mjj_min
            write(*,*) ''
            write(*,*) 'lepton cuts:'
            write(*,*) 'ptl_min    = ',ptl_min 
            write(*,*) 'ptl1_min    = ',ptl1_min 
            write(*,*) 'ptl2_min    = ',ptl2_min 
            write(*,*) 'etal_max   = ',etal_max 
            write(*,*) ''
            write(*,*) 'mll_min    = ',mll_min
            write(*,*) 'mll_max    = ',mll_max
            write(*,*) 'phill_max    = ',phill_max
            write(*,*) 'Rll_min    = ',Rll_min
            write(*,*) 'ptmiss_min    = ',ptmiss_min
            write(*,*) ''
            write(*,*) 'extra cuts:'
            write(*,*) 'Rjl_min  = ',Rjl_min
            write(*,*) 'rap_gap  = ',rap_gap
            write(*,*) 'rap_sign  = ',rap_sign
            write(*,*) ''
            write(*,*) 'extra CJV settings'
            write(*,*) 'pt_veto = ',pt_veto
            write(*,*) 'rap_veto = ',rap_veto
            write(*,*) ''
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
      if(nleptons.lt.2.or.nnu.lt.2) then
         write(*,*) ' crazy event, missing leptons'
         return
      endif
c sort by pt
      call sortbypt(nleptons,ileptons)
      call sortbypt(nnu,inu)         
c
cccccccccc
c
c this is needed/meaningful for computing m(34),m(56) in e-mu mode only:
c find electrons and muons and neutrinos
      nel  = 0
      nenu = 0
      nmu  = 0
      nmnu = 0
      do ihep=1,nhep
         if(isthep(ihep).eq.1) then
            if(abs(idhep(ihep)).eq.11) then
               nel=nel+1
               iel(nel)=ihep
            elseif(abs(idhep(ihep)).eq.13) then
               nmu=nmu+1
               imu(nmu)=ihep
            elseif(abs(idhep(ihep)).eq.12) then
               nenu=nenu+1
               ienu(nenu)=ihep
            elseif(abs(idhep(ihep)).eq.14) then
               nmnu=nmnu+1
               imnu(nmnu)=ihep
            endif
         endif
      enddo
      
c sort by pt
      call sortbypt(nel,iel)
      call sortbypt(nenu,ienu)
      call sortbypt(nmu,imu)
      call sortbypt(nmnu,imnu)

      pen = phep(:4,iel(1))+phep(:4,ienu(1))
      pmn = phep(:4,imu(1))+phep(:4,imnu(1))
c
ccccccccccc      
c
c initialize:
	do j = 1,maxjet
	   kt(j) = 0d0
	enddo   

      mjets=10
      call buildjets(mjets,kt,eta,rap,phi,pj,ptrel)

      if (mjets.lt.2) return

      if (kt(1).lt.pt_jet_min) return
      if (kt(2).lt.pt_jet_min) return

c xsec with incl. cuts:
      diag = 1
      call pwhgfill(diag,0.5d0,dsig)

c************************************

      etael1=etafromp(phep(1,ileptons(1)))
      etael2=etafromp(phep(1,ileptons(2)))

      ptel1=ptfromp(phep(1,ileptons(1)))
      ptel2=ptfromp(phep(1,ileptons(2)))

      pel12 = phep(:4,ileptons(1))+phep(:4,ileptons(2))

      pmiss_sum1 = 0d0
      pmiss_sum2 = 0d0
      do il = 1,nnu ! neutrinos
         pmiss_sum1 = pmiss_sum1 + phep(1,inu(il))
         pmiss_sum2 = pmiss_sum2 + phep(2,inu(il))
      enddo   
      ptmiss = sqrt((pmiss_sum1)**2+
     1              (pmiss_sum2)**2)

      pl_sum(1:4) = pel12(1:4)
      do il = 1,nnu ! neutrinos
         pl_sum(1:4) = pl_sum(1:4) + phep(1:4,inu(il))
      enddo     

      m4l = invmass(pl_sum)
      mll = invmass(pel12)

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

c identify 2 hardest leptons: 
      ptl_1 = 0d0
      ptl_2 = 0d0
      ltag1 = 0
      ltag2 = 0
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
      do j = 1,nleptons
           if ( ptl(j).gt.ptl_2 .and. ptl(j).gt.ptl_min 
     &          .and. abs(etal(j)).lt.etal_max
     &          .and. j.ne.ltag1 
     &         ) then
              ltag2 = j    
              ptl_2  = ptl(j)
           end if         
      end do
      if (ltag1.eq.0.or.ltag2.eq.0) passcuts_vbf = .false.
      if (ltag1.eq.0.or.ltag2.eq.0) return 

      pmiss_sum1 = 0d0
      pmiss_sum2 = 0d0
      do il = 1,nnu !neutrinos
         pmiss_sum1 = pmiss_sum1 + phep(1,inu(il))
         pmiss_sum2 = pmiss_sum2 + phep(2,inu(il))
      enddo  
      ptmiss = sqrt((pmiss_sum1)**2+
     1              (pmiss_sum2)**2)

      pel12 = phep(:4,ileptons(ltag1))+phep(:4,ileptons(ltag2))

      pl_sum(1:4) = pel12(1:4)
      do il = 1,nnu ! neutrinos
         pl_sum(1:4) = pl_sum(1:4) + phep(1:4,inu(il))
      enddo       

      phill = fphi(phep(:4,ileptons(ltag1)),phep(:4,ileptons(ltag2)))
c
cccccccccccccccccccccccccccccccccccccc
c
c apply VBF cuts:
c
      if(mjets.lt.2) passcuts_vbf = .false.      
      if (itag1.eq.0.or.itag2.eq.0) passcuts_vbf = .false. 
      if (ltag1.eq.0.or.ltag2.eq.0) passcuts_vbf = .false.

      if (kt(itag1).lt.pt_tag_min) passcuts_vbf = .false.   
      if (kt(itag2).lt.pt_tag_min) passcuts_vbf = .false. 
c
      if(abs(rap(itag1)).gt.yj_max) passcuts_vbf = .false. 
      if(abs(rap(itag2)).gt.yj_max) passcuts_vbf = .false.
      if(abs(rap(itag1)-rap(itag2)).lt.yjj_min) passcuts_vbf = .false.
      if (rap_sign.and.(rap(itag1)*rap(itag2).ge.0d0)) 
     &     passcuts_vbf = .false.
c      
      Rjj_tmp = r(pj(1:4,itag1),pj(1:4,itag2))
      if (Rjj_tmp.lt.Rjj_min) passcuts_vbf = .false.
      
      if(ptmiss.lt.ptmiss_min) passcuts_vbf = .false. 
     
      etael1=etal(ltag1)
      if(abs(etael1).gt.etal_max) passcuts_vbf = .false. 
      etael2=etal(ltag2)
      if(abs(etael2).gt.etal_max)  passcuts_vbf = .false. 

      ptel1=ptl(ltag1)
      if(ptel1.lt.ptl1_min)  passcuts_vbf = .false. 
      ptel2=ptl(ltag2)
      if(ptel2.lt.ptl2_min)  passcuts_vbf = .false. 


      pel12 = phep(:4,ileptons(ltag1))+phep(:4,ileptons(ltag2))
      mll = invmass(pel12)
      if(mll.lt.mll_min) passcuts_vbf = .false.
      if(mll.gt.mll_max) passcuts_vbf = .false.

      if(phill.gt.phill_max) passcuts_vbf = .false.

c separation of tagging jets from leading charged leptons:
      Rjl_tmp1 = 1d10
      do ij = 1,mjets
         if(ij.eq.itag1.or.ij.eq.itag2) then
         do il = 1,nleptons
            if(il.eq.ltag1.or.il.eq.ltag2) then
               Rjl_tmp0 = r(pj(1:4,ij),phep(1:4,ileptons(il)))
               Rjl_tmp1 = min(Rjl_tmp0,Rjl_tmp1)
            endif   
         enddo !il
         endif
      enddo !ij
      if (Rjl_tmp1.lt.Rjl_min) passcuts_vbf = .false.

      Rll_tmp = r(phep(1:4,ileptons(ltag1)),
     &            phep(1:4,ileptons(ltag2)))
      if (Rll_tmp.lt.Rll_min) passcuts_vbf = .false.

      pj12 = pj(:4,itag1)+pj(:4,itag2)
      mjj = invmass(pj12)
      if(mjj.lt.mjj_min)  passcuts_vbf = .false. 

      if (rap_gap.and.
     &     (min(rap(itag1),rap(itag2)).ge.min(etael1,etael2)))
     & passcuts_vbf = .false.

      if (rap_gap.and.
     &     (max(rap(itag1),rap(itag2)).le.max(etael1,etael2)))
     & passcuts_vbf = .false.

c*************************

      if (passcuts_vbf) then

      diag = 2  
      call pwhgfill(diag,0.5d0,dsig)

c     three jet inclusive.
      diag=3
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
      diag=4
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
      diag=5

c     PT of jet 1
      diag=6
      call pwhgfill(diag,kt(itag1),dsig/binsize(diag))
c     PT of jet 2
      diag=7
      call pwhgfill(diag,kt(itag2),dsig/binsize(diag))
c     Y of jet 1
      diag=8
      call pwhgfill(diag,rap(itag1),dsig/binsize(diag))
c     Y of jet 2
      diag=9
      call pwhgfill(diag,rap(itag2),dsig/binsize(diag))

C     M(2l2v) 
      diag=10
      call pwhgfill(diag,m4l,dsig/binsize(diag))

c     Y jet1 - Y jet 2
      diag=11
      pj12 = pj(:4,itag1)+pj(:4,itag2)
      y12 = rap(itag1)-rap(itag2)
      call pwhgfill(diag,y12,dsig/binsize(diag))

C     Pt lept1     
      diag=12
      call pwhgfill(diag,ptel1,dsig/binsize(diag))
C     Pt lept2     
      diag=13
      call pwhgfill(diag,ptel2,dsig/binsize(diag))

C     Pt miss
      diag=14
      call pwhgfill(diag,ptmiss,dsig/binsize(diag))

C     ETA lept1
      diag=15
      call pwhgfill(diag,etael1,dsig/binsize(diag))
C     eta lept2     
      diag=16
      call pwhgfill(diag,etael2,dsig/binsize(diag))

C     ETA(l1)- ETA(l2) 
      diag=17
      etael12 = etael1-etael2

      call pwhgfill(diag,etael12,dsig/binsize(diag))

C     M(l1l2) 
      diag=18
      pel12 = phep(:4,ileptons(ltag1))+phep(:4,ileptons(ltag2))
      mll = invmass(pel12)
      call pwhgfill(diag,mll,dsig/binsize(diag))

C     MT_WW
      diag=19
      Etmisstilde = sqrt(ptmiss**2+mll**2)
      Etll=sqrt(ptfromp(pel12)**2
     1     + mll**2)
      mtww = sqrt((Etll+Etmisstilde)**2
     1     -ptfromp(pl_sum)**2)
      call pwhgfill(diag,mtww,dsig/binsize(diag))

c     Phi(l1l2) 
      diag=20 
      phill=fphi(phep(1,ileptons(ltag1)),phep(1,ileptons(ltag2)))
      call pwhgfill(diag,phill,dsig/binsize(diag))


c     pt of third jet
      diag=21
      if(itag3.ne.0) call pwhgfill(diag,kt(itag3),dsig/binsize(diag))

c     y of third jet, kt3>pt_jet_min
      diag=22
      if(itag3.ne.0 .and. kt(itag3) > pt_jet_min ) 
     .     call pwhgfill(diag,rap(itag3),dsig/binsize(diag))
c     y of third jet, kt3>pt_cut
      diag=23
      if(itag3.ne.0 .and. kt(itag3) > pt_cut ) 
     .     call pwhgfill(diag,rap(itag3),dsig/binsize(diag))

C     M_j1j2
      diag=24
      mjj = invmass(pj12)
      if(mjets.ge.2) call pwhgfill(diag,mjj,dsig/binsize(diag))

c     Phi(j1j2) 
      diag=25
      phijj=fphi(pj(1:4,itag1),pj(1:4,itag2))
      if(mjets.ge.2) call pwhgfill(diag,phijj,dsig/binsize(diag))

C     M(34) 
      diag=26
      m34 = invmass(pen)
      call pwhgfill(diag,m34,dsig/binsize(diag))

C     M(56) 
      diag=27
      m56 = invmass(pmn)
      call pwhgfill(diag,m56,dsig/binsize(diag))     

      if(itag3.ne.0 .and. kt(itag3) > pt_jet_min ) then 
          ystar     = rap(itag3)-(rap(itag1)+rap(itag2))/2d0
c         y* 
          diag=28
          call pwhgfill(diag,ystar,dsig/binsize(diag))
      endif

      if(itag3.ne.0 .and. kt(itag3) > pt_cut ) then 
          ystar     = rap(itag3)-(rap(itag1)+rap(itag2))/2d0
c         y* 
          diag=29
          call pwhgfill(diag,ystar,dsig/binsize(diag))
      endif

c     (y1+y2)/2
      diag=30
      rap_av = (rap(itag1)+rap(itag2))/2d0
      call pwhgfill(diag,rap_av,dsig/binsize(diag))
c
ccccccccccccccccccccccccccccccccccccccccccccc
c
      passcuts_cjv = .true.

c identify 2 hardest veto jets:
      pt_v1 = 0d0
      pt_v2 = 0d0
      iveto1 = 0
      iveto2 = 0

      nveto = 0

      if (rap_veto) then 
c hardest veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_v1.and.kt(j).gt.pt_veto .and.(
     &          rap(j).gt.min(rap(itag1),rap(itag2)).and.
     &          rap(j).lt.max(rap(itag1),rap(itag2)))
     &        ) then              
              iveto1 = j
              pt_v1 = kt(j)
           end if
      end do

c 2nd hardest veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_v2.and.kt(j).gt.pt_veto .and.(
     &          rap(j).gt.min(rap(itag1),rap(itag2)).and.
     &          rap(j).lt.max(rap(itag1),rap(itag2))).and. 
     &          j.ne.iveto1 ) then
              iveto2 = j
              pt_v2 = kt(j)
           end if
      end do

c count veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_veto .and.
     &          rap(j).gt.min(rap(itag1),rap(itag2)).and.
     &          rap(j).lt.max(rap(itag1),rap(itag2))) then
              nveto = nveto+1
           end if
      end do


      else ! no rap-constraint on veto jet 
c hardest veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_v1.and.kt(j).gt.pt_veto
     &        ) then              
              iveto1 = j
              pt_v1 = kt(j)
           end if
      end do

c 2nd hardest veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_v2.and.kt(j).gt.pt_veto.and. 
     &          j.ne.iveto1 ) then
              iveto2 = j
              pt_v2 = kt(j)
           end if
      end do

c count veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_veto) then
              nveto = nveto+1
           end if
      end do

      endif !identification of veto jets

      if (nveto.ne.0) passcuts_cjv = .false.       

c
cccccccccccccccccccccccc
c
c histograms with VBF and extra CJV cuts:
c
      if (passcuts_cjv) then

      diag = 31  
      call pwhgfill(diag,0.5d0,dsig)

      endif !passcuts_cjv



      endif !passcuts_vbf

      end
      
      function etafromp(p)
      implicit none
      real * 8 p(4),etafromp,pp
      pp=sqrt(p(1)**2+p(2)**2+p(3)**2)
      if (pp-abs(p(3)) .lt. 1d-13) then 
         etafromp = 100d0 
      else
         etafromp=log((pp+p(3))/(pp-p(3)))/2
      endif
      end
      
      function ptfromp(p)
      implicit none
      real * 8 p(4),ptfromp
      ptfromp=sqrt(p(1)**2+p(2)**2)
      end

      function invmass(p)
      implicit none
      real * 8 p(1:4),invmass
      invmass=sqrt(abs(p(4)**2-p(1)**2-p(2)**2-p(3)**2))
      end


c---- calculate the jets separation between p1 and p2
      double precision function r(p1,p2)
      implicit none
      double precision p1(4),p2(4),r2,dely,delphi,e1,e2
      
      e1=dsqrt(p1(1)**2+p1(2)**2+p1(3)**2)
      e2=dsqrt(p2(1)**2+p2(2)**2+p2(3)**2)
      
      dely = (e1+p1(3))*(e2-p2(3))/
     .     ((e2+p2(3))*(e1-p1(3)))
      dely = 0.5d0*dlog(dely)
      
      r2= (p1(1)*p2(1)+p1(2)*p2(2))
     .     /dsqrt((p1(1)**2+p1(2)**2)*(p2(1)**2+p2(2)**2))
      if (r2 .gt. +0.9999999D0) r2=+1d0
      if (r2 .lt. -0.9999999D0) r2=-1d0
      delphi=dacos(r2)

      r=dsqrt(dely**2+delphi**2)
      
      end


C---  calculate azimuthal angle between vectors 
      double precision function fphi(p1,p2)
      implicit none
      double precision p1(4),p2(4)
      double precision pi
      parameter(pi=3.14159265358979d0)
   
      fphi=p1(1)*p2(1)+p1(2)*p2(2)
      fphi=fphi/dsqrt(p1(1)**2+p1(2)**2)
      fphi=fphi/dsqrt(p2(1)**2+p2(2)**2)
      if     (fphi .gt. +0.9999999D0) then
        fphi=0d0
      elseif (fphi .lt. -0.9999999D0) then
        fphi=pi
      else
        fphi=dacos(fphi)
      endif

      end
      
      function idigit(k,l)
      implicit none
      integer idigit,k,l
      idigit=abs(mod(l,10**k)/10**(k-1))
      end

      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,ptrel)
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

c this is for kt-algorithm:
c      call fastjetktwhich(ptrack,ntracks,ptmin,R,
c     #     pjet,njets,jetvec) 
c
c anti-kT (c.f. ZZ code):
      palg=-1
c      r=rr
c      ptmin=1d0 
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)

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
      if(tmp.gt.1d-4) then
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

      logical function bson(j,jb)
      implicit none
      integer j,jb
      include   'hepevt.h'
      integer jcurr
      logical bhadrstable
      jcurr=j
c     This only happens in parton level analysis
      if(abs(idhep(jcurr)).eq.5) then
         bson=.true.
         jb=jcurr
         return
      endif
 1    continue
      bson=.false.
      if(bhadrstable(idhep(jcurr))) then
         bson=.true.
         jb=jcurr
         return
      endif
      jcurr=jmohep(1,jcurr)
      if(idhep(jcurr).eq.0) then
         bson=.false.
         return
      endif
      goto 1
      end

      logical function bhadr(idhep)
      implicit none
      integer idhep
      integer i1,i2,idigit
      i1=idigit(1,idhep)
      if(i1.eq.1) then
c         is a bottomed meson
         i2=idigit(3,idhep)
      elseif(i1.eq.2) then
c is a bottomed barion
         i2=idigit(4,idhep)
      endif
      if(i2.eq.5) then
         bhadr=.true.
      else
         bhadr=.false.
      endif
      end

      logical function bhadrstable(idhep)
      implicit none
      integer idhep
      integer i1,i2,i3,idigit
      i1=idigit(1,idhep)
      if(i1.eq.1) then
c         is a bottomed meson
         i2=idigit(3,idhep)
         i3=idhep/1000
      elseif(i1.eq.2) then
c is a bottomed barion
         i2=idigit(4,idhep)
         i3=idhep/10000
      else
         bhadrstable=.false.
c         write(133,*) idhep,' false '
         return
      endif
      if(i2.eq.5.and.i3.eq.0) then
         bhadrstable=.true.
c         write(133,*) idhep,' true '
      else
         bhadrstable=.false.
c         write(133,*) idhep,' false '
      endif
      end

      subroutine incbhadrons(jb)
      implicit none
      include   'hepevt.h'      
      integer jb
      integer barray(1000),bnum
      common/cbarray/barray,bnum
      integer k
      if(jb.eq.-1) then
         bnum=0
         return
      endif
      do k=1,bnum
         if(jb.eq.barray(k)) return
         if(jmohep(1,barray(k)).eq.jb) return
         if(jmohep(1,jb).eq.barray(k)) then
            barray(k)=jb
            return
         endif
      enddo
      bnum=bnum+1
      barray(bnum)=jb
      end

      subroutine sortbypt(n,iarr)
      implicit none
      integer n,iarr(n)
      include 'hepevt.h'
      integer j,k
      real * 8 tmp,pt(nmxhep)
      logical touched
      do j=1,n
         pt(j)=sqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
      enddo
c bubble sort
      touched=.true.
      do while(touched)
         touched=.false.
         do j=1,n-1
            if(pt(j).lt.pt(j+1)) then
               k=iarr(j)
               iarr(j)=iarr(j+1)
               iarr(j+1)=k
               tmp=pt(j)
               pt(j)=pt(j+1)
               pt(j+1)=tmp
               touched=.true.
            endif
         enddo
      enddo
      end


      subroutine computeptrel(ptracks,ntracks,rapjets,ktjets,phijets,
     1     njets,jetvec,ptrel)
      implicit none
      integer ntracks,njets,jetvec(ntracks)
      real * 8 ptracks(4,ntracks),rapjets(njets),
     1     ktjets(njets),phijets(njets),ptrel(njets)
      integer j,i
      real * 8 yj,kj1,kj2,y,pt(3)
      do j=1,njets
         ptrel(j)=0
      enddo
      do i=1,ntracks
         j=jetvec(i)
         if(j.gt.0.and.j.le.njets) then
c Track i belongs to jet j
            yj=rapjets(j)
            kj1=ktjets(j)*cos(phijets(j))
            kj2=ktjets(j)*sin(phijets(j))
c rapidity of track i
            y=0.5d0*log((ptracks(4,i)+ptracks(3,i))
     1                 /(ptracks(4,i)-ptracks(3,i)))
c rapidity of track i in frame where the jet has zero rapidity
            y=y-yj
c find momentum of track i in frame where the jet has zero rapidity
            pt(1)=ptracks(1,i)
            pt(2)=ptracks(2,i)
            pt(3)=sqrt(pt(1)**2+pt(2)**2)*sinh(y)
c pt rel is sum of the ptrack momentum projection ortogonal to the jet
c momentum in the frame where the jet has zero rapidity
            ptrel(j)=sqrt(((pt(1)*kj2-pt(2)*kj1)**2+
     1                     (         -pt(3)*kj2)**2+
     2                     (pt(3)*kj1          )**2)/
     3                     (kj1**2+kj2**2)) + ptrel(j)
         endif
      enddo
      end
