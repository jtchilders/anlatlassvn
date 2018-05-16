c  The next subroutines opens some histograms and prepares them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist_slm
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
      call pwhgbookup(diag,'sig(inc cuts)','LIN',binsize(diag),0d0,1d0)

C -----------------------------------------------------
C HISTOGRAMS WITH VBF CUTS:  

      diag = 2
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(all VBF cuts)','LIN',
     &                binsize(diag),0d0,1d0)
      
c     cross section sanity check
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
      binsize(diag) = 20d0
      call pwhgbookup(diag,'M_reconstructed(2l2j)- VBF CUTS','LOG',
     .     binsize(diag),0d0,600d0)

      diag = 11
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y j1j2(tag)- VBF CUTS','LIN',
     .     binsize(diag),-6d0,6d0)

      diag = 12
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT decay jet1 - VBF CUTS','LOG',
     .     binsize(diag),0d0,520d0)

      diag = 13
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT decay jet2 - VBF CUTS','LOG',
     .     binsize(diag),0d0,520d0)

      diag = 14
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT LEPT- VBF CUTS','LOG',
     .     binsize(diag),0d0,440d0)

      diag = 15
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT MISS- VBF CUTS','LOG',
     .     binsize(diag),0d0,600d0)

      diag = 16
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Eta_lept- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 17
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Eta_dec_jet1- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag = 18
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Eta_dec_jet2- VBF CUTS','LIN',
     .     binsize(diag),-4d0,4d0)


      diag = 19
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Pt J3- VBF CUTS','LOG',
     .     binsize(diag),0d0,200d0)

      diag = 20
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y J3- VBF CUTS, kt3>pt_tag_min','LOG',
     .     binsize(diag),-5d0,5d0)

      diag = 21
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y J3- VBF CUTS, kt3>pt_cut','LOG',
     .     binsize(diag),-5d0,5d0)

      diag = 22
      binsize(diag) = 50d0
      call pwhgbookup(diag,'M_j1j2(tag) - VBF CUTS','LIN',
     .     binsize(diag),0d0,3000d0)

      diag = 23
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Phi(j1j2)- VBF CUTS','LIN',
     .     binsize(diag),0d0,3.2d0)

      diag = 24
      binsize(diag) = 2d0
      call pwhgbookup(diag,'M(dec_jets)- VBF CUTS','Lin',
     .     binsize(diag),40d0,120d0)
  
      diag = 25
      binsize(diag) = 2d0
      call pwhgbookup(diag,'m(34)- VBF CUTS','LIN',
     .     binsize(diag),40d0,120d0)

      diag = 26
      binsize(diag) = 2d0
      call pwhgbookup(diag,'m(56)- VBF CUTS','LIN',
     .     binsize(diag),60d0,100d0)
 
      diag = 27
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y*- VBF CUTS, kt3>pt_jet_min','LIN',
     .     binsize(diag),-4d0,4d0)
      
      diag = 28
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y*- VBF CUTS, kt3>pt_cut','LIN',
     .     binsize(diag),-4d0,4d0)
   
c=================================
c XSEC WITH VBF CUTS and CJV CUTS: 

      diag = 29
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(CJV cuts)','LIN',
     &                binsize(diag),0d0,1d0)
  
ccccccc
c
      end 
      
c ========================

      subroutine analysis_slm(dsig0)
      implicit none
      real * 8 dsig0,dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'nlegborn.h'
      include  'LesHouches.h'
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
      real * 8  mll,phill
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
      real*8 ystar,rap_av

      integer iveto1,iveto2,nveto
      real*8 pt_v1,pt_v2
      real*8 pt_veto
      logical rap_veto
      logical passcuts_cjv
      real*8 ptl1_min,ptl2_min,mll_min,mll_max,phill_max

      real * 8 powheginput
      external powheginput

      real*8 mjjdec_diff,mjj_pair(maxjet,maxjet),mjj_diff(maxjet,maxjet)
      real*8 pjj(4,maxjet,maxjet)
      integer i,jdec1,jdec2
      integer pdec1,pdec2

      real*8   kt_dec(2),eta_dec(2),rap_dec(2),
     1         phi_dec(2),pj_dec(4,2),ptrel_dec(2)
      real*8   kt_j(maxjet),eta_j(maxjet),rap_j(maxjet),
     1         phi_j(maxjet),pj_j(4,maxjet),ptrel_j(maxjet)
      integer ii,mjets_j
      
      real*8 mdec,mdec_jj
      real*8 eta_dec_max,Rjd_min,Rjd_tmp0,Rjd_tmp1
      real*8 rdec(1:4),rmdec

      logical mdec_window
      real*8 mdec_min,mdec_max

      complex*16 sol1,sol2
      real*8 arg
      real*8 plep(1:4),pneu(1:4),rneu(1:4)
      real*8 rap_veto_eta
      real*8 pw1(1:4),pw2(1:4)
c
c
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
      
      yjj_min = 0d0
      rap_sign = .false.
      rap_gap  = .false.

      Rjj_min = 0d0
      mjj_tag_min = 0d0

      Rjl_min = 0d0
      Rjd_min = 0d0

      mdec_window = .false.
      mdec_min = 0d0
      mdec_max = 1d10

      pt_veto = 1d10
      rap_veto = .false.
      rap_veto_eta = 0d0
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c VBF cuts (can be modified by user):
C
        ptl_min    = 30d0
        ptmiss_min = 30d0

        pt_tag_min = 25d0  !tag jets
        pt_dec_min = 25d0  !W-decay jets
        ptmin = 1d0        !parameter for jet algorithm
        pt_cut = 10d0      !parameter for y3 distributions
        pt_dec_cut = 25d0  !parameter for MW reconstruction
        yj_max  = 4.5d0

        yjj_min  = 3.0d0
        rap_sign = .true.
        rap_gap  = .true.

        Rjj_min = 0.4d0 
        mjj_tag_min = 600d0

        Rjl_min = 0.3d0
        Rjd_min = 0.3d0

        mdec_window = .true.
        mdec_min = 71d0
        mdec_max = 91d0
C
        passcuts_vbf = .true.
        passcuts_cjv = .true.

C CJV cuts:
        pt_veto = 25d0
        rap_veto = .false.
        rap_veto_eta = 3.2d0
C
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
         
            if(((abs(vdecaymodeWm).eq.101.and.abs(vdecaymodeWp).eq.11).or. 
     $        (abs(vdecaymodeWm).eq.103.and.abs(vdecaymodeWp).eq.11).or.
     $        (abs(vdecaymodeWm).eq.107.and.abs(vdecaymodeWp).eq.11).or.  
     $        (abs(vdecaymodeWm).eq.101.and.abs(vdecaymodeWp).eq.13).or.
     $        (abs(vdecaymodeWm).eq.103.and.abs(vdecaymodeWp).eq.13).or.
     $        (abs(vdecaymodeWm).eq.107.and.abs(vdecaymodeWp).eq.13)
     $        ).and.vdecaymodeWp*vdecaymodeWm.lt.0) then

              continue
            else
             write(*,*) '**************************************'
             write(*,*) ' template analysis works only for W+  '
             write(*,*) ' bosons decaying to electrons or muons'
             write(*,*) ' and W- decaying hadronically         '
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
         write(*,*) '                ANALYSIS CUTS                     '
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
            write(*,*) 'W-decay jets pt_dec_min = ',pt_dec_min
            write(*,*) 'W-decay jets pt_dec_cut = ',pt_dec_cut
            write(*,*) 'yj_max  = ',yj_max
            write(*,*) 'yjj_min = ',yjj_min
            write(*,*) 'Rjj_min = ',Rjj_min
            write(*,*) 'mjj_tag_min = ',mjj_tag_min
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
            if (mdec_window) then
               write(*,*) 'mdec_window: W decay jets must fulfill'
               write(*,*) mdec_min,'< Mjj_dec < ',mdec_max
               write(*,*) ''
            endif   
            write(*,*) ''
            write(*,*) 'extra CJV settings'
            write(*,*) 'pt_veto = ',pt_veto
            write(*,*) 'rap_veto = ',rap_veto
            write(*,*) 'rap_veto_eta = ',rap_veto_eta
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
        kt_dec(:) = 0d0

      mjets_j=10
      call buildjets(mjets_j,kt_j,eta_j,rap_j,phi_j,pj_j,ptrel_j)

cccccccccccccccccccccccccccccccccccccccc
c
c now start analysis:

      if (mjets_j.lt.4) return

c require at least 4 hard jets:
      if (kt_j(4).lt.pt_dec_min) return 

c compute m(jet_i,jet_j) to find W->jj pair:
      mjjdec_diff = 1d10
      pj(:,:) = 0d0
      mjj_pair(:,:) = 0d0
      mjj_diff(:,:) = 2d10
      jdec1 = 0
      jdec2 = 0
      do i = 1,4
         do j = i+1,4
            pjj(:,i,j) = pj_j(:,i)+pj_j(:,j)
            mjj_pair(i,j) = invmass(pjj(1,i,j))

            mjj_diff(i,j) = abs(ph_wmass-mjj_pair(i,j))

            if (mjj_diff(i,j).lt.mjjdec_diff) then
               mjjdec_diff = mjj_diff(i,j)
               jdec1 = i
               jdec2 = j
            endif !mjj
         enddo !j
      enddo  !i 

      if (mdec_window) then   
      if ((jdec1.ne.0).and.(jdec2.ne.0)) then 
         if ( (mjj_pair(jdec1,jdec2).lt.mdec_min).or.
     &        (mjj_pair(jdec1,jdec2).gt.mdec_max) ) then 
            return
         endif   
      endif   
      endif     
     
      ii = 0
      do i = 1,mjets_j
         if (i.eq.jdec1) then
            kt_dec(1)    = kt_j(jdec1)
            eta_dec(1)   = eta_j(jdec1)
            rap_dec(1)   = rap_j(jdec1)
            phi_dec(1)   = phi_j(jdec1)
            ptrel_dec(1) = ptrel_j(jdec1)

            pj_dec(1:4,1) = pj_j(1:4,jdec1)
         elseif (i.eq.jdec2) then
            kt_dec(2)    = kt_j(jdec2)
            eta_dec(2)   = eta_j(jdec2)
            rap_dec(2)   = rap_j(jdec2)
            phi_dec(2)   = phi_j(jdec2)
            ptrel_dec(2) = ptrel_j(jdec2)

            pj_dec(1:4,2) = pj_j(1:4,jdec2)
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

      if (mjets.ne.(mjets_j-2)) then
         write(*,*) ' crazy event, wrong number of jets'
         write(*,*) ' number of non-decay jets =',mjets
         write(*,*) ' full number of jets =',mjets_j
         return
      endif 

cccccccccccccc

      if (kt(1).lt.pt_tag_min) return
      if (kt(2).lt.pt_tag_min) return

      diag = 1
      call pwhgfill(diag,0.5d0,dsig)

c************************************

	do ihep = 1,nhep
	if (idhep(ihep).eq.24) then
	 pw1(1:4) = phep(1:4,ihep)
	endif 
	if (idhep(ihep).eq.-24) then
	 pw2(1:4) = phep(1:4,ihep)	 
	endif 
	enddo   

      etael1=etafromp(phep(1,ileptons(1)))
      ptel1=ptfromp(phep(1,ileptons(1)))

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
     &            pj_dec(1:4,1)+pj_dec(1:4,2)

      rmdec = invmass(rdec)

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

      if (kt_dec(1).lt.pt_dec_min)  passcuts_vbf = .false.
      if (kt_dec(2).lt.pt_dec_min)  passcuts_vbf = .false.

      if(abs(eta_dec(1)).gt.eta_dec_max) passcuts_vbf = .false. 
      if(abs(eta_dec(2)).gt.eta_dec_max) passcuts_vbf = .false. 

      if(abs(rap(itag1)).gt.yj_max) passcuts_vbf = .false. 
      if(abs(rap(itag2)).gt.yj_max) passcuts_vbf = .false.
      if(abs(rap(itag1)-rap(itag2)).lt.yjj_min) passcuts_vbf = .false.
      if (rap_sign.and.(rap(itag1)*rap(itag2).ge.0d0)) 
     &     passcuts_vbf = .false.

      if(abs(rap_dec(1)).gt.yj_max) passcuts_vbf = .false. 
      if(abs(rap_dec(2)).gt.yj_max) passcuts_vbf = .false.
       
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

c separation of tagging jets from decay jets:
      Rjd_tmp1 = 1d10
      do ij = 1,mjets
         if(ij.eq.itag1.or.ij.eq.itag2) then
         do il = 1,2
               Rjd_tmp0 = r(pj(1:4,ij),pj_dec(1:4,il))
               Rjd_tmp1 = min(Rjd_tmp0,Rjd_tmp1)
         enddo !il
         endif
      enddo !ij
      if (Rjd_tmp1.lt.Rjd_min) passcuts_vbf = .false. 

      pj12 = pj(:4,itag1)+pj(:4,itag2)
      mjj = invmass(pj12)
      if(mjj.lt.mjj_tag_min)  passcuts_vbf = .false. 

      if (rap_gap.and.(min(rap(itag1),rap(itag2))
     &            .ge.min(etael1,eta_dec(1),eta_dec(2))))
     & passcuts_vbf = .false. 

      if (rap_gap.and.(max(rap(itag1),rap(itag2))
     &            .le.max(etael1,eta_dec(1),eta_dec(2))))
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

C     M_dec(reconstruced)
      diag=10
      call pwhgfill(diag,rmdec,dsig/binsize(diag))

c     Y jet1 - Y jet 2
      diag=11
      pj12 = pj(:4,itag1)+pj(:4,itag2)
      y12 = rap(itag1)-rap(itag2)
      call pwhgfill(diag,y12,dsig/binsize(diag))

C     Pt dec jet1    
      diag=12
      call pwhgfill(diag,kt_dec(1),dsig/binsize(diag))
C     Pt dec jet2   
      diag=13
      call pwhgfill(diag,kt_dec(2),dsig/binsize(diag))

C     Pt lept     
      diag=14
      call pwhgfill(diag,ptel1,dsig/binsize(diag))

C     Pt miss
      diag=15
      call pwhgfill(diag,ptmiss,dsig/binsize(diag))

C     ETA lept
      diag=16
      call pwhgfill(diag,etael1,dsig/binsize(diag))

C     eta dec jet1    
      diag=17
      call pwhgfill(diag,eta_dec(1),dsig/binsize(diag))
C     eta dec jet2   
      diag=18
      call pwhgfill(diag,eta_dec(2),dsig/binsize(diag))


c     pt of third jet
      diag=19
      if(itag3.ne.0) call pwhgfill(diag,kt(itag3),dsig/binsize(diag))

c     y of third jet, kt3>pt_tag_min
      diag=20
      if(itag3.ne.0 .and. kt(itag3) > pt_tag_min ) 
     .     call pwhgfill(diag,rap(itag3),dsig/binsize(diag))
c     y of third jet, kt3>pt_cut
      diag=21
      if(itag3.ne.0 .and. kt(itag3) > pt_cut ) 
     .     call pwhgfill(diag,rap(itag3),dsig/binsize(diag))

C     M_j1j2
      diag=22
      mjj = invmass(pj12)
      if(mjets.ge.2) call pwhgfill(diag,mjj,dsig/binsize(diag))

c     Phi(j1j2) 
      diag=23 
      phijj=fphi(pj(1:4,itag1),pj(1:4,itag2))
      if(mjets.ge.2) call pwhgfill(diag,phijj,dsig/binsize(diag))
     
C     Minv(dec_jets) 
      diag=24
      mdec_jj = mjj_pair(jdec1,jdec2)
      call pwhgfill(diag,mdec_jj,dsig/binsize(diag))

C     M(34) 
      diag=25
      m34 = invmass(pw1)
      call pwhgfill(diag,m34,dsig/binsize(diag))

C     M(56) 
      diag=26
      m56 = invmass(pw2)
      call pwhgfill(diag,m56,dsig/binsize(diag))     

      if(itag3.ne.0 .and. kt(itag3) > pt_tag_min ) then 
          ystar     = rap(itag3)-(rap(itag1)+rap(itag2))/2d0
c         y* 
          diag=27
          call pwhgfill(diag,ystar,dsig/binsize(diag))
      endif

      if(itag3.ne.0 .and. kt(itag3) > pt_cut ) then 
          ystar     = rap(itag3)-(rap(itag1)+rap(itag2))/2d0
c         y* 
          diag=28
          call pwhgfill(diag,ystar,dsig/binsize(diag))
      endif
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


      else ! no rap-gap-constraint on veto jet 
c hardest veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_v1.and.kt(j).gt.pt_veto
     &        .and. abs(rap(j)).lt.rap_veto_eta) then              
              iveto1 = j
              pt_v1 = kt(j)
           end if
      end do

c 2nd hardest veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_v2.and.kt(j).gt.pt_veto.and. 
     &          j.ne.iveto1
     &        .and. abs(rap(j)).lt.rap_veto_eta ) then
              iveto2 = j
              pt_v2 = kt(j)
           end if
      end do

c count veto jet:
      do j = 1, mjets
           if ( kt(j).gt.pt_veto
     &        .and. abs(rap(j)).lt.rap_veto_eta) then
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

      diag = 29
      call pwhgfill(diag,0.5d0,dsig)

      endif !passcuts_cjv

      endif !passcuts_vbf


      end

