c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include '../pwhg_book.h'
      integer diag
      real * 8 binsize(100)
      common/pwhghistcommon/binsize
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
            
      call pwhginihist

      diag=0

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig','LIN',binsize(diag),0d0,10d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(jetcuts)','LIN',binsize(diag),0d0,1d0)


C     -- HISTOGRAMS WITH LEPTON CUTS 
      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'sig(all cuts)','LIN',binsize(diag),0d0,1d0)

      
c     total cross section sanity check
      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'2 jet incl.','LIN',binsize(diag),25d0,75d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'2 jet excl.','LIN',binsize(diag),25d0,75d0)

      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'PT jet 1','LOG',binsize(diag),0d0,900d0)

      diag=diag+1
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT jet 2','LOG',binsize(diag),0d0,520d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta jet 1','LIN',binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'Eta jet 1 asym','LIN',binsize(diag),0d0,4d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta jet 2','LIN',binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'Eta jet 2 asym','LIN',binsize(diag),0d0,4d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta_j1j2','LIN',binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 200d0
      call pwhgbookup(diag,'HT,TOT','LOG',binsize(diag),0d0,3000d0)

      diag=diag+1
      binsize(diag) = 200d0
      call pwhgbookup(diag,'HT,TOT 3j','LOG',binsize(diag),0d0,3000d0)

      diag=diag+1
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT LEPT','LOG',binsize(diag),0d0,450d0)

      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'PT MISS','LOG',binsize(diag),0d0,600d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta_lept','LIN',binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta_l1l2','LIN',binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 60d0
      call pwhgbookup(diag,'M(l1l2)','LOG',binsize(diag),0d0,1000d0)

      diag=diag+1
      binsize(diag) = 100d0
      call pwhgbookup(diag,'MT_WW','LOG',binsize(diag),0d0,1500d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'DelR(j1lep)','LIN',binsize(diag),0d0,6d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'DelR(j2lept)','LIN',binsize(diag),0d0,6d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Phi(l1l2)','LIN',binsize(diag),0d0,3.2d0)

      diag=diag+1
      binsize(diag) = 5d0
      call pwhgbookup(diag,'Pt J3','LOG',binsize(diag),0d0,500d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y J3','LOG',binsize(diag),-5d0,5d0)

      diag=diag+1
      binsize(diag) = 5d0
      call pwhgbookup(diag,'Y J3 asym','LOG',binsize(diag),0d0,5d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Ptrel J1','LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Ptrel J2','LOG',binsize(diag),0d0,200d0)


C     -- SAME HISTOGRAMS WITHOUT LEPTON CUTS 
      
c     total cross section sanity check
      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'2 jet incl.- NO CUTS','LIN',
     .     binsize(diag),25d0,75d0)

      diag=diag+1
      binsize(diag) = 10d0
      call pwhgbookup(diag,'2 jet excl.- NO CUTS','LIN',
     .     binsize(diag),25d0,75d0)

      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'PT jet 1- NO CUTS','LOG',
     .     binsize(diag),0d0,900d0)

      diag=diag+1
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT jet 2- NO CUTS','LOG',
     .     binsize(diag),0d0,520d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta jet 1- NO CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'Eta jet 1 asym- NO CUTS','LIN',
     .     binsize(diag),0d0,4d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta jet 2- NO CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'Eta jet 2 asym- NO CUTS','LIN',
     .     binsize(diag),0d0,4d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta_j1j2- NO CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 200d0
      call pwhgbookup(diag,'HT,TOT- NO CUTS','LOG',
     .     binsize(diag),0d0,3000d0)

      diag=diag+1
      binsize(diag) = 200d0
      call pwhgbookup(diag,'HT,TOT 3j- NO CUTS','LOG',
     .     binsize(diag),0d0,3000d0)

      diag=diag+1
      binsize(diag) = 20d0
      call pwhgbookup(diag,'PT LEPT- NO CUTS','LOG',
     .     binsize(diag),0d0,450d0)

      diag=diag+1
      binsize(diag) = 40d0
      call pwhgbookup(diag,'PT MISS- NO CUTS','LOG',
     .     binsize(diag),0d0,600d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta_lept- NO CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'Eta_l1l2- NO CUTS','LIN',
     .     binsize(diag),-4d0,4d0)

      diag=diag+1
      binsize(diag) = 60d0
      call pwhgbookup(diag,'M(l1l2)- NO CUTS','LOG',
     .     binsize(diag),0d0,1000d0)

      diag=diag+1
      binsize(diag) = 100d0
      call pwhgbookup(diag,'MT_WW- NO CUTS','LOG',
     .     binsize(diag),0d0,1500d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'DelR(j1lep)- NO CUTS','LIN',
     .     binsize(diag),0d0,6d0)

      diag=diag+1
      binsize(diag) = 0.4d0
      call pwhgbookup(diag,'DelR(j2lept)- NO CUTS','LIN',
     .     binsize(diag),0d0,6d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Phi(l1l2)- NO CUTS','LIN',
     .     binsize(diag),0d0,3.2d0)

      diag=diag+1
      binsize(diag) = 5d0
      call pwhgbookup(diag,'Pt J3- NO CUTS','LOG',
     .     binsize(diag),0d0,500d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Y J3- NO CUTS','LOG',
     .     binsize(diag),-5d0,5d0)

      diag=diag+1
      binsize(diag) = 5d0
      call pwhgbookup(diag,'Y J3 asym- NO CUTS','LOG',
     .     binsize(diag),0d0,5d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Ptrel J1- NO CUTS','LOG',
     .     binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Ptrel J2- NO CUTS','LOG',
     .     binsize(diag),0d0,200d0)

      end 

      
      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0,dsig
      include 'hepevt.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h' 
      include 'pwhg_rad.h' 
      include  'LesHouches.h'
      logical ini
      data ini/.true./
      save ini
c     binsize
      integer diag
      real * 8 binsize(100)
      common/pwhghistcommon/binsize
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer nelectrons,nnue,nmuons,nnumu
      integer ielectrons(100),inue(100),imuons(100),inumu(100)
      integer ntaus,nnutau
      integer itaus(100),inutau(100)

      integer   maxjet
      parameter (maxjet=2048)
      real * 8  kt(maxjet),eta(maxjet),rap(maxjet),
     1    phi(maxjet),pj(4,maxjet),ptrel(maxjet)
      real * 8 ptel1,ptel2,etael1,etael2,ptmiss,httot
      real * 8 pmiss_sum1,pmiss_sum2
      real * 8 pj12(4),eta12, pel12(4),etael12
      real * 8 mll, mtww, rj1l, rj2l, phill, Etmisstilde, Etll,invmass 
      real * 8 r,fphi
      integer ihep,j,mjets
      real * 8 etafromp,ptfromp
      logical passcuts 
      integer ikinreg,iuborn,ialr,lept1,lept2,il
      common/iargs/ikinreg,iuborn,ialr
c      if(ikinreg.gt.0) then
c         if(rad_kinreg.ne.ikinreg) return
c      endif
c      if(iuborn.gt.0) then
c         if(rad_ubornidx.ne.iuborn) return
c      endif
c      if(ialr.gt.0) then
c         if(rad_realalr.ne.ialr) return
c      endif
      passcuts = .true. 
c from pico to femto      
      dsig=dsig0*1000
      diag = 0 
      if (ini) then
         write (*,*)
         write (*,*) '********************************************'
         if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
            write (*,*) '           NLO ANALYSIS CALLED        '
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS CALLED     '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS CALLED     '
         endif
         write (*,*) '********************************************'
         ini=.false.
      endif
         
         
      

c find electrons and muons
      nelectrons=0
      nnue=0
      nmuons=0
      nnumu=0
      ntaus=0
      nnutau=0

      do ihep=1,nhep
         if(isthep(ihep).eq.1) then
            if(abs(idhep(ihep)).eq.11) then
               nelectrons=nelectrons+1
               ielectrons(nelectrons)=ihep
            elseif(abs(idhep(ihep)).eq.12) then
               nnue=nnue+1
               inue(nnue)=ihep
            elseif(abs(idhep(ihep)).eq.13) then
               nmuons=nmuons+1
               imuons(nmuons)=ihep
            elseif(abs(idhep(ihep)).eq.14) then
               nnumu=nnumu+1
               inumu(nnumu)=ihep
            elseif(abs(idhep(ihep)).eq.15) then
               ntaus=ntaus+1
               itaus(ntaus)=ihep
            elseif(abs(idhep(ihep)).eq.16) then
               nnutau=nnutau+1
               inutau(nnutau)=ihep

            endif
         endif
      enddo
      if(nelectrons.gt.100.or.nnue.gt.100
     1     .or.nmuons.gt.100.or.nnumu.gt.100 
     1     .or.ntaus.gt.100.or.nnutau.gt.100) then
         write(*,*) ' crazy event, too many leptons'
         return
      endif
      if(nelectrons+nmuons+ntaus.lt.2 .or. nnue+nnumu+nnutau.lt.2) then
         write(*,*) ' crazy event, missing leptons'
         return
      endif
c sort by pt
      call sortbypt(nelectrons,ielectrons)
      call sortbypt(nnue,inue)
      call sortbypt(nmuons,imuons)
      call sortbypt(nnumu,inumu)
      call sortbypt(ntaus,itaus)
      call sortbypt(nnutau,inutau)


      if (nelectrons >= 2) then 
         lept1 = ielectrons(1)
         lept2 = ielectrons(2)
      elseif (nmuons >= 2) then 
         lept1 = imuons(1)
         lept2 = imuons(2)
      elseif (ntaus >= 2) then 
         lept1 = itaus(1)
         lept2 = itaus(2)
      elseif (nmuons == 1 .and. nelectrons == 1) then 
         lept1 = ielectrons(1)
         lept2 = imuons(1)
      elseif (ntaus == 1 .and. nelectrons == 1) then 
         lept1 = ielectrons(1)
         lept2 = itaus(1)
      elseif (ntaus == 1 .and. nmuons == 1) then 
         lept1 = itaus(1)
         lept2 = imuons(1)
      else
         stop 'can not find two charged leptons' 
      endif 



      mjets=10
      call buildjets(mjets,kt,eta,rap,phi,pj,ptrel)

      diag = diag+1
c avoid configurations with light partons in underlying Born configuration
c with too low pt, or too low invariant mass
      call pwhgfill(diag,0.5d0,dsig)


c********** CUTS
C      if(mjets.lt.2) return

      diag = diag+1
      if (kt(2) .ge.30d0) call pwhgfill(diag,0.5d0,dsig)

      pmiss_sum1 = 0d0
      pmiss_sum2 = 0d0
      do il = 1,nnue !electron neutrinos
         pmiss_sum1 = pmiss_sum1 + phep(1,inue(il))
         pmiss_sum2 = pmiss_sum2 + phep(2,inue(il))
      enddo  
      do il = 1,nnumu !muon neutrinos
         pmiss_sum1 = pmiss_sum1 + phep(1,inumu(il))
         pmiss_sum2 = pmiss_sum2 + phep(2,inumu(il))
      enddo   
      do il = 1,nnutau !tau neutrinos
         pmiss_sum1 = pmiss_sum1 + phep(1,inutau(il))
         pmiss_sum2 = pmiss_sum2 + phep(2,inutau(il))
      enddo   

      ptmiss = sqrt((pmiss_sum1)**2+
     1              (pmiss_sum2)**2)

C      if(ptmiss.lt.30) return
      if(ptmiss.lt.30) passcuts = .false. 
      
      etael1=etafromp(phep(1,lept1))
C      if(abs(etael1).gt.2.4) return
      if(abs(etael1).gt.2.4d0) passcuts = .false. 

      etael2=etafromp(phep(1,lept2))
C      if(abs(etael2).gt.2.4) return
      if(abs(etael2).gt.2.4d0)  passcuts = .false. 

      ptel1=ptfromp(phep(1,lept1))
C      if(ptel1.lt.20) return
      if(ptel1.lt.20)  passcuts = .false. 

      ptel2=ptfromp(phep(1,lept2))
C      if(ptel2.lt.20) return
      if(ptel2.lt.20)  passcuts = .false. 
      
      diag = diag+1
      if (kt(2) .ge. 30d0) call pwhgfill(diag,0.5d0,dsig)

C     -- first old histograms with lepton cuts 
      if (passcuts) then 


C      diag=0
c     two jet inclusive and two jet exclusive.
c     It is exclusive if the third jet is below the cut
c     inclusive otherwise
      diag=diag+1
      do j=1,3
         if(mjets.lt.j) kt(j)=0
      enddo
      do j=30,70,10
         if(kt(1).gt.j.and.kt(2).gt.j) then
c inclusive
            call pwhgfill(diag,dble(j),dsig)
            if(kt(3).lt.j) then
c exclusive
               call pwhgfill(diag+1,dble(j),dsig)
            endif
         endif
      enddo
      diag=diag+1

c     PT of jet 1
      diag=diag+1
      call pwhgfill(diag,kt(1),dsig/binsize(diag))
c     PT of jet 2
      diag=diag+1
      call pwhgfill(diag,kt(2),dsig/binsize(diag))
c     Eta of jet 1
      diag=diag+1
      call pwhgfill(diag,eta(1),dsig/binsize(diag))
c     Eta of jet 1 asym
      diag=diag+1
      call pwhgfill(diag,abs(eta(1)),sign(dsig/binsize(diag),eta(1)))
c     Eta of jet 2
      diag=diag+1
      call pwhgfill(diag,eta(2),dsig/binsize(diag))
c     Eta of jet 2 asym
      diag=diag+1
      call pwhgfill(diag,abs(eta(2)),sign(dsig/binsize(diag),eta(2)))
c     Eta jet1 - Eta jet 2
      diag=diag+1
      pj12 = pj(:4,1)+pj(:4,2)
      eta12 = etafromp(pj12)
      call pwhgfill(diag,eta12,dsig/binsize(diag))
c     H_T tot
      httot=ptel1+ptel2+ptmiss
      do j=1,mjets
         httot=httot+kt(j)
      enddo
      diag=diag+1
      call pwhgfill(diag,httot,dsig/binsize(diag))
c     H_T tot, 3j
      httot=ptel1+ptel2+ptmiss
      do j=1,min(mjets,3)
         httot=httot+kt(j)
      enddo
      diag=diag+1
      call pwhgfill(diag,httot,dsig/binsize(diag))

C     Pt lept     
      diag=diag+1
      call pwhgfill(diag,ptel1,dsig/binsize(diag)/2d0)
      call pwhgfill(diag,ptel2,dsig/binsize(diag)/2d0)

C     Pt miss
      diag=diag+1
      call pwhgfill(diag,ptmiss,dsig/binsize(diag))

C     Eta lept
      diag=diag+1
      call pwhgfill(diag,etael1,dsig/binsize(diag)/2d0)
      call pwhgfill(diag,etael2,dsig/binsize(diag)/2d0)

C     eta(l1)- eta(l2) 
      diag=diag+1
      pel12 = phep(:4,lept1)+phep(:4,lept2)
      etael12 = etafromp(pel12)
      call pwhgfill(diag,etael12,dsig/binsize(diag))

C     M(l1l2) 
      diag=diag+1
      mll = invmass(pel12)
      call pwhgfill(diag,mll,dsig/binsize(diag))

C     MT_WW
      diag=diag+1
      Etmisstilde = sqrt(ptmiss**2+mll**2)
      Etll=sqrt(ptfromp(pel12)**2
     1     + mll**2)
      mtww = sqrt((Etll+Etmisstilde)**2
     1     -ptfromp(pel12
     2     +phep(1:4,inue(1))+phep(1:4,inue(2)))**2)
      call pwhgfill(diag,mtww,dsig/binsize(diag))

C     Delta R(j1 lept) 
      diag=diag+1
      if (mjets .ge. 2) then 
      Rj1l = r(pj(1:4,1),phep(1:4,lept1))
      call pwhgfill(diag,Rj1l,dsig/binsize(diag)/2d0)
      Rj1l = r(pj(1:4,1),phep(1:4,lept2))
      call pwhgfill(diag,Rj1l,dsig/binsize(diag)/2d0)
      endif

C     Delta R(j2 lept) 
      diag=diag+1
      if (mjets .ge. 2) then 
      Rj2l = r(pj(1:4,2),phep(1:4,lept1))
      call pwhgfill(diag,Rj2l,dsig/binsize(diag)/2d0)
      Rj2l = r(pj(1:4,2),phep(1:4,lept2))
      call pwhgfill(diag,Rj2l,dsig/binsize(diag)/2d0)
      endif

c     Phi(l1l2) 
      diag=diag+1 
      phill=fphi(phep(1,lept1),phep(1,lept2))
      call pwhgfill(diag,phill,dsig/binsize(diag))

c     pt of third jet
      diag=diag+1
      if(mjets.ge.3) call pwhgfill(diag,kt(3),dsig/binsize(diag))
c     y of jet 3
      diag=diag+1
      if(mjets.ge.3 .and. kt(3) > 30d0 ) 
     .     call pwhgfill(diag,rap(3),dsig/binsize(diag))
c     y of jet 3 asym
      diag=diag+1
      if(mjets.ge.3 .and. kt(3) > 30d0 ) 
     .call pwhgfill(diag,abs(rap(3)),sign(dsig/binsize(diag),rap(3)))
c     pt-rel of first jet
      diag=diag+1
      if(mjets.ge.1) call pwhgfill(diag,ptrel(1),dsig/binsize(diag))
c     pt-rel of second jet
      diag=diag+1
      if(mjets.ge.2) call pwhgfill(diag,ptrel(2),dsig/binsize(diag))
      else
         diag = diag+25
      endif

C     now same histagrams WITHOUT lepton cuts 

C      diag=0
c     two jet inclusive and two jet exclusive.
c     It is exclusive if the third jet is below the cut
c     inclusive otherwise
      diag=diag+1
      do j=1,3
         if(mjets.lt.j) kt(j)=0
      enddo
      do j=30,70,10
         if(kt(1).gt.j.and.kt(2).gt.j) then
c inclusive
            call pwhgfill(diag,dble(j),dsig)
            if(kt(3).lt.j) then
c exclusive
               call pwhgfill(diag+1,dble(j),dsig)
            endif
         endif
      enddo
      diag=diag+1

c     PT of jet 1
      diag=diag+1
      call pwhgfill(diag,kt(1),dsig/binsize(diag))
c     PT of jet 2
      diag=diag+1
      call pwhgfill(diag,kt(2),dsig/binsize(diag))
c     Eta of jet 1
      diag=diag+1
      call pwhgfill(diag,eta(1),dsig/binsize(diag))
c     Eta of jet 1 asym
      diag=diag+1
      call pwhgfill(diag,abs(eta(1)),sign(dsig/binsize(diag),eta(1)))
c     Eta of jet 2
      diag=diag+1
      call pwhgfill(diag,eta(2),dsig/binsize(diag))
c     Eta of jet 2 asym
      diag=diag+1
      call pwhgfill(diag,abs(eta(2)),sign(dsig/binsize(diag),eta(2)))
c     Eta jet1 - Eta jet 2
      diag=diag+1
      pj12 = pj(:4,1)+pj(:4,2)
      eta12 = etafromp(pj12)
      call pwhgfill(diag,eta12,dsig/binsize(diag))
c     H_T tot
      httot=ptel1+ptel2+ptmiss
      do j=1,mjets
         httot=httot+kt(j)
      enddo
      diag=diag+1
      call pwhgfill(diag,httot,dsig/binsize(diag))

c     H_T tot, 3j
      httot=ptel1+ptel2+ptmiss
      do j=1,min(mjets,3)
         httot=httot+kt(j)
      enddo
      diag=diag+1
      call pwhgfill(diag,httot,dsig/binsize(diag))

C     Pt lept     
      diag=diag+1
      call pwhgfill(diag,ptel1,dsig/binsize(diag)/2d0)
      call pwhgfill(diag,ptel2,dsig/binsize(diag)/2d0)

C     Pt miss
      diag=diag+1
      call pwhgfill(diag,ptmiss,dsig/binsize(diag))

C     Eta lept
      diag=diag+1
      call pwhgfill(diag,etael1,dsig/binsize(diag)/2d0)
      call pwhgfill(diag,etael2,dsig/binsize(diag)/2d0)

C     eta(l1)- eta(l2) 
      diag=diag+1
      pel12 = phep(:4,lept1)+phep(:4,lept2)
      etael12 = etafromp(pel12)
      call pwhgfill(diag,etael12,dsig/binsize(diag))

C     M(l1l2) 
      diag=diag+1
      mll = invmass(pel12)
      call pwhgfill(diag,mll,dsig/binsize(diag))

C     MT_WW
      diag=diag+1
      Etmisstilde = sqrt(ptmiss**2+mll**2)
      Etll=sqrt(ptfromp(pel12)**2
     1     + mll**2)
      mtww = sqrt((Etll+Etmisstilde)**2
     1     -ptfromp(pel12
     2     +phep(1:4,inue(1))+phep(1:4,inue(2)))**2)
      call pwhgfill(diag,mtww,dsig/binsize(diag))

C     Delta R(j1 lept) 
      diag=diag+1
      if (mjets .ge. 2) then 
      Rj1l = r(pj(1:4,1),phep(1:4,lept1))
      call pwhgfill(diag,Rj1l,dsig/binsize(diag)/2d0)
      Rj1l = r(pj(1:4,1),phep(1:4,lept2))
      call pwhgfill(diag,Rj1l,dsig/binsize(diag)/2d0)
      endif

C     Delta R(j2 lept) 
      diag=diag+1
      if (mjets .ge. 2) then 
      Rj2l = r(pj(1:4,2),phep(1:4,lept1))
      call pwhgfill(diag,Rj2l,dsig/binsize(diag)/2d0)
      Rj2l = r(pj(1:4,2),phep(1:4,lept2))
      call pwhgfill(diag,Rj2l,dsig/binsize(diag)/2d0)
      endif

c     Phi(l1l2) 
      diag=diag+1 
      phill=fphi(phep(1,lept1),phep(1,lept2))
      call pwhgfill(diag,phill,dsig/binsize(diag))

c     pt of third jet
      diag=diag+1
      if(mjets.ge.3) call pwhgfill(diag,kt(3),dsig/binsize(diag))
c     y of third jet
      diag=diag+1
      if(mjets.ge.3 .and. kt(3) > 30d0 ) 
     .     call pwhgfill(diag,rap(3),dsig/binsize(diag))
c     y of jet 3 asym
      diag=diag+1
      if(mjets.ge.3 .and. kt(3) > 30d0 ) 
     .  call pwhgfill(diag,abs(rap(3)),sign(dsig/binsize(diag),rap(3)))
c     pt-rel of first jet
      diag=diag+1
      if(mjets.ge.1) call pwhgfill(diag,ptrel(1),dsig/binsize(diag))
c     pt-rel of second jet
      diag=diag+1
      if(mjets.ge.2) call pwhgfill(diag,ptrel(2),dsig/binsize(diag))


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
     1     pjet(4,mjets),ptrel(mjets)
      include   'hepevt.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu
      real * 8 r,palg,ptmin,pp,tmp
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
C     R = 0.7   radius parameter
C     f = 0.75  overlapping fraction
      palg=-1
      r=0.4d0
      ptmin=30
      ptmin=0.1d0
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
