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

      call pwhginihist

c     total cross section sanity check
      diag=1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'total','LOG',binsize(diag),0d0,3d0)

      diag=diag+1
      binsize(diag) = 0.5d0
      call pwhgbookup(diag,'mt(W)','LOG',binsize(diag),50d0,100d0)

      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'m(W)','LOG',binsize(diag),20d0,100d0)

      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'pt(lep)','LOG',binsize(diag)
     $   ,25d0,55d0)


      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'Et_miss','LOG',binsize(diag),0d0
     $     ,100d0)

      diag=diag+1
      binsize(diag) = 0.1d0
      call pwhgbookup(diag,'Delta_phi','LOG'
     $     ,binsize(diag),0d0,3.2d0)

      diag=diag+1
      binsize(diag) = 0.1d0
      call pwhgbookup(diag,'eta(lep)','LOG'
     $     ,binsize(diag),-4d0,4d0)
  
      diag=diag+1
      binsize(diag) = 0.1d0
      call pwhgbookup(diag,'y(W)','LOG',binsize(diag),-3d0,3d0)

    
      diag=diag+1
      binsize(diag) = 0.25d0
      call pwhgbookup(diag,'pt(W) zoom','LOG',binsize(diag),0d0
     $     ,25d0)
 
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'pt(W)','LOG',binsize(diag),0d0,400d0)

      end 
c***********************************************************************
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h' 
      include 'pwhg_flg.h'  !WZGRAD EDIT
      include 'pwhg_wzgrad.h'  !WZGRAD EDIT
      integer ikini
      common/phspinfo/ikini
      real *8 p_neutrino(0:3),p_lepton(0:3),pcm(0:3),p_ll(0:3)
      real *8 dR,p1(4),p2(4),p3(4),sp1(4),sp2(4),sp3(4)!WZGRAD EDIT
      real *8 pt_lepton,pt_neutrino,eta_lepton,eta_neutrino,
     $     delphi,mt_v,mvv,ptv,yv
      integer flag1,flag2 !WZGRAD EDIT
      logical flg_insigreal
      common/caloflags/flg_insigreal
      common/leptonidf/flag1,flag2
      logical flg_saverand
      common/pwhg_flg_saverand/flg_saverand
      integer ihep,mu
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
      integer vdecaytemp,vdecay2temp,idvecbos
      save vdecaytemp,vdecay2temp,idvecbos
      integer maxnumlep
      parameter (maxnumlep=10)
      integer lepvec(maxnumlep),nuvec(maxnumlep)
      integer ilep,inu,lep,nu,nlep,nnu,jlep,jnu,i
      real *8 mV2ref,mV2
      real *8 Wmass,Wwidth,Wmass2low,Wmass2high
      logical foundlep
      real*8 powheginput,lcut,ncut,etlcut !WZGRAD EDIT
      external powheginput !WZGRAD EDIT

      if (ini) then!===================================================
         write (*,*)
         write (*,*) '********************************************'
         if(whcprg.eq.'NLO') then
            write (*,*) '           NLO ANALYSIS CALLED        '
         elseif(WHCPRG.eq.'LHE   ') then
            write (*,*) '           LHE ANALYSIS CALLED        '
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS CALLED     '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS CALLED     '
         endif
         write (*,*) '********************************************'
         write (*,*)
!        id of the charged decay product of the W
         vdecaytemp=lprup(1)-10000
         if(vdecaytemp.lt.0) then 
            vdecay2temp=-vdecaytemp+1 
            idvecbos=24
         elseif(vdecaytemp.gt.0) then
            vdecay2temp=-(vdecaytemp+1)
            idvecbos=-24
         else
            write(*,*) 'Error in decay mode in pwhg_analysis'
            call exit(1)
         endif
         if(abs(vdecaytemp).eq.11.or.abs(vdecaytemp).eq.13
     $        .or.abs(vdecaytemp).eq.15) then
            continue
         else
            write(*,*) '**************************************'
            write(*,*) ' Analysis works only for e, mu or tau decays'
            write(*,*) '                 STOP     '
            write(*,*) '**************************************'
            call exit(1)
         endif
         write (*,*)
         ini=.false.
      endif !==========================================================
      
      diag=0
      
      do mu=0,3
      pcm(mu)=0d0
      enddo

      if(whcprg.eq.'NLO    '.or.whcprg.eq.'LHE   ')then!================
c     find W decay products by trivial analysis
         do ihep=1,nhep
            if(isthep(ihep).eq.1) then
               if(idhep(ihep).eq.vdecay2temp) then
                  p_neutrino(0)=phep(4,ihep)
                  do mu=1,3
                     p_neutrino(mu)=phep(mu,ihep)
                  enddo
               elseif(idhep(ihep).eq.vdecaytemp) then
                  p_lepton(0)=phep(4,ihep)
                  do mu=1,3
                     p_lepton(mu)=phep(mu,ihep)
                  enddo
               endif
            endif
         enddo
      elseif ((WHCPRG.eq.'HERWIG').or.(WHCPRG.eq.'PYTHIA'))then        
c     find W decay products looking for their invariant mass
      Wmass = 80.398d0
      Wwidth = 2.141d0
      Wmass2low = (Wmass-30d0*Wwidth)**2
      Wmass2high = (Wmass+30d0*Wwidth)**2
      nlep=0
      nnu=0

         do i=1,maxnumlep
         lepvec(i) = 0
         nuvec(i) = 0
         enddo
         do ihep=1,nhep
            if(isthep(ihep).eq.1) then
c           Scan over final state particles and record the entries
               if(idhep(ihep).eq.vdecay2temp) then 
c              with a neutrino
               nnu=nnu+1
               nuvec(nnu)=ihep
               elseif(idhep(ihep).eq.vdecaytemp) then
c              with a lepton
               nlep=nlep+1
               lepvec(nlep)=ihep
               endif
            endif
         enddo

         if(nlep.eq.0.or.nnu.eq.0) then
         write(*,*)" not enough leptons! drop event"
         call exit(1)
         endif   

C     Now reconstruct the invariant masses from lepton-neutrino candidates
      foundlep = .false.
      mV2ref = 1d30
         do lep=1,nlep
            do nu=1,nnu
            ilep=lepvec(lep)
            inu=nuvec(nu)
            mV2 = (phep(4,ilep)+phep(4,inu))**2
     $          -(phep(1,ilep)+phep(1,inu))**2
     $          -(phep(2,ilep)+phep(2,inu))**2
     $          -(phep(3,ilep)+phep(3,inu))**2          
               if((Wmass2low.lt.mV2).and.(mV2.lt.Wmass2high))then
c                  if (foundlep) then
c                     write(*,*) 
c     $                    'two lepton couples satisfy W mass window'
c                     write(*,*) 'event dropped!'
c                     return
c                  endif
               foundlep=.true.
                  if(abs(mV2-Wmass**2).lt.abs(mV2ref-Wmass**2))then
                  mV2ref = mV2
                  jlep = ilep
                  jnu = inu               
                  endif
               endif
            enddo
         enddo

         if(.not.foundlep) then
         write(*,*)" No leptons in the mass window! drop event"
         return
         endif   
c     Assign lepton and neutrino momenta          
      p_neutrino(0)=phep(4,jnu)
      p_lepton(0)=phep(4,jlep)
         do mu=1,3
         p_neutrino(mu)=phep(mu,jnu)
         p_lepton(mu)=phep(mu,jlep)
         enddo
      else
      write(*,*) 'Not yet implemented analysis'
      call exit(1)
      endif 
!===========================================================
!     perform smearing + recomb - if user set calo = 1 
!     =====================================================
      if(powheginput('calo').eq.1)then 
c calculate smeared momenta (with possible recombination)
         if ((WHCPRG.eq.'HERWIG').or.(WHCPRG.eq.'PYTHIA').or.
     $        (WHCPRG.eq.'LHE   '))then 
            flg_saverand=.true.
            flg_insigreal=.false.
         endif
         call momentumprep_anl(p_lepton,p_neutrino,p1,p2,p3)
         call smear(p1,p2,p3,sp1,sp2,sp3)
c only allow for recombination when NLO histograms are filled
c and only for electrons
         if(whcprg.eq.'NLO    '.and.
     $        powheginput('vdecaymode').eq.1) then
         if(powheginput('wgrad2').eq.1.and.
     $        flg_insigreal.and.ikini.eq.1)then
            flag1=0
            call PStest(p_lepton)
c if the photon is hard and non-collinear
c recombination if dR<0.1 in electron case
            if(flag1.eq.0)call delR(sp1,sp2,sp3,dR)
         endif
         endif
         p_lepton(0)=sp1(4)
         p_lepton(1)=sp1(1)
         p_lepton(2)=sp1(2)
         p_lepton(3)=sp1(3)
         p_neutrino(0)=sp2(4)
         p_neutrino(1)=sp2(1)
         p_neutrino(2)=sp2(2)
         p_neutrino(3)=sp2(3)
      endif
!     =====================================================
!     calculate observables
!     =====================================================
c     neutrino transverse momentum
      pt_neutrino=sqrt(p_neutrino(1)**2 + p_neutrino(2)**2)
      call get_pseudorap(p_neutrino,eta_neutrino)
c     lepton transverse momentum
      pt_lepton=sqrt(p_lepton(1)**2 + p_lepton(2)**2)
      call get_pseudorap(p_lepton,eta_lepton)
c     invariant mass of the lepton-neutrino system
      do mu=0,3
         p_ll(mu)=p_lepton(mu)+p_neutrino(mu)
      enddo
      call getinvmass(p_ll,mvv)
c     azimuthal separation between lepton and neutrino
      delphi = dabs(atan2(p_lepton(2),p_lepton(1)) - 
     $        atan2(p_neutrino(2),p_neutrino(1)))
      delphi=min(delphi,2d0*pi-delphi)
c     transverse mass of the lepton-neutrino system
      mt_v=sqrt(2d0*pt_lepton*pt_neutrino*(1d0-dcos(delphi)))
c     rapidity of the lepton-neutrino system
      call getrapidity(p_ll,yv)

      ptv=sqrt((p_lepton(1)+p_neutrino(1))**2
     $     + (p_lepton(2)+p_neutrino(2))**2)
!     =====================================================
!     perform bare cuts - if user set bare =  1
!     =====================================================
      if(powheginput('bare').eq.1)then
         lcut=powheginput('cut1')
         ncut=powheginput('cut2')
         etlcut=powheginput('cut3')
         if(pt_lepton.lt.lcut.or.pt_neutrino.lt.ncut.or.
     &        dabs(eta_lepton).gt.etlcut)then
            goto 99
         endif
      endif
!     =====================================================
!     fill histograms
!     =====================================================
c     total sigma 
      diag=diag+1
      call pwhgfill(diag,1.5d0,dsig/binsize(diag))

c     transverse mass of the lepton-neutrino system
      diag=diag+1
      call pwhgfill(diag,mt_v,dsig/binsize(diag))

c     invariant mass of the lepton-neutrino system
      diag=diag+1
      call pwhgfill(diag,mvv,dsig/binsize(diag))

c     pt(l)
      diag=diag+1
      call pwhgfill(diag,pt_lepton,dsig/binsize(diag))

c     Et_miss
      diag=diag+1
      call pwhgfill(diag,pt_neutrino,dsig/binsize(diag))

c     azimuthal separation betwen lepton and neutrino
      diag=diag+1
      call pwhgfill(diag,delphi,dsig/binsize(diag))

c     eta(l)
      diag=diag+1
      call pwhgfill(diag,eta_lepton,dsig/binsize(diag))

c     y(W)
      diag=diag+1
      call pwhgfill(diag,yv,dsig/binsize(diag))

c     pt(W), zoom
      diag=diag+1
      call pwhgfill(diag,ptv,dsig/binsize(diag))

c     pt(W)
      diag=diag+1
      call pwhgfill(diag,ptv,dsig/binsize(diag))
  99  continue
      end
c***********************************************************************      
      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(0:3),y
      y=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end
c***********************************************************************      
      subroutine getinvmass(p,m)
      implicit none
      real * 8 p(0:3),m
      m=sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end
c***********************************************************************      
      subroutine get_pseudorap(p,eta)
      implicit none
      real*8 p(0:3),eta,pt,th
      real *8 tiny
      parameter (tiny=1.d-5)

      pt=sqrt(p(1)**2+p(2)**2)
      if(pt.lt.tiny.and.abs(p(3)).lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      elseif(pt.lt.tiny) then   !: added this elseif
         eta=sign(1.d0,p(3))*1.d8
      else
         th=atan2(pt,p(3))
         eta=-log(tan(th/2.d0))
      endif
      
      end
c***********************************************************************      
      subroutine PStest(p_lepton)
      implicit none
      
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_wzgrad.h'
      real*8 dotpr,sinv(3,5),p_lepton(0:3)
      external dotpr

      integer flag1,flag2
      common/leptonidf/flag1,flag2

      real*8 epmin,ephoton,smin
!     =====================================================
!     reset flags - they were already used, don't want to
!     carry over their old values
!     =====================================================
      flag1=0
      flag2=0
!     =====================================================
!     define invariants ito the new momenta
!     =====================================================
      sinv(1,2)=2d0*dotpr(kn_preal(0,1),kn_preal(0,2))
      sinv(1,5)=2d0*dotpr(kn_preal(0,1),kn_preal(0,5))
      sinv(2,5)=2d0*dotpr(kn_preal(0,2),kn_preal(0,5))
      sinv(3,5)=2d0*dotpr(p_lepton(0),kn_preal(0,5))

      epmin=deltas*dsqrt(sinv(1,2))/2d0
      if(test(3).eq.2.and.test(4).eq.2)epmin=deltas*mw/2d0
      ephoton=kn_cmpreal(0,5) !photon energy in COM frame!
      smin=deltac*dsqrt(sinv(1,2))*ephoton
!     cut 1 ==> is photon soft?
      if(ephoton.lt.epmin)then
      flag1=1
      endif
!     cut2 ==> quark-photon angle too small?
      if(dabs(sinv(1,5)).lt.smin.or.dabs(sinv(2,5)).lt.smin)then
      flag1=1 
      endif
!     cut 3 ==> lepton-photon angle too small?
      if(collcut.eq.1)then
      smin=smin*(1d0-2d0*ephoton/dsqrt(sinv(1,2)))
         if(dabs(sinv(3,5)).lt.smin)then
         flag1=1 
         endif
      endif
      end
c***********************************************************************
      subroutine momentumprep_anl(p_lepton,p_neutrino,p1,p2,p3)
c.....
c     this subroutine prepares the momenta for smearing - only for analysis
c     input: p_lepton(0:3),p_neutrino(0:3),kn_preal(0:3,5) 
c     output: p1,p2,p3 
c.....
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      real*8 p_lepton(0:3),p_neutrino(0:3),p1(4),p2(4),p3(4)
      logical flg_insigreal
      common/caloflags/flg_insigreal
c     lepton
      p1(1) = p_lepton(1)
      p1(2) = p_lepton(2)
      p1(3) = p_lepton(3)
      p1(4) = p_lepton(0)
c     neutrino
      p2(1) = p_neutrino(1)
      p2(2) = p_neutrino(2)
      p2(3) = p_neutrino(3)
      p2(4) = p_neutrino(0)
c     photon
         if(flg_insigreal)then
         p3(1) = kn_preal(1,5)
         p3(2) = kn_preal(2,5)
         p3(3) = kn_preal(3,5)
         p3(4) = kn_preal(0,5)
         else
         p3(1) = 0d0
         p3(2) = 0d0
         p3(3) = 0d0
         p3(4) = 0d0
         endif
      end 
c***********************************************************************

