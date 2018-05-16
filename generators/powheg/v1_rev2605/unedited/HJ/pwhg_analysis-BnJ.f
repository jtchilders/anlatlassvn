c Analysis used for arXiv:1212.4504.

c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include     'LesHouches.h'
      include     'pwhg_math.h'
      include     'pwhg_bookhist-new.h'
      integer      maxjet
      parameter   (maxjet=2)
      integer      nptmin
      parameter   (nptmin=4)
      character*1  cnums(9)
      data         cnums/'1','2','3','4','5','6','7','8','9'/
      character*4  cptcuts(nptmin)
      data         cptcuts/'-001','-020','-040','-080'/
      real*8       ptcuts(nptmin)
      data         ptcuts/1d0,20d0,40d0,80d0/
      character*4  cycuts(5)
      data         cycuts/'-002','-004','-008','-012','-020'/
      real*8       ycuts(5)
      data         ycuts/0.2d0,0.4d0,0.8d0,1.2d0,2.0d0/
      character*4  calgs(3)
      data         calgs/'--kt','--ca','-akt'/
      real*8       jetalgs(3)
      data         jetalgs/1d0,0d0,-1d0/
      common/infohist/ptcuts,ycuts,jetalgs,cnums,cptcuts,cycuts,calgs
      save  /infohist/
      character*1  pr
      common/pwhgprocess/pr
      integer      ixx,jxx,kxx
      real*8       dy,dpt,dr,dptzm,dptzm_alt
      character*4  cptcut,cycut,calg
      character*1  ctmp1
      character*8  suffix
      character*12 suffix12

      call inihists

C - REMEMBER to set process: H/W/Z
      pr='H'

C - Set bin widths
      dy=0.25d0
      dpt=10d0
      dr=0.1d0
      dptzm=2d0
      dptzm_alt=1d0

C     - Total cross-section 
      call bookupeqbins('sigmatot', 1d0,0d0,1d0)

C - Inclusive boson:
      call bookupeqbins('B-y-inc',          dy,-5d0,5d0)
      call bookupeqbins('B-eta-inc',        dy,-5d0,5d0)
      call bookupeqbins('B-pt-inc',        dpt,0d0,400d0)
      call bookupeqbins('B-ptzm-inc',dptzm_alt,0d0,100d0)
      call bookupeqbins('B-m-inc',       dptzm,0d0,400d0)

C - Inclusive boson pT spectrum again but in bins of rapidity:
      call bookupeqbins('B-ptzm-y-0-1',  dptzm_alt,0d0,100d0)
      call bookupeqbins('B-ptzm-y-1-2',  dptzm_alt,0d0,100d0)
      call bookupeqbins('B-ptzm-y-2-3',  dptzm_alt,0d0,100d0)
      call bookupeqbins('B-ptzm-y-3-inf',dptzm_alt,0d0,100d0)

C - Inclusive lepton
      call bookupeqbins('l-y-inc',          dy,-5d0,5d0)
      call bookupeqbins('l-eta-inc',        dy,-5d0,5d0)
      call bookupeqbins('l-pt-inc',        dpt,0d0,400d0)
      call bookupeqbins('l-ptzm-inc',dptzm_alt,0d0,100d0)

C - Inclusive anti-lepton
      call bookupeqbins('al-y-inc',          dy,-5d0,5d0)
      call bookupeqbins('al-eta-inc',        dy,-5d0,5d0)
      call bookupeqbins('al-pt-inc',        dpt,0d0,400d0)
      call bookupeqbins('al-ptzm-inc',dptzm_alt,0d0,100d0)

C - Loop over jet algorithms:
      do ixx=1,3
        calg=calgs(ixx)

C - Differential jet rates - only the kT algorithm as anti kT gives
C - very different results changing from NLO -> LHE -> shower.
        if(ixx.eq.1) then
           call bookupeqbins('y01'//calg,0.030d0,0.3d0,3.0d0)
           call bookupeqbins('y12'//calg,0.030d0,0.3d0,3.0d0)
           call bookupeqbins('y23'//calg,0.025d0,0.3d0,2.5d0)
           call bookupeqbins('y34'//calg,0.020d0,0.3d0,2.0d0)

           call bookupeqbins('y12-y01gt50sq'//calg,0.030d0,0.3d0,3.0d0)
           call bookupeqbins('y12-y01gt100sq'//calg,0.030d0,0.3d0,3.0d0)
        endif

C - Jet pT spectra:
        do kxx=1,3
           ctmp1=cnums(kxx)
C - "Inclusively"
           call bookupeqbins('j'//ctmp1//'-pt'//calg,dpt,0d0,400d0)
C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
           do jxx=2,nptmin
              cptcut=cptcuts(jxx)
              suffix=cptcut//calg
C - in >=1 jet events
              call bookupeqbins('Nj>=1-j'//ctmp1//'-pt'//suffix,
     $                           dpt,0d0,400d0)
C - in >=2 jet events
              call bookupeqbins('Nj>=2-j'//ctmp1//'-pt'//suffix,
     $                           dpt,0d0,400d0)
           enddo
C - End loop over jet pT cuts
        enddo

C - Zoomed jet pT spectra:
        do kxx=1,3
           ctmp1=cnums(kxx)
C - "Inclusively"
           call bookupeqbins('j'//ctmp1//'-ptzm'//calg,
     $                       dptzm_alt,0d0,100d0)
C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
           do jxx=2,nptmin
              cptcut=cptcuts(jxx)
              suffix=cptcut//calg
C - in >=1 jet events
              call bookupeqbins('Nj>=1-j'//ctmp1//'-ptzm'//suffix,
     $                           dpt,0d0,400d0)
C - in >=2 jet events
              call bookupeqbins('Nj>=2-j'//ctmp1//'-ptzm'//suffix,
     $                           dpt,0d0,400d0)
           enddo
C - End loop over jet pT cuts
        enddo

C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
        do jxx=2,2
           cptcut=cptcuts(jxx)
           suffix=cptcut//calg
C - Inclusive jet multiplicties
           call bookupeqbins('Njet-inc'//suffix,1d0,-0.5d0,5.5d0)
C - Exclusive jet multiplicties
           call bookupeqbins('Njet-exc'//suffix,1d0,-0.5d0,5.5d0)

C - Boson+>=1-jet inclusive obsverables:
           call bookupeqbins('B-y-1j-inc'//suffix,dy,-5d0,5d0)
           call bookupeqbins('B-pt-1j-inc'//suffix,dpt,0d0,400d0)
           call bookupeqbins('B-ptzm-1j-inc'//suffix,
     $                        dptzm_alt,0d0,100d0)

C - Boson+>=2-jet inclusive obsverables:
           call bookupeqbins('B-y-2j-inc'//suffix,dy,-5d0,5d0)
           call bookupeqbins('B-pt-2j-inc'//suffix,dpt,0d0,400d0)
           call bookupeqbins('B-ptzm-2j-inc'//suffix,
     $                        dptzm_alt,0d0,100d0)

C - Boson+0-jet exclusive obsverables:
           call bookupeqbins('B-y-0j-exc'//suffix,dy,-5d0,5d0)
           call bookupeqbins('B-pt-0j-exc'//suffix,dpt,0d0,400d0)
           call bookupeqbins('B-ptzm-0j-exc'//suffix,
     $                        dptzm_alt,0d0,100d0)

C - Boson+1-jet exclusive obsverables:
           call bookupeqbins('B-y-1j-exc'//suffix,dy,-5d0,5d0)
           call bookupeqbins('B-pt-1j-exc'//suffix,dpt,0d0,400d0)
           call bookupeqbins('B-ptzm-1j-exc'//suffix,
     $                        dptzm_alt,0d0,100d0)

C - Boson+2-jet exclusive obsverables:
           call bookupeqbins('B-y-2j-exc'//suffix,dy,-5d0,5d0)
           call bookupeqbins('B-pt-2j-exc'//suffix,dpt,0d0,400d0)
           call bookupeqbins('B-ptzm-2j-exc'//suffix,
     $                        dptzm_alt,0d0,100d0)

C - Delta-y, -eta, -phi and -R of the boson w.r.t. 
C - the 1st-jet in events with at least 1 jet
           call bookupeqbins('Nj>=1-Bj1-dy'//suffix,dy,-5d0,5d0)
           call bookupeqbins('Nj>=1-Bj1-deta'//suffix,dy,-5d0,5d0)
           call bookupeqbins('Nj>=1-Bj1-delphi'//suffix,pi/50,0d0,pi)
           call bookupeqbins('Nj>=1-Bj1-dr'//suffix,dr,0d0,8d0)  

C - y and m for jets 1 and 2 separately:
           do kxx=1,2
             ctmp1=cnums(kxx)
             call bookupeqbins('j'//ctmp1//'-y'//suffix,dy,-5d0,5d0)
             call bookupeqbins('j'//ctmp1//'-m'//suffix,dptzm,1d0,101d0)
           enddo

C - y and m for jet 1 in >=2 jet events:
           do kxx=1,1
             ctmp1=cnums(kxx)
             call bookupeqbins('Nj>=2-j'//ctmp1//'-y'//suffix,
     $                          dy,-5d0,5d0)
             call bookupeqbins('Nj>=2-j'//ctmp1//'-m'//suffix,
     $                          dptzm,1d0,101d0)
           enddo

C - Delta-y, -eta, -phi and -R of the 1st jet w.r.t. the 2nd-jet
           call bookupeqbins('j1j2-dy'//suffix,dy,-5d0,5d0)
           call bookupeqbins('j1j2-deta'//suffix,dy,-5d0,5d0)
           call bookupeqbins('j1j2-delphi'//suffix,pi/50,0d0,pi)
           call bookupeqbins('j1j2-dr'//suffix,dr,0d0,8d0)

C - End loop over jet pT cuts:
        enddo

C - Booking GPS's evil pT spectra ~ the volcano region:
        do jxx=1,5
           suffix=cycuts(jxx)//calg
           call bookupeqbins('j1-ptzm-dyj1'//suffix,1d0,1d0,101d0)
        enddo
        do jxx=1,5
           suffix=cycuts(jxx)//calg
           call bookupeqbins('j1-ptzm-dyBj1'//suffix,1d0,1d0,101d0)
        enddo



C - End loop over jet algorithms:
      enddo

      call bookupeqbins('vbf-pt3',5d0,0d0,100d0)
      call bookupeqbins('vbf-pt3-veto',5d0,0d0,100d0)
      call bookupeqbins('vbf-y3-020',0.1d0,-5d0,5d0)
      call bookupeqbins('vbf-y3-030',0.1d0,-5d0,5d0)
      call bookupeqbins('vbf-y3-040',0.1d0,-5d0,5d0)

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C
C - Relative differences in pT spectra between jet algorithms - C
C - (done in pwhgfinalopshist):                               - C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C

C - Loop over CA and AkT jet algorithms:
      do ixx=2,3
        suffix=calgs(ixx)//calgs(1)

C - Jet pT spectra:
        do kxx=1,3
          ctmp1=cnums(kxx)
C - "Inclusively"
          call bookupeqbins('j'//ctmp1//'-pt'//suffix,dpt,0d0,400d0)
C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
          do jxx=2,nptmin
             cptcut=cptcuts(jxx)
             suffix12=cptcut//suffix
C - in >=1 jet events
             call bookupeqbins('Nj>=1-j'//ctmp1//'-pt'//suffix12,
     $                          dpt,0d0,400d0)
C - in >=2 jet events
             call bookupeqbins('Nj>=2-j'//ctmp1//'-pt'//suffix12,
     $                          dpt,0d0,400d0)
          enddo
C - End loop over jet pT cuts
        enddo

C - Zoomed jet pT spectra:
        do kxx=1,3
          ctmp1=cnums(kxx)
C - "Inclusively"
          call bookupeqbins('j'//ctmp1//'-ptzm'//suffix,
     $                      dptzm_alt,0d0,100d0)
C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
          do jxx=2,nptmin
             cptcut=cptcuts(jxx)
             suffix12=cptcut//suffix
C - in >=1 jet events
             call bookupeqbins('Nj>=1-j'//ctmp1//'-ptzm'//suffix12,
     $                          dpt,0d0,400d0)
C - in >=2 jet events
             call bookupeqbins('Nj>=2-j'//ctmp1//'-ptzm'//suffix12,
     $                          dpt,0d0,400d0)
          enddo
C - End loop over jet pT cuts
        enddo

C - End loop over CA and AkT algorithms
      enddo

      write(*,*) 'booked ',jhist,' histograms !'
      end
     
      subroutine analysis(dsig0)
      implicit  none
      real * 8     dsig(7),dsig0
      include     'hepevt.h'
      include     'nlegborn.h'
      include     'pwhg_flst.h'
      include     'pwhg_kn.h'
      include     'pwhg_math.h' 
      include     'pwhg_rad.h' 
      include     'pwhg_flg.h'
      include     'LesHouches.h'
      include     'pwhg_weights.h'
      integer      maxjet
      parameter   (maxjet=2048)
      integer      nptmin
      parameter   (nptmin=4)
      character*1  cnums(9)
      character*4  cptcuts(nptmin)
      real*8       ptcuts(nptmin)
      character*4  cycuts(5)
      real*8       ycuts(5)
      character*4  calgs(3)
      real*8       jetalgs(3)
      common/infohist/ptcuts,ycuts,jetalgs,cnums,cptcuts,cycuts,calgs
      save  /infohist/
      character*1  pr
      common/pwhgprocess/pr
      integer      ixx,jxx,kxx,lxx
      character*4  cptcut,cycut,calg
      character*1  ctmp1
      character*8  suffix
      character*12 suffix12
C - We need to tell to this analysis file which program is running it
      character*6 WHCPRG
      common/cWHCPRG/WHCPRG
      data    WHCPRG/'NLO   '/
      integer      mjets,njets
      real*8       ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1             phij(maxjet),pj(4,maxjet),jetRadius,ptrel(4)
      real*8       pCandidate(4)
      real*8       yijs(4)
      real*8       pB(4),pal(4),pl(4)
      real*8       y,eta,pt,m
      integer      ihep
      real*8       powheginput
      external     powheginput
      real*8       yB,etaB,ptB,mB,ymin,ymax
      logical vbfcuts 
      integer      iMaxPt
      integer      nBosons

      logical ini
      data ini/.true./
      save ini
      if(ini) then
         if(WHCPRG.eq.'NLO'.or.WHCPRG.eq.'LHE') then
            weights_num=0
         endif
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         ini=.false.
      endif

      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
      endif

      if(sum(abs(dsig)).eq.0) return

      nBosons=0



C - Find the boson:
      nBosons = 0
      if(pr.eq.'H') then
C     For HERWIG or PYTHIA, the last H boson is kept,
C     which is likely the Higgs before decaying, if it decays,
c     or the final Higgs, if undecayed
         do ihep=1,nhep
            if(idhep(ihep).eq.25) then
               pB = phep(1:4,ihep)
               nBosons = 1
            endif
         enddo
         if(nBosons.eq.0) then
            write(*,*) 'pwhg_analysis-BnJ.f: Fatal error.'
            write(*,*) ' Higgs boson not found'
            call exit(-1)
         endif
      else if(pr.eq.'Z') then
         if(WHCPRG.eq.'NLO   ') then
            pB =phep(1:4,3)+phep(1:4,4)
            pl =phep(1:4,3)
            pal=phep(1:4,4)
         elseif(WHCPRG.eq.'LHE   ') then
            pB =phep(1:4,4)+phep(1:4,5)
            pl =phep(1:4,4)
            pal=phep(1:4,5)
         else
           do ihep=1,nhep
             if(idhep(ihep).eq.23) then
                pCandidate=phep(1:4,ihep)
             endif
           enddo
           do ihep=1,nhep
             if(isthep(ihep).eq.1.and.
     1          abs(idhep(ihep)).ge.11.and.abs(idhep(ihep)).le.16
     2         ) then
               if(idhep(jmohep(1,ihep)).eq.23) then
                 pB=phep(1:4,jmohep(1,ihep))
                 if(idhep(ihep).gt.0) then
                    pl =phep(1:4,ihep)
                 else
                    pal=phep(1:4,ihep)
                 endif
               elseif(idhep(jmohep(1,jmohep(1,ihep))).eq.23) then
                 pB=phep(1:4,jmohep(1,jmohep(1,ihep)))
                 if(idhep(ihep).gt.0) then
                    pl =phep(1:4,ihep)
                 else
                    pal=phep(1:4,ihep)
                 endif
               endif
             endif
           enddo
           nBosons = 0
           if(abs(
     $        ((pCandidate(4)-pCandidate(3))
     $        *(pCandidate(4)+pCandidate(3))
     $        -pCandidate(1)**2-pCandidate(2)**2)
     $       -((pB(4)-pB(3))*(pB(4)+pB(3))
     $        -pB(1)**2-pB(2)**2)).lt.0.5d0) then
              nBosons = 1
           endif
           if(nBosons.eq.0) then
              write(*,*) 'pwhg_analysis-BnJ.f: Fatal error.'
              write(*,*) 'Z boson not found'
              call exit(-1)
           endif
         endif
      else if(pr.eq.'W') then
         if(WHCPRG.eq.'NLO   ') then
            pB =phep(1:4,3)+phep(1:4,4)
            pl =phep(1:4,3)
            pal=phep(1:4,4)
         elseif(WHCPRG.eq.'LHE   ') then
            pB =phep(1:4,4)+phep(1:4,5)
            pl =phep(1:4,4)
            pal=phep(1:4,5)
         else
           do ihep=1,nhep
             if(abs(idhep(ihep)).eq.24) then
                pCandidate=phep(1:4,ihep)
             endif
           enddo
           do ihep=1,nhep
             if(isthep(ihep).eq.1.and.
     1          abs(idhep(ihep)).ge.11.and.abs(idhep(ihep)).le.16
     2         ) then
               if(abs(idhep(jmohep(1,ihep))).eq.24) then
                 pB=phep(1:4,jmohep(1,ihep))
                 if(idhep(ihep).gt.0) then
                    pl =phep(1:4,ihep)
                 else
                    pal=phep(1:4,ihep)
                 endif
               elseif(abs(idhep(jmohep(1,jmohep(1,ihep)))).eq.24) then
                 pB=phep(1:4,jmohep(1,jmohep(1,ihep)))
                 if(idhep(ihep).gt.0) then
                    pl =phep(1:4,ihep)
                 else
                    pal=phep(1:4,ihep)
                 endif
               endif
             endif
           enddo
           nBosons = 0
           if(abs(
     $        ((pCandidate(4)-pCandidate(3))
     $        *(pCandidate(4)+pCandidate(3))
     $        -pCandidate(1)**2-pCandidate(2)**2)
     $       -((pB(4)-pB(3))*(pB(4)+pB(3))
     $        -pB(1)**2-pB(2)**2)).lt.0.5d0) then
              nBosons = 1
           endif
           if(nBosons.eq.0) then
              write(*,*) 'pwhg_analysis-BnJ.f: Fatal error.'
              write(*,*) 'W boson not found'
              call exit(-1)
           endif
         endif
      endif

C - Boson kinematics fully inclusive w.r.t radiation:
      call getyetaptmass(pB,y,eta,pt,m)

C     - Fill total cross-section 
      call filld('sigmatot',0.5d0,dsig)

C - Inclusive boson:
      call filld('B-y-inc',    y,dsig)
      call filld('B-eta-inc',eta,dsig)
      call filld('B-pt-inc',  pt,dsig)
      call filld('B-ptzm-inc',pt,dsig)
      call filld('B-m-inc',    m,dsig)

C - Inclusive boson pT spectrum again but in bins of rapidity:
      if(abs(y).lt.1d0) then
         call filld('B-ptzm-y-0-1',  pt,dsig)
      elseif(abs(y).lt.2d0) then
         call filld('B-ptzm-y-1-2',  pt,dsig)
      elseif(abs(y).lt.3d0) then
         call filld('B-ptzm-y-2-3',  pt,dsig)
      else
         call filld('B-ptzm-y-3-inf',pt,dsig)
      endif

C - Lepton kinematics fully inclusive w.r.t radiation:
      call getyetaptmass(pl,y,eta,pt,m)

C - Inclusive lepton:
      call filld('l-y-inc',    y,dsig)
      call filld('l-eta-inc',eta,dsig)
      call filld('l-pt-inc',  pt,dsig)
      call filld('l-ptzm-inc',pt,dsig)

C - Anti-Lepton kinematics fully inclusive w.r.t radiation:
      call getyetaptmass(pal,y,eta,pt,m)

C - Inclusive anti-lepton:
      call filld('al-y-inc',    y,dsig)
      call filld('al-eta-inc',eta,dsig)
      call filld('al-pt-inc',  pt,dsig)
      call filld('al-ptzm-inc',pt,dsig)

C - Re-do Boson kinematics fully inclusive w.r.t radiation (just in case):
      call getyetaptmass(pB,y,eta,pt,m)

C - Loop over jet algorithms:
      do ixx=1,3
        calg=calgs(ixx)

C - Erase any jet building output leftover from previous jet building:
        mjets=0
        njets=0
        do jxx=1,maxjet
           ktj(jxx)  = 0d0
           etaj(jxx) = 0d0
           rapj(jxx) = 0d0
           phij(jxx) = 0d0
           do kxx=1,4
              pj(kxx,jxx) = 0d0
           enddo
        enddo
        do jxx=1,4
           ptrel(jxx) = 0d0
        enddo
        do jxx=1,4
           yijs(jxx)  = 0d0
        enddo

C - Build jets:
        jetRadius=0.5d0
        call buildjets(1,jetalgs(ixx),jetRadius,1d0,mjets,
     $                 ktj,etaj,rapj,phij,ptrel,pj,yijs)



        mjets=min(mjets,3)

c VBF plots
        if(ixx.eq.3.and.mjets.ge.2) then
           vbfcuts = rapj(1)*rapj(2).lt.0
           vbfcuts = vbfcuts .and. ktj(1).gt.20.and.ktj(2).gt.20
           vbfcuts = vbfcuts .and. abs(rapj(1)).lt.4.5
     1          .and.abs(rapj(2)).lt.4.5
           vbfcuts = vbfcuts .and.
     1         sqrt( (pj(4,1)+pj(4,2))**2- (pj(1,1)+pj(1,2))**2
     2          - (pj(2,1)+pj(2,2))**2- (pj(3,1)+pj(3,2))**2).gt.600d0
           vbfcuts = vbfcuts .and. abs(rapj(1)-rapj(2)) .gt. 4d0 
           if(vbfcuts.and.mjets.ge.3) then
              call filld('vbf-pt3',ktj(3),dsig)
              if (ktj(3).gt.20d0) call filld('vbf-y3-020',rapj(3),dsig)
              if (ktj(3).gt.30d0) call filld('vbf-y3-030',rapj(3),dsig)
              if (ktj(3).gt.40d0) call filld('vbf-y3-040',rapj(3),dsig)
              ymin = min(rapj(1),rapj(2))
              ymax = max(rapj(1),rapj(2))
              if (ymin .lt. rapj(3) .and. rapj(3) .lt. ymax) 
     1            call filld('vbf-pt3-veto',ktj(3),dsig)
           endif
        endif


           
C - Differential jet rates.
        if(ixx.eq.1) then
          if(yijs(1).gt.0)
     $      call filld('y01'//calg,log10(sqrt(yijs(1))),dsig)
          if(yijs(2).gt.0) then 
             call filld('y12'//calg,log10(sqrt(yijs(2))),dsig)
             if (yijs(1) .gt. 50d0**2) 
     $     call filld('y12-y01gt50sq'//calg,log10(sqrt(yijs(2))),dsig)
             if (yijs(1) .gt. 100d0**2) 
     $     call filld('y12-y01gt100sq'//calg,log10(sqrt(yijs(2))),dsig)
          endif
          if(yijs(3).gt.0)
     $      call filld('y23'//calg,log10(sqrt(yijs(3))),dsig)
          if(yijs(4).gt.0)
     $      call filld('y34'//calg,log10(sqrt(yijs(4))),dsig)
        endif

C - Jet pT spectra:
        do kxx=1,mjets
           ctmp1=cnums(kxx)
           call getyetaptmass(pj(:,kxx),y,eta,pt,m)
C - "Inclusively"
           call filld('j'//ctmp1//'-pt'//calg,pt,dsig)
           call filld('j'//ctmp1//'-ptzm'//calg,pt,dsig)
C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
           do jxx=2,nptmin
              cptcut=cptcuts(jxx)
              suffix=cptcut//calg
C - Count number of jets above cut threshold (njets):
              njets=0
              do lxx=1,mjets
                 if (ktj(lxx).gt.ptcuts(jxx)) njets=njets+1
              enddo
C     - in >=1 jet events
              if(njets.ge.1) then 
                 call filld('Nj>=1-j'//ctmp1//'-pt'//suffix,pt,dsig)
                 call filld('Nj>=1-j'//ctmp1//'-ptzm'//suffix,pt,dsig)
              endif
C - in >=2 jet events
              if(njets.ge.2) then  
                 call filld('Nj>=2-j'//ctmp1//'-pt'//suffix,pt,dsig)
                 call filld('Nj>=2-j'//ctmp1//'-ptzm'//suffix,pt,dsig)
              endif
           enddo
C - End loop over jet pT cuts
        enddo

C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
        do jxx=2,2
          cptcut=cptcuts(jxx)
          suffix=cptcut//calg
C - Count number of jets above cut threshold (njets):
          njets=0
          do lxx=1,mjets
            if (ktj(lxx).gt.ptcuts(jxx)) njets=njets+1
          enddo
C - Plot inclusive jet multiplicties at current jet pT threshold:
          if(njets.ge.0) call filld('Njet-inc'//suffix,0d0,dsig)
          if(njets.ge.1) call filld('Njet-inc'//suffix,1d0,dsig)
          if(njets.ge.2) call filld('Njet-inc'//suffix,2d0,dsig)
          if(njets.ge.3) call filld('Njet-inc'//suffix,3d0,dsig)
          if(njets.ge.4) call filld('Njet-inc'//suffix,4d0,dsig)
          if(njets.ge.5) call filld('Njet-inc'//suffix,5d0,dsig)
C - Plot exclusive jet multiplicties at current jet pT threshold:
          call filld('Njet-exc'//suffix,njets*1d0,dsig)

C - Compute y, eta, pt, and m for the boson.
          call getyetaptmass(pB,y,eta,pt,m)

C - Boson+>=1-jet inclusive obsverables:
          if(njets.ge.1) then
            call filld('B-y-1j-inc'//suffix,    y,dsig)
            call filld('B-pt-1j-inc'//suffix,  pt,dsig)
            call filld('B-ptzm-1j-inc'//suffix,pt,dsig)
          endif

C - Boson+>=2-jet inclusive obsverables:
          if(njets.ge.2) then
            call filld('B-y-2j-inc'//suffix,    y,dsig)
            call filld('B-pt-2j-inc'//suffix,  pt,dsig)
            call filld('B-ptzm-2j-inc'//suffix,pt,dsig)
          endif

C - Boson+0-jet exclusive obsverables:
          if(njets.eq.0) then
            call filld('B-y-0j-exc'//suffix,    y,dsig)
            call filld('B-pt-0j-exc'//suffix,  pt,dsig)
            call filld('B-ptzm-0j-exc'//suffix,pt,dsig)
          endif

C - Boson+1-jet exclusive obsverables:
          if(njets.eq.1) then
            call filld('B-y-1j-exc'//suffix,    y,dsig)
            call filld('B-pt-1j-exc'//suffix,  pt,dsig)
            call filld('B-ptzm-1j-exc'//suffix,pt,dsig)
          endif

C - Boson+2-jet exclusive obsverables:
          if(njets.eq.2) then
            call filld('B-y-2j-exc'//suffix,    y,dsig)
            call filld('B-pt-2j-exc'//suffix,  pt,dsig)
            call filld('B-ptzm-2j-exc'//suffix,pt,dsig)
          endif
         
C - Delta-y, -eta, -phi and -R of the boson w.r.t.
C - the 1st-jet in events with at least 1 jet
          if(njets.ge.1) then
            call deltaplot(pB,pj(:,1),dsig,'Nj>=1-Bj1',suffix)
          endif

C - y and m for jets 1 and 2 separately:
          do kxx=1,min(2,njets)
             ctmp1=cnums(kxx)
             call getyetaptmass(pj(:,kxx),y,eta,pt,m)
             call filld('j'//ctmp1//'-y'//suffix,y,dsig)
             call filld('j'//ctmp1//'-m'//suffix,m,dsig)
          enddo

C - y and m for jet 1 in >=2 jet events:
          if(njets.ge.2) then
             do kxx=1,1
                ctmp1=cnums(kxx)
                call getyetaptmass(pj(:,kxx),y,eta,pt,m)
                call filld('Nj>=2-j'//ctmp1//'-y'//suffix,y,dsig)
                call filld('Nj>=2-j'//ctmp1//'-m'//suffix,m,dsig)
             enddo
          endif

C - Delta-y, -eta, -phi and -R of the 1st jet w.r.t. the 2nd-jet
          if(njets.ge.2) then
            call deltaplot(pj(:,1),pj(:,2),dsig,'j1j2',suffix)
          endif

C - End loop over jet pT cuts:
        enddo

C - Compute y, eta, pt, and m for the boson.
        call getyetaptmass(pB,yB,etaB,ptB,mB)

C - GPS's evil pT spectra ~ the volcano region:
        do jxx=1,5 ! <--- loop over jet delta y windows
           suffix=cycuts(jxx)//calg

           iMaxPt=-1
           do kxx=1,mjets
              call getyetaptmass(pj(:,kxx),y,eta,pt,m)
              if(ktj(kxx).gt.ptcuts(1).and.abs(y).le.ycuts(jxx)) then
                 if(iMaxPt.eq.-1) iMaxPt=kxx
                 if(ktj(kxx).gt.ktj(iMaxPt)) iMaxPt=kxx
              endif
           enddo
           if(iMaxPt.gt.0) then
              call filld('j1-ptzm-dyj1'//suffix,ktj(iMaxPt),dsig)
           endif
           
           iMaxPt=-1
           do kxx=1,mjets
              call getyetaptmass(pj(:,kxx),y,eta,pt,m)
              if(ktj(kxx).gt.ptcuts(1).and.abs(yB-y).le.ycuts(jxx)) then
                 if(iMaxPt.eq.-1) iMaxPt=kxx
                 if(ktj(kxx).gt.ktj(iMaxPt)) iMaxPt=kxx
              endif
           enddo
           if(iMaxPt.gt.0) then
              call filld('j1-ptzm-dyBj1'//suffix,ktj(iMaxPt),dsig)
           endif
           
        enddo

C - End loop over jet algorithms:
      enddo

      end


      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig(*)
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-delphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end


      subroutine buildjets(iflag,palg,rr,ptmin,mjets,kt,eta,rap,phi,
     $     ptrel,pjet,yijs)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
      integer   iflag,mjets
      real * 8  rr,ptmin,kt(*),eta(*),rap(*),
     1          phi(*),ptrel(3),pjet(4,*), pB(4) 
      include  'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      real * 8 ptsq_tot, ptot(4) 
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,jb,i
      real * 8  r,palg,pp,tmp
      logical   islept
      external  islept
      real * 8  vec(3),pjetin(0:3),pjetout(0:3),beta,
     $          ptrackin(0:3),ptrackout(0:3)
      real * 8  get_ptrel
      external  get_ptrel
      real * 8  yijs(4)
      character*1 pr
      common/pwhgprocess/pr
C - We need to tell to this analysis file which program is running it
      character*6 WHCPRG
      common/cWHCPRG/WHCPRG
      data    WHCPRG/'NLO   '/
      logical debug
      parameter (debug=.false.)

C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,maxjet
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
      enddo
      do mu=1,4
         pB(mu)=0d0
         ptot(mu)=0d0
      enddo
      if(iflag.eq.1) then
C - Extract final state particles to feed to jet finder
         do j=1,nhep
C - (all but the boson).
            if(WHCPRG.eq.'NLO   ') then
               if(pr.eq.'H') then
                  if(j.eq.3) then
                     pB = phep(1:4,j) 
                     cycle
                  endif
               elseif(pr.eq.'Z') then
                  if(j.eq.3.or.j.eq.4) then
                     pB=pB+phep(1:4,j) 
                     cycle
                  endif
               elseif(pr.eq.'W') then
                  if(j.eq.3.or.j.eq.4) then
                     pB=pB+phep(1:4,j) 
                     cycle
                  endif
               endif
            elseif(WHCPRG.eq.'LHE   ') then
               if(pr.eq.'H') then
                  if(j.eq.3) then
                     pB=phep(1:4,j) 
                     cycle
                  endif
               elseif(pr.eq.'Z') then
                  if(j.eq.4.or.j.eq.5) then
                     pB=pB+phep(1:4,j) 
                  endif
                  if(j.eq.3.or.j.eq.4.or.j.eq.5) cycle
               elseif(pr.eq.'W') then
                  if(j.eq.4.or.j.eq.5) then
                     pB=pB+phep(1:4,j) 
                  endif
                  if(j.eq.3.or.j.eq.4.or.j.eq.5) cycle
               endif
            elseif ((WHCPRG.eq.'HERWIG'.and.(isthep(j).eq.147.or.
     $                                       isthep(j).eq.148.or.
     $                                       isthep(j).eq.149.or.
     $                                       isthep(j).eq.155.or.
     $                                       isthep(j).eq.190))
     $           .or.
     $              (WHCPRG.eq.'PYTHIA'.and.isthep(j).eq.1)) then
               if(pr.eq.'H') then
                  if(idhep(j).eq.25) then
                     pB=phep(1:4,j)
                     cycle
                  endif
               elseif(pr.eq.'Z') then
                  if(idhep(jmohep(1,j)).eq.23) then
                     pB=phep(1:4,jmohep(1,j))
                     cycle
                  elseif(jmohep(1,jmohep(1,j)).ne.0) then
                     if(idhep(jmohep(1,jmohep(1,j))).eq.23) then
                        pB=phep(1:4,jmohep(1,jmohep(1,j)))
                        cycle
                     endif
                  endif
               elseif(pr.eq.'W') then
                  if(abs(idhep(jmohep(1,j))).eq.24) then
                     pB=phep(1:4,jmohep(1,j))
                     cycle
                  elseif(jmohep(1,jmohep(1,j)).ne.0) then
                     if(abs(idhep(jmohep(1,jmohep(1,j)))).eq.24) then
                        pB=phep(1:4,jmohep(1,jmohep(1,j)))
                        cycle
                     endif
                  endif
               endif
            else 
               cycle
            endif
            if(ntracks.eq.maxtrack) then
               write(*,*) 'analyze: need to increase maxtrack!'
               write(*,*) 'ntracks: ',ntracks
               stop
            endif
            ntracks=ntracks+1
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,j)
               ptot(mu)=ptot(mu)+phep(mu,j)
            enddo
            itrackhep(ntracks)=j
         enddo

         do mu=1,4
            ptot(mu)=ptot(mu)+pB(mu)
         enddo
         ptsq_tot = ptot(1)**2+ptot(2)**2
         if (ptsq_tot .ge. 1d-3) then 
            write(*,*) '----> WARNING: ptot(1:2) /= (0,0)', ptot 
            write(*,*) '----> WARNING: ptsq_tot          ', ptsq_tot
         endif

      else
C     this will not work for Higgs at the moment... 
         do j=1,nup
            if (istup(j).eq.1.and..not.islept(idup(j))) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=pup(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      endif
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------------------------------------------------------------- C
C - rr = radius parameter (R=0.4 / 0.7 typical for ATLAS - anti kT)   - C
C - palg=1 is standard kt, 0 is Cambridge Aachen, -1 is antikt        - C
C --------------------------------------------------------------------- C
      r=rr
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                    jetvec,yijs)
      mjets=njets
      if(njets.eq.0) return
      if(debug) then
c     check consistency
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
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         call getyetaptmass(pjet(:,j),rap(j),eta(j),kt(j),tmp)
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo

c     loop over the hardest 3 jets
      do j=1,min(njets,3)
         do mu=1,3
            pjetin(mu) = pjet(mu,j)
         enddo
         pjetin(0) = pjet(4,j)         
         vec(1)=0d0
         vec(2)=0d0
         vec(3)=1d0
         beta = -pjet(3,j)/pjet(4,j)
         call mboost(1,vec,beta,pjetin,pjetout)         
         ptrel(j) = 0
         do i=1,ntracks
            if (jetvec(i).eq.j) then
               do mu=1,3
                  ptrackin(mu) = ptrack(mu,i)
               enddo
               ptrackin(0) = ptrack(4,i)
               call mboost(1,vec,beta,ptrackin,ptrackout) 
               ptrel(j) = ptrel(j) + get_ptrel(ptrackout,pjetout)
            endif
         enddo
      enddo
      end

      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.16) then
         islept = .true.
      else
         islept = .false.
      endif
      end

      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end

      subroutine pwhgfinalopshist
      implicit none
      include     'pwhg_bookhist-new.h'
      integer      maxjet
      parameter   (maxjet=2048)
      integer      nptmin
      parameter   (nptmin=4)
      character*1  cnums(9)
      character*4  cptcuts(nptmin)
      real*8       ptcuts(nptmin)
      character*4  cycuts(5)
      real*8       ycuts(5)
      character*4  calgs(3)
      real*8       jetalgs(3)
      common/infohist/ptcuts,ycuts,jetalgs,cnums,cptcuts,cycuts,calgs
      save  /infohist/
      character*4  cptcut,cycut,calg
      character*1  ctmp
      integer      a_idx,b_idx,d_idx
      integer      ixx,jxx,kxx
      integer      indexhist
      character*1  ctmp1
      character*8  suffix
      character*12 suffix12

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C
C - Relative differences in pT spectra between jet algorithms - C
C - (done in pwhgfinalopshist):                               - C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C

C - Loop over CA and AkT jet algorithms:
      do ixx=2,3
        suffix=calgs(ixx)//calgs(1)

C - Jet pT spectra:
        do kxx=1,3
          ctmp1=cnums(kxx)
C - "Inclusively"
          a_idx=indexhist('j'//ctmp1//'-pt'//calgs(1))
          b_idx=indexhist('j'//ctmp1//'-pt'//calgs(ixx))
          d_idx=indexhist('j'//ctmp1//'-pt'//suffix)
          call get_diff_hists(a_idx,b_idx,d_idx)

C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
          do jxx=2,nptmin
          cptcut=cptcuts(jxx)
          suffix12=cptcut//suffix
C - in >=1 jet events
          a_idx=indexhist('Nj>=1-j'//ctmp1//'-pt'//cptcut//calgs(1))
          b_idx=indexhist('Nj>=1-j'//ctmp1//'-pt'//cptcut//calgs(ixx))
          d_idx=indexhist('Nj>=1-j'//ctmp1//'-pt'//suffix12)
          call get_diff_hists(a_idx,b_idx,d_idx)
C - in >=2 jet events
          a_idx=indexhist('Nj>=2-j'//ctmp1//'-pt'//cptcut//calgs(1))
          b_idx=indexhist('Nj>=2-j'//ctmp1//'-pt'//cptcut//calgs(ixx))
          d_idx=indexhist('Nj>=2-j'//ctmp1//'-pt'//suffix12)
          call get_diff_hists(a_idx,b_idx,d_idx)
          enddo
C - End loop over jet pT cuts

        enddo

C - Zoomed jet pT spectra:
        do kxx=1,3
          ctmp1=cnums(kxx)
C - "Inclusively"
          a_idx=indexhist('j'//ctmp1//'-ptzm'//calgs(1))
          b_idx=indexhist('j'//ctmp1//'-ptzm'//calgs(ixx))
          d_idx=indexhist('j'//ctmp1//'-ptzm'//suffix)
          call get_diff_hists(a_idx,b_idx,d_idx)

C - Loop over jet pT cuts (skipping 1st one which goes down to 1 GeV jets):
          do jxx=2,nptmin
          cptcut=cptcuts(jxx)
          suffix12=cptcut//suffix
C - in >=1 jet events
          a_idx=indexhist('Nj>=1-j'//ctmp1//'-ptzm'//cptcut//calgs(1))
          b_idx=indexhist('Nj>=1-j'//ctmp1//'-ptzm'//cptcut//calgs(ixx))
          d_idx=indexhist('Nj>=1-j'//ctmp1//'-ptzm'//suffix12)
          call get_diff_hists(a_idx,b_idx,d_idx)
C - in >=2 jet events
          a_idx=indexhist('Nj>=2-j'//ctmp1//'-ptzm'//cptcut//calgs(1))
          b_idx=indexhist('Nj>=2-j'//ctmp1//'-ptzm'//cptcut//calgs(ixx))
          d_idx=indexhist('Nj>=2-j'//ctmp1//'-ptzm'//suffix12)
          call get_diff_hists(a_idx,b_idx,d_idx)
          enddo
C - End loop over jet pT cuts
        enddo

C - End loop over CA and AkT algorithms
      enddo

      end

      subroutine get_diff_hists(a_idx,b_idx,d_idx)
      implicit none
      include 'pwhg_bookhist-new.h'
      real*8   a,ea,b,eb
      integer  ixx,a_idx,b_idx,d_idx
      
      do ixx=1,nbins(a_idx)
         a=yhistarr2(ixx,a_idx)
         ea=errhistarr2(ixx,a_idx)
         b=yhistarr2(ixx,b_idx)
         eb=errhistarr2(ixx,b_idx)
         if((a+b).gt.0d0) then         ! Guard against division by zero.
            yhistarr2(ixx,d_idx)=(a-b)/(a+b)
            errhistarr2(ixx,d_idx)=2*sqrt((b*ea)**2+(a*eb)**2)
     1                              /(a+b)**2
         else
            yhistarr2(ixx,d_idx)=0d0
            errhistarr2(ixx,d_idx)=0d0
         endif
      enddo

      end


