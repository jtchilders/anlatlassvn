c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'

      call inihists

c     total cross section sanity check
      call bookupeqbins('total',1d0,0d0,1d0)

      call bookupeqbins('HTTOT',20d0,0d0,500d0)

c     TM doing m34 and m56 regardless
      call bookupeqbins('massZ34',5d0,0d0,300d0)
      call bookupeqbins('massZ56',5d0,0d0,300d0)
      call bookupeqbins('massZ34b',2d0,20d0,150d0)
      call bookupeqbins('massZ56b',2d0,20d0,150d0)

      call bookupeqbins('massZ1',2d0,0d0,500d0)
      call bookupeqbins('mt_Z1',10d0,0d0,500d0)
      call bookupeqbins('pt_Z1',10d0,0d0,500d0)
      call bookupeqbins('y_Z1',0.2d0,-4d0,4d0)
      call bookupeqbins('delphi_Z1',Pi/10,0d0,Pi)

      call bookupeqbins('massZ2',2d0,0d0,500d0)
      call bookupeqbins('mt_Z2',10d0,0d0,500d0)
      call bookupeqbins('pt_Z2',10d0,0d0,500d0)
      call bookupeqbins('y_Z2',0.2d0,-4d0,4d0)
      call bookupeqbins('delphi_Z2',Pi/10,0d0,Pi)

      call bookupeqbins('ZZmass',2d0,0d0,500d0)
      call bookupeqbins('pt_ZZ',10d0,0d0,500d0)
      call bookupeqbins('eta_ZZ',0.2d0,-4d0,4d0)
      call bookupeqbins('y_ZZ',0.2d0,-4d0,4d0)

      call bookupeqbins('pt(l+,1)',10d0,0d0,500d0)
      call bookupeqbins('pt(l+,2)',10d0,0d0,500d0)
      call bookupeqbins('pt(l-,1)',10d0,0d0,500d0)
      call bookupeqbins('pt(l-,2)',10d0,0d0,500d0)

      call bookupeqbins('eta(l+,1)',0.2d0,-4d0,4d0)
      call bookupeqbins('eta(l+,2)',0.2d0,-4d0,4d0)
      call bookupeqbins('eta(l-,1)',0.2d0,-4d0,4d0)
      call bookupeqbins('eta(l-,2)',0.2d0,-4d0,4d0)

      call bookupeqbins('deltall3',0.01d0,0d0,1d0)

      call bookupeqbins('ch_y_asyml',0.2d0,-4d0,4d0)
      call bookupeqbins('ch_pt_asyml',5d0,0d0,500d0)

C     -- variable with jets 
      call bookupeqbins('dpt_jnocut',10d0,0d0,500d0)

      call bookupeqbins('eta_j020',0.2d0,-4d0,4d0)
      call bookupeqbins('dy_jZZ020', 0.2d0,-4d0,4d0)
      call bookupeqbins('deta_jZZ020',0.2d0,-4d0,4d0)

      call bookupeqbins('eta_j040',0.2d0,-4d0,4d0)
      call bookupeqbins('dy_jZZ040', 0.2d0,-4d0,4d0)
      call bookupeqbins('deta_jZZ040',0.2d0,-4d0,4d0)

      call bookupeqbins('eta_j100',0.2d0,-4d0,4d0)
      call bookupeqbins('dy_jZZ100', 0.2d0,-4d0,4d0)
      call bookupeqbins('deta_jZZ100',0.2d0,-4d0,4d0)

C     -- additions... 
      call bookupeqbins('massZ1b',2d0,0d0,200d0)
      call bookupeqbins('massZ2b',2d0,0d0,200d0)


      end
     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      real *8 pt_lplus,pt_lminus,eta_lplus,eta_lminus,
     $delphi,mt_v,mv,ptv,yv
      real *8 pt_muplus,pt_muminus,eta_muplus,eta_muminus
      logical passcuts,passjetcuts
      integer ihep,mu
      logical ini
      data ini/.true./
      save ini
c     binsize
c     we need to tell to this analysis file which program is running it
      integer   maxjet,mjets
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer vdecaytemp1,vdecaytemp2
      save vdecaytemp1,vdecaytemp2

      integer nlm1, nlm2, nlp1,nlp2
      integer ilm1(100),ilm2(100),ilp1(100),ilp2(100)
      real * 8 plm1(4),plp1(4),plm2(4),plp2(4)
      real * 8 httot,y3,y4,y5,y6,y,eta,mass,ptl3,ptl4,ptl5,ptl6,maxptl,
     1     etal3,etal4,etal5,etal6
      real * 8 dy,deta,dphi,dr
      real * 8 m34,m36,m45,m56,y34,y36,y45,y56,pt34,pt36,pt45,pt56,
     1     delphi34,delphi36,delphi45,delphi56
      real * 8 Zmassdistr(4), ptdistr(4),ydistr(4),delphidistr(4)
      real * 8 DZmass,massZ1,ptZ1,yZ1,delphiZ1, mtZ1
      real * 8 massZ2,ptZ2,yZ2,delphiZ2, mtZ2
      real * 8 pt_l1minus, pt_l1plus  
      real * 8 pt_l2minus, pt_l2plus  
      real * 8 fb3,fb4,fb5,fb6
      real * 8 yZZ, etaZZ,ptZZ,ZZmass 
      real * 8 pt_j,eta_j, dy_jZZ,deta_jZZ, jetcut(3), deltall3
      character * 12 dyjZZstring, detajZZstring,etajstring
      character * 8 ptjstring
      character * 3 jetcutstring(3)
      integer j, iZmd,ijc,izmdmax,bestiZmd, iZ1, iZ2 
      real * 8 powheginput
      external powheginput
      real * 8 zmass
      parameter (zmass=91.1876d0)
      real * 8 mllmin
      save mllmin

      if(dsig.eq.0) return

      if (ini) then
         vdecaytemp1=powheginput("vdecaymodeZ1")
         vdecaytemp2=powheginput("vdecaymodeZ2")
         mllmin=powheginput('mllmin')
         write (*,*)
         write (*,*) '********************************************'
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO analysis'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE analysis'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         write(*,*) '*****************************'

         if((vdecaytemp1.eq.11.or.vdecaytemp1.eq.13
     $        .or.vdecaytemp1.eq.15) .and. 
     $        (vdecaytemp2.eq.11.or.vdecaytemp2.eq.13
     $        .or.vdecaytemp2.eq.15)) then
            continue
         else
            write(*,*) '**************************************'
            write(*,*) 'template analysis works only for e, mu or tau'
            write(*,*) '                 STOP     '
            write(*,*) '**************************************'
            call exit(1)
         endif
         ini=.false.
      endif


c     find Z decay products
      nlm1=0                   ! number of l-
      nlm2=0                   ! number of l-
      nlp1=0                   ! number of l+
      nlp2=0                   ! number of l+
      do ihep=1,nhep
         if(isthep(ihep).eq.1) then
            if(idhep(ihep).eq.vdecaytemp1) then
               nlm1 = nlm1+1
               ilm1(nlm1) = ihep
            elseif(idhep(ihep).eq.-vdecaytemp1) then
               nlp1 = nlp1+1
               ilp1(nlp1) = ihep
            elseif(idhep(ihep).eq.vdecaytemp2) then
               nlm2 = nlm2+1
               ilm2(nlm2) = ihep
            elseif(idhep(ihep).eq.-vdecaytemp2) then
               nlp2 = nlp2+1
               ilp2(nlp2) = ihep
            endif
         endif
      enddo

      if (nlm1 .ge. 100 .or. nlm2 .ge. 100 .or. nlp1 .ge. 100 .or.
     .     nlp2 .ge. 100) then
         write(*,*) 'crazy event, too many leptons'
         return
      endif
      
      if (vdecaytemp1 .eq. vdecaytemp2) then
c - we won't be able to tell the difference between the leptons, so they'll
c - all be in nlm1 and nlp1.
         if(nlm1.lt.2 .or. nlp1.lt.2) then 
            write(*,*) ' crazy event, missing leptons',
     .           nlm1,nlm2,nlp1,nlp2
            return
         endif
      else
         if (nlm1.lt.1 .or. nlm2.lt.1 .or. nlp1.lt.1 .or. nlp2.lt.1)then
            write(*,*) ' crazy event, missing leptons',
     .           nlm1,nlm2,nlp1,nlp2
            return
         endif
      endif


c     sort by pt 
      if (vdecaytemp1 .eq. vdecaytemp2) then
c can't tell where the leptons are from, so use two hardest l- and l+
c should be fine, since we consider all combinations of l-l+         
         call sortbypt(nlm1,ilm1(1:nlm1))
         call sortbypt(nlp1,ilp1(1:nlp1))

         plm1(1:4)=phep(1:4,ilm1(1))
         plm2(1:4)=phep(1:4,ilm1(2))
         plp1(1:4)=phep(1:4,ilp1(1))
         plp2(1:4)=phep(1:4,ilp1(2))
      else
c now we can tell the difference between the lepton pairs         
         call sortbypt(nlm1,ilm1(1:nlm1))
         call sortbypt(nlm2,ilm2(1:nlm2))
         call sortbypt(nlp1,ilp1(1:nlp1))
         call sortbypt(nlp2,ilp2(1:nlp2))

         plm1(1:4)=phep(1:4,ilm1(1))
         plm2(1:4)=phep(1:4,ilm2(1))
         plp1(1:4)=phep(1:4,ilp1(1))
         plp2(1:4)=phep(1:4,ilp2(1))
      endif
C PN&GZ change R:0.7->0.6 and kt -> antikt 7/6/2011, ptmin 20 -> 1 
      rr=0.6d0 
      call buildjets(1,mjets,rr,ktj,etaj,rapj,phij,pj)

      call getyetaptmass(plm1+plp1,y34,eta,pt34,m34)
      call getyetaptmass(plm2+plp2,y56,eta,pt56,m56)      

      if(vdecaytemp1.eq.vdecaytemp2) then
         call getyetaptmass(plp1+plm2,y45,eta,pt45,m45)
         call getyetaptmass(plm1+plp2,y36,eta,pt36,m36)
      endif
      
      call getyetaptmass(plm1,y3,etal3,ptl3,mass)
      call getyetaptmass(plp1,y4,etal4,ptl4,mass)
      call getyetaptmass(plm2,y5,etal5,ptl5,mass)
      call getyetaptmass(plp2,y6,etal6,ptl6,mass)

      maxptl=max(ptl3,ptl4,ptl5,ptl6)

c -- lepton cuts: highest pt_l > 20 GeV, all other pt_l > 10 GeV
c -- |eta| < 2.5
c -- independent of lepton flavour      
      if (maxptl.lt.20d0) return

C     PN&GZ corrected but on pt_l -- (max -> min) 7/6/2011
      if(maxptl.eq.ptl3) then
         if(min(ptl4,ptl5,ptl6).lt.10) return
      elseif(maxptl.eq.ptl4) then
         if(min(ptl3,ptl5,ptl6).lt.10) return
      elseif(maxptl.eq.ptl5) then
         if(min(ptl3,ptl4,ptl6).lt.10) return
      elseif(maxptl.eq.ptl6) then
         if(min(ptl3,ptl4,ptl5).lt.10) return
      endif

      if (abs(etal3) .gt. 2.5 .or. abs(etal4) .gt. 2.5 .or. 
     .     abs(etal5) .gt. 2.5 .or. abs(etal6) .gt. 2.5) return

      httot=ptl3+ptl4+ptl5+ptl6
      do j=1,mjets
         if (ktj(j) > 20d0) httot=httot+ktj(j)
      enddo


c     azimuthal separation between e^+e^-

      call getdydetadphidr(plm1,plp1,dy,deta,delphi34,dr)
      call getdydetadphidr(plm2,plp2,dy,deta,delphi56,dr)
      if(vdecaytemp1.eq.vdecaytemp2) then         
         call getdydetadphidr(plm1,plp2,dy,deta,delphi36,dr)
         call getdydetadphidr(plp1,plm2,dy,deta,delphi45,dr)
      endif
      
c --  find the lepton pair with the inv mass closest to M_Z,
      Zmassdistr = (/m34,m56,m45,m36/)
      ptdistr = (/pt34,pt56,pt45,pt36/)
      ydistr = (/y34,y56,y45,y36/)
      delphidistr = (/delphi34,delphi56,delphi45,delphi36/)
      
      DZmass = 1d6
      if(vdecaytemp1.eq.vdecaytemp2) then
         izmdmax = 4
      else
         izmdmax = 2
      endif

      do iZmd = 1,izmdmax
         if (abs(Zmassdistr(iZmd) - zmass) .le. DZmass) then

            DZmass=abs(Zmassdistr(iZmd) -zmass)
            iZ1=iZmd
            
            massZ1 = Zmassdistr(iZmd)
            ptZ1 = ptdistr(iZmd)
            yZ1 = ydistr(iZmd)
            delphiZ1 = delphidistr(iZmd)
         endif
      enddo
 
     
      if (iZ1 .eq. 1) then
         pt_l1minus = ptl3
         pt_l1plus=ptl4
         pt_l2minus = ptl5
         pt_l2plus=ptl6
         massZ2 = m56
         ptZ2 = pt56
         yZ2 = y56
         delphiZ2 = delphi56
      elseif (iZ1 .eq. 2) then
         pt_l1minus = ptl5
         pt_l1plus=ptl6
         pt_l2minus = ptl3
         pt_l2plus=ptl4
         massZ2 = m34
         ptZ2 = pt34
         yZ2 = y34
         delphiZ2 = delphi34
      elseif (iZ1 .eq. 3) then
         pt_l1minus = ptl5
         pt_l1plus=ptl4
         pt_l2minus = ptl3
         pt_l2plus=ptl6
         massZ2 = m36
         ptZ2 = pt36
         yZ2 = y36
         delphiZ2 = delphi36
      elseif (iZ1 .eq. 4) then
         pt_l1minus = ptl3
         pt_l1plus=ptl6
         pt_l2minus = ptl5
         pt_l2plus=ptl4
         massZ2 = m45
         ptZ2 = pt45
         yZ2 = y45
         delphiZ2 = delphi45
      endif

c     transverse mass of e^+e^- system
      mtZ1=sqrt(2*pt_l1plus*pt_l1minus*(1d0-dcos(delphiZ1)))
      mtZ2=sqrt(2*pt_l2plus*pt_l2minus*(1d0-dcos(delphiZ2)))
      
      call filld('total',0.5d0,dsig)
      call filld('HTTOT',httot,dsig)

c     here plot the z masses without finding Z_1 and Z_2

      call filld('massZ34',m34,dsig)
      call filld('massZ56',m56,dsig)
      call filld('massZ34b',m34,dsig)
      call filld('massZ56b',m56,dsig)
      
      call filld('massZ1',massZ1,dsig)
      call filld('massZ1b',massZ1,dsig)
      call filld('pt_Z1', ptZ1,dsig)
      call filld('y_Z1', yZ1,dsig)
      call filld('delphi_Z1', delphiZ1,dsig)
      call filld('mt_Z1', mtZ1,dsig)

      call filld('massZ2',massZ2,dsig)
      call filld('massZ2b',massZ2,dsig)
      call filld('pt_Z2', ptZ2,dsig)
      call filld('y_Z2', yZ2,dsig)
      call filld('delphi_Z2', delphiZ2,dsig)
      call filld('mt_Z2', mtZ2,dsig)
      
      call getyetaptmass(plm1+plp1+plm2+plp2,yZZ,etaZZ,ptZZ,ZZmass)

      call filld('ZZmass',ZZmass,dsig)
      call filld('pt_ZZ', ptZZ, dsig)
      call filld('eta_ZZ',etaZZ,dsig)
      call filld('y_ZZ',yZZ,dsig)
      
c --  these are mostly for checking the calculation      
c     pt(l+) (all)
      call filld('pt(l+,1)',ptl4,dsig)
      call filld('pt(l+,2)',ptl6,dsig)
c     pt(l-) (all)
      call filld('pt(l-,1)',ptl3,dsig)
      call filld('pt(l-,2)',ptl5,dsig)

c     eta(l+) (all)
         call filld('eta(l+,1)',etal4,dsig)
         call filld('eta(l+,2)',etal6,dsig)
c     eta(l-) (all)
         call filld('eta(l-,1)',etal3,dsig)
         call filld('eta(l-,2)',etal5,dsig)


c     y charge asymmetry 
      call filld('ch_y_asyml', etal3,-dsig/4d0)
      call filld('ch_y_asyml', etal4, dsig/4d0)
      call filld('ch_y_asyml', etal5,-dsig/4d0)
      call filld('ch_y_asyml', etal6, dsig/4d0)

c     pt charge asymmetry 
      call filld('ch_pt_asyml', ptl3,-dsig/4d0)
      call filld('ch_pt_asyml', ptl4, dsig/4d0)
      call filld('ch_pt_asyml', ptl5,-dsig/4d0)
      call filld('ch_pt_asyml', ptl6, dsig/4d0)

      deltall3 = (  ((plm1(1)-plm2(1))**2
     .     +(plm1(2)-plm2(2))**2
     .     +(plm1(3)-plm2(3))**2   )/
     .      (plm1(1)**2+plm2(1)**2
     .     + plm1(2)**2+plm2(2)**2
     .     + plm1(3)**2+plm2(3)**2 )   )**(1.5d0) 
      
      call filld('deltall3',deltall3,dsig/2d0)  
      deltall3 = ((plp1(1)-plp2(1))**2
     .     +(plp1(2)-plp2(2))**2
     .     +(plp1(3)-plp2(3))**2)**(1.5d0) 
      call filld('deltall3',deltall3,dsig/2d0)  
         
         
c     --- jet plots      
         if(mjets.ge.1) then
            
            passjetcuts=.true.

            call filld('dpt_jnocut', ktj(1),dsig)

c -- jet cuts: pt >  input jet cut (20, 40, 100 GeV), |eta| < 3.5           
            jetcut = (/20d0, 40d0, 100d0/)
            jetcutstring = (/'020', '040', '100'/)
            do ijc = 1,3
               
               if ((ktj(1) .le. jetcut(ijc)) .or. 
     .              (abs(etaj(1)) .ge. 3.5d0)) then
                  passjetcuts = .false.
               endif


               if (passjetcuts) then

                  etajstring = 'eta_j'//(jetcutstring(ijc))
                  dyjZZstring = 'dy_jZZ'//(jetcutstring(ijc))
                  detajZZstring = 'deta_jZZ'//(jetcutstring(ijc))
                  
                  call getdydetadphidr(plm1+plp1+plm2+plp2,pj(1:4,1),
     .                 dy_jZZ,deta_jZZ,dphi,dr)

                  call filld(etajstring,etaj(1),dsig)
                  call filld(dyjZZstring, dy_jZZ, dsig)
                  call filld(detajZZstring, deta_jZZ, dsig)
               endif
            enddo
         endif
            
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

      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(0:3),y
      y=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end

      subroutine getinvmass(p,m)
      implicit none
      real * 8 p(0:3),m
      m=sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end

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



      subroutine buildjets(iflag,mjets,rr,kt,eta,rap,phi,pjet)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
      integer iflag,mjets
      real * 8  rr,kt(*),eta(*),rap(*),
     1     phi(*),pjet(4,*)
      include   'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,jb
      real * 8 r,palg,ptmin,pp,tmp
      logical islept
      external islept
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
      if(iflag.eq.1) then
C     - Extract final state particles to feed to jet finder
         do j=1,nhep
            if (isthep(j).eq.1.and..not.islept(idhep(j))) then
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
      else
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
C - Inclusive jet pT and Y spectra are to be compared to CDF data:    - C    
C --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
C     f = 0.75  overlapping fraction
c palg=1 is standard kt, -1 is antikt
      palg=-1
      r=rr
      ptmin=1d0 
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)
      mjets=njets
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

      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.15) then
         islept = .true.
      else
         islept = .false.
      endif
      end
