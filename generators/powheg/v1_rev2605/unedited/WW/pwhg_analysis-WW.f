c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  bookup  : opens histograms
c  filld   : fills histograms with data

c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include '../pwhg_book.h'

      call inihists

c$$$ total cross section sanity check:
      call bookupeqbins('total',1d0,0d0,1d0)
      call bookupeqbins('HTTOT',10d0,0d0,500d0)
      call bookupeqbins('eta_lplus',0.2d0,-4d0,4d0)
      call bookupeqbins('pt_lplus',10d0,0d0,500d0)
      call bookupeqbins('eta_lminus',0.2d0,-4d0,4d0)
      call bookupeqbins('pt_lminus',10d0,0d0,500d0)
      call bookupeqbins('eta_ll',0.2d0,-4d0,4d0)
      call bookupeqbins('y_ll',0.2d0,-4d0,4d0)
      call bookupeqbins('m_ll',10d0,0d0,500d0)
      call bookupeqbins('pt_ll',10d0,0d0,500d0)
      call bookupeqbins('pt_miss',10d0,0d0,500d0)
      call bookupeqbins('phi_ll',0.314d0,0d0,3.14d0)
      call bookupeqbins('mT_WW',10d0,0d0,500d0)
      call bookupeqbins('FB_asym_lm',1d0,-3d0,3d0)
      call bookupeqbins('FB_asym_lp',1d0,-3d0,3d0)
      call bookupeqbins('pt_jnocut', 10d0,0d0,500d0)
      call bookupeqbins('eta_j',0.2d0,-4d0,4d0)
      call bookupeqbins('pt_j',10d0,0d0,500d0)
      call bookupeqbins('R_lpj',0.1d0,0d0,4d0)
      call bookupeqbins('R_lmj',0.1d0,0d0,4d0)
      call bookupeqbins('pt_leadlept', 10d0,0d0,500d0)
      call bookupeqbins('pt_alllept', 10d0,0d0,500d0)

      end


      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      include 'cvecbos.h'

      integer ihep,mu
      logical ini
      data ini/.true./
      save ini
      integer   maxjet,mjets
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer vdecaytemp1,vdecaytemp2,vdecaytemp3,vdecaytemp4
      integer nlminus, nlplus, nnu,nnubar
      integer ilminus(100),ilplus(100),inu(100),inubar(100)
      real * 8 plplus(4),plminus(4),pnu(4),pnubar(4),pmiss(4)
      real *8 pt_lplus,pt_lminus,eta_lplus,eta_lminus,pt_lead,pt_alllept
      real * 8 eta_nu, eta_nubar, pt_nu, pt_nubar
      real * 8 y_lplus,y_lminus,y_nu,y_nubar, eta_34,pt_34,mass_34,y_34
      real * 8 pt_345,eta_345,y_345,mass_345,dr_45,dphi_45
      real * 8 eta_56,pt_56,dr_56,dphi_56,mass_56,dr_57,dr_67,y_56,dr_47
      real * 8 mass_567,etadiff
      real * 8 httot,ht,y,eta,pt,mass,dy,dr,deta,dphi
      real * 8 ET_miss,ET_ll,dphi_llmiss,MT_WWsq,MT_WW,pt_miss,pxnu,pynu
      real * 8 pt_ll,eta_ll,y_ll,deta_ll,dphi_ll,mass_ll
      real * 8 pt_j, eta_j, dr_jlmi, dr_jlpl
      real * 8 fblminus, fblplus
      real * 8 powheginput
      external powheginput

      integer i3,i4,i5,i6,j,nu,ijc

      if (ini) then
         vdecaymodeWp=powheginput("vdecaymodeWp")
         vdecaymodeWm=powheginput("vdecaymodeWm")
         write (*,*)
         write (*,*) '********************************************'
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO analysis'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE analysis'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
            write(*,*) 'not implemented analysis'
            write(*,*) 'no plots will be present at the end of the run'
          elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
            write(*,*) 'not implemented analysis'
            write(*,*) 'no plots will be present at the end of the run'
         endif
         write(*,*) '*****************************'

         if((abs(vdecaymodeWp).eq.11.or.abs(vdecaymodeWp).eq.13
     $        .or.abs(vdecaymodeWp).eq.15) .and. 
     $      (abs(vdecaymodeWm).eq.11.or.abs(vdecaymodeWm).eq.13
     $        .or.abs(vdecaymodeWm).eq.15)) then
            continue
         else
            write(*,*) '**************************************'
            write(*,*) ' template analysis works only for e, mu and tau'
            write(*,*) '                 STOP     '
            write(*,*) '**************************************'
            call exit(1)
         endif
         ini=.false.
      endif

	vdecaytemp1=vdecaymodeWm   ! l-
	vdecaytemp3=vdecaymodeWp   ! l+
         
      vdecaytemp2 =-sign(1,vdecaytemp1)*(abs(vdecaytemp1)+1)   ! nubar
      vdecaytemp4 =-sign(1,vdecaytemp3)*(abs(vdecaytemp3)+1)   ! nu

c     find W decay products
         nlminus = 0           ! number of l-
         nlplus = 0            ! number of l+
         nnu = 0               ! number of nu
         nnubar=0              ! number of nubar

c     conventions: i3=l-,i5=l+, i4,i6 neutrinos 
         do ihep=1,nhep
            if(isthep(ihep).eq.1) then
               if(idhep(ihep).eq.vdecaytemp1) then
                  nlminus = nlminus+1
                  ilminus(nlminus)=ihep
               elseif(idhep(ihep).eq.vdecaytemp2) then
                  nnubar = nnubar+1
                  inubar(nnubar)=ihep
               elseif(idhep(ihep).eq.vdecaytemp3) then
                  nlplus = nlplus+1
                  ilplus(nlplus)=ihep
               elseif(idhep(ihep).eq.vdecaytemp4) then
                  nnu = nnu+1
                  inu(nnu)=ihep
               endif
            endif
         enddo

         if (nlminus .ge. 100 .or. nlplus .ge. 100 .or. nnu .ge. 100
     .        .or. nnubar .ge. 100) then
            write(*,*) 'crazy event, too many leptons'
            return
         endif

         if(nlminus.lt.1 .or. nlplus.lt.1 .or. 
     .        nnu.lt.1 .or. nnubar.lt.1) then 
            write(*,*) ' crazy event, missing leptons',
     .           nlminus,nlplus,nnu,nnubar
            return
         endif
         
         call sortbypt(nlminus,ilminus(1:nlminus))
         call sortbypt(nlplus,ilplus(1:nlplus))
         call sortbypt(nnu,inu(1:nnu))
         call sortbypt(nnubar,inubar(1:nnubar))
         
         
         plminus(1:4)=phep(1:4,ilminus(1))
         plplus(1:4)=phep(1:4,ilplus(1))
         pnu(1:4)=phep(1:4,inu(1))
         pnubar(1:4)=phep(1:4,inubar(1))
         
c     --  Leptons identified: W^+ -> (l+(plplus)+nu(pnu))
c     --                      W^- -> (l-(plminus)+nubar(pnubar))
         
         rr=0.6d0
         call buildjets(1,mjets,rr,ktj,etaj,rapj,phij,pj)
         
         call getyetaptmass(plminus,y_lminus,eta_lminus,pt_lminus,mass)
         call getyetaptmass(pnu,y_nu,eta_nu,pt_nu,mass)
         call getyetaptmass(plplus,y_lplus,eta_lplus,pt_lplus,mass)
         call getyetaptmass(pnubar,y_nubar,eta_nubar,pt_nubar,mass)
         
c --- pt > 20 GeV, |eta| < 2.5, pt_miss >  20 GeV      
      if ((pt_lplus .le. 20d0) .or. (pt_lminus .le. 20d0)) return

      if ((abs(eta_lplus) .gt. 2.5d0) .or. 
     .     (abs(eta_lminus) .gt. 2.5d0)) then 
         return
      endif

      call getyetaptmass(pnu+pnubar,y,eta,pt_miss,mass)

      if (pt_miss .le. 20d0) return

      httot=pt_lminus+pt_lplus+pt_miss
      ht = pt_nu +pt_nubar+pt_lplus + pt_lminus 

      do j=1,mjets
c --  need to apply an additional cut here, as there are jets with pt < 20 GeV
         if (ktj(j) .ge. 20d0) then
            httot=httot+ktj(j)
            ht = ht + ktj(j)
         endif
      enddo

c --  find which lepton is the hardest
      if (pt_lminus .ge. pt_lplus) then
         pt_lead = pt_lminus
      else
         pt_lead = pt_lplus
      endif

      call getyetaptmass(plminus+plplus,y_ll,eta_ll,pt_ll,
     &     mass_ll)

      call getdydetadphidr(plminus,plplus,dy,deta_ll,
     .     dphi_ll,dr)

      pt_alllept = pt_ll + pt_miss
      
c --  transverse mass of the WW system
      ET_ll = sqrt(pt_ll**2 + mass_ll**2)
      ET_miss=sqrt(pt_miss**2 + mass_ll**2)
      call getdydetadphidr(plminus+plplus,pnu+pnubar,dy,deta,
     .     dphi_llmiss,dr)
      MT_WWsq = 2d0*(ET_ll*ET_miss + mass_ll**2 - 
     .     pt_miss*pt_ll*dcos(dphi_llmiss))
      MT_WW = sqrt(MT_WWsq)
      
      call filld('HTTOT',httot,dsig)
      call filld('total', 0.0d0, dsig)
c     -- transverse momenta and rapidity of leptons
      call filld('eta_lplus',eta_lplus,dsig)
      call filld('pt_lplus',pt_lplus,dsig)
      call filld('eta_lminus',eta_lminus,dsig)
      call filld('pt_lminus',pt_lminus,dsig)
      call filld('pt_leadlept', pt_lead,dsig)
      call filld('pt_alllept', pt_alllept,dsig)
      call filld('eta_ll',eta_ll,dsig)
      call filld('y_ll',y_ll,dsig)
      call filld('m_ll',mass_ll,dsig)
      call filld('pt_ll', pt_ll,dsig)
      call filld('pt_miss',pt_miss,dsig)
      call filld('phi_ll',dphi_ll,dsig)
c     -- transverse momenta and masses of Ws
      call filld('mT_WW',MT_WW,dsig)


c -- forward-backward asymm of charged leptons, as a check
c -- RR: this needs to be looked at.
      fblminus = sign(1d0,eta_lminus)*1.5d0
      fblplus = sign(1d0,eta_lplus)*1.5d0
      
      call filld('FB_asym_lm', fblminus,dsig)
      call filld('FB_asym_lp', fblplus,dsig)


c --- jet distributions

      if(mjets.ge.1) then
         
c     call getyetaptmass(pj(:,1),y,eta_j,pt_j,mass)
         call filld ('pt_jnocut',ktj(1),dsig)

c     -- jet cuts: pt > 20, |eta| < 3.5           
         if ((ktj(1) .le. 20d0) .or. 
     .        (abs(etaj(1)) .ge. 3.5d0)) then
            return
         endif
         call getdydetadphidr(plminus,pj(:,1),dy,deta,
     .        dphi,dr_jlmi)
         call getdydetadphidr(plplus,pj(:,1),dy,deta,
     .        dphi,dr_jlpl)

         call filld('eta_j',etaj(1),dsig)
         call filld('pt_j',ktj(1),dsig)
         call filld('R_lpj',dr_jlpl,dsig)
         call filld('R_lmj',dr_jlmi,dsig)

      endif

      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(*),y,eta,pt,mass,pv
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      eta=0.5d0*log((pv+p(3))/(pv-p(3)))
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
c      ptmin=20
      ptmin=0.1d0
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


      subroutine boostx(p_in,pt,ptt,p_out)
      implicit none
c--- Boost input vector p_in to output vector p_out using the same
c--- transformation as required to boost massive vector pt to ptt
      double precision p_in(4),pt(4),ptt(4),p_out(4),
     . p_tmp(4),beta(3),mass,gam,bdotp
      integer j
    
      mass=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2  
      if (mass .lt. 0d0) then
        write(6,*) 'mass**2 .lt. 0 in boostx.f, mass**2=',mass,pt
        stop
      endif
      mass=dsqrt(mass)

c--- boost to the rest frame of pt
      gam=pt(4)/mass

      bdotp=0d0
      do j=1,3
        beta(j)=-pt(j)/pt(4)
        bdotp=bdotp+beta(j)*p_in(j)
      enddo
      p_tmp(4)=gam*(p_in(4)+bdotp)
      do j=1,3
        p_tmp(j)=p_in(j)+gam*beta(j)/(1d0+gam)*(p_in(4)+p_tmp(4))
      enddo     

c--- boost from rest frame of pt to frame in which pt is identical
c--- with ptt, thus completing the transformation          
      gam=ptt(4)/mass

      bdotp=0d0
      do j=1,3
        beta(j)=+ptt(j)/ptt(4)
        bdotp=bdotp+beta(j)*p_tmp(j)
      enddo
      p_out(4)=gam*(p_tmp(4)+bdotp)
      do j=1,3
        p_out(j)=p_tmp(j)+gam*beta(j)/(1d0+gam)*(p_out(4)+p_tmp(4))
      enddo

      return
      end
      

      
