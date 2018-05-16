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
      integer j,k
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/

      call inihists

      call bookupeqbins('total',1d0,-0.5d0,0.5d0)

c      call bookupeqbins('Nphot',1d0,-0.5d0,5.5d0)

      call bookupeqbins('V_mt',0.5d0,50d0,100d0)
      call bookupeqbins('V_pt',0.25d0,0d0,25d0)
      call bookupeqbins('V_pt2',1d0,0d0,350d0)

      call bookupeqbins('V_m',2.4d0,60d0,120d0)

      call bookupeqbins('forward',2.4d0,60d0,120d0)
      call bookupeqbins('backward',2.4d0,60d0,120d0)

      call bookupeqbins('X_m',0.006d0,0.6d0,1.2d0)

      call bookupeqbins('l1_y',0.2d0,-4d0,4d0)
      call bookupeqbins('l1_eta',0.2d0,-2.5d0,2.5d0)
      call bookupeqbins('l1_pt',1.4d0,25d0,60d0)

      call bookupeqbins('l2_y',0.2d0,-4d0,4d0)
      call bookupeqbins('l2_eta',0.2d0,-2.5d0,2.5d0)
      call bookupeqbins('l2_pt',1.4d0,25d0,60d0)

      call bookupeqbins('X_p',0.008d0,0.4d0,1.2d0)

      call bookupeqbins('delta_phi',0.1d0,0d0,3.2d0)

      end
     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      integer   maxphot,nphot
      parameter (maxphot=2048)
      real * 8 pg(4,maxphot)
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/
      save cnum
      integer j,k
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      real * 8 pz(4),pl1(4),pl2(4)
      real * 8 y1,eta1,ptl1,m1
      real * 8 y2,eta2,ptl2,m2
      real * 8 pt,m,mtv
      real * 8 dy,deta,delphi,dr
      real * 8 getpt,getdphi,getmass,geteta
      external getpt,getdphi,getmass,geteta
      integer ihep
      real * 8 powheginput,dotp
      external powheginput,dotp
      integer vdecaytemp
      logical accepted
      real*8 mw,mz
      real*8 cs
      real*8 cstar
      external cstar

      if(dsig.eq.0) return

      mw=80.385d0
      mz=91.1876d0

      pz = (/0,0,0,0/)
      pl1= (/0,0,0,0/)
      pl2= (/0,0,0,0/)
      nphot = 0

      vdecaytemp=lprup(1)-10000 ! Z decay product, with positive id

      do ihep=1,nhep
c p_Z = p_l1 + p_l2
         if( idhep(ihep).eq.vdecaytemp  ) then
             if (phep(4,ihep).gt.pl1(4)) pl1 = phep(1:4,ihep)
         endif
         if( idhep(ihep).eq.-vdecaytemp ) then
             if (phep(4,ihep).gt.pl2(4)) pl2 = phep(1:4,ihep)
         endif
         pz = pl1 + pl2
         if( idhep(ihep).eq.22 ) then 
             if (phep(4,ihep).gt.10d0)then
                 nphot = nphot + 1
                 pg(1:4,nphot) = phep(1:4,ihep)
             endif
         endif
      enddo


      call getyetaptmass(pl1,y1,eta1,ptl1,m1)
      call getyetaptmass(pl2,y2,eta2,ptl2,m2)
      delphi = getdphi(pl1,pl2)
      pt=getpt(pz)
      m=getmass(pz)
      mtv = sqrt(2*ptl1*ptl2*(1d0-cos(delphi)))

      cs = cstar(pl1,pl2)

      call filld('total',0d0,dsig)

c      call filld('Nphot',dble(nphot),dsig)

c lepton 1
      call filld('l1_y',    y1, dsig)
      call filld('l1_eta',eta1, dsig)
      call filld('l1_pt', ptl1, dsig)

c lepton 2
      call filld('l2_y',    y2, dsig)
      call filld('l2_eta',eta2, dsig)
      call filld('l2_pt', ptl2, dsig)


c azimuthal separation betwen lepton and neutrino
      call filld('delta_phi',delphi,dsig)

c W
      call filld('V_pt',pt, dsig)
      call filld('V_pt2',pt, dsig)
      call filld('V_m',  m, dsig)
c transverse mass of the lepton-neutrino system
      call filld('V_mt',mtv,dsig)

      call filld('X_m',mtv/mz,dsig)
      call filld('X_p',ptl1*2d0/mz,dsig)

      if (cs.lt.0d0) then
          call filld('backward', m, dsig)
      else
          call filld('forward', m, dsig)
      endif


      end

      subroutine yetaptmassplot(p,dsig,prefix)
      implicit none
      real * 8 p(4),dsig
      character *(*) prefix
      real * 8 y,eta,pt,m
      call getyetaptmass(p,y,eta,pt,m)
      call filld(prefix//'_y',y,dsig)
      call filld(prefix//'_eta',eta,dsig)
      call filld(prefix//'_pt',pt,dsig)
      call filld(prefix//'_m',m,dsig)
      end

      subroutine deltaplot(p1,p2,dsig,prefix)
      implicit none
      real * 8 p1(4),p2(4),dsig
      character *(*) prefix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'_dy',dy,dsig)
      call filld(prefix//'_deta',deta,dsig)
      call filld(prefix//'_delphi',delphi,dsig)
      call filld(prefix//'_dr',dr,dsig)
      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real * 8 gety,getpt,geteta,getmass
      external gety,getpt,geteta,getmass
      y  = gety(p)
      pt = getpt(p)
      eta = geteta(p)
      mass = getmass(p)
      end


      function gety(p)
      implicit none
      real * 8 gety,p(4)
      gety=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function getpt(p)
      implicit none
      real * 8 getpt,p(4)
      getpt = sqrt(p(1)**2+p(2)**2)
      end

      function getmass(p)
      implicit none
      real * 8 getmass,p(4)
      getmass=sqrt(abs(p(4)**2-p(3)**2-p(2)**2-p(1)**2))
      end

      function geteta(p)
      implicit none
      real * 8 geteta,p(4),pv
      real * 8 tiny
      parameter (tiny=1.d-5)
      pv = sqrt(p(1)**2+p(2)**2+p(3)**2)
      if(pv.lt.tiny)then
         geteta=sign(1.d0,p(3))*1.d8
      else
         geteta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      end



      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 getdy,getdeta,getdphi,getdr
      external getdy,getdeta,getdphi,getdr
      dy=getdy(p1,p2)
      deta=getdeta(p1,p2)
      dphi=getdphi(p1,p2)
      dr=getdr(deta,dphi)
      end

      function getdy(p1,p2)
      implicit none
      real*8 p1(*),p2(*),getdy
      real*8 y1,y2
      real*8 gety
      external gety
      y1 = gety(p1)
      y2 = gety(p2)
      getdy = y1-y2
      end

      function getdeta(p1,p2)
      implicit none
      real*8 p1(*),p2(*),getdeta
      real*8 eta1,eta2
      real*8 geteta
      external geteta
      eta1 = geteta(p1)
      eta2 = geteta(p2)
      getdeta = eta1-eta2
      end

      function getdphi(p1,p2)
      implicit none
      include 'pwhg_math.h' 
      real*8 p1(*),p2(*),getdphi
      real*8 phi1,phi2
      real*8 geteta
      external geteta
      phi1=atan2(p1(2),p1(1))
      phi2=atan2(p2(2),p2(1))
      getdphi=abs(phi1-phi2)
      getdphi=min(getdphi,2d0*pi-getdphi)
      end

      function getdr(deta,dphi)
      implicit none
      real*8 getdr,deta,dphi 
      getdr=sqrt(deta**2+dphi**2)
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

*
**
*
      real*8 function cstar(p1,p2)
      implicit none
      real*8 p1(4),p2(4),psum(4)
*
      real*8 dotp
      external dotp
*
      real*8 rq2,cs1p,cs2p,cs1m,cs2m,qmass,pt2,sig
      integer k
*
      psum = p1 + p2
      rq2 = sqrt(2.d0)
! Collins - Soper momenta for particle 1 and 2 
      cs1p = (p1(4) + p1(3))/rq2
      cs2p = (p2(4) + p2(3))/rq2
      cs1m = (p1(4) - p1(3))/rq2
      cs2m = (p2(4) - p2(3))/rq2
      qmass = sqrt(psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2)
      pt2 = psum(1)**2 + psum(2)**2
      cstar = 2.d0/qmass/sqrt(qmass**2 + pt2)*(cs1p*cs2m - cs1m*cs2p)
! for a ppbar should end here
c      if (hadr1.eq.hadr2) then
         sig = 1.d0
         if (psum(3).ne.0.d0) sig = abs(psum(3))/psum(3)
         cstar = cstar * sig
c      endif
      return
      end
