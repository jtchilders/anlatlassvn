      subroutine fillPhyspars(esq,Gf,alphaem,
     & hmass, wmass, zmass, tmass, bmass
     & hwidth,wwidth,zwidth,twidth)
      implicit none
      include 'PhysPars.h' 
      double precision esq,Gf,alphaem,
     & hmass, wmass, zmass, tmass, bmass
     & hwidth,wwidth,zwidth,twidth
      double precision pi
      pi=2d0*asin(1d0)
      ph_unit_e=sqrt(esq)
      ph_CKM(:,:)=0d0
      ph_CKM(1,1)=1d0
      ph_CKM(2,2)=1d0
      ph_CKM(3,3)=1d0
      ph_gfermi=Gf
      ph_alphaem=esq/(4d0*pi)
      ph_Zmass=zmass
      ph_Zwidth=zwidth
      ph_Wmass=wmass
      ph_Wwidth=wwidth
      ph_cthw=ph_Wmass/ph_Zmass
      ph_sthw2=1d0-ph_cthw**2
      ph_sthw=sqrt(ph_sthw2)
      ph_Zmass2=ph_Zmass**2
      ph_Zmass2low=0d0
      ph_Zmass2high=0d0
      ph_Wmass2=ph_Wmass**2
      ph_Wmass2low=0d0
      ph_Wmass2high=0d0
      ph_ZmZw=ph_Zmass*ph_Zwidth
      ph_WmWw=ph_Wmass*ph_Wwidth
      ph_gfermi=Gf
      ph_alphaem=esq/(4d0*pi)

      write(6,*) 'fillPhyspars: ph_Zmass,ph_Zwidth,ph_Wmass,ph_Wwidth',
     & ph_Zmass,ph_Zwidth,ph_Wmass,ph_Wwidth
      write(6,*) 'fillPhyspars: ph_gfermi,1d0/ph_alphaem',
     & ph_gfermi,1d0/ph_alphaem

      return
      end
