      subroutine fillMCFMCommon(eesq,gfermi,
     &     mh, mw, mz, tmass,   bmass,
     &     width_h,width_w,width_z,width_t,
     &     st_alpha,st_muren2,st_mufac2)
      implicit none
      double precision
     &     eesq,gfermi,mh, mw, mz, tmass,   bmass,
     &     width_h,width_w,width_z,width_t,
     &     st_alpha,st_muren2,st_mufac2
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'scale.f'
      include 'scheme.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'lc.f'

      hmass=mh
      wmass=mw
      zmass=mz
      md=0d0
      mu=0d0
      ms=0d0
      mc=0d0
      mb=bmass
      mt=tmass
      hwidth=width_h
      wwidth=width_w
      zwidth=width_z
      twidth=width_t
      write(6,*) 'md=        ',md
      write(6,*) 'mu=        ',mu
      write(6,*) 'ms=        ',ms
      write(6,*) 'mc=        ',mc
      write(6,*) 'mb=        ',mb
      write(6,*) 'mt=        ',mt
      write(6,*) 'wmass=     ',wmass
      write(6,*) 'wwidth=    ',wwidth
      write(6,*) 'zmass=     ',zmass
      write(6,*) 'zwidth=    ',zwidth
      Gf=gfermi
      gwsq=Gf/sqrt(2d0)*(8d0*wmass**2)
      gw=sqrt(gwsq)
      esq=eesq
      xw=1d0-(wmass/zmass)**2
      vevsq=(2d0*wmass/gw)**2
      write(6,*) 'Gf=        ',Gf
      write(6,*) 'gw=         ',gw
      write(6,*) 'gwsq=       ',gwsq
      write(6,*) 'xw=         ',xw
      write(6,*) 'esq=       ',esq
      write(6,*) 'vevsq=     ',vevsq

      call couplz(xw)

      epinv=0d0
      epinv2=0d0

      colourchoice=0

      scheme='msbr'

      write(6,*) 
      return
      end


      block data wsalam1
      implicit none
      include 'constants.f'
      include 'ewcharge.f'
      data Q(-5)/+0.333333333333333d0/
      data Q(-4)/-0.666666666666667d0/
      data Q(-3)/+0.333333333333333d0/
      data Q(-2)/-0.666666666666667d0/
      data Q(-1)/+0.333333333333333d0/
      data Q(0)/+0d0/
      data Q(+1)/-0.333333333333333d0/
      data Q(+2)/+0.666666666666667d0/
      data Q(+3)/-0.333333333333333d0/
      data Q(+4)/+0.666666666666667d0/
      data Q(+5)/-0.333333333333333d0/
      data tau/1d0,-1d0,1d0,-1d0,1d0,0d0,-1d0,1d0,-1d0,1d0,-1d0/
      end 

************************************************************************
*     CKM matrix entries                                               *
************************************************************************
      block data block_ckm
      implicit none
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      data  Vud  ,  Vus  ,  Vub  ,
     .      Vcd  ,  Vcs  ,  Vcb
     .   /1d0,0d0,0d0,
     .    0d0,1d0,0d0/
      end
