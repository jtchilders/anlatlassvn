      real * 8 function getdeltaew(ampl1loop)
      implicit none
      include 'PhysPars.h'
      include 'pwhg_math.h'
      include 'Flags.h'
      integer ewbin
      complex *  16 ampl1loop
      real * 8 deltaew
      external deltaew

      if (flg_fast_ew.eq.1) then
         ewbin = nint(mh)
         getdeltaew = cached_ew_corr(ewbin)
      else
         getdeltaew = deltaew(ampl1loop)
      end if

      end function getdeltaew


      real * 8 function deltaew(ampl1loop)
      implicit none
      include 'PhysPars.h'
      include 'pwhg_math.h'
      real * 8 tmp
      complex * 16 ampl1loop,ewampl
      complex * 16 amplew

!      write(*,*) 'mh', mh
!      write(*,*) 'ampl1loop', ampl1loop
c$$$      write(*,*) '2*ampl1loop', 2d0*ampl1loop
            
      ewampl = amplew()
c$$$      tmp =   abs(ampl1loop)**2+
c$$$     -     2d0*dble(ampl1loop) *
c$$$     -     dble(ewampl)  +
c$$$     -     2d0*dimag(ampl1loop) *
c$$$     -     dimag(ewampl)

      tmp =
     -     2d0*dble(ampl1loop) *
     -     dble(ewampl)  +
     -     2d0*dimag(ampl1loop) *
     -     dimag(ewampl)

!      deltaewfull = (tmp+abs(ewampl)**2)/abs(ampl1loop)**2
      deltaew = tmp/abs(ampl1loop)**2
!      write(*,*) 'ewampl******',ewampl
!      write(*,*)  'deltaew',deltaew
!      write(*,*) 'deltaewfull',deltaewfull
!      write(*,*) 'tmp',tmp
!      write(*,*) 'tmpfull',tmpfull

      end function deltaew


c     Used from here.
      complex * 16 function amplew()
      implicit none
      include 'PhysPars.h'
      include 'pwhg_math.h'
      include 'Flags.h'
      complex * 16 amp_glu_fus_2lew_lf,cmw,norm
      external amp_glu_fus_2lew_lf
c     alpha(mz)

!      write(*,*) 'ph_Zmass',  ph_Zmass
!      write(*,*) 'ph_Wmass',  ph_Wmass
!      write(*,*) 'ph_alphaem', ph_alphaem
!      write(*,*) 'ph_alphammz', ph_alphaemmz
!      write(*,*) 'ph_Zwidth', ph_Zwidth
!      write(*,*) 'pw_Wwidth', ph_Wwidth
!      write(*,*) 'ph_sthw2', ph_sthw2

      cmw = ph_Wmass-dcmplx(0d0,1d0)*ph_Wwidth/2d0

      norm = 1/4d0/pi/ph_sthw2*
     $     cmw*cmw/mh/mh

c      write(*,*) 'pi', pi

c     MSSM
c        Light Higgs
      if (ih.eq.1) then
c              sin(beta-alpha) = sin(beta)cos(alpha)-cos(beta)sin(alpha)
        norm = ph_alphaem*(ph_sb*cos(ph_alpha)-ph_cb*sin(ph_alpha))
     $*norm
c        Heavy Higgs
      else
c              cos(beta-alpha) = cos(beta)cos(alpha)+sin(beta)sin(alpha)
        norm = ph_alphaem*(ph_cb*cos(ph_alpha)+ph_sb*sin(ph_alpha))
     $*norm
      endif

!      write(*,*) 'fac', ph_sb*cos(ph_alpha)-ph_cb*sin(ph_alpha)

      amplew = norm*amp_glu_fus_2lew_lf()

      end function amplew

      complex * 16 function amp_glu_fus_2lew_lf()
      implicit none
      include 'PhysPars.h'
      complex * 16 A1, cmw,cmz,zh,wh
      external A1

c$$$      cmw = ph_Wmass-dcmplx(0d0,1d0)*ph_Wwidth/2d0
c$$$      cmz = ph_Zmass-dcmplx(0d0,1d0)*ph_Zwidth/2d0
c$$$
c$$$      zh=cmz*cmz/mh/mh
c$$$      wh=cmw*cmw/mh/mh
c$$$
c$$$      amp_glu_fus_2lew_lf=
c$$$     $ 2d0/ph_cthw2/ph_cthw2*(5d0/4d0-7d0/3d0*ph_sthw2+22d0/9d0*
c$$$     $ ph_sthw2*ph_sthw2)*A1(zh)+4d0*A1(wh)


c      write(*,*) 'pw_Zwidth', ph_Zwidth

      cmw = ph_Wmass-dcmplx(0d0,1d0)*ph_Wwidth/2d0
      cmz = ph_Zmass-dcmplx(0d0,1d0)*ph_Zwidth/2d0
      ph_cthw2 = 1d0-ph_sthw2
      zh=cmz*cmz/mh/mh
      wh=cmw*cmw/mh/mh

c$$$      write(*,*) 'cmw', cmw
c$$$      write(*,*) 'cmz', cmz
c$$$      write(*,*) 's2', ph_sthw2
c$$$      write(*,*) 'c2', ph_cthw2
c$$$      write(*,*) 'zh', zh
c$$$      write(*,*) 'wh', wh

      amp_glu_fus_2lew_lf=
     -  2d0/ph_cthw2/ph_cthw2*(5d0/4d0-7d0/3d0*ph_sthw2+22d0/9d0*
     -  ph_sthw2*ph_sthw2)*A1(zh)+
     -  4d0*A1(wh)

c      write(*,*) '***** INNER AMPLS ******'
c      write(*,*) wh, A1(wh)
c      write(*,*) zh, A1(zh)
c      write(*,*) '*****  ******'
c      stop
      end function amp_glu_fus_2lew_lf

      double complex function A1(x)
      implicit none

      double precision IEPS
      double complex H1,H2,H3,x
      external H1,H2,H3

      IEPS = 1d-9

      A1 = -4d0 -
     -  2d0* (-1d0+x)* H1(-1,-1d0/x) -
     -       2d0* x* H2(0,-1,-1d0/x) -
     -  2d0* (-1d0+3d0*x)* H3(0,0,-1,-1d0/x) -
     -  2d0* (-1d0+2d0* x)* H3(0,-2,-2,-1d0/x) + 
     -  3d0* (-1d0+2d0* x)* H3(-2,-2,-1,-1d0/x) +
c     -   dcmplx(0d0,1d0)*sqrt(4d0*x-dcmplx(0d0,IEPS)-1d0)*
     -   sqrt(1d0-4d0*x)*
     -      ( -2d0*  H1(-2,-1d0/x)  -
     -              3d0* (-1d0+2d0* x)* H3(-4,-2,-1,-1d0/x) +
     -        2d0* (-1d0+2d0* x)* H3(-2, 0,-1,-1d0/x) +
     -        2d0* (-1d0+2d0* x)* H3(-2,-2,-2,-1d0/x)   )

c       write(*,*) 
c       write(*,*) x
c       write(*,*) 'h1 1',H1(-1,-1d0/x)
c       write(*,*) 'h2 1',H2(0,-1,-1d0/x)
c       write(*,*) 'h3 2',H3(0,0,-1,-1d0/x)
c       write(*,*) 'h3 3',H3(0,-2,-2,-1d0/x)
c       write(*,*) 'h3 4',H3(-2,-2,-1,-1d0/x)
c       write(*,*) 'h1 2',H1(-2,-1d0/x)
c       write(*,*) 'h3 5',H3(-4,-2,-1,-1d0/x)
c       write(*,*) 'h3 6',H3(-2, 0,-1,-1d0/x)
c       write(*,*) 'h3 7',H3(-2,-2,-2,-1d0/x)
c       write(*,*) 
c       stop

      end function A1
