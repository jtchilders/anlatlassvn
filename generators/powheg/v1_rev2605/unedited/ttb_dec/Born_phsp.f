      subroutine born_phsp(xxborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
C      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'constants.f'
      include 'masses.f'
      include 'zerowidth.f'
      include 'process.f'
      include 'phasemin.f'
      include 'mxdim.f'
c      include 'cvecbos.h' 
c      include 'vvsettings.f' 
      real * 8 xxborn(ndiminteg-3)
      real * 8 xborn(mxdim)
      real * 8 xjac,tau,beta,vec(3)

      double precision q(mxpart,4),sqrts
      double precision wt6
      real * 8 xx(2)
      common/x1x2/xx
      integer i,k
      real * 8 powheginput
      external powheginput
      logical debug,ini
      common/energy/sqrts
      data ini/.true./
      data debug/.false./

      double precision lntaum,ymax,ycm

      integer ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10,ii11,ii12,ii13
      common/inputprocind/ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,
     1     ii9,ii10,ii11,ii12,ii13

      xborn(ndiminteg-2:mxdim)=0
      xborn(1:ndiminteg-3)=xxborn
C----set up parameters needed by gen6
      sqrts = sqrt(kn_sbeams)

      zerowidth=.false.
      taumin=4d0*mt**2/sqrts**2
      case='tt_bbl'
C----endset up parameters needed by gen6
      
c      write(6,*) 'xborn',xborn
c      write(6,*) 'mt',mt
c      write(6,*) 'mb',mb
c      write(6,*) 'sqrts',sqrts
      call gen6(xborn,q,wt6,*99)
      if (debug)  write(*,*) 'Entering Born_phsp: ndiminteg', ndiminteg

      kn_jacborn = wt6


      kn_pborn(0,ii1)= -q(1,4)
      kn_pborn(1:3,ii1)= -q(1,1:3)

      kn_pborn(0,ii2)= -q(2,4)
      kn_pborn(1:3,ii2)= -q(2,1:3)

      kn_pborn(0,ii3)= q(3,4)+q(4,4)+q(5,4)
      kn_pborn(1:3,ii3)= q(3,1:3)+q(4,1:3)+q(5,1:3)

      kn_pborn(0,ii4)= q(6,4)+q(7,4)+q(8,4)
      kn_pborn(1:3,ii4)= q(6,1:3)+q(7,1:3)+q(8,1:3)

      kn_pborn(0,ii5)= q(3,4)+q(4,4)
      kn_pborn(1:3,ii5)= q(3,1:3)+q(4,1:3)

      kn_pborn(0,ii6)= q(7,4)+q(8,4)
      kn_pborn(1:3,ii6)= q(7,1:3)+q(8,1:3)

      kn_pborn(0,ii7)= q(4,4)
      kn_pborn(1:3,ii7)= q(4,1:3)

      kn_pborn(0,ii8)= q(3,4)
      kn_pborn(1:3,ii8)= q(3,1:3)

      kn_pborn(0,ii9)= q(7,4)
      kn_pborn(1:3,ii9)= q(7,1:3)

      kn_pborn(0,ii10)= q(8,4)
      kn_pborn(1:3,ii10)= q(8,1:3)

      kn_pborn(0,ii11)= q(5,4)
      kn_pborn(1:3,ii11)= q(5,1:3)

      kn_pborn(0,ii12)= q(6,4)
      kn_pborn(1:3,ii12)= q(6,1:3)

      kn_xb1 = xx(1)
      kn_xb2 = xx(2)
      beta=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=-1      
      call mboost(nlegborn,vec,beta,kn_pborn,kn_cmpborn)
      kn_sborn=(2*kn_cmpborn(0,1))**2


      return
 99   continue
      kn_jacborn=0d0
      return
      end



      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real * 8 fact,pt2supp,powheginput,pt
      save ini,pt2supp,pt     
      if (ini) then
         pt = powheginput("#ptsupp")         
         ini = .false.
         pt2supp = pt**2
      endif
      if (pt.lt.0) then
         fact=1d0
      else         
         fact=1d0
c CAVEAT!!!   No suppression
c         pt2=kn_pborn(1,5)**2+kn_pborn(2,5)**2
c         fact=pt2/(pt2+pt2supp)         
c      if (pt2.gt.10) then
c         fact = 1d0
c      else
c         fact = 0d0
c      endif
      endif
      end


      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 muf,mur
      logical ini
      data ini/.true./
      real *8 dotp,powheginput
      external dotp,powheginput
      logical fixedscale
      save fixedscale,ini
      if (ini) then
         if(powheginput("#fixedscale").lt.0) then
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization'
            write(*,*) '    Scales mur=muf=sqrt(mt^2+pttop^2)'
            write(*,*) '*************************************'
            fixedscale=.false.
         else
c            write(*,*) '*************************************'
c            write(*,*) '    Factorization and renormalization'
c            write(*,*) '    Scales mur=muf=mt)               '
c            write(*,*) '*************************************'
c--- for comparison with MCFM
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization'
            write(*,*) '    Scales mur=muf=mZ)               '
            write(*,*) '*************************************'
c--- for comparison with MCFM
            fixedscale=.true.
         endif
         ini=.false.
      endif
      if(fixedscale) then
         muf = ph_tmass
         mur = ph_tmass
c         muf = ph_zmass
c         mur = ph_zmass
      else
         muf=sqrt(kn_cmpborn(0,3)**2-kn_cmpborn(3,3)**2)
         mur=muf
      endif
      end

