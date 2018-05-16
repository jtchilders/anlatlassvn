c -*- Fortran -*-
      real * 8 physpar_ml(3)
      real * 8 physpar_mq(6)
      real * 8 physpar_aem
      common/pwhg_physpar/physpar_ml,physpar_mq,physpar_aem
      save /pwhg_physpar/

      complex*16 mw,mz,mh,mh2,mw2,mz2
      common/masses/mw,mz,mh,mw2,mz2,mh2
*
      real*8 me,mm,mtl,me2,mm2,mtl2,
     +                     mu,md,mu2,md2,
     +                     mc,ms,mc2,ms2,
     +                     mt,mb,mt2,mb2
      common/fermionmasses/me,mm,mtl,me2,mm2,mtl2,
     +                     mu,md,mu2,md2,
     +                     mc,ms,mc2,ms2,
     +                     mt,mb,mt2,mb2
*
      complex*16 gn(0:1),gl(0:1),gu(0:1),gd(0:1)
      common/vectorassial/gn,gl,gu,gd
*
      real*8 qu,qd,ql,qq
      common/charges/qu,qd,ql,qq
*
      complex*16 sw,cw,sw2,cw2,sw4,cw4,alpha,el2,alsu4pi,el2_scheme
      common/couplings/sw,sw2,sw4,cw,cw2,cw4,alpha,el2,alsu4pi,
     +                 el2_scheme
*
      complex*16 i3q,i3l,gq(0:1)
      common/updown/i3q,i3l,gq
*
      complex*16 a(0:1,0:1),chiz
      common/helic/a,chiz
*
