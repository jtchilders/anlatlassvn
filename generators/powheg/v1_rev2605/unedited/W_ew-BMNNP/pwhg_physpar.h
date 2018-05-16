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
      complex*16 gnm,gnp,glm,glp,gum,gup,gdm,gdp
      common/vectorassial/gnm,gnp,glm,glp,gum,gup,gdm,gdp
*
      real*8 qu,qd,ql,qnu
      common/charges/qu,qd,qnu,ql
*
      complex*16 sw,cw,sw2,cw2,sw4,cw4,alpha,el2,alsu4pi,g2
      common/couplings/sw,sw2,sw4,cw,cw2,cw4,alpha,el2,alsu4pi,g2

