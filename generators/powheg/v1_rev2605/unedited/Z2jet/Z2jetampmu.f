c--- File written by FORM program Z2jet.frm on Sat Jul 28 09:41:52 CDT 2012
      subroutine Z2jetampmu(p1,p2,p3,p4,p5,p6,f56,f65,fqed)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex zab(mxpart,4,mxpart)
      common/zamub/zab
      integer p1,p2,p3,p4,p5,p6,mu
      double precision s34,s56,s134,s234,s15,s25,s16,s26
      double complex iza,izb,f56(4,2,2,2),f65(4,2,2,2),fqed(4,2,2,2)
C     first index of f is mu
C     second index of f  is helicity gluon 6,
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
      s134=s(p1,p3)+s(p1,p4)+s(p3,p4)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s15=s(p1,p5)
      s16=s(p1,p6)
      s25=s(p2,p5)
      s26=s(p2,p6)
      do mu=1,4
      f56(mu,1,1,1)= + zab(p2,mu,p1)*s34**(-1) * (  - za(p1,p3)*za(p5,
     &    p6)*zb(p1,p4)*zb(p1,p5)*izb(p1,p6)*s134**(-1)*s56**(-1) + za(
     &    p2,p3)*za(p5,p6)*zb(p1,p5)*zb(p2,p4)*izb(p1,p6)*s234**(-1)*
     &    s56**(-1) )
      f56(mu,1,1,1) = f56(mu,1,1,1) + zab(p2,mu,p2)*s34**(-1) * ( za(p2
     &    ,p6)*za(p3,p4)*zb(p1,p4)**2*izb(p1,p6)*s25**(-1)*s134**(-1) )
      f56(mu,1,1,1) = f56(mu,1,1,1) + zab(p2,mu,p4)*s34**(-1) * ( za(p3
     &    ,p4)*za(p5,p6)*zb(p1,p4)*zb(p1,p5)*izb(p1,p6)*s134**(-1)*
     &    s56**(-1) )
      f56(mu,1,1,1) = f56(mu,1,1,1) + zab(p2,mu,p5)*s34**(-1) * ( za(p3
     &    ,p4)*za(p5,p6)*zb(p1,p4)**2*izb(p1,p6)*s25**(-1)*s134**(-1) )
      f56(mu,1,1,1) = f56(mu,1,1,1) + zab(p3,mu,p1)*s34**(-1) * ( za(p2
     &    ,p3)*za(p5,p6)*zb(p1,p5)*zb(p3,p4)*izb(p1,p6)*s234**(-1)*
     &    s56**(-1) )
      f56(mu,1,1,1) = f56(mu,1,1,1) + zab(p6,mu,p1)*s34**(-1) * ( za(p1
     &    ,p3)*za(p2,p6)*zb(p1,p4)*s134**(-1)*s56**(-1) - za(p2,p3)*za(
     &    p2,p6)*zb(p2,p4)*s234**(-1)*s56**(-1) - za(p2,p3)*za(p3,p6)*
     &    zb(p3,p4)*s234**(-1)*s56**(-1) - za(p2,p6)*za(p3,p4)*zb(p1,p4
     &    )*zb(p4,p6)*izb(p1,p6)*s134**(-1)*s56**(-1) )
      f56(mu,1,1,1) = f56(mu,1,1,1) + zab(p6,mu,p6)*s34**(-1) * (  - 
     &    za(p2,p6)*za(p3,p4)*zb(p1,p4)**2*izb(p1,p6)*s134**(-1)*
     &    s56**(-1) )

      f56(mu,2,1,1)= + zab(p2,mu,p1)*s34**(-1) * (  - za(p1,p5)*za(p2,
     &    p3)*zb(p1,p6)*zb(p2,p4)*iza(p5,p6)*s234**(-1)*s16**(-1) )
      f56(mu,2,1,1) = f56(mu,2,1,1) + zab(p2,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p2,p5)*zb(p1,p4)*zb(p1,p6)*iza(p5,p6)*s25**(-1)*
     &    s134**(-1) - za(p1,p5)*za(p2,p3)*zb(p1,p4)*zb(p1,p6)*iza(p5,
     &    p6)*s16**(-1)*s25**(-1) - za(p2,p3)*zb(p1,p6)*zb(p4,p6)*
     &    s16**(-1)*s25**(-1) + za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)
     &    *iza(p5,p6)*s25**(-1)*s134**(-1) )
      f56(mu,2,1,1) = f56(mu,2,1,1) + zab(p2,mu,p5)*s34**(-1) * ( za(p1
     &    ,p5)*za(p3,p5)*zb(p1,p4)*zb(p1,p6)*iza(p5,p6)*s16**(-1)*
     &    s25**(-1) + za(p3,p5)*zb(p1,p6)*zb(p4,p6)*s16**(-1)*s25**(-1)
     &     )
      f56(mu,2,1,1) = f56(mu,2,1,1) + zab(p2,mu,p6)*s34**(-1) * ( za(p2
     &    ,p3)*zb(p1,p6)*zb(p2,p4)*s234**(-1)*s16**(-1) )
      f56(mu,2,1,1) = f56(mu,2,1,1) + zab(p3,mu,p1)*s34**(-1) * (  - 
     &    za(p1,p5)*za(p2,p3)*zb(p1,p6)*zb(p3,p4)*iza(p5,p6)*s234**(-1)
     &    *s16**(-1) )
      f56(mu,2,1,1) = f56(mu,2,1,1) + zab(p3,mu,p6)*s34**(-1) * ( za(p2
     &    ,p3)*zb(p1,p6)*zb(p3,p4)*s234**(-1)*s16**(-1) )
      f56(mu,2,1,1) = f56(mu,2,1,1) + zab(p5,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p2,p6)*zb(p1,p4)*zb(p1,p6)*iza(p5,p6)*s134**(-1)
     &    *s56**(-1) + za(p2,p3)*za(p2,p6)*zb(p1,p6)*zb(p2,p4)*iza(p5,
     &    p6)*s234**(-1)*s56**(-1) + za(p2,p3)*za(p3,p6)*zb(p1,p6)*zb(
     &    p3,p4)*iza(p5,p6)*s234**(-1)*s56**(-1) + za(p2,p6)*za(p3,p4)*
     &    zb(p1,p4)*zb(p4,p6)*iza(p5,p6)*s134**(-1)*s56**(-1) )
      f56(mu,2,1,1) = f56(mu,2,1,1) + zab(p6,mu,p6)*s34**(-1) * ( za(p1
     &    ,p3)*za(p2,p5)*zb(p1,p4)*zb(p1,p6)*iza(p5,p6)*s134**(-1)*
     &    s56**(-1) - za(p2,p3)*za(p2,p5)*zb(p1,p6)*zb(p2,p4)*iza(p5,p6
     &    )*s234**(-1)*s56**(-1) - za(p2,p3)*za(p3,p5)*zb(p1,p6)*zb(p3,
     &    p4)*iza(p5,p6)*s234**(-1)*s56**(-1) - za(p2,p5)*za(p3,p4)*zb(
     &    p1,p4)*zb(p4,p6)*iza(p5,p6)*s134**(-1)*s56**(-1) )

      f65(mu,1,1,1)= + zab(p1,mu,p1)*s34**(-1) * (  - za(p2,p3)*za(p2,
     &    p6)*zb(p1,p2)*zb(p1,p4)*izb(p1,p6)*s15**(-1)*s26**(-1) + za(
     &    p2,p6)*za(p3,p6)*zb(p1,p4)*s15**(-1)*s26**(-1) )
      f65(mu,1,1,1) = f65(mu,1,1,1) + zab(p2,mu,p1)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p2,p6)*zb(p1,p2)*zb(p1,p4)*izb(p1,p6)*s134**(-1)
     &    *s26**(-1) + za(p1,p3)*za(p5,p6)*zb(p1,p4)*zb(p1,p5)*izb(p1,
     &    p6)*s134**(-1)*s56**(-1) - za(p2,p3)*za(p5,p6)*zb(p1,p5)*zb(
     &    p2,p4)*izb(p1,p6)*s234**(-1)*s56**(-1) )
      f65(mu,1,1,1) = f65(mu,1,1,1) + zab(p2,mu,p4)*s34**(-1) * ( za(p2
     &    ,p6)*za(p3,p4)*zb(p1,p2)*zb(p1,p4)*izb(p1,p6)*s134**(-1)*
     &    s26**(-1) - za(p3,p4)*za(p5,p6)*zb(p1,p4)*zb(p1,p5)*izb(p1,p6
     &    )*s134**(-1)*s56**(-1) )
      f65(mu,1,1,1) = f65(mu,1,1,1) + zab(p3,mu,p1)*s34**(-1) * (  - 
     &    za(p2,p3)*za(p5,p6)*zb(p1,p5)*zb(p3,p4)*izb(p1,p6)*s234**(-1)
     &    *s56**(-1) )
      f65(mu,1,1,1) = f65(mu,1,1,1) + zab(p5,mu,p1)*s34**(-1) * ( za(p2
     &    ,p3)*za(p2,p6)*zb(p1,p2)*zb(p4,p5)*izb(p1,p6)*s15**(-1)*
     &    s26**(-1) - za(p2,p3)*za(p2,p6)*zb(p1,p5)*zb(p2,p4)*izb(p1,p6
     &    )*s234**(-1)*s15**(-1) - za(p2,p3)*za(p3,p6)*zb(p1,p5)*zb(p3,
     &    p4)*izb(p1,p6)*s234**(-1)*s15**(-1) - za(p2,p6)*za(p3,p6)*zb(
     &    p4,p5)*s15**(-1)*s26**(-1) )
      f65(mu,1,1,1) = f65(mu,1,1,1) + zab(p6,mu,p1)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p2,p6)*zb(p1,p4)*s134**(-1)*s56**(-1) - za(p1,p3
     &    )*za(p2,p6)*zb(p1,p4)*s134**(-1)*s26**(-1) + za(p2,p3)*za(p2,
     &    p6)*zb(p2,p4)*s234**(-1)*s56**(-1) + za(p2,p3)*za(p3,p6)*zb(
     &    p3,p4)*s234**(-1)*s56**(-1) + za(p2,p6)*za(p3,p4)*zb(p1,p4)*
     &    zb(p4,p6)*izb(p1,p6)*s134**(-1)*s56**(-1) )
      f65(mu,1,1,1) = f65(mu,1,1,1) + zab(p6,mu,p4)*s34**(-1) * ( za(p2
     &    ,p6)*za(p3,p4)*zb(p1,p4)*s134**(-1)*s26**(-1) )
      f65(mu,1,1,1) = f65(mu,1,1,1) + zab(p6,mu,p6)*s34**(-1) * ( za(p2
     &    ,p6)*za(p3,p4)*zb(p1,p4)**2*izb(p1,p6)*s134**(-1)*s56**(-1) )

      f65(mu,2,1,1)= + zab(p1,mu,p1)*s34**(-1) * (  - za(p2,p3)*za(p2,
     &    p5)*zb(p1,p4)*zb(p2,p6)*iza(p5,p6)*s15**(-1)*s26**(-1) - za(
     &    p2,p3)*za(p2,p5)*zb(p1,p6)*zb(p2,p4)*iza(p5,p6)*s234**(-1)*
     &    s15**(-1) - za(p2,p3)*za(p3,p5)*zb(p1,p6)*zb(p3,p4)*iza(p5,p6
     &    )*s234**(-1)*s15**(-1) )
      f65(mu,2,1,1) = f65(mu,2,1,1) + zab(p2,mu,p1)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p2,p5)*zb(p1,p4)*zb(p2,p6)*iza(p5,p6)*s134**(-1)
     &    *s26**(-1) )
      f65(mu,2,1,1) = f65(mu,2,1,1) + zab(p2,mu,p4)*s34**(-1) * ( za(p2
     &    ,p5)*za(p3,p4)*zb(p1,p4)*zb(p2,p6)*iza(p5,p6)*s134**(-1)*
     &    s26**(-1) )
      f65(mu,2,1,1) = f65(mu,2,1,1) + zab(p5,mu,p1)*s34**(-1) * (  - 
     &    za(p2,p3)*za(p2,p5)*zb(p2,p4)*zb(p5,p6)*iza(p5,p6)*s234**(-1)
     &    *s15**(-1) + za(p2,p3)*za(p2,p5)*zb(p2,p6)*zb(p4,p5)*iza(p5,
     &    p6)*s15**(-1)*s26**(-1) - za(p2,p3)*za(p3,p5)*zb(p3,p4)*zb(p5
     &    ,p6)*iza(p5,p6)*s234**(-1)*s15**(-1) )
      f65(mu,2,1,1) = f65(mu,2,1,1) + zab(p5,mu,p6)*s34**(-1) * ( za(p1
     &    ,p3)*za(p2,p6)*zb(p1,p4)*zb(p1,p6)*iza(p5,p6)*s134**(-1)*
     &    s56**(-1) - za(p2,p3)*za(p2,p6)*zb(p1,p6)*zb(p2,p4)*iza(p5,p6
     &    )*s234**(-1)*s56**(-1) - za(p2,p3)*za(p3,p6)*zb(p1,p6)*zb(p3,
     &    p4)*iza(p5,p6)*s234**(-1)*s56**(-1) - za(p2,p6)*za(p3,p4)*zb(
     &    p1,p4)*zb(p4,p6)*iza(p5,p6)*s134**(-1)*s56**(-1) )
      f65(mu,2,1,1) = f65(mu,2,1,1) + zab(p6,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p2,p5)*zb(p1,p4)*zb(p1,p6)*iza(p5,p6)*s134**(-1)
     &    *s56**(-1) + za(p2,p3)*za(p2,p5)*zb(p1,p6)*zb(p2,p4)*iza(p5,
     &    p6)*s234**(-1)*s56**(-1) + za(p2,p3)*za(p3,p5)*zb(p1,p6)*zb(
     &    p3,p4)*iza(p5,p6)*s234**(-1)*s56**(-1) + za(p2,p5)*za(p3,p4)*
     &    zb(p1,p4)*zb(p4,p6)*iza(p5,p6)*s134**(-1)*s56**(-1) )

      f56(mu,1,1,2)= + zab(p2,mu,p1)*s34**(-1) * (  - za(p1,p4)*za(p5,
     &    p6)*zb(p1,p3)*zb(p1,p5)*izb(p1,p6)*s134**(-1)*s56**(-1) + za(
     &    p2,p4)*za(p5,p6)*zb(p1,p5)*zb(p2,p3)*izb(p1,p6)*s234**(-1)*
     &    s56**(-1) )
      f56(mu,1,1,2) = f56(mu,1,1,2) + zab(p2,mu,p2)*s34**(-1) * (  - 
     &    za(p2,p6)*za(p3,p4)*zb(p1,p3)**2*izb(p1,p6)*s25**(-1)*
     &    s134**(-1) )
      f56(mu,1,1,2) = f56(mu,1,1,2) + zab(p2,mu,p3)*s34**(-1) * (  - 
     &    za(p3,p4)*za(p5,p6)*zb(p1,p3)*zb(p1,p5)*izb(p1,p6)*s134**(-1)
     &    *s56**(-1) )
      f56(mu,1,1,2) = f56(mu,1,1,2) + zab(p2,mu,p5)*s34**(-1) * (  - 
     &    za(p3,p4)*za(p5,p6)*zb(p1,p3)**2*izb(p1,p6)*s25**(-1)*
     &    s134**(-1) )
      f56(mu,1,1,2) = f56(mu,1,1,2) + zab(p4,mu,p1)*s34**(-1) * (  - 
     &    za(p2,p4)*za(p5,p6)*zb(p1,p5)*zb(p3,p4)*izb(p1,p6)*s234**(-1)
     &    *s56**(-1) )
      f56(mu,1,1,2) = f56(mu,1,1,2) + zab(p6,mu,p1)*s34**(-1) * ( za(p1
     &    ,p4)*za(p2,p6)*zb(p1,p3)*s134**(-1)*s56**(-1) - za(p2,p4)*za(
     &    p2,p6)*zb(p2,p3)*s234**(-1)*s56**(-1) + za(p2,p4)*za(p4,p6)*
     &    zb(p3,p4)*s234**(-1)*s56**(-1) + za(p2,p6)*za(p3,p4)*zb(p1,p3
     &    )*zb(p3,p6)*izb(p1,p6)*s134**(-1)*s56**(-1) )
      f56(mu,1,1,2) = f56(mu,1,1,2) + zab(p6,mu,p6)*s34**(-1) * ( za(p2
     &    ,p6)*za(p3,p4)*zb(p1,p3)**2*izb(p1,p6)*s134**(-1)*s56**(-1) )

      f56(mu,2,1,2)= + zab(p2,mu,p1)*s34**(-1) * (  - za(p1,p5)*za(p2,
     &    p4)*zb(p1,p6)*zb(p2,p3)*iza(p5,p6)*s234**(-1)*s16**(-1) )
      f56(mu,2,1,2) = f56(mu,2,1,2) + zab(p2,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p5)*zb(p1,p3)*zb(p1,p6)*iza(p5,p6)*s25**(-1)*
     &    s134**(-1) - za(p1,p5)*za(p2,p4)*zb(p1,p3)*zb(p1,p6)*iza(p5,
     &    p6)*s16**(-1)*s25**(-1) - za(p2,p4)*zb(p1,p6)*zb(p3,p6)*
     &    s16**(-1)*s25**(-1) - za(p2,p5)*za(p3,p4)*zb(p1,p3)*zb(p3,p6)
     &    *iza(p5,p6)*s25**(-1)*s134**(-1) )
      f56(mu,2,1,2) = f56(mu,2,1,2) + zab(p2,mu,p5)*s34**(-1) * ( za(p1
     &    ,p5)*za(p4,p5)*zb(p1,p3)*zb(p1,p6)*iza(p5,p6)*s16**(-1)*
     &    s25**(-1) + za(p4,p5)*zb(p1,p6)*zb(p3,p6)*s16**(-1)*s25**(-1)
     &     )
      f56(mu,2,1,2) = f56(mu,2,1,2) + zab(p2,mu,p6)*s34**(-1) * ( za(p2
     &    ,p4)*zb(p1,p6)*zb(p2,p3)*s234**(-1)*s16**(-1) )
      f56(mu,2,1,2) = f56(mu,2,1,2) + zab(p4,mu,p1)*s34**(-1) * ( za(p1
     &    ,p5)*za(p2,p4)*zb(p1,p6)*zb(p3,p4)*iza(p5,p6)*s234**(-1)*
     &    s16**(-1) )
      f56(mu,2,1,2) = f56(mu,2,1,2) + zab(p4,mu,p6)*s34**(-1) * (  - 
     &    za(p2,p4)*zb(p1,p6)*zb(p3,p4)*s234**(-1)*s16**(-1) )
      f56(mu,2,1,2) = f56(mu,2,1,2) + zab(p5,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p6)*zb(p1,p3)*zb(p1,p6)*iza(p5,p6)*s134**(-1)
     &    *s56**(-1) + za(p2,p4)*za(p2,p6)*zb(p1,p6)*zb(p2,p3)*iza(p5,
     &    p6)*s234**(-1)*s56**(-1) - za(p2,p4)*za(p4,p6)*zb(p1,p6)*zb(
     &    p3,p4)*iza(p5,p6)*s234**(-1)*s56**(-1) - za(p2,p6)*za(p3,p4)*
     &    zb(p1,p3)*zb(p3,p6)*iza(p5,p6)*s134**(-1)*s56**(-1) )
      f56(mu,2,1,2) = f56(mu,2,1,2) + zab(p6,mu,p6)*s34**(-1) * ( za(p1
     &    ,p4)*za(p2,p5)*zb(p1,p3)*zb(p1,p6)*iza(p5,p6)*s134**(-1)*
     &    s56**(-1) - za(p2,p4)*za(p2,p5)*zb(p1,p6)*zb(p2,p3)*iza(p5,p6
     &    )*s234**(-1)*s56**(-1) + za(p2,p4)*za(p4,p5)*zb(p1,p6)*zb(p3,
     &    p4)*iza(p5,p6)*s234**(-1)*s56**(-1) + za(p2,p5)*za(p3,p4)*zb(
     &    p1,p3)*zb(p3,p6)*iza(p5,p6)*s134**(-1)*s56**(-1) )

      f65(mu,1,1,2)= + zab(p1,mu,p1)*s34**(-1) * (  - za(p2,p4)*za(p2,
     &    p6)*zb(p1,p2)*zb(p1,p3)*izb(p1,p6)*s15**(-1)*s26**(-1) + za(
     &    p2,p6)*za(p4,p6)*zb(p1,p3)*s15**(-1)*s26**(-1) )
      f65(mu,1,1,2) = f65(mu,1,1,2) + zab(p2,mu,p1)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p6)*zb(p1,p2)*zb(p1,p3)*izb(p1,p6)*s134**(-1)
     &    *s26**(-1) + za(p1,p4)*za(p5,p6)*zb(p1,p3)*zb(p1,p5)*izb(p1,
     &    p6)*s134**(-1)*s56**(-1) - za(p2,p4)*za(p5,p6)*zb(p1,p5)*zb(
     &    p2,p3)*izb(p1,p6)*s234**(-1)*s56**(-1) )
      f65(mu,1,1,2) = f65(mu,1,1,2) + zab(p2,mu,p3)*s34**(-1) * (  - 
     &    za(p2,p6)*za(p3,p4)*zb(p1,p2)*zb(p1,p3)*izb(p1,p6)*s134**(-1)
     &    *s26**(-1) + za(p3,p4)*za(p5,p6)*zb(p1,p3)*zb(p1,p5)*izb(p1,
     &    p6)*s134**(-1)*s56**(-1) )
      f65(mu,1,1,2) = f65(mu,1,1,2) + zab(p4,mu,p1)*s34**(-1) * ( za(p2
     &    ,p4)*za(p5,p6)*zb(p1,p5)*zb(p3,p4)*izb(p1,p6)*s234**(-1)*
     &    s56**(-1) )
      f65(mu,1,1,2) = f65(mu,1,1,2) + zab(p5,mu,p1)*s34**(-1) * ( za(p2
     &    ,p4)*za(p2,p6)*zb(p1,p2)*zb(p3,p5)*izb(p1,p6)*s15**(-1)*
     &    s26**(-1) - za(p2,p4)*za(p2,p6)*zb(p1,p5)*zb(p2,p3)*izb(p1,p6
     &    )*s234**(-1)*s15**(-1) + za(p2,p4)*za(p4,p6)*zb(p1,p5)*zb(p3,
     &    p4)*izb(p1,p6)*s234**(-1)*s15**(-1) - za(p2,p6)*za(p4,p6)*zb(
     &    p3,p5)*s15**(-1)*s26**(-1) )
      f65(mu,1,1,2) = f65(mu,1,1,2) + zab(p6,mu,p1)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p6)*zb(p1,p3)*s134**(-1)*s56**(-1) - za(p1,p4
     &    )*za(p2,p6)*zb(p1,p3)*s134**(-1)*s26**(-1) + za(p2,p4)*za(p2,
     &    p6)*zb(p2,p3)*s234**(-1)*s56**(-1) - za(p2,p4)*za(p4,p6)*zb(
     &    p3,p4)*s234**(-1)*s56**(-1) - za(p2,p6)*za(p3,p4)*zb(p1,p3)*
     &    zb(p3,p6)*izb(p1,p6)*s134**(-1)*s56**(-1) )
      f65(mu,1,1,2) = f65(mu,1,1,2) + zab(p6,mu,p3)*s34**(-1) * (  - 
     &    za(p2,p6)*za(p3,p4)*zb(p1,p3)*s134**(-1)*s26**(-1) )
      f65(mu,1,1,2) = f65(mu,1,1,2) + zab(p6,mu,p6)*s34**(-1) * (  - 
     &    za(p2,p6)*za(p3,p4)*zb(p1,p3)**2*izb(p1,p6)*s134**(-1)*
     &    s56**(-1) )

      f65(mu,2,1,2)= + zab(p1,mu,p1)*s34**(-1) * (  - za(p2,p4)*za(p2,
     &    p5)*zb(p1,p3)*zb(p2,p6)*iza(p5,p6)*s15**(-1)*s26**(-1) - za(
     &    p2,p4)*za(p2,p5)*zb(p1,p6)*zb(p2,p3)*iza(p5,p6)*s234**(-1)*
     &    s15**(-1) + za(p2,p4)*za(p4,p5)*zb(p1,p6)*zb(p3,p4)*iza(p5,p6
     &    )*s234**(-1)*s15**(-1) )
      f65(mu,2,1,2) = f65(mu,2,1,2) + zab(p2,mu,p1)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p5)*zb(p1,p3)*zb(p2,p6)*iza(p5,p6)*s134**(-1)
     &    *s26**(-1) )
      f65(mu,2,1,2) = f65(mu,2,1,2) + zab(p2,mu,p3)*s34**(-1) * (  - 
     &    za(p2,p5)*za(p3,p4)*zb(p1,p3)*zb(p2,p6)*iza(p5,p6)*s134**(-1)
     &    *s26**(-1) )
      f65(mu,2,1,2) = f65(mu,2,1,2) + zab(p5,mu,p1)*s34**(-1) * (  - 
     &    za(p2,p4)*za(p2,p5)*zb(p2,p3)*zb(p5,p6)*iza(p5,p6)*s234**(-1)
     &    *s15**(-1) + za(p2,p4)*za(p2,p5)*zb(p2,p6)*zb(p3,p5)*iza(p5,
     &    p6)*s15**(-1)*s26**(-1) + za(p2,p4)*za(p4,p5)*zb(p3,p4)*zb(p5
     &    ,p6)*iza(p5,p6)*s234**(-1)*s15**(-1) )
      f65(mu,2,1,2) = f65(mu,2,1,2) + zab(p5,mu,p6)*s34**(-1) * ( za(p1
     &    ,p4)*za(p2,p6)*zb(p1,p3)*zb(p1,p6)*iza(p5,p6)*s134**(-1)*
     &    s56**(-1) - za(p2,p4)*za(p2,p6)*zb(p1,p6)*zb(p2,p3)*iza(p5,p6
     &    )*s234**(-1)*s56**(-1) + za(p2,p4)*za(p4,p6)*zb(p1,p6)*zb(p3,
     &    p4)*iza(p5,p6)*s234**(-1)*s56**(-1) + za(p2,p6)*za(p3,p4)*zb(
     &    p1,p3)*zb(p3,p6)*iza(p5,p6)*s134**(-1)*s56**(-1) )
      f65(mu,2,1,2) = f65(mu,2,1,2) + zab(p6,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p5)*zb(p1,p3)*zb(p1,p6)*iza(p5,p6)*s134**(-1)
     &    *s56**(-1) + za(p2,p4)*za(p2,p5)*zb(p1,p6)*zb(p2,p3)*iza(p5,
     &    p6)*s234**(-1)*s56**(-1) - za(p2,p4)*za(p4,p5)*zb(p1,p6)*zb(
     &    p3,p4)*iza(p5,p6)*s234**(-1)*s56**(-1) - za(p2,p5)*za(p3,p4)*
     &    zb(p1,p3)*zb(p3,p6)*iza(p5,p6)*s134**(-1)*s56**(-1) )

      f56(mu,1,2,1)= + zab(p1,mu,p2)*s34**(-1) * ( za(p1,p6)*za(p2,p3)*
     &    zb(p1,p5)*zb(p2,p4)*izb(p5,p6)*s234**(-1)*s16**(-1) )
      f56(mu,1,2,1) = f56(mu,1,2,1) + zab(p1,mu,p4)*s34**(-1) * (  - 
     &    za(p1,p6)*za(p3,p4)*zb(p1,p5)*zb(p2,p4)*izb(p5,p6)*s234**(-1)
     &    *s16**(-1) )
      f56(mu,1,2,1) = f56(mu,1,2,1) + zab(p2,mu,p2)*s34**(-1) * ( za(p1
     &    ,p3)*za(p1,p6)*zb(p1,p4)*zb(p2,p5)*izb(p5,p6)*s25**(-1)*
     &    s134**(-1) + za(p1,p3)*za(p1,p6)*zb(p1,p5)*zb(p2,p4)*izb(p5,
     &    p6)*s16**(-1)*s25**(-1) + za(p1,p3)*za(p3,p6)*zb(p2,p5)*zb(p3
     &    ,p4)*izb(p5,p6)*s25**(-1)*s134**(-1) + za(p1,p6)*za(p3,p6)*
     &    zb(p2,p4)*s16**(-1)*s25**(-1) )
      f56(mu,1,2,1) = f56(mu,1,2,1) + zab(p5,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p1,p6)*zb(p1,p5)*zb(p4,p5)*izb(p5,p6)*s16**(-1)*
     &    s25**(-1) - za(p1,p6)*za(p3,p6)*zb(p4,p5)*s16**(-1)*s25**(-1)
     &     )
      f56(mu,1,2,1) = f56(mu,1,2,1) + zab(p6,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p6)*za(p2,p3)*zb(p2,p4)*s234**(-1)*s16**(-1) )
      f56(mu,1,2,1) = f56(mu,1,2,1) + zab(p6,mu,p4)*s34**(-1) * ( za(p1
     &    ,p6)*za(p3,p4)*zb(p2,p4)*s234**(-1)*s16**(-1) )
      f56(mu,1,2,1) = f56(mu,1,2,1) + zab(p6,mu,p5)*s34**(-1) * ( za(p1
     &    ,p3)*za(p1,p6)*zb(p1,p4)*zb(p2,p6)*izb(p5,p6)*s134**(-1)*
     &    s56**(-1) + za(p1,p3)*za(p3,p6)*zb(p2,p6)*zb(p3,p4)*izb(p5,p6
     &    )*s134**(-1)*s56**(-1) - za(p1,p6)*za(p2,p3)*zb(p2,p4)*zb(p2,
     &    p6)*izb(p5,p6)*s234**(-1)*s56**(-1) + za(p1,p6)*za(p3,p4)*zb(
     &    p2,p4)*zb(p4,p6)*izb(p5,p6)*s234**(-1)*s56**(-1) )
      f56(mu,1,2,1) = f56(mu,1,2,1) + zab(p6,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p1,p6)*zb(p1,p4)*zb(p2,p5)*izb(p5,p6)*s134**(-1)
     &    *s56**(-1) - za(p1,p3)*za(p3,p6)*zb(p2,p5)*zb(p3,p4)*izb(p5,
     &    p6)*s134**(-1)*s56**(-1) + za(p1,p6)*za(p2,p3)*zb(p2,p4)*zb(
     &    p2,p5)*izb(p5,p6)*s234**(-1)*s56**(-1) - za(p1,p6)*za(p3,p4)*
     &    zb(p2,p4)*zb(p4,p5)*izb(p5,p6)*s234**(-1)*s56**(-1) )

      f56(mu,2,2,1)= + zab(p1,mu,p2)*s34**(-1) * (  - za(p1,p2)*za(p2,
     &    p3)*zb(p1,p6)*zb(p2,p4)*iza(p2,p6)*s234**(-1)*s16**(-1) + za(
     &    p1,p3)*za(p2,p5)*zb(p1,p4)*zb(p5,p6)*iza(p2,p6)*s134**(-1)*
     &    s56**(-1) - za(p2,p3)*za(p2,p5)*zb(p2,p4)*zb(p5,p6)*iza(p2,p6
     &    )*s234**(-1)*s56**(-1) )
      f56(mu,2,2,1) = f56(mu,2,2,1) + zab(p1,mu,p4)*s34**(-1) * ( za(p1
     &    ,p2)*za(p3,p4)*zb(p1,p6)*zb(p2,p4)*iza(p2,p6)*s234**(-1)*
     &    s16**(-1) + za(p2,p5)*za(p3,p4)*zb(p2,p4)*zb(p5,p6)*iza(p2,p6
     &    )*s234**(-1)*s56**(-1) )
      f56(mu,2,2,1) = f56(mu,2,2,1) + zab(p2,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p1,p3)*zb(p1,p4)*zb(p2,p6)*iza(p2,p6)*s25**(-1)*
     &    s134**(-1) - za(p1,p2)*za(p1,p3)*zb(p1,p6)*zb(p2,p4)*iza(p2,
     &    p6)*s16**(-1)*s25**(-1) + za(p1,p3)*za(p2,p3)*zb(p2,p6)*zb(p3
     &    ,p4)*iza(p2,p6)*s25**(-1)*s134**(-1) )
      f56(mu,2,2,1) = f56(mu,2,2,1) + zab(p2,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p1,p6)*zb(p1,p4)*zb(p2,p6)*iza(p2,p6)*s134**(-1)
     &    *s56**(-1) - za(p1,p3)*za(p3,p6)*zb(p2,p6)*zb(p3,p4)*iza(p2,
     &    p6)*s134**(-1)*s56**(-1) + za(p1,p6)*za(p2,p3)*zb(p2,p4)*zb(
     &    p2,p6)*iza(p2,p6)*s234**(-1)*s56**(-1) - za(p1,p6)*za(p3,p4)*
     &    zb(p2,p4)*zb(p4,p6)*iza(p2,p6)*s234**(-1)*s56**(-1) )
      f56(mu,2,2,1) = f56(mu,2,2,1) + zab(p3,mu,p2)*s34**(-1) * ( za(p1
     &    ,p3)*za(p2,p5)*zb(p3,p4)*zb(p5,p6)*iza(p2,p6)*s134**(-1)*
     &    s56**(-1) )
      f56(mu,2,2,1) = f56(mu,2,2,1) + zab(p5,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p1,p3)*zb(p1,p4)*zb(p5,p6)*iza(p2,p6)*s25**(-1)*
     &    s134**(-1) + za(p1,p2)*za(p1,p3)*zb(p1,p6)*zb(p4,p5)*iza(p2,
     &    p6)*s16**(-1)*s25**(-1) + za(p1,p3)*za(p2,p3)*zb(p3,p4)*zb(p5
     &    ,p6)*iza(p2,p6)*s25**(-1)*s134**(-1) )
      f56(mu,2,2,1) = f56(mu,2,2,1) + zab(p6,mu,p6)*s34**(-1) * ( za(p1
     &    ,p2)*za(p1,p3)*zb(p1,p4)*zb(p2,p6)*iza(p2,p6)*s134**(-1)*
     &    s56**(-1) - za(p1,p2)*za(p2,p3)*zb(p2,p4)*zb(p2,p6)*iza(p2,p6
     &    )*s234**(-1)*s56**(-1) + za(p1,p2)*za(p3,p4)*zb(p2,p4)*zb(p4,
     &    p6)*iza(p2,p6)*s234**(-1)*s56**(-1) - za(p1,p3)*za(p2,p3)*zb(
     &    p2,p6)*zb(p3,p4)*iza(p2,p6)*s134**(-1)*s56**(-1) )

      f65(mu,1,2,1)= + zab(p1,mu,p1)*s34**(-1) * ( za(p1,p3)*za(p2,p6)*
     &    zb(p2,p4)*zb(p2,p5)*izb(p5,p6)*s15**(-1)*s26**(-1) + za(p1,p6
     &    )*za(p2,p3)*zb(p2,p4)*zb(p2,p5)*izb(p5,p6)*s234**(-1)*
     &    s15**(-1) - za(p1,p6)*za(p3,p4)*zb(p2,p4)*zb(p4,p5)*izb(p5,p6
     &    )*s234**(-1)*s15**(-1) )
      f65(mu,1,2,1) = f65(mu,1,2,1) + zab(p1,mu,p2)*s34**(-1) * ( za(p1
     &    ,p3)*za(p2,p6)*zb(p1,p4)*zb(p2,p5)*izb(p5,p6)*s134**(-1)*
     &    s26**(-1) )
      f65(mu,1,2,1) = f65(mu,1,2,1) + zab(p1,mu,p5)*s34**(-1) * ( za(p2
     &    ,p3)*za(p5,p6)*zb(p2,p4)*zb(p2,p5)*izb(p5,p6)*s234**(-1)*
     &    s15**(-1) - za(p2,p6)*za(p3,p5)*zb(p2,p4)*zb(p2,p5)*izb(p5,p6
     &    )*s15**(-1)*s26**(-1) - za(p3,p4)*za(p5,p6)*zb(p2,p4)*zb(p4,
     &    p5)*izb(p5,p6)*s234**(-1)*s15**(-1) )
      f65(mu,1,2,1) = f65(mu,1,2,1) + zab(p3,mu,p2)*s34**(-1) * ( za(p1
     &    ,p3)*za(p2,p6)*zb(p2,p5)*zb(p3,p4)*izb(p5,p6)*s134**(-1)*
     &    s26**(-1) )
      f65(mu,1,2,1) = f65(mu,1,2,1) + zab(p6,mu,p5)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p1,p6)*zb(p1,p4)*zb(p2,p6)*izb(p5,p6)*s134**(-1)
     &    *s56**(-1) - za(p1,p3)*za(p3,p6)*zb(p2,p6)*zb(p3,p4)*izb(p5,
     &    p6)*s134**(-1)*s56**(-1) + za(p1,p6)*za(p2,p3)*zb(p2,p4)*zb(
     &    p2,p6)*izb(p5,p6)*s234**(-1)*s56**(-1) - za(p1,p6)*za(p3,p4)*
     &    zb(p2,p4)*zb(p4,p6)*izb(p5,p6)*s234**(-1)*s56**(-1) )
      f65(mu,1,2,1) = f65(mu,1,2,1) + zab(p6,mu,p6)*s34**(-1) * ( za(p1
     &    ,p3)*za(p1,p6)*zb(p1,p4)*zb(p2,p5)*izb(p5,p6)*s134**(-1)*
     &    s56**(-1) + za(p1,p3)*za(p3,p6)*zb(p2,p5)*zb(p3,p4)*izb(p5,p6
     &    )*s134**(-1)*s56**(-1) - za(p1,p6)*za(p2,p3)*zb(p2,p4)*zb(p2,
     &    p5)*izb(p5,p6)*s234**(-1)*s56**(-1) + za(p1,p6)*za(p3,p4)*zb(
     &    p2,p4)*zb(p4,p5)*izb(p5,p6)*s234**(-1)*s56**(-1) )

      f65(mu,2,2,1)= + zab(p1,mu,p1)*s34**(-1) * (  - za(p1,p2)*za(p2,
     &    p3)*zb(p2,p4)*zb(p2,p6)*iza(p2,p6)*s234**(-1)*s15**(-1) + za(
     &    p1,p2)*za(p3,p4)*zb(p2,p4)*zb(p4,p6)*iza(p2,p6)*s234**(-1)*
     &    s15**(-1) - za(p1,p3)*zb(p2,p6)*zb(p4,p6)*s15**(-1)*s26**(-1)
     &     )
      f65(mu,2,2,1) = f65(mu,2,2,1) + zab(p1,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p2,p5)*zb(p1,p4)*zb(p5,p6)*iza(p2,p6)*s134**(-1)
     &    *s56**(-1) + za(p2,p3)*za(p2,p5)*zb(p2,p4)*zb(p5,p6)*iza(p2,
     &    p6)*s234**(-1)*s56**(-1) )
      f65(mu,2,2,1) = f65(mu,2,2,1) + zab(p1,mu,p4)*s34**(-1) * (  - 
     &    za(p2,p5)*za(p3,p4)*zb(p2,p4)*zb(p5,p6)*iza(p2,p6)*s234**(-1)
     &    *s56**(-1) )
      f65(mu,2,2,1) = f65(mu,2,2,1) + zab(p1,mu,p5)*s34**(-1) * ( za(p2
     &    ,p3)*za(p2,p5)*zb(p2,p4)*zb(p2,p6)*iza(p2,p6)*s234**(-1)*
     &    s15**(-1) - za(p2,p5)*za(p3,p4)*zb(p2,p4)*zb(p4,p6)*iza(p2,p6
     &    )*s234**(-1)*s15**(-1) + za(p3,p5)*zb(p2,p6)*zb(p4,p6)*
     &    s15**(-1)*s26**(-1) )
      f65(mu,2,2,1) = f65(mu,2,2,1) + zab(p1,mu,p6)*s34**(-1) * ( za(p1
     &    ,p3)*zb(p1,p4)*zb(p2,p6)*s134**(-1)*s26**(-1) )
      f65(mu,2,2,1) = f65(mu,2,2,1) + zab(p2,mu,p6)*s34**(-1) * ( za(p1
     &    ,p3)*za(p1,p6)*zb(p1,p4)*zb(p2,p6)*iza(p2,p6)*s134**(-1)*
     &    s56**(-1) + za(p1,p3)*za(p3,p6)*zb(p2,p6)*zb(p3,p4)*iza(p2,p6
     &    )*s134**(-1)*s56**(-1) - za(p1,p6)*za(p2,p3)*zb(p2,p4)*zb(p2,
     &    p6)*iza(p2,p6)*s234**(-1)*s56**(-1) + za(p1,p6)*za(p3,p4)*zb(
     &    p2,p4)*zb(p4,p6)*iza(p2,p6)*s234**(-1)*s56**(-1) )
      f65(mu,2,2,1) = f65(mu,2,2,1) + zab(p3,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p3)*za(p2,p5)*zb(p3,p4)*zb(p5,p6)*iza(p2,p6)*s134**(-1)
     &    *s56**(-1) )
      f65(mu,2,2,1) = f65(mu,2,2,1) + zab(p3,mu,p6)*s34**(-1) * ( za(p1
     &    ,p3)*zb(p2,p6)*zb(p3,p4)*s134**(-1)*s26**(-1) )
      f65(mu,2,2,1) = f65(mu,2,2,1) + zab(p6,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p1,p3)*zb(p1,p4)*zb(p2,p6)*iza(p2,p6)*s134**(-1)
     &    *s56**(-1) + za(p1,p2)*za(p2,p3)*zb(p2,p4)*zb(p2,p6)*iza(p2,
     &    p6)*s234**(-1)*s56**(-1) - za(p1,p2)*za(p3,p4)*zb(p2,p4)*zb(
     &    p4,p6)*iza(p2,p6)*s234**(-1)*s56**(-1) + za(p1,p3)*za(p2,p3)*
     &    zb(p2,p6)*zb(p3,p4)*iza(p2,p6)*s134**(-1)*s56**(-1) )

      f56(mu,1,2,2)= + zab(p1,mu,p2)*s34**(-1) * ( za(p1,p6)*za(p2,p4)*
     &    zb(p1,p5)*zb(p2,p3)*izb(p5,p6)*s234**(-1)*s16**(-1) )
      f56(mu,1,2,2) = f56(mu,1,2,2) + zab(p1,mu,p3)*s34**(-1) * ( za(p1
     &    ,p6)*za(p3,p4)*zb(p1,p5)*zb(p2,p3)*izb(p5,p6)*s234**(-1)*
     &    s16**(-1) )
      f56(mu,1,2,2) = f56(mu,1,2,2) + zab(p2,mu,p2)*s34**(-1) * ( za(p1
     &    ,p4)*za(p1,p6)*zb(p1,p3)*zb(p2,p5)*izb(p5,p6)*s25**(-1)*
     &    s134**(-1) + za(p1,p4)*za(p1,p6)*zb(p1,p5)*zb(p2,p3)*izb(p5,
     &    p6)*s16**(-1)*s25**(-1) - za(p1,p4)*za(p4,p6)*zb(p2,p5)*zb(p3
     &    ,p4)*izb(p5,p6)*s25**(-1)*s134**(-1) + za(p1,p6)*za(p4,p6)*
     &    zb(p2,p3)*s16**(-1)*s25**(-1) )
      f56(mu,1,2,2) = f56(mu,1,2,2) + zab(p5,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p1,p6)*zb(p1,p5)*zb(p3,p5)*izb(p5,p6)*s16**(-1)*
     &    s25**(-1) - za(p1,p6)*za(p4,p6)*zb(p3,p5)*s16**(-1)*s25**(-1)
     &     )
      f56(mu,1,2,2) = f56(mu,1,2,2) + zab(p6,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p6)*za(p2,p4)*zb(p2,p3)*s234**(-1)*s16**(-1) )
      f56(mu,1,2,2) = f56(mu,1,2,2) + zab(p6,mu,p3)*s34**(-1) * (  - 
     &    za(p1,p6)*za(p3,p4)*zb(p2,p3)*s234**(-1)*s16**(-1) )
      f56(mu,1,2,2) = f56(mu,1,2,2) + zab(p6,mu,p5)*s34**(-1) * ( za(p1
     &    ,p4)*za(p1,p6)*zb(p1,p3)*zb(p2,p6)*izb(p5,p6)*s134**(-1)*
     &    s56**(-1) - za(p1,p4)*za(p4,p6)*zb(p2,p6)*zb(p3,p4)*izb(p5,p6
     &    )*s134**(-1)*s56**(-1) - za(p1,p6)*za(p2,p4)*zb(p2,p3)*zb(p2,
     &    p6)*izb(p5,p6)*s234**(-1)*s56**(-1) - za(p1,p6)*za(p3,p4)*zb(
     &    p2,p3)*zb(p3,p6)*izb(p5,p6)*s234**(-1)*s56**(-1) )
      f56(mu,1,2,2) = f56(mu,1,2,2) + zab(p6,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p1,p6)*zb(p1,p3)*zb(p2,p5)*izb(p5,p6)*s134**(-1)
     &    *s56**(-1) + za(p1,p4)*za(p4,p6)*zb(p2,p5)*zb(p3,p4)*izb(p5,
     &    p6)*s134**(-1)*s56**(-1) + za(p1,p6)*za(p2,p4)*zb(p2,p3)*zb(
     &    p2,p5)*izb(p5,p6)*s234**(-1)*s56**(-1) + za(p1,p6)*za(p3,p4)*
     &    zb(p2,p3)*zb(p3,p5)*izb(p5,p6)*s234**(-1)*s56**(-1) )

      f56(mu,2,2,2)= + zab(p1,mu,p2)*s34**(-1) * (  - za(p1,p2)*za(p2,
     &    p4)*zb(p1,p6)*zb(p2,p3)*iza(p2,p6)*s234**(-1)*s16**(-1) + za(
     &    p1,p4)*za(p2,p5)*zb(p1,p3)*zb(p5,p6)*iza(p2,p6)*s134**(-1)*
     &    s56**(-1) - za(p2,p4)*za(p2,p5)*zb(p2,p3)*zb(p5,p6)*iza(p2,p6
     &    )*s234**(-1)*s56**(-1) )
      f56(mu,2,2,2) = f56(mu,2,2,2) + zab(p1,mu,p3)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p3,p4)*zb(p1,p6)*zb(p2,p3)*iza(p2,p6)*s234**(-1)
     &    *s16**(-1) - za(p2,p5)*za(p3,p4)*zb(p2,p3)*zb(p5,p6)*iza(p2,
     &    p6)*s234**(-1)*s56**(-1) )
      f56(mu,2,2,2) = f56(mu,2,2,2) + zab(p2,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p1,p4)*zb(p1,p3)*zb(p2,p6)*iza(p2,p6)*s25**(-1)*
     &    s134**(-1) - za(p1,p2)*za(p1,p4)*zb(p1,p6)*zb(p2,p3)*iza(p2,
     &    p6)*s16**(-1)*s25**(-1) - za(p1,p4)*za(p2,p4)*zb(p2,p6)*zb(p3
     &    ,p4)*iza(p2,p6)*s25**(-1)*s134**(-1) )
      f56(mu,2,2,2) = f56(mu,2,2,2) + zab(p2,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p1,p6)*zb(p1,p3)*zb(p2,p6)*iza(p2,p6)*s134**(-1)
     &    *s56**(-1) + za(p1,p4)*za(p4,p6)*zb(p2,p6)*zb(p3,p4)*iza(p2,
     &    p6)*s134**(-1)*s56**(-1) + za(p1,p6)*za(p2,p4)*zb(p2,p3)*zb(
     &    p2,p6)*iza(p2,p6)*s234**(-1)*s56**(-1) + za(p1,p6)*za(p3,p4)*
     &    zb(p2,p3)*zb(p3,p6)*iza(p2,p6)*s234**(-1)*s56**(-1) )
      f56(mu,2,2,2) = f56(mu,2,2,2) + zab(p4,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p5)*zb(p3,p4)*zb(p5,p6)*iza(p2,p6)*s134**(-1)
     &    *s56**(-1) )
      f56(mu,2,2,2) = f56(mu,2,2,2) + zab(p5,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p1,p4)*zb(p1,p3)*zb(p5,p6)*iza(p2,p6)*s25**(-1)*
     &    s134**(-1) + za(p1,p2)*za(p1,p4)*zb(p1,p6)*zb(p3,p5)*iza(p2,
     &    p6)*s16**(-1)*s25**(-1) - za(p1,p4)*za(p2,p4)*zb(p3,p4)*zb(p5
     &    ,p6)*iza(p2,p6)*s25**(-1)*s134**(-1) )
      f56(mu,2,2,2) = f56(mu,2,2,2) + zab(p6,mu,p6)*s34**(-1) * ( za(p1
     &    ,p2)*za(p1,p4)*zb(p1,p3)*zb(p2,p6)*iza(p2,p6)*s134**(-1)*
     &    s56**(-1) - za(p1,p2)*za(p2,p4)*zb(p2,p3)*zb(p2,p6)*iza(p2,p6
     &    )*s234**(-1)*s56**(-1) - za(p1,p2)*za(p3,p4)*zb(p2,p3)*zb(p3,
     &    p6)*iza(p2,p6)*s234**(-1)*s56**(-1) + za(p1,p4)*za(p2,p4)*zb(
     &    p2,p6)*zb(p3,p4)*iza(p2,p6)*s134**(-1)*s56**(-1) )

      f65(mu,1,2,2)= + zab(p1,mu,p1)*s34**(-1) * ( za(p1,p4)*za(p2,p6)*
     &    zb(p2,p3)*zb(p2,p5)*izb(p5,p6)*s15**(-1)*s26**(-1) + za(p1,p6
     &    )*za(p2,p4)*zb(p2,p3)*zb(p2,p5)*izb(p5,p6)*s234**(-1)*
     &    s15**(-1) + za(p1,p6)*za(p3,p4)*zb(p2,p3)*zb(p3,p5)*izb(p5,p6
     &    )*s234**(-1)*s15**(-1) )
      f65(mu,1,2,2) = f65(mu,1,2,2) + zab(p1,mu,p2)*s34**(-1) * ( za(p1
     &    ,p4)*za(p2,p6)*zb(p1,p3)*zb(p2,p5)*izb(p5,p6)*s134**(-1)*
     &    s26**(-1) )
      f65(mu,1,2,2) = f65(mu,1,2,2) + zab(p1,mu,p5)*s34**(-1) * ( za(p2
     &    ,p4)*za(p5,p6)*zb(p2,p3)*zb(p2,p5)*izb(p5,p6)*s234**(-1)*
     &    s15**(-1) - za(p2,p6)*za(p4,p5)*zb(p2,p3)*zb(p2,p5)*izb(p5,p6
     &    )*s15**(-1)*s26**(-1) + za(p3,p4)*za(p5,p6)*zb(p2,p3)*zb(p3,
     &    p5)*izb(p5,p6)*s234**(-1)*s15**(-1) )
      f65(mu,1,2,2) = f65(mu,1,2,2) + zab(p4,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p6)*zb(p2,p5)*zb(p3,p4)*izb(p5,p6)*s134**(-1)
     &    *s26**(-1) )
      f65(mu,1,2,2) = f65(mu,1,2,2) + zab(p6,mu,p5)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p1,p6)*zb(p1,p3)*zb(p2,p6)*izb(p5,p6)*s134**(-1)
     &    *s56**(-1) + za(p1,p4)*za(p4,p6)*zb(p2,p6)*zb(p3,p4)*izb(p5,
     &    p6)*s134**(-1)*s56**(-1) + za(p1,p6)*za(p2,p4)*zb(p2,p3)*zb(
     &    p2,p6)*izb(p5,p6)*s234**(-1)*s56**(-1) + za(p1,p6)*za(p3,p4)*
     &    zb(p2,p3)*zb(p3,p6)*izb(p5,p6)*s234**(-1)*s56**(-1) )
      f65(mu,1,2,2) = f65(mu,1,2,2) + zab(p6,mu,p6)*s34**(-1) * ( za(p1
     &    ,p4)*za(p1,p6)*zb(p1,p3)*zb(p2,p5)*izb(p5,p6)*s134**(-1)*
     &    s56**(-1) - za(p1,p4)*za(p4,p6)*zb(p2,p5)*zb(p3,p4)*izb(p5,p6
     &    )*s134**(-1)*s56**(-1) - za(p1,p6)*za(p2,p4)*zb(p2,p3)*zb(p2,
     &    p5)*izb(p5,p6)*s234**(-1)*s56**(-1) - za(p1,p6)*za(p3,p4)*zb(
     &    p2,p3)*zb(p3,p5)*izb(p5,p6)*s234**(-1)*s56**(-1) )

      f65(mu,2,2,2)= + zab(p1,mu,p1)*s34**(-1) * (  - za(p1,p2)*za(p2,
     &    p4)*zb(p2,p3)*zb(p2,p6)*iza(p2,p6)*s234**(-1)*s15**(-1) - za(
     &    p1,p2)*za(p3,p4)*zb(p2,p3)*zb(p3,p6)*iza(p2,p6)*s234**(-1)*
     &    s15**(-1) - za(p1,p4)*zb(p2,p6)*zb(p3,p6)*s15**(-1)*s26**(-1)
     &     )
      f65(mu,2,2,2) = f65(mu,2,2,2) + zab(p1,mu,p2)*s34**(-1) * (  - 
     &    za(p1,p4)*za(p2,p5)*zb(p1,p3)*zb(p5,p6)*iza(p2,p6)*s134**(-1)
     &    *s56**(-1) + za(p2,p4)*za(p2,p5)*zb(p2,p3)*zb(p5,p6)*iza(p2,
     &    p6)*s234**(-1)*s56**(-1) )
      f65(mu,2,2,2) = f65(mu,2,2,2) + zab(p1,mu,p3)*s34**(-1) * ( za(p2
     &    ,p5)*za(p3,p4)*zb(p2,p3)*zb(p5,p6)*iza(p2,p6)*s234**(-1)*
     &    s56**(-1) )
      f65(mu,2,2,2) = f65(mu,2,2,2) + zab(p1,mu,p5)*s34**(-1) * ( za(p2
     &    ,p4)*za(p2,p5)*zb(p2,p3)*zb(p2,p6)*iza(p2,p6)*s234**(-1)*
     &    s15**(-1) + za(p2,p5)*za(p3,p4)*zb(p2,p3)*zb(p3,p6)*iza(p2,p6
     &    )*s234**(-1)*s15**(-1) + za(p4,p5)*zb(p2,p6)*zb(p3,p6)*
     &    s15**(-1)*s26**(-1) )
      f65(mu,2,2,2) = f65(mu,2,2,2) + zab(p1,mu,p6)*s34**(-1) * ( za(p1
     &    ,p4)*zb(p1,p3)*zb(p2,p6)*s134**(-1)*s26**(-1) )
      f65(mu,2,2,2) = f65(mu,2,2,2) + zab(p2,mu,p6)*s34**(-1) * ( za(p1
     &    ,p4)*za(p1,p6)*zb(p1,p3)*zb(p2,p6)*iza(p2,p6)*s134**(-1)*
     &    s56**(-1) - za(p1,p4)*za(p4,p6)*zb(p2,p6)*zb(p3,p4)*iza(p2,p6
     &    )*s134**(-1)*s56**(-1) - za(p1,p6)*za(p2,p4)*zb(p2,p3)*zb(p2,
     &    p6)*iza(p2,p6)*s234**(-1)*s56**(-1) - za(p1,p6)*za(p3,p4)*zb(
     &    p2,p3)*zb(p3,p6)*iza(p2,p6)*s234**(-1)*s56**(-1) )
      f65(mu,2,2,2) = f65(mu,2,2,2) + zab(p4,mu,p2)*s34**(-1) * ( za(p1
     &    ,p4)*za(p2,p5)*zb(p3,p4)*zb(p5,p6)*iza(p2,p6)*s134**(-1)*
     &    s56**(-1) )
      f65(mu,2,2,2) = f65(mu,2,2,2) + zab(p4,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p4)*zb(p2,p6)*zb(p3,p4)*s134**(-1)*s26**(-1) )
      f65(mu,2,2,2) = f65(mu,2,2,2) + zab(p6,mu,p6)*s34**(-1) * (  - 
     &    za(p1,p2)*za(p1,p4)*zb(p1,p3)*zb(p2,p6)*iza(p2,p6)*s134**(-1)
     &    *s56**(-1) + za(p1,p2)*za(p2,p4)*zb(p2,p3)*zb(p2,p6)*iza(p2,
     &    p6)*s234**(-1)*s56**(-1) + za(p1,p2)*za(p3,p4)*zb(p2,p3)*zb(
     &    p3,p6)*iza(p2,p6)*s234**(-1)*s56**(-1) - za(p1,p4)*za(p2,p4)*
     &    zb(p2,p6)*zb(p3,p4)*iza(p2,p6)*s134**(-1)*s56**(-1) )

      enddo
      fqed(:,:,:,:)=f56(:,:,:,:)+f65(:,:,:,:)
      return
      end
