      subroutine setBmunu(p,Bmunu)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     April, 2012.                                                     *
*----Matrix element squared for ttbar production                       *
*----averaged over initial colours and spins                           *
*    with one gluon line uncontracted                                  *
*     g(-p1)+g(-p2)--> t(p3,p4,p5)+tb(p6,p7,p8)                        *
*                                                                      *
*     g(-p1) +g(-p2)= nu(p3)+e+(p4)+b(p5)+bbar(p6)+e-(p7)+nubar(p8)    *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'process.f'
      include 'masses.f'
      include 'zprods_com.f'
      double precision p(mxpart,4),s,mt2,fac,q(mxpart,4),
     & s1t,s2t,s12,c1,c2,Bmunu(4,4,2)
      integer j,k,nu
      double complex prop
C-----statement function
      s(j,k)=2d0
     & *(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
C-----end statement function

      s1t=s(1,3)+s(1,4)+s(1,5)
      s2t=s(2,3)+s(2,4)+s(2,5)
      s12=s(1,2)
      mt2=mt**2
      c1=mt2/(s(3,4)+s(4,5))
      c2=mt2/(s(6,7)+s(7,8))
      prop=dcmplx(s(3,4)-wmass**2,wmass*wwidth)
     .    *dcmplx(s(7,8)-wmass**2,wmass*wwidth)
     .    *dcmplx(zip,mt*twidth)**2
      fac=avegg*V*xn/4d0*gwsq**4*gsq**2/abs(prop)**2*s(5,3)*s(6,8)
      fac=fac*2d0 
C extra factor of two from ggn.frm
      
C--include factor for hadronic decays
c      if ((case .eq. 'tt_bbh') .or. (case .eq. 'tt_hdk')) fac=2d0*xn*fac


c--- pt=p3+p4+p5, ptb=p6+p7+p8
      do nu=1,4
      q(1,nu)=p(1,nu)    
      q(2,nu)=p(2,nu)    
      q(3,nu)=p(3,nu)+p(4,nu)*(1d0-c1)+p(5,nu)
      q(4,nu)=p(4,nu)    
      q(5,nu)=p(6,nu)+p(7,nu)*(1d0-c2)+p(8,nu)
      q(6,nu)=p(7,nu)    
      enddo

      call spinoru(6,q,za,zb)
      call Bmunusrc(q,s1t,s2t,s12,c1,c2,Bmunu)
      
      Bmunu(:,:,:)=fac*Bmunu(:,:,:)
      return
      end
