      subroutine pwhg_p_p4(p,p4)
      implicit none
      real * 8 p(0:3),p4(1:4)
      p4(4)=p(0)      
      p4(1)=p(1)
      p4(2)=p(2)
      p4(3)=p(3)        
      end

      subroutine pwhg_getpt(p,pt)
      implicit none
      real * 8 p(4),pt
      pt=sqrt(p(1)**2+p(2)**2)
      end      

      function pwhg_pt(p)
      implicit none      
      real * 8 pwhg_pt,p(0:3)
      real * 8 p4(4),pt
      call pwhg_p_p4(p,p4)
      call pwhg_getpt(p4,pt)
      pwhg_pt=pt
      end      


      subroutine pwhg_getinvmass(p,m)
      implicit none
      real * 8 p(4),m
      real * 8 arg 
      arg = p(4)**2-p(1)**2-p(2)**2-p(3)**2
c     SIGN(A,B) returns the value of A with the sign of B 
      m=sign(1.d0,arg)*sqrt(abs(arg))
      end

      function pwhg_invmass(p)
      implicit none
      real * 8 pwhg_invmass,p(0:3)
      real * 8 p4(4),m
      call pwhg_p_p4(p,p4)
      call pwhg_getinvmass(p4,m)
      pwhg_invmass=m      
      end

      subroutine pwhg_getrapidity(p,y)
      implicit none
      real * 8 p(4),y
      real * 8 x
      if (p(4).le.0d0) then
c     should NEVER enter here
         write(*,*) 'WARNING!! pwhg_getrapidity called with energy <= 0'
         y=sign(1.d0,p(3))*1.d8
      else
         x=p(3)/p(4)
         if (x.ge.1d0.or.x.le.-1d0) then
            y = sign(1d0,x)*1d8
         else
            y=0.5*log((1+x)/(1-x))
         endif
      endif         
      end      


      function pwhg_rapidity(p)
      implicit none
      real * 8 pwhg_rapidity,p(0:3)
      real * 8 p4(4),y
      call pwhg_p_p4(p,p4)
      call pwhg_getrapidity(p4,y)
      pwhg_rapidity=y
      end


      subroutine pwhg_getpseudorapidity(p,eta)
      implicit none
      real * 8 p(4),eta
      real * 8 mod, costh      
      mod = sqrt(p(1)**2+p(2)**2+p(3)**2)
      if (mod.eq.0d0) then
c     should NEVER enter here
         write(*,*) 
     $        'WARNING!! pwhg_getpseudorapidity called with vector 0'
         eta=sign(1.d0,p(3))*1.d8
      else
         costh = p(3)/mod
         if (costh.ge.1d0.or.costh.le.-1d0) then
            eta = sign(1d0,costh)*1d8
         else
            eta=0.5*log((1+costh)/(1-costh))
         endif
      endif
      end

      function pwhg_pseudorapidity(p)
      implicit none
      real * 8 pwhg_pseudorapidity,p(0:3)
      real * 8 p4(4),eta
      call pwhg_p_p4(p,p4)
      call pwhg_getpseudorapidity(p4,eta)
      pwhg_pseudorapidity=eta
      end


      subroutine pwhg_getazi(p,azi)
      implicit none
      real * 8 p(4),azi
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      azi=atan2(p(2),p(1))
      return
c$$$
c$$$      if (p(1).ne.0d0) then
c$$$         azi = atan(p(2)/p(1))
c$$$      else
c$$$         azi = sign(pi/2,p(2))
c$$$      return
c$$$      endif      
c$$$      if (p(1).lt.0d0) then
c$$$         if (azi.gt.0d0) then               
c$$$            azi = azi - pi
c$$$         else
c$$$            azi = azi + pi
c$$$         endif
c$$$      endif   
      end

      function pwhg_azi(p)
      implicit none
      real * 8 p(0:3),pwhg_azi
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      real * 8 p4(4),azi
      call pwhg_p_p4(p,p4)
      call pwhg_getazi(p4,azi)
      pwhg_azi=azi
      end


c     compute the azimuthal angle between two momenta
      subroutine pwhg_getdelta_azi(p1,p2,delphi)
      implicit none
      real * 8 p1(4),p2(4),delphi
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      real * 8 phi1, phi2
      call pwhg_getazi(p1,phi1)
      call pwhg_getazi(p2,phi2)
      delphi = abs(phi1-phi2)
      if (delphi.gt.pi) then
         delphi = 2*pi-delphi
      endif
      if (delphi.lt.0 .or. delphi.gt.pi) then
         write(*,*) 'WARNING!! Problem in pwhg_getdelta_azi!!'
      endif
      end



c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 in azi and pseudorapidity
      subroutine pwhg_getR_phieta(p1,p2,R)
      implicit none
      real * 8 p1(4),p2(4),R
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      real * 8 eta1,eta2
      real * 8 delphi
      call pwhg_getdelta_azi(p1,p2,delphi)
      call pwhg_getpseudorapidity(p1,eta1)
      call pwhg_getpseudorapidity(p2,eta2)
      R = sqrt( (eta1-eta2)**2 + delphi**2 )
      end


c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 in azi and pseudorapidity
      function pwhg_R_phieta(p1,p2)
      implicit none
      real * 8 p1(0:3),p2(0:3),pwhg_R_phieta
      real * 8 p14(4),p24(4),R
      call pwhg_p_p4(p1,p14)
      call pwhg_p_p4(p2,p24)
      call pwhg_getR_phieta(p14,p24,R)
      pwhg_R_phieta=R
      end



c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 in azi and pseudorapidity
      subroutine pwhg_getR_phiy(p1,p2,R)
      implicit none
      real * 8 p1(4),p2(4),R
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      real * 8 y1,y2
      real * 8 delphi
      call pwhg_getdelta_azi(p1,p2,delphi)
      call pwhg_getrapidity(p1,y1)
      call pwhg_getrapidity(p2,y2)
      R = sqrt( (y1-y2)**2 + delphi**2 )
      end


c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 in azi and pseudorapidity
      function pwhg_R_phiy(p1,p2)
      implicit none
      real * 8 p1(0:3),p2(0:3),pwhg_R_phiy
      real * 8 p14(4),p24(4),R
      call pwhg_p_p4(p1,p14)
      call pwhg_p_p4(p2,p24)
      call pwhg_getR_phiy(p14,p24,R)
      pwhg_R_phiy=R
      end



      function pwhg_ptrel(pin,pjet)
      implicit none
      real * 8 pwhg_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      pwhg_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end



      subroutine pwhg_getrapidity_old(p,y)
      implicit none
      real * 8 p(4),y
      real * 8 tiny,xplus,xminus
      parameter (tiny=1.d-10)
      xplus=p(4)+p(3)
      xminus=p(4)-p(3)
      if(xplus.gt.tiny.and.xminus.gt.tiny) then
        if((xplus/xminus).gt.tiny) then
          y=0.5d0*log(xplus/xminus)
        else
c     SIGN(A,B) returns the value of A with the sign of B 
          y=sign(1.d0,p(3))*1.d8
        endif
      else
        y=sign(1.d0,p(3))*1.d8
      endif
      end      

