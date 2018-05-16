C$Modified: Wed Jul 14 10:34:34 2004 by puwer $
C*************************************************************************
      subroutine SoftConf(isoft,delta,pout)
      implicit none
      double precision pout(4,-2:4) 
      double precision delta
      double precision p(4,100) ! momenta
      double precision w        ! weight 
      double precision masses(100) ! masses
      double precision koff(4),kg(4),kqb(4),kgh(4),kqbh(4),vtmp
      double precision mt,moff,sqrts,dotp
      integer mu,isoft,err
      mt = 175.d0
      moff = mt+delta
      masses(1) = mt
      masses(2) = moff
      masses(3) = 0.d0
*     lw = 0
      sqrts = 500.d0

      call rambo(3, sqrts, masses,p, w)

      do mu=1,4
         pout(mu,1) = p(mu,1)
         koff(mu)   = p(mu,2)
         if (isoft.eq.4) then
            pout(mu,3) = p(mu,3)
         else
            pout(mu,4) = p(mu,3)
         endif
      enddo

      masses(1) = mt      
      masses(2) = 0.d0

      call rambo(2, moff, masses, p, w)
      do mu=1,4
         kqbh(mu) = p(mu,1)
         kgh(mu) = p(mu,2)
      enddo

      call boost(koff,kqbh,kqb,err)
      call boost(koff,kgh,kg,err)


      do mu=1,4
         pout(mu,2) = kqb(mu)
         if (isoft.eq.4) then
            pout(mu,4) = kg(mu)
         else
            pout(mu,3) = kg(mu)
         endif
      enddo
      masses(1) = 0.d0
      masses(2) = 0.d0 
      call rambo(2, sqrts, masses,p, w)
      do mu=1,4
         pout(mu,-1) = p(mu,1)
         pout(mu,-2) = p(mu,2)
      enddo
      end
