C     $Modified: Wed Jul 14 15:54:12 2004 by puwer $
      parameter (nf = 5.d0, 
     &     g_q = 3.d0/2.d0*CF,
     &     g_g = 11.d0/6.d0 * CA - 2.d0/3.d0 * Tr * Nf, 
     &     K_g = (67.d0/18.d0-pi**2/6.d0)*CA-10.d0/9.d0*TR*nf,
     &     K_q = (7.d0/2.d0 - pi**2/6.d0)*CF,
     &     Tgq = CA,
     &     Tqq = CF)
      theta(z) = (1.d0 + dsign(1.d0,z))/2.d0
      pqqreg(z) = -CF * (1.d0+z)
      pgqreg(z) = Tr*(z*z+(1.d0-z)**2)
      pqgreg(z) = CF * (1.d0+(1.d0-z)**2)/z
      pggreg(z) = 2.d0 * CA * ((1.d0-z)/z-1.d0+z*(1.d0-z))
      
      pggnreg(z) = 2.d0 * Ca /(1.d0 -z)

      psqq(z) = CF*(1.d0 -z)
      psqg(z) = CF*z
      psgq(z) = 2*Tr * z *(1.d0-z)
      psgg(z) = 0.d0

      Kgq_q_reg(z,sv,mj) = 0.d0
      Kgq_q_p(z,sv,mj) = 0.d0

      Kqq_q_reg(z,sv,mj) = -2.d0*dlog(2.d0-z)/(1.d0-z)


      Kqg_q_reg(z,sv,mj) = 2.d0*CF/CA*mj**2/z/sv
     &     * dlog(mj**2/((1.d0-z)*sv+mj**2))

      Kgg_q_reg(z,sv,mj) = Kqq_q_reg(z,sv,mj)
     &     + CA/CF*Kqg_q_reg(z,sv,mj)


      Kgg_g_reg(z,sv,mf) = 0.d0 
      Kgg_g_p(z,sv,mf) = -g_g/Ca/(1.d0-z) 


      Kbargg_reg(z) = pggreg(z)*dlog((1.d0-z)/z) + psgg(z)
      Kbargg_p(z) = CA * 2.d0/(1.d0-z)*dlog((1.d0-z)/z)
