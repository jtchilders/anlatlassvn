
      subroutine bottomshifts(mt,mtrun,mqs,mts,mbs,ma,Ab,At,mg,mu,
     $     tb,q,v,alphas,asq,mw,mz,mbmb,mbmssm,mode,mbshift,msusybl)
      
      implicit none

      real*8 mt,mqs,mts,mbs,ma,Ab,At,mg,mu,tb,q,v,
     $     alphas,asq,mw,mz,mbmb,mbmssm,mbshift,msusybl
      real*8 T1,T2,B1,B2,sst,cst,ssb,csb,mstop2(2),msbot2(2)
      real*8 mt2,mg2,mu2,ma2,sb,cb,v2,runb
      real*8 db,eb,dhb,dbp,ebp,dhbp,dv,dmsbl,dmsblp
      real*8 gs,httree,hbtree,hbrunsm,hbrun,hbos,mbrun,mtrun,vrun,Pi

      integer mode

      pi = 3.1415926535898d0
      
      mt2=mt**2
      mg2=mg**2
      mu2=mu**2
      ma2=ma**2
      v2 = v**2

      sb = dsin(datan(tb))
      cb = dcos(datan(tb))

      gs = sqrt(4d0*pi*asq)

c
c     SM DRbar running mass at q      
c

      mbrun = runb(q,mbmb,mt,alphas,mz) 
c      if(q.eq.mt) mbrun = 2.67576012224198d0 ! CHECK FEYNHIGGS

c
c     stop masses and top Yukawa
c

      call diagonalize(mt,mw,mz,mqs,mts,At,mu,tb,1,0,mstop2,sst,cst)

      T1 = mstop2(1)
      T2 = mstop2(2)
      
c
c     sbottom masses and bottom Yukawa using the SM running mass
c

      call diagonalize(mbrun,mw,mz,mqs,mbs,Ab,mu,tb,2,0,msbot2,ssb,csb)

      B1 = msbot2(1)
      B2 = msbot2(2)

c$$$      write(*,*) 'normalTB'
c$$$      write(*,*) sqrt(T1),sqrt(T2),cst,sst
c$$$      write(*,*) sqrt(B1),sqrt(B2),csb,ssb

c
c     compute epsilon and epsilon'
c     
      
      httree = dsqrt(2d0)*mtrun/v/sb
      hbtree = dsqrt(2d0)*mbrun/v/cb
      
      eb = -2d0*mu*mg/(B1-B2)*tb*(B1/(B1-mg2)*Log(B1/mg2)-
     $     B2/(B2-mg2)*Log(B2/mg2))
      eb  = gs**2/12d0/Pi**2*eb
      
      ebp = -At*mu/(T1-T2)*tb*
     $     (T1/(T1-mu2)*Log(T1/mu2) - T2/(T2-mu2)*Log(T2/mu2))
      ebp = httree**2/16d0/pi**2*ebp
      
c
c     compute the remaining strong corrections
c     
      
      call shiftbotstr2(hbtree,mg,B1,B2,Ab,q**2,gs,db,dhb)
      
      call shiftmbsstr2(mt2,mg,T1,T2,B1,B2,sst,cst,
     $     ssb,csb,gs,dmsbl)
      
c     
c     squark masses in the large tanB limit:
c     

      call diagonalize(mt,mw,mz,mqs,mts,At,mu,tb,1,1,mstop2,sst,cst)

      T1 = mstop2(1)
      T2 = mstop2(2)

      call diagonalize(mbrun,mw,mz,mqs,mbs,Ab,mu,tb,2,1,msbot2,ssb,csb)

      B1 = msbot2(1)
      B2 = msbot2(2)

c$$$      write(*,*) 'largeTB'
c$$$      write(*,*) sqrt(T1),sqrt(T2),cst,sst
c$$$      write(*,*) sqrt(B1),sqrt(B2),csb,ssb

c
c     compute the remaining Yukawa corrections (in the large tanB limit)
c

      call shiftbotyuk2(mt2,hbtree,ma2,T1,T2,B1,B2,sst,cst,ssb,csb,
     $     mu,Ab,v2,q**2,dbp,dhbp,dv)

      call shiftmbsyuk2(mt2,hbtree,ma2,T1,T2,B1,B2,sst,cst,ssb,csb,
     $     Ab,mu,v2,q**2,dmsblp)
      
c
c     compute mbshift and msusybl
c

      write(*,*) 'bottomshifts',eb,ebp,db,dbp,dhb/hbtree,dhbp/hbtree,
     $     dmsbl,dmsblp,dv/v

      if(mode.eq.1) then        ! strong corrections only

         hbrunsm = sqrt(2d0)*mbrun/v/cb
         
         hbrun = hbrunsm*(1d0+db)/abs(1d0+eb)
         hbos = hbrun - dhb

         mbshift = hbos*v*cb/sqrt(2d0)
         mbmssm = hbrun*v*cb/sqrt(2d0)

         msusybl = sqrt(mqs**2+dmsbl)       

      elseif(mode.eq.2) then        ! strong + resummation O(at)

         eb = eb+ebp            

         hbrunsm = sqrt(2d0)*mbrun/v/cb
         
         hbrun = hbrunsm*(1d0+db)/abs(1d0+eb)
         hbos = hbrun - dhb

         mbshift = hbos*v*cb/sqrt(2d0)
         mbmssm = hbrun*v*cb/sqrt(2d0)

         msusybl = sqrt(mqs**2+dmsbl)       

      elseif(mode.eq.3) then        ! full Yukawa corrections

         eb = eb+ebp            
         db = db+dbp
         dhb = dhb+dhbp
         dmsbl = dmsbl+dmsblp   

         vrun = v+dv
         hbrunsm = sqrt(2d0)*mbrun/vrun/cb

         hbrun = hbrunsm*(1d0+db)/abs(1d0+eb)
         hbos = hbrun - dhb

         mbshift = hbos*v*cb/sqrt(2d0)
c        Fixed by Pietro on 04/06/2012 v->vrun
         mbmssm = hbrun*vrun*cb/sqrt(2d0)

         msusybl = sqrt(mqs**2+dmsbl)       

      else

         write(*,*) 'ERROR in bottomshifts - mode out of range'
         stop

      endif

      end

*
*******************************************************************
*

      subroutine shiftbotstr2(hb,mg,B1,B2,A,q,gs,db,dhb)
      
      implicit none      
      real*8 hb,g,mg,B1,B2,A,q,gs,db,dhb
      real*8 CF,Pi,c1

      pi = 3.1415926535898d0
      CF = 4d0/3d0
      g = mg**2

      db = 3d0/2d0-Log(g/q)+1d0/2d0*(
     $     B1/(g-B1)*(1d0-((2d0*g-B1)/(g-B1)
     $     -4d0*mg*A/(B1-B2))*Log(g/B1))+
     $     B2/(g-B2)*(1d0-((2d0*g-B2)/(g-B2)
     $     +4d0*mg*A/(B1-B2))*Log(g/B2)))

      c1 = -4d0 - 2d0*Log(g/q)
     $     + 2d0*B1/(B1-B2)*(2d0*Log(B1/q)-(1-g/B1)**2*Log(Abs(1-B1/g)))
     $     - 2d0*B2/(B1-B2)*(2d0*Log(B2/q)-(1-g/B2)**2*Log(Abs(1-B2/g)))

      dhb = hb*c1

      db  = gs**2*CF/16d0/Pi**2*db
      dhb = gs**2*CF/16d0/Pi**2*dhb 
      
      return
      end

*
*******************************************************************
*

      subroutine shiftbotyuk2(t,hb,A0,T1,T2,B1,B2,st,ct,sb,cb,
     $     mu,Ab,vv,q,db,dhb,dv)

      implicit none
      real*8 t,A0,T1,T2,B1,B2,st,ct,sb,cb,mu,vv,Ab,q,db
      real*8 mu2,ct2,st2,cb2,sb2,ht,hb,s2t,c2t,s2b,c2b
      real*8 v,dv,dmu,dhb,ffW
      real*8 Pi,Nc,lhiggs,hhiggs,higgsino
      real*8 Pi11b1,Pi22b2,Pi12b1,Pi12b2,fPi11b,fPi12b,myB1,myB0
      mu2 = mu**2

      pi = 3.1415926535898d0
      Nc = 3d0

      c2t = ct**2-st**2
      s2t = 2d0*ct*st
      c2b = cb**2-sb**2
      s2b = 2d0*cb*sb
      ct2 = ct**2
      st2 = st**2
      cb2 = cb**2
      sb2 = sb**2

      v = dsqrt(vv)

      ht = dsqrt(2d0/vv)*Sqrt(t) ! in the limit

      lhiggs = -ht**2/4d0*(5d0-6d0*Log(t/q))

      hhiggs = 2d0*ht**2*(t-A0+A0*Log(A0/q)-t*Log(t/q))/(t-A0)
     $     +hb**2*t*(t-A0+t*Log(A0/t))/2d0/(t-A0)**2
     $     +3d0*hb**2/4d0*(1-2d0*Log(A0/q))

      higgsino = hb**2*(myB1(0d0,mu2,B1,q)+myB1(0d0,mu2,B2,q))
     $     +(hb**2*ct2+ht**2*st2)*myB1(0d0,mu2,T1,q)
     $     +(hb**2*st2+ht**2*ct2)*myB1(0d0,mu2,T2,q)
     $     -2d0*ht**2*mu2/(T1-T2)*
     $     (myB0(0d0,mu2,T1,q)-myB0(0d0,mu2,T2,q))
      
      db = (lhiggs + hhiggs + higgsino)/32d0/pi**2
      
      dv = ht**2*Nc*v/16d0/pi**2*((2d0*Log(t/q)-1d0)/4d0
     $     -(ct2*cb2*ffW(T1,B1)+ct2*sb2*ffW(T1,B2)
     $     +st2*cb2*ffW(T2,B1)+st2*sb2*ffW(T2,B2))/t)

      dmu = 0d0

      Pi11b1 = fPi11b(B1,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2
      Pi22b2 = fPi11b(B2,t,hb,A0,T1,T2,B2,B1,s2t,c2t,-s2b,-c2b,q,
     $     mu,vv,Ab)/16d0/pi**2
      Pi12b1 = fPi12b(B1,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2
      Pi12b2 = fPi12b(B2,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2

      dhb = (hb*(c2b/s2b*(Pi12b1+Pi12b2)/(B1-B2)
     $     +(Pi11b1-Pi22b2)/(B1-B2) - dmu/mu - dv/v))

      return
      end

*     
***********************************************************************
*

      subroutine shiftmbsstr2(t,mg,T1,T2,B1,B2,stt,ctt,stb,ctb,gs,dmsbl)
      
      implicit none

      real*8 t,g,mt,mg,T1,T2,B1,B2,stt,ctt,stb,ctb,dmsbl
      real*8 Pi,gs,CF,msdr
      real*8 myB0,dmt1,dmt2,dmb1,dmb2,dmt,dtht,dthb
      real*8 pit11,pit22,pitt,pib11,pib22,dmQt,dmQb
      real*8 s2t,c2t,s2b,c2b,q

      pi = 3.1415926535898d0
      msdr = -5d0

      q = t*1d10                ! irrelevant!
      g = mg**2
      mt = dsqrt(t)

      s2t = 2d0*ctt*stt
      c2t = ctt**2 - stt**2

      s2b = 2d0*ctb*stb
      c2b = ctb**2 - stb**2

      CF = 4d0/3d0
      
c     top sector

      pitt =                    
     $     mt*(3*Log(t/q) + msdr + .5d0*(2*g/t*(Log(g/q)-1)
     $      -T1/t*(Log(T1/q)-1) - T2/t*(Log(T2/q)-1)
     $      +(g+t-T1 - 2*s2t*mt*mg)/t*myB0(t,g,T1,q)
     $      +(g+t-T2 + 2*s2t*mt*mg)/t*myB0(t,g,T2,q)))

      dmt = gs**2*CF/16/Pi**2*pitt

      pit11 =
     $     T1*(3*Log(T1/q) - 7 - c2t**2*(Log(T1/q)-1)
     $      -s2t**2*T2/T1*(Log(T2/q)-1) + 2*(
     $      g/T1*(Log(g/q)-1) + t/T1*(Log(t/q)-1)
     $      +(T1-g-t + 2*s2t*mt*mg)/T1*myB0(T1,t,g,q)))
      
      pit22 =
     $     T2*(3*Log(T2/q) - 7 - c2t**2*(Log(T2/q)-1)
     $      -s2t**2*T1/T2*(Log(T1/q)-1) + 2*(
     $      g/T2*(Log(g/q)-1) + t/T2*(Log(t/q)-1)
     $      +(T2-g-t - 2*s2t*mt*mg)/T2*myB0(T2,t,g,q)))
      
      dmt1 = gs**2*CF/16/Pi**2*pit11
      dmt2 = gs**2*CF/16/Pi**2*pit22

      dtht = (4d0*mt*mg*c2t*(myB0(T1,t,g,q)+myB0(T2,t,g,q)) +
     $     2d0*c2t*s2t*(T2*(1d0-Log(T2/q))-T1*(1d0-Log(T1/q))))/
     $     2d0/(T1-T2)
      dtht = gs**2*CF/16/Pi**2*dtht
      
      dmQt = dmt1*ctt**2+dmt2*stt**2-2d0*(T1-T2)*stt*ctt*dtht -
     $     2d0*mt*dmt

c     bottom sector

      pib11 =                     
     $     B1*(3*Log(B1/q) - 3 - c2b**2*(Log(B1/q)-1)
     $     -s2b**2*B2/B1*(Log(B2/q)-1)-6*g/B1
     $     -2*(1-2*g/B1)*Log(g/q)-2*(1-g/B1)**2*Log(Abs(1-B1/g)))

      pib22 =                     
     $     B2*(3*Log(B2/q) - 3 - c2b**2*(Log(B2/q)-1)
     $     -s2b**2*B1/B2*(Log(B1/q)-1)-6*g/B2
     $     -2*(1-2*g/B2)*Log(g/q)-2*(1-g/B2)**2*Log(Abs(1-B2/g)))

      dmb1 = gs**2*CF/16/Pi**2*pib11
      dmb2 = gs**2*CF/16/Pi**2*pib22

      dthb = c2b*s2b*(B2*(1d0-Log(B2/q))-B1*(1d0-Log(B1/q)))/(B1-B2)  
      dthb = gs**2*CF/16/Pi**2*dthb

      dmQb = dmb1*ctb**2+dmb2*stb**2-2d0*(B1-B2)*stb*ctb*dthb
      
c     final combination

      dmsbl = dmQt-dmQb
      
      return 
      end

*     
***********************************************************************
*

      subroutine shiftmbsyuk2(t,hb,A0,T1,T2,B1,B2,st,ct,sb,cb,
     $     Ab,mu,vv,qq,dmsbl)
      
      implicit none

      real*8 t,T1,T2,B1,B2,st,ct,sb,cb,mu,Ab,vv,qq,A0,dmsbl,Pi
      real*8 myB0,dmt1,dmt2,dmb1,dmb2,dmt,dtht,dthb,dmQt,dmQb
      real*8 s2t,c2t,s2b,c2b,mt,mu2,hb,ht,q
      real*8 Pi12t1,Pi12t2,fPi11t,fPi12t,Pi12b1,Pi12b2,fPi11b,fPi12b
      real*8 myB1

      pi = 3.1415926535898d0

      q = qq!*1d10

      mt = dsqrt(t)
      mu2=mu**2

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      s2b = 2d0*cb*sb
      c2b = cb**2 - sb**2

      ht = dsqrt(2d0/vv)*Sqrt(t) ! in the limit
      
c     top sector

      dmt = mt/2d0*(ht**2
     $     *(2d0*myB1(t,t,0d0,q)+myB1(t,mu2,T1,q)+myB1(t,mu2,T2,q))
     $     +ht**2*myB1(t,0d0,0d0,q)+hb**2*myB1(t,0d0,A0,q)
     $     +(ht**2*cb**2+hb**2*sb**2)*myB1(t,mu2,B1,q)
     $     +(ht**2*sb**2+hb**2*cb**2)*myB1(t,mu2,B2,q)
     $     -ht*hb*s2b*mu/mt*(myB0(t,mu2,B1,q)-myB0(t,mu2,B2,q)))
     $     /16d0/pi**2
      
      dmt1 = fPi11t(T1,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2

      dmt2 = fPi11t(T2,t,hb,A0,T2,T1,B1,B2,-s2t,-c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2

      Pi12t1 = fPi12t(T1,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2
      Pi12t2 = fPi12t(T2,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2

      dtht = (Pi12t1+Pi12t2)/2d0/(T1-T2)
      
      dmQt = dmt1*ct**2+dmt2*st**2-2d0*(T1-T2)*st*ct*dtht -
     $     2d0*mt*dmt

c     bottom sector

      dmb1 = fPi11b(B1,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2

      dmb2 = fPi11b(B2,t,hb,A0,T1,T2,B2,B1,s2t,c2t,-s2b,-c2b,q,
     $     mu,vv,Ab)/16d0/pi**2

      Pi12b1 = fPi12b(B1,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2
      Pi12b2 = fPi12b(B2,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)/16d0/pi**2

      dthb = (Pi12b1+Pi12b2)/2d0/(B1-B2)

      dmQb = dmb1*cb**2+dmb2*sb**2-2d0*(B1-B2)*sb*cb*dthb
      
c     final combination

      dmsbl = dmQt-dmQb

      return 
      end

*
*********************************************************************
*

      function fPi11t(p,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)

      implicit none

      real*8 fPi11t,p,t,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv
      real*8 ht,hb,mt,ct2,st2,cb2,sb2,At,Ab,mu2
      real*8 Nc
      real*8 myG,myB0,myB1,myAA
      external myG,myB0,myB1,myAA

      Nc=3d0
      mt = dsqrt(t)
      ht = dsqrt(2d0/vv)*mt

      mu2 = mu**2
      
      At = (T1-T2)*s2t/2d0/mt    

      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0
      cb2 = (1d0+c2b)/2d0
      sb2 = (1d0-c2b)/2d0

      fPi11t = (ht**2*(myG(p,t,mu2,q)
     $     +st2*cb2*myAA(B1,q)+st2*sb2*myAA(B2,q)
     $     +(c2t**2-(Nc-1d0)/2d0*s2t**2)*myAA(T2,q)
     $     +(Nc+1d0)/2d0*s2t**2*myAA(T1,q)
     $     +((2d0*mt+s2t*At)**2*myB0(p,T1,0d0,q)
     $     + mu**2*s2t**2*myB0(p,T1,A0,q)
     $     +(1d0+c2t**2)*At**2*myB0(p,T2,0d0,q)
     $     +(1d0+c2t**2)*mu**2*myB0(p,T2,A0,q))/2d0)
     $     +(ht**2*st2+hb**2*ct2)*myG(p,0d0,mu2,q)
     $     +hb**2*ct2*(myAA(A0,q)+sb2*myAA(B1,q)+cb2*myAA(B2,q))
     $     +(hb**2*sb2*ct2*mu**2
     $     +ht**2*cb2*(ct2*mt**2+s2t*mt*At+st2*At**2)
     $     -ht*hb*s2b*mu*(s2t*At/2d0+mt*ct2))*myB0(p,B1,0d0,q)
     $     +(hb**2*sb2*(st2*mt**2+s2t*mt*Ab+ct2*Ab**2)
     $     +ht**2*mu**2*st2*cb2	
     $     -ht*hb*(s2b*s2t*Ab*mu/2d0+s2b*mt*st2*mu))*myB0(p,B1,A0,q)
     $     +(hb**2*cb2*ct2*mu**2
     $     +ht**2*sb2*(ct2*mt**2+s2t*mt*At+st2*At**2)
     $     +ht*hb*s2b*mu*(s2t*At/2d0+mt*ct2))*myB0(p,B2,0d0,q)
     $     +(hb**2*cb2*(st2*mt**2+s2t*mt*Ab+ct2*Ab**2)
     $     +ht**2*mu**2*st2*sb2	
     $     +ht*hb*s2b*mu*(s2t*Ab/2d0+mt*st2))*myB0(p,B2,A0,q))

      return 
      end

*
*********************************************************************
*

      function fPi12t(p,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)

      implicit none

      real*8 fPi12t,p,t,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv
      real*8 ht,hb,mt,cb2,sb2,At,Ab,mu2
      real*8 Nc
      real*8 myG,myB0,myB1,myAA
      external myG,myB0,myB1,myAA

      Nc=3d0
      mt = dsqrt(t)
      ht = dsqrt(2d0/vv)*mt

      mu2 = mu**2
      
      At = (T1-T2)*s2t/2d0/mt    

      cb2 = (1d0+c2b)/2d0
      sb2 = (1d0-c2b)/2d0

      fPi12t = (ht**2/2d0*((Nc+1d0)*c2t*s2t*(myAA(T1,q)-myAA(T2,q))
     $     +s2t*cb2*myAA(B1,q)+s2t*sb2*myAA(B2,q)
     $     +c2t*At*(2d0*mt+s2t*At)*myB0(p,T1,0d0,q)
     $     +c2t*s2t*mu**2*myB0(p,T1,A0,q)
     $     +c2t*At*(2d0*mt-s2t*At)*myB0(p,T2,0d0,q)
     $     -c2t*s2t*mu**2*myB0(p,T2,A0,q))
     $     +(ht**2-hb**2)/2d0*s2t*myG(p,0d0,mu2,q)
     $     -hb**2/2d0*s2t*(myAA(A0,q)+sb2*myAA(B1,q)+cb2*myAA(B2,q))
     $     +((ht**2*cb2*(s2t*At**2-s2t*mt**2+2d0*c2t*mt*At)
     $     -hb**2*s2t*sb2*mu**2
     $     -ht*hb*s2b*mu*(c2t*At-s2t*mt))*myB0(p,B1,0d0,q)
     $     +(ht**2*cb2*s2t*mu**2
     $     -hb**2*(s2t*sb2*Ab**2-s2t*sb2*mt**2-2d0*c2t*sb2*mt*Ab)
     $     -ht*hb*s2b*mu*(c2t*Ab+s2t*mt))*myB0(p,B1,A0,q))/2d0
     $     +((ht**2*sb2*(s2t*At**2-s2t*mt**2+2d0*c2t*mt*At)
     $     -hb**2*s2t*cb2*mu**2
     $     +ht*hb*s2b*mu*(c2t*At-s2t*mt))*myB0(p,B2,0d0,q)
     $     +(ht**2*sb2*s2t*mu**2
     $     -hb**2*(s2t*cb2*Ab**2-s2t*cb2*mt**2-2d0*c2t*cb2*mt*Ab)
     $     +ht*hb*s2b*mu*(c2t*Ab+s2t*mt))*myB0(p,B2,A0,q))/2d0)

      return 
      end

*
*********************************************************************
*

      function fPi11b(p,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)

      implicit none

      real*8 fPi11b,p,t,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv
      real*8 ht,hb,mt,ct2,st2,cb2,sb2,At,Ab,mu2
      real*8 Nc
      real*8 myG,myB0,myB1,myAA
      external myG,myB0,myB1,myAA

      Nc=3d0
      mt = dsqrt(t)
      ht = dsqrt(2d0/vv)*mt

      mu2 = mu**2
      
      At = (T1-T2)*s2t/2d0/mt    

      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0
      cb2 = (1d0+c2b)/2d0
      sb2 = (1d0-c2b)/2d0

      fPi11b = (hb**2*(myG(p,0d0,mu2,q)+(1d0+sb2)*myAA(A0,q)
     $     +sb2*ct2*myAA(T1,q)+sb2*st2*myAA(T2,q)
     $     +(c2b**2-(Nc-1d0)/2d0*s2b**2)*myAA(B2,q)
     $     +(Nc+1d0)/2d0*s2b**2*myAA(B1,q)
     $     +(s2b**2*mu**2*myB0(p,B1,0d0,q)
     $     +s2b**2*Ab**2*myB0(p,B1,A0,q)
     $     +(1d0+c2b**2)*mu**2*myB0(p,B2,0d0,q)
     $     +(1d0+c2b**2)*Ab**2*myB0(p,B2,A0,q))/2d0)
     $     +(hb**2*sb2+ht**2*cb2)*myG(p,t,mu2,q)
     $     +2d0*ht*hb*mt*mu*s2b*myB0(p,t,mu2,q)
     $     +ht**2*cb2*(st2*myAA(T1,q)+ct2*myAA(T2,q))
     $     +(hb**2*sb2*ct2*mu**2
     $     +ht**2*cb2*(ct2*mt**2+s2t*mt*At+st2*At**2)
     $     -ht*hb*s2b*mu*(s2t*At/2d0+mt*ct2))*myB0(p,T1,0d0,q)
     $     +(hb**2*sb2*(st2*mt**2+s2t*mt*Ab+ct2*Ab**2)
     $     +ht**2*mu**2*st2*cb2	
     $     -ht*hb*(s2b*s2t*Ab*mu/2d0+s2b*mt*st2*mu))*myB0(p,T1,A0,q)
     $     +(hb**2*sb2*st2*mu**2
     $     +ht**2*cb2*(st2*mt**2-s2t*mt*At+ct2*At**2)
     $     -ht*hb*s2b*mu*(-s2t*At/2d0+mt*st2))*myB0(p,T2,0d0,q)
     $     +(hb**2*sb2*(ct2*mt**2-s2t*mt*Ab+st2*Ab**2)
     $     +ht**2*mu**2*ct2*cb2	
     $     -ht*hb*(-s2b*s2t*Ab*mu/2d0+s2b*mt*ct2*mu))*myB0(p,T2,A0,q))

      return 
      end

*
*********************************************************************
*

      function fPi12b(p,t,hb,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,
     $     mu,vv,Ab)

      implicit none

      real*8 fPi12b,p,t,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv
      real*8 ht,hb,mt,ct2,st2,At,Ab,mu2
      real*8 Nc
      real*8 myG,myB0,myB1,myAA
      external myG,myB0,myB1,myAA

      Nc=3d0
      mt = dsqrt(t)
      ht = dsqrt(2d0/vv)*mt

      mu2 = mu**2
      
      At = (T1-T2)*s2t/2d0/mt    

      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0

      fPi12b = (hb**2/2d0*((Nc+1d0)*c2b*s2b*(myAA(B1,q)-myAA(B2,q))
     $     +s2b*(myAA(A0,q)+ct2*myAA(T1,q)+st2*myAA(T2,q))
     $     +c2b*s2b*mu**2*(myB0(p,B1,0d0,q)-myB0(p,B2,0d0,q))
     $     +c2b*s2b*Ab**2*(myB0(p,B1,A0,q)-myB0(p,B2,A0,q)))
     $     +(hb**2-ht**2)/2d0*s2b*myG(p,t,mu2,q)
     $     +2d0*ht*hb*mt*mu*c2b*myB0(p,t,mu2,q)
     $     -ht**2/2d0*s2b*(st2*myAA(T1,q)+ct2*myAA(T2,q))
     $     +((hb**2*s2b*ct2*mu**2
     $     -ht**2*s2b*(ct2*mt**2+s2t*mt*At+st2*At**2)
     $     -ht*hb*(c2b*s2t*mu*At+2d0*c2b*ct2*mu*mt))*myB0(p,T1,0d0,q)
     $     +(hb**2*s2b*(st2*mt**2+s2t*mt*Ab+ct2*Ab**2)
     $     -ht**2*s2b*st2*mu**2
     $     -ht*hb*(c2b*s2t*mu*Ab+2d0*c2b*st2*mu*mt))*myB0(p,T1,A0,q))
     $     /2d0
     $     +((hb**2*s2b*st2*mu**2
     $     -ht**2*s2b*(st2*mt**2-s2t*mt*At+ct2*At**2)
     $     -ht*hb*(-c2b*s2t*mu*At+2d0*c2b*st2*mu*mt))*myB0(p,T2,0d0,q)
     $     +(hb**2*s2b*(ct2*mt**2-s2t*mt*Ab+st2*Ab**2)
     $     -ht**2*s2b*ct2*mu**2
     $     -ht*hb*(-c2b*s2t*mu*Ab+2d0*c2b*ct2*mu*mt))*myB0(p,T2,A0,q))
     $     /2d0)

      return 
      end

*
***********************************************************************
*

      real*8 function myAA(m,q)      
      real*8 m,q

      if(m.ne.0d0) then
         myAA = m*(1d0-Log(m/q))
      else
         myAA = 0d0
      endif

      return
      end

c$$$*
c$$$***********************************************************************
c$$$*
c$$$
c$$$
c$$$      real*8 function myB0(q,m1,m2,mu2) 
c$$$
c$$$c     from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104.
c$$$      
c$$$      real*8 q,m1,m2,Omega,mu2
c$$$
c$$$      if(q.eq.0d0) then
c$$$
c$$$         if(m1.eq.0d0.and.m2.ne.0d0) then
c$$$            myB0 = 1d0-Log(m2/mu2)
c$$$         elseif(m1.ne.0d0.and.m2.eq.0d0) then
c$$$            myB0 = 1d0-Log(m1/mu2)
c$$$         elseif(m1.eq.m2) then
c$$$            myB0 = -Log(m1/mu2)
c$$$         else
c$$$            myB0 = 1d0 - Log(m2/mu2) + m1/(m1-m2)*Log(m2/m1)
c$$$         endif
c$$$         
c$$$      else
c$$$
c$$$         if(m1.eq.0d0.and.m2.ne.0d0) then
c$$$            
c$$$            if(m2.ne.q) then
c$$$               myB0 = -(Log(m2/mu2)-2-(m2/q-1d0)*Log(dabs(1d0-q/m2))) 
c$$$            else 
c$$$               myB0 = -(Log(m2/mu2) - 2)
c$$$            endif
c$$$            
c$$$         elseif(m2.eq.0d0.and.m1.ne.0d0) then
c$$$            
c$$$            if(m1.ne.q) then
c$$$               myB0 = -(Log(m1/mu2)-2-(m1/q-1d0)*Log(dabs(1d0-q/m1))) 
c$$$            else
c$$$               myB0 = -(Log(m1/mu2) - 2)
c$$$            endif
c$$$            
c$$$         elseif(m2.eq.0d0.and.m1.eq.0d0) then
c$$$            
c$$$            myB0 = -(Log(q/mu2) - 2) ! cut the imaginary part (I Pi)
c$$$            
c$$$         else
c$$$            
c$$$            myB0 = -( dlog(q/mu2)-2.d0 + 
c$$$     1           1.d0/2.d0*( 1.d0 + (m1/q-m2/q))*dlog(m1/q) +
c$$$     2           1.d0/2.d0*( 1.d0 - (m1/q-m2/q))*dlog(m2/q) +
c$$$     3           2.d0*Omega(m1/q,m2/q))
c$$$            
c$$$         endif
c$$$         
c$$$      endif
c$$$
c$$$      return
c$$$      end
c$$$      
c$$$c     function Omega(a,b) contained in myB0
c$$$      real*8 function Omega(a,b)
c$$$      real*8 a,b,cbig
c$$$      Cbig = (a+b)/2.d0 - (a-b)**2.d0/4.d0 -1.d0/4.d0
c$$$      if(Cbig.gt.0.d0) then
c$$$         Omega = dsqrt(Cbig)*
c$$$     1        (datan((1.d0 + a - b)/(2.d0*dsqrt(Cbig))) +
c$$$     2        datan((1.d0 - a + b)/(2.d0*dsqrt(Cbig))) )
c$$$      elseif(Cbig.lt.0d0) then
c$$$         Cbig = - Cbig
c$$$         Omega = 1.d0/2.d0*dsqrt(Cbig)*
c$$$     1        dlog((a/2.d0 +b/2.d0 -1.d0/2.d0 -dsqrt(Cbig))/
c$$$     2        (a/2.d0 + b/2.d0 -1.d0/2.d0 + dsqrt(Cbig)))
c$$$      else
c$$$         Omega = 0         
c$$$      endif
c$$$
c$$$      return
c$$$      end

*
**********************************************************************
*
      
      real*8 function myB1(p,m1,m2,q)

      implicit none

      real*8 p,m1,m2,q
      real*8 myAA,myB0
      
      if(p.eq.0d0) then
         myB1 = (1d0-Log(m2/q)+m1**2/(m1-m2)**2*Log(m2/m1)
     $        +(m1+m2)/(m1-m2)/2d0)/2d0
      else
         myB1 = (myAA(m2,q)-myAA(m1,q)+(p+m1-m2)*myB0(p,m1,m2,q))/2d0/p
      endif
      
      return
      end

*
**********************************************************************
*
      
      real*8 function myG(p,m1,m2,q)

      implicit none

      real*8 p,m1,m2,q
      real*8 myAA,myB0

      myG = (p-m1-m2)*myB0(p,m1,m2,q) - myAA(m1,q) - myAA(m2,q)

      return
      end
      
*
*********************************************************************
*

      function ffW(x,y)

      implicit none
      real*8 ffw,x,y

      ffw = (x+y)/4d0 -x*y/2d0/(x-y)*Log(x/y)

      return
      end

c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      double precision function asf(x,tc,asc,zm)

c     compute the running (SM, MSbar) alpha_s
      
      implicit double precision (a-h,o-z)

      pi = 4d0*atan(1d0)

      fn=5d0
      
      b0=11d0-2*fn/3d0
      b1=102d0-38d0*fn/3d0
      vvv=1-b0*asc/(2d0*pi)*log(zm/x)
      
      if(x.le.tc) then          ! five flavors

         asf=asc/vvv*(1-b1/b0*asc/(4d0*pi*vvv)*log(vvv))
         
      else                      

         vvv=1-b0*asc/(2d0*pi)*log(zm/tc) ! first evolve up to q=mt
         
         ast=asc/vvv*(1-b1/b0*asc/(4*pi*vvv)*log(vvv))
         
         b0t=b0-2d0/3d0         ! six flavours
         b1t=b1-38d0/3d0
         vvv=1-b0t*ast/(2*pi)*log(tc/x) !     now evolve up to the scale >mt

         asf=ast/vvv*(1-b1t/b0t*ast/(4*pi*vvv)*log(vvv))

      endif
      
      return
      end

c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      double precision function runb(x,mbmb,mt,asmz,mz)

c     compute the running (SM, DRbar) bottom mass
 
      implicit none
      
      double precision x,mbmb,mt,asmz,mz
      double precision pi,b0,b1,g0,g1,asx,asb,asf,ast,rrr,corr,mbmt
      integer nf

      pi=4*atan(1d0)
      
      if(x.le.mt) then

         nf = 5                 ! evolve from mb to x with nf=5
         b0 = 11-2*nf/3d0
         b1 = 102-38*nf/3d0

         g0 = 8d0
         g1 = 404d0/3d0-40*nf/9d0
      
         asx=asf(x,mt,asmz,mz)
         asb=asf(mbmb,mt,asmz,mz)

         rrr=mbmb*(asx/asb)**(g0/(2*b0))      
         corr=1+asb*g0/(4*pi*2*b0)*(-b1/b0+g1/g0)*(asx/asb-1)
      
         runb=rrr*corr

      else

         nf = 5                 ! first evolve from mb to mt with nf=5
         b0 = 11-2*nf/3d0
         b1 = 102-38*nf/3d0

         g0 = 8d0
         g1 = 404d0/3d0-40*nf/9d0
      
         ast=asf(mt,mt,asmz,mz)
         asb=asf(mbmb,mt,asmz,mz)

         rrr=mbmb*(ast/asb)**(g0/(2*b0))      
         corr=1+asb*g0/(4*pi*2*b0)*(-b1/b0+g1/g0)*(ast/asb-1)
      
         mbmt=rrr*corr
         
         nf = 6                 ! then evolve from mt to x with nf=6
         b0 = 11-2*nf/3d0
         b1 = 102-38*nf/3d0

         g0 = 8d0
         g1 = 404d0/3d0-40*nf/9d0
      
         asx=asf(x,mt,asmz,mz)

         rrr=mbmt*(asx/ast)**(g0/(2*b0))      
         corr=1+ast*g0/(4*pi*2*b0)*(-b1/b0+g1/g0)*(asx/ast-1)
      
         runb=rrr*corr

      endif

      runb = runb*(1d0-asx/pi/3d0-(asx/pi)**2*29d0/72d0) ! shift to DRbar
c      runb = runb*(1d0-asx/pi/3d0) ! 1-loop only (as in FeynHiggs)

      return
      end

c$$$*
c$$$**********************************************************************
c$$$*
c$$$
c$$$      subroutine topshift(mt,mstl,mstr,At,mg,mu,tanb,qq,as,mw,mz,mtrun)
c$$$
c$$$      implicit none
c$$$      
c$$$      real*8 mt,mstl,mstr,At,mg,mu,tanb,qq,as,mw,mz,mtrun
c$$$      real*8 T1,T2,mstop2(2),myB0,st,ct,s2t,pi,dmt,g,t,q,msdr,CF
c$$$
c$$$      pi = 3.1415926535898d0
c$$$      msdr = -5d0
c$$$      CF = 4d0/3d0
c$$$
c$$$      g = mg**2
c$$$      t = mt**2
c$$$      q = qq**2
c$$$
c$$$      call diagonalize(mt,mw,mz,mstl,mstr,At,mu,tanb,1,0,mstop2,
c$$$     $     st,ct)
c$$$      
c$$$      T1 = mstop2(1)
c$$$      T2 = mstop2(2)
c$$$      s2t = 2*st*ct
c$$$
c$$$      dmt =                     ! eq. (B2) of DSZ
c$$$     $     CF*mt*(3*Log(t/q) + msdr + .5d0*(2*g/t*(Log(g/q)-1)
c$$$     $     -T1/t*(Log(T1/q)-1) - T2/t*(Log(T2/q)-1)
c$$$     $     +(g+t-T1 - 2*s2t*mg*mt)/t*myB0(t,g,T1,q)
c$$$     $     +(g+t-T2 + 2*s2t*mg*mt)/t*myB0(t,g,T2,q)))
c$$$      
c$$$      dmt = dmt*as/4d0/pi
c$$$
c$$$      mtrun = mt + dmt
c$$$
c$$$      return
c$$$      end
c$$$
