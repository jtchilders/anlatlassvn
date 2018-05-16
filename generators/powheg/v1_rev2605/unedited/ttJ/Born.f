      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nlegs,mu,j,k
      parameter (nlegs=nlegborn)
      double precision p(0:3,nlegs),bornjk(nlegs,nlegs)
      integer bflav(nlegs)
      double precision bmunu(0:3,0:3,nlegs),born
      double precision amp2
      logical ini
      data ini/.true./
      save ini
cccccccccccccccccccc
      include 'pwhg_st.h'
      logical debug
      parameter(debug=.false.)
      logical compare_with_DUW_paper
      parameter(compare_with_DUW_paper=.false.)
      double precision kr_mad(0:3,nlegs),amp2MNR,tmp,amp2alt
      integer ileg
      double precision tiny
      parameter (tiny=1d-3)
      logical alt
      parameter (alt=.false.)
c
      if(ini) then
c     Set MADGRAPH parameters at first call
         call set_madgraph_parameters
      endif
      ini=.false.     
      amp2 = 0d0
      if (abs(bflav(3)).le.5.or.abs(bflav(4)).le.5) then
         write(*,*) 'setborn: ERROR in flavor assignement'
         write(*,*) bflav(1),' ',bflav(2),'->',bflav(3),' ',bflav(4),' '
     $        ,bflav(5)
         call exit(1)
      endif

ccccccccccccccccccc
      if(compare_with_DUW_paper) then
         call check_DUW_paper_born
      endif
cccccccccccccccccccccc

c     Reset MADGRAPH parameters at each call, only those 
c     that depends on the kinematics
      call mad_setparam
      
      do ileg=1,nlegs
         do mu=0,3
            kr_mad(mu,ileg)=p(mu,ileg)
         enddo
      enddo

c     to avoid bugs in HELAS, restore exact masslessness of  incoming partons 
      kr_mad(0,1)=dabs(kr_mad(3,1))
      kr_mad(0,2)=dabs(kr_mad(3,2))

      call compborn(kr_mad,bflav,amp2,bornjk,bmunu)
      
      if(amp2.eq.0d0) then
        write(*,*) 'WARNING setborn: returning 0 amplitude!'
        write(*,*) bflav(1),' ',bflav(2),'->',bflav(3),' ',bflav(4),' '
     $        ,bflav(5)
      endif

      if(amp2.lt.0d0) then
        write(*,*) 'ERROR setborn: returning squared amplitude < 0!'
        write(*,*) bflav(1),' ',bflav(2),'->',bflav(3),' ',bflav(4),' '
     $        ,bflav(5)
        call exit(0);
      endif
      
      

      born=amp2
c      print *,((bornjk,j=1,nlegborn),k=1,nlegborn)
      do j=1,nlegs
         do k=j+1,nlegs
            if(abs(bflav(k)).le.6) then
               if (bornjk(k,j).ne.bornjk(j,k)) then
                  write(*,*) 'Error: bornjk NOT symmetric'
                  call exit(1)
               endif
            endif
         enddo
      enddo


c compare with MNR and alternative code
      if (debug) then
         if(ini) then 
            if(alt) then
               write(*,*) "BORN: PERFORMING CHECK AGAINST ALT. CODE"
            else
               write(*,*) "BORN: PERFORMING CHECK AGAINST MNR CODE"
            endif
         endif
         
         if(.not.alt) then
            call MNR_ttbarj(p,bflav,amp2MNR)
         else
            call alt_ttqqg(p,bflav,amp2alt)
            amp2MNR=amp2alt
         endif

         tmp=abs(amp2MNR/amp2 -1d0)
         
         if  (tmp.gt.tiny) then 
            write(*,*) bflav,'\n BORN: MUST BE EQUAL =====> ', amp2MNR,
     $           amp2,'\n RATIO: ', amp2MNR/amp2
            
            write(*,*) "P:",((kr_mad(mu,ileg),mu=0,3),"\n",ileg=1
     $           ,nlegs)
            stop
         endif
      
         tmp=(abs(amp2MNR-amp2)/amp2MNR)
         
         if  (tmp.gt.tiny) then 
            write(*,*) bflav,'\n BORN: % DIFFERENCE MUST BE 0 =====> ',
     $           abs(amp2MNR-amp2)/amp2MNR
            stop
         endif
      endif
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC    CHECKS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_DUW_paper_born
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nlegs,mu,nu,j,k
      parameter (nlegs=nlegborn)
      integer bflav(nlegs)
      double precision born,bornjk(nlegs,nlegs),bmunu(0:3,0:3,nlegs)
      include 'pwhg_st.h'
      double precision kr_mad(0:3,nlegs),amp2MNR,tmp,amp2alt
      double precision uwer_res,amp2
      integer ileg
      double precision tiny
      parameter (tiny=1d-3)
c colored output
      character*1 red(5)
      character*1 reset(4) 
      data red /' ','[','3','1','m'/ 
      data reset /' ','[','0','m'/ 
      red(1) = char(27)     ! Escape character (ASCII 27).
      reset(1) = char(27)

      print *,red,"###########################################",reset
      print *,red,"###########################################",reset
      print *,red," CHECK BORN ROUTINES AGAINST DUW PAPER ",reset

      print *,red," AND AGAINST MNR AND ALTERNATIVE CODES ",reset

      print *,red,"###########################################",reset 
      print *,red,"###########################################",reset
      write(*,*) "P:"
c     NUMBERS to compare with P.Uwer  paper
      kr_mad(0,1)=500d0
      kr_mad(1,1)=0d0
      kr_mad(2,1)=0d0
      kr_mad(3,1)=500d0
c
      kr_mad(0,2)=500d0
      kr_mad(1,2)=0d0
      kr_mad(2,2)=0d0
      kr_mad(3,2)=-500d0
c
      kr_mad(0,3)=458.5331753852783d0
      kr_mad(1,3)=207.0255169909440d0
      kr_mad(2,3)=0d0
      kr_mad(3,3)=370.29327328961670
c                                   
      kr_mad(0,4)=206.6000026080000d0
      kr_mad(1,4)=-10.656936772525890d0
      kr_mad(2,4)=42.52372780926147d0
      kr_mad(3,4)=-102.3998210421085d0
c                                   
      kr_mad(0,5)=334.8668220067217d0
      kr_mad(1,5)=-196.3685802184181d0
      kr_mad(2,5)=-42.52372780926147d0
      kr_mad(3,5)=-267.8934522475083d0

      write(*,*) ((kr_mad(mu,ileg),mu=0,3),'\n',ileg=1,nlegs)

      call checkmomzeronew(nlegs,kr_mad)
      do ileg=1,nlegs
         call checkmassnew(ileg,kr_mad)
      enddo

      bflav(1)=0
      bflav(2)=0
      bflav(3)=6
      bflav(4)=-6
      bflav(5)=0

      call compborn(kr_mad,bflav,amp2,bornjk,bmunu)
      
      call MNR_ttbarj(kr_mad,bflav,amp2MNR)
      
      amp2alt=0

      uwer_res=0.6566843362709775d-03

      print *,'as =',st_alpha
      print *,bflav
      print *,' UWER, MADGRAPH, MNR, alt'
      print *, uwer_res,amp2/(4d0*pi*st_alpha)**3,amp2MNR /(4d0*pi
     $     *st_alpha)**3, amp2alt /(4d0*pi*st_alpha)**3
      print *,red,' RATIO MUST BE 1 =====> ',amp2/(4d0*pi*st_alpha)**3
     $     /uwer_res,amp2MNR/(4d0*pi*st_alpha)**3/uwer_res,amp2alt/(4d0
     $     *pi*st_alpha)**3/uwer_res,reset

      bflav(1)=1
      bflav(2)=-1
      bflav(3)=6
      bflav(4)=-6
      bflav(5)=0

      call compborn(kr_mad,bflav,amp2,bornjk,bmunu)

      call MNR_ttbarj(kr_mad,bflav,amp2MNR)

      call alt_ttqqg(kr_mad,bflav,amp2alt)

      uwer_res=0.5790368001550938d-04

      print *,'as =',st_alpha
      print *,bflav
      print *,' UWER, MADGRAPH, MNR, alt'
      print *, uwer_res,amp2/(4d0*pi*st_alpha)**3,amp2MNR /(4d0*pi
     $     *st_alpha)**3, amp2alt /(4d0*pi*st_alpha)**3
      print *,red,' RATIO MUST BE 1 =====> ',amp2/(4d0*pi*st_alpha)**3
     $     /uwer_res,amp2MNR/(4d0*pi*st_alpha)**3/uwer_res,amp2alt/(4d0
     $     *pi*st_alpha)**3/uwer_res,reset
     
  
      bflav(1)=1
      bflav(2)=0
      bflav(3)=6
      bflav(4)=-6
      bflav(5)=1

      call compborn(kr_mad,bflav,amp2,bornjk,bmunu)
      
      call MNR_ttbarj(kr_mad,bflav,amp2MNR)

      call alt_ttqqg(kr_mad,bflav,amp2alt)
      
      uwer_res=0.1607845322071585d-04

      print *,'as =',st_alpha
      print *,bflav
      print *,' UWER, MADGRAPH, MNR, alt'
      print *, uwer_res,amp2/(4d0*pi*st_alpha)**3,amp2MNR /(4d0*pi
     $     *st_alpha)**3, amp2alt /(4d0*pi*st_alpha)**3
      print *,red,' RATIO MUST BE 1 =====> ',amp2/(4d0*pi*st_alpha)**3
     $     /uwer_res,amp2MNR/(4d0*pi*st_alpha)**3/uwer_res,amp2alt/(4d0
     $     *pi*st_alpha)**3/uwer_res,reset

      bflav(1)=0
      bflav(2)=-1
      bflav(3)=6
      bflav(4)=-6
      bflav(5)=-1

      call compborn(kr_mad,bflav,amp2,bornjk,bmunu)

      call MNR_ttbarj(kr_mad,bflav,amp2MNR)

      call alt_ttqqg(kr_mad,bflav,amp2alt)
      
      uwer_res=0.2603527972645620d-03

      print *,'as =',st_alpha
      print *,bflav
      print *,' UWER, MADGRAPH, MNR, alt'
      print *, uwer_res,amp2/(4d0*pi*st_alpha)**3,amp2MNR /(4d0*pi
     $     *st_alpha)**3, amp2alt /(4d0*pi*st_alpha)**3
      print *,red,' RATIO MUST BE 1 =====> ',amp2/(4d0*pi*st_alpha)**3
     $     /uwer_res,amp2MNR/(4d0*pi*st_alpha)**3/uwer_res,amp2alt/(4d0
     $     *pi*st_alpha)**3/uwer_res,reset

      stop "COMPARISON WITH P.UWER NUMBERS FINISHED. PROGRAM STOPPED"
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC    MNR ROUTINES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine MNR_ttbarj(p,flav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      double precision p(0:3,nlegborn),amp2
      integer flav(nlegborn)
      character * 2 prc
      common/process/prc
      double precision s,q1q,xm2,tk,uk,q2q,w1h,w2h,x,y,cth2
      integer ifl
      double precision dotp,fppx
      external dotp,fppx
c compute mnr invariants
      s=2*dotp(p(0,1),p(0,2))
      q1q=-2*dotp(p(0,1),p(0,3))
      xm2=kn_masses(3)**2
      tk=-2*dotp(p(0,1),p(0,5))
      uk=-2*dotp(p(0,2),p(0,5))
      q2q=-2*dotp(p(0,2),p(0,4))
c     w1=-q1+q2-tk, w2=q1-q2-uk, w1(w2) = w1h(w2h) *s*(1-x)/2
      w1h=2*(q2q-q1q-tk)/(-tk-uk)
      w2h=2*(q1q-q2q-uk)/(-tk-uk)
      x=1+(tk+uk)/s
      y=(tk-uk)/(s*(1-x))
c     cth2 is not used by fgg,fqg,fqq, unless we are very near the collinear limit
      cth2=0
      if(flav(1).eq.0.and.flav(2).eq.0) then
         prc='gg'
         ifl=0
      elseif(flav(1).gt.0.and.flav(2).lt.0) then
         prc='qq'
         ifl=1
      elseif(flav(1).lt.0.and.flav(2).gt.0) then
         prc='qq'
         ifl=-1
      elseif(flav(1).gt.0.and.flav(2).eq.0) then
         prc='qg'
         ifl=1
      elseif(flav(1).lt.0.and.flav(2).eq.0) then
         prc='qg'
         ifl=-1
      elseif(flav(1).eq.0.and.flav(2).gt.0) then
         prc='gq'
         ifl=-1
      elseif(flav(1).eq.0.and.flav(2).lt.0) then
         prc='gq'
         ifl=+1
      endif
      amp2= fppx(ifl,s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)/(4*tk*uk)
      amp2=amp2 * (4*pi*st_alpha)**3 * 2 * s
      end

      function fppx(ifl,s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
c front-end to fpp in hvqcrsx.f that deals with the cases:
c prc=gq,
c ifl=-1, i.e. light quark -> light antiquark
c Thus
c prc ifl
c gg  0     -> gluon gluon
c qq  1     -> quark antiquark
c qq -1     -> antiquark quark
c qg  1     -> quark gluon
c qg -1     -> antiquark gluon
c gq -1     -> gluon quark
c gq  1     -> gluon antiquark
      implicit none
      double precision fppx,s,x,y,xm2,q1q,q2q,w1h,w2h,cth2
      integer ifl
      character * 2 prc
      common/process/prc
      double precision fpp
      external fpp
      if(prc.eq.'gg') then
         if(ifl.eq.0) then
            continue
         else
            goto 998
         endif
      elseif(prc.eq.'qq') then
         if(ifl.eq.1) then
            continue
         elseif(ifl.eq.-1) then
            call exchk1k2(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
         else
            goto 998
         endif
      elseif(prc.eq.'qg') then
         if(ifl.eq.1) then
            continue
         elseif(ifl.eq.-1) then
            call exchk1k2(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
         else
            goto 998
         endif
      elseif(prc.eq.'gq') then
         prc='qg'
         call exchp1p2(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
         if(ifl.eq.-1) then
            continue
         elseif(ifl.eq.1) then
            call exchk1k2(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
         else
            goto 998
         endif
      else
         goto 998
      endif
      fppx=fpp(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
      goto 999
 998  write(*,*) 'fppx called with inconsistent prc,jfl=',prc,ifl
      call exit(-1)
 999  end



      subroutine exchk1k2(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
c exchange of momenta of heavy quarks
      implicit none
      double precision s,x,y,xm2,q1q,q2q,w1h,w2h,cth2
      double precision tk,uk,tmp
      tk=-0.5d0*s*(1-x)*(1-y)
      uk=-0.5d0*s*(1-x)*(1+y)
c exchange
      q1q=-s-tk-q1q
      q2q=-s-uk-q2q
      tmp=w1h
      w1h=w2h
      w2h=tmp
      cth2=-cth2
      end

      subroutine exchp1p2(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
c exchange of momenta of incoming partons
      implicit none
      double precision s,x,y,xm2,q1q,q2q,w1h,w2h,cth2
      double precision tk,uk,tmp
      tk=-0.5d0*s*(1-x)*(1-y)
      uk=-0.5d0*s*(1-x)*(1+y)
      tmp=q1q
      q1q=-s-uk-q2q
      q2q=-s-tk-tmp
      y=-y
c What happens to cth2?  theta2-> 2 pi - theta2, no change
      end




c NR06:  Normalization: (NB: 0 < th2 < pi)
c NR06:   d sigma = fgg/qq/qg * f(x1)*f(x2) * gs^6/(16 pi^2 64 pi^2)
c NR06:   1/(1-x)_rho (1/(1-y)_+ + 1/(1-y)_-) betax/s dcosth1 dth2 dy dx dx1 dx2
c NR06:  NB: fqg -> q(p1)+g(p2)
      function fpp(s0,x,y,xm20,q1q0,q2q0,w1h,w2h,cth2)
      implicit double precision (a-h,o-z)
      character * 2 prc
      common/process/prc
      s=1
      xm2=xm20/s0
      q1q=q1q0/s0
      q2q=q2q0/s0
      if(prc.eq.'gg') then
         fpp = fgg(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
      elseif(prc.eq.'qg') then
         fpp = fqg(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
      elseif(prc.eq.'qq') then
         fpp = fqq(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
      else
         write(*,*) 'FPP: non existent subprocess',prc
         stop
      endif
      end


       FUNCTION FGG(S,X,Y,XM2,q1Q,q2Q,W1H,W2H,CTH2)
c
c      d sigma_gg (f) = N g^6 / (64 pi^2 s) beta_x d costh1 d th2 dy dx
c                    1/(1-x)_rho ( 1/(1-y)_+ + 1/(1+y)_+ ) fgg
c
c           N = 1 / (4 pi)^2
c
      implicit double precision(A-H,O-Z)
      double precision N
       TINY = 0
       V = 8
       N = 3
       N2 = N*N
       CF = 4.d0/3.d0
       CA = 3
       TK = - (1-X)*(1-Y)*S/2
       UK = - (1-X)*(1+Y)*S/2
       IF(1-X.LT.TINY)THEN
          q2C=-S-q2Q
          P13 = -q1Q/2
          P23 = -q2C/2
          P12 = S/2
          BORN = 1/(2*V*N)*(V/(P13*P23)-2*N**2/P12**2)*
     #           (P13**2+P23**2+2*XM2*P12-(XM2*P12)**2/(P13*P23))
          XXX  = N**2/(4*V)*( (P12+2*XM2)/P23-(P12+2*XM2)/P13
     #           +XM2**2/P13**2-XM2**2/P23**2+2*(P23-P13)/P12 )
          YYY  = 1/(4*V*N**2)*
     #           (P13**2+P23**2+2*XM2*P12-(XM2*P12)**2/(P13*P23))
     #          *(1/(P13*P23)+2*N**2/P12**2)
C         --------------------------------------------------------------
C         Fattori iconali moltiplicati per 4*tk*uk
C         P1.K = -TK/2
C         P2.K = -UK/2
C         P3.K = W1/2
C         P4.K = W2/2
C
          P14 = P23
          P24 = P13
C
          q13 = 16*(1+Y)*P13/W1H
          q14 = 16*(1+Y)*P14/W2H
          q23 = 16*(1-Y)*P23/W1H
          q24 = 16*(1-Y)*P24/W2H
          q12 = 16*P12
          q33 = 16*(1-Y)*(1+Y)* XM2/W1H**2
          q44 = 16*(1-Y)*(1+Y)* XM2/W2H**2
          q34 = 16*(1-Y)*(1+Y)* (S/2-XM2)/(W1H*W2H)
          SUM = BORN*( CF*(q13+q14+q23+q24-2*q12-q33-q44) +2*N*q12 )
     #        +        XXX *(q14+q23-q13-q24)
     #        +        YYY *(2*q12+2*q34-q13-q24-q14-q23)
          FGG = 1/(2*S)*SUM
       ELSEIF(1-Y.LT.TINY)THEN
          q2C = -S-UK-q2Q
          P13 = - X*q1Q/2
          P23 = -q2C/2
          P12 = S*X/2
          BX = SQRT(1-4*XM2/(S*X))
          CTH1 = (P23-P13)/P12/BX
          AZIdEP =
     #    -512*BX**2*(CTH1-1)*(CTH1+1)*V*(2*V+BX**2*CTH1**2*N2-N2)
     #    *(X-1)**2*XM2/((BX*CTH1-1)**2*(BX*CTH1+1)**2*X**2)
          AZIdEP = AZIdEP * (2*CTH2**2-1)/2 /(4*64)
          BORN = 1/(2*V*N)*(V/(P13*P23)-2*N**2/P12**2)*
     #           (P13**2+P23**2+2*XM2*P12-(XM2*P12)**2/(P13*P23))
          SUM = - BORN*(8*UK)*2*CA*(X/(1-X)+(1-X)/X+X*(1-X))
     #    + AZIdEP
          FGG = 1/(2*X*S)*SUM
       ELSEIF(1+Y.LT.TINY)THEN
          q2C = -S-UK-q2Q
          P13 = -q1Q/2
          P23 = -X*q2C/2
          P12 = S*X/2
          BX = SQRT(1-4*XM2/(S*X))
          CTH1 = (P23-P13)/P12/BX
          AZIdEP =
     #    -512*BX**2*(CTH1-1)*(CTH1+1)*V*(2*V+BX**2*CTH1**2*N2-N2)
     #    *(X-1)**2*XM2/((BX*CTH1-1)**2*(BX*CTH1+1)**2*X**2)
          AZIdEP = AZIdEP * (2*CTH2**2-1)/2 /(4*64)
          BORN = 1/(2*V*N)*(V/(P13*P23)-2*N**2/P12**2)*
     #           (P13**2+P23**2+2*XM2*P12-(XM2*P12)**2/(P13*P23))
          SUM = - BORN*(8*TK)*2*CA*(X/(1-X)+(1-X)/X+X*(1-X))
     #          + AZIdEP
          FGG = 1/(2*X*S)*SUM
       ELSE
       S2 = S+TK+UK
       q1C=-S-TK-q1Q
       q2C=-S-UK-q2Q
       W1 =q2Q-q1Q-TK
       W2 =q1Q-q2Q-UK
C
       P12 = (S2 - 2*XM2)/2
       P13 = W1/2
       P14 = q1Q/2
       P15 = q2C/2
       P23 = W2/2
       P24 = q1C/2
       P25 = q2Q/2
       P34 = TK/2
       P35 = UK/2
       P45 = S/2
      ANS = FBB1(XM2,P12,P13,P14,P15,P23,P24,P25,P34,P35,P45)
      ANS = 4*TK*UK*ANS/(2*S*4*64)
      FGG = ANS
      ENDIF
      RETURN
      END


       FUNCTION FQQ(S,X,Y,XM2,q1Q,q2Q,W1H,W2H,CTH2)
c
c      q(p1) + q_barra(p2) --> Q(k1) + Q_barra(k2) + g(k)
c
c      d sigma_qq (f) = N g^6 / (64 pi^2 s) beta_x d costh1 d th2 dy dx
c                       1/(1-x)_rho ( 1/(1-y)_+ + 1/(1+y)_+ ) FQQ
c
c           N = 1 / (4 pi)^2
c
      implicit double precision(A-H,O-Z)
      double precision N
       TINY = 0
       V = 8
       N = 3
       N2 = N*N
       CF = 4.d0/3.d0
       CA = 3
       TK = - (1-X)*(1-Y)*S/2
       UK = - (1-X)*(1+Y)*S/2
       IF(1-X.LT.TINY)THEN
          q2C=-S-q2Q
          P13 = -q1Q/2
          P23 = -q2C/2
          P12 = S/2
          BORN = V/(2*N2)*( (P13**2+P23**2)/P12**2 + XM2/P12 )
C
C         Fattori iconali moltiplicati per 4*tk*uk
C         P1.K = -TK/2
C         P2.K = -UK/2
C         P3.K = W1/2
C         P4.K = W2/2
C
          P14 = P23
          P24 = P13
C
          q13 = 16*(1+Y)*P13/W1H
          q14 = 16*(1+Y)*P14/W2H
          q23 = 16*(1-Y)*P23/W1H
          q24 = 16*(1-Y)*P24/W2H
          q12 = 16*P12
          q33 = 16*(1-Y)*(1+Y)* XM2/W1H**2
          q44 = 16*(1-Y)*(1+Y)* XM2/W2H**2
          q34 = 16*(1-Y)*(1+Y)* (S/2-XM2)/(W1H*W2H)
          SUM = BORN*( CF*(2*q13+2*q24-q33-q44)
     #        +           (2*q14+2*q23-q13-q24-q12-q34)/N )
          FQQ = 1/(2*S)*SUM
       ELSEIF(1-Y.LT.TINY)THEN
          q2C = -S-UK-q2Q
          P13 = - X*q1Q/2
          P23 = -q2C/2
          P12 = S*X/2
          BORN = V/(2*N2)*( (P13**2+P23**2)/P12**2 + XM2/P12 )
          SUM = - BORN*(8*UK)*CF*(1+X**2)/(1-X)
          FQQ = 1/(2*X*S)*SUM
       ELSEIF(1+Y.LT.TINY)THEN
          q2C = -S-UK-q2Q
          P13 = -q1Q/2
          P23 = -X*q2C/2
          P12 = S*X/2
          BORN = V/(2*N2)*( (P13**2+P23**2)/P12**2 + XM2/P12 )
          SUM = - BORN*(8*TK)*CF*(1+X**2)/(1-X)
          FQQ = 1/(2*X*S)*SUM
       ELSE
       S2 = S+TK+UK
       q1C=-S-TK-q1Q
       q2C=-S-UK-q2Q
       W1 =q2Q-q1Q-TK
       W2 =q1Q-q2Q-UK
C
       P12 = S/2
       P13 = q1Q/2
       P14 = q1C/2
       P15 = TK/2
       P23 = q2C/2
       P24 = q2Q/2
       P25 = UK/2
       P34 = (S2-2*XM2)/2
       P35 = W1/2
       P45 = W2/2
      ans = q4G1M(XM2,P34,P24,P14,P45,P23,P13,P35,P12,P25,P15)
      ANS = 4*TK*UK*ANS/(2*S*4*9)
      FQQ = ANS
      ENDIF
      RETURN
      END

       FUNCTION FQG(S,X,Y,XM2,q1Q,q2Q,W1H,W2H,CTH2)
c
c      q(p1) + g(p2) --> Q(k1) + Q_barra(k2) + g(k)
c
c      d sigma_qg (f) = N g^6 / (64 pi^2 s) beta_x d costh1 d th2 dy dx
c                       1/(1-x)_rho ( 1/(1-y)_+ + 1/(1+y)_+ ) FQG
c
c           N = 1 / (4 pi)^2
c
      implicit double precision(A-H,O-Z)
      double precision N
       TINY = 0
       V = 8
       N = 3
       N2 = N*N
       TF = .5d0
       CF = 4.d0/3.d0
       CA = 3
       TK = - (1-X)*(1-Y)*S/2
       UK = - (1-X)*(1+Y)*S/2
       IF(1-X.LT.TINY)THEN
          FQG = 0
       ELSEIF(1-Y.LT.TINY)THEN
          q2C = -S-UK-q2Q
          P13 = - X*q1Q/2
          P23 = -q2C/2
          P12 = S*X/2
          BX = SQRT(1-4*XM2/(S*X))
          CTH1 = (P23-P13)/P12/BX
          AZIdEP =
     #    -512*BX**2*(CTH1-1)*(CTH1+1)*V*(2*V+BX**2*CTH1**2*N2-N2)
     #    *(X-1)**2*XM2/((BX*CTH1-1)**2*(BX*CTH1+1)**2*X**2)
          AZIdEP = AZIdEP * (2*CTH2**2-1)/2 /(4*64)
          AZIdEP = AZIdEP * CF / CA
          BORN = 1/(2*V*N)*(V/(P13*P23)-2*N**2/P12**2)*
     #           (P13**2+P23**2+2*XM2*P12-(XM2*P12)**2/(P13*P23))
          SUM = - BORN*(8*UK)*CF*(1+(1-X)**2)/X
     #    + AZIdEP
          FQG = 1/(2*X*S)*SUM
       ELSEIF(1+Y.LT.TINY)THEN
          q2C = -S-UK-q2Q
          P13 = -q1Q/2
          P23 = -X*q2C/2
          P12 = S*X/2
          BORN = V/(2*N2)*( (P13**2+P23**2)/P12**2 + XM2/P12 )
          SUM = - BORN*(8*TK)*TF*(X**2+(1-X)**2)
          FQG = 1/(2*X*S)*SUM
       ELSE
       S2 = S+TK+UK
       q1C=-S-TK-q1Q
       q2C=-S-UK-q2Q
       W1 =q2Q-q1Q-TK
       W2 =q1Q-q2Q-UK
C
       P12 = S/2
       P13 = q1Q/2
       P14 = q1C/2
       P15 = TK/2
       P23 = q2C/2
       P24 = q2Q/2
       P25 = UK/2
       P34 = (S2-2*XM2)/2
       P35 = W1/2
       P45 = W2/2
      ANS = q4G1M(XM2,P34,P45,P14,P24,P35,P13,P23,P15,P25,P12)
      ANS = - 4*TK*UK*ANS/(2*S*4*3*8)
      FQG = ANS
      ENDIF
      RETURN
      END



c
c       This function subroutine, q4G1M, calculates the
c       invariant matrix element squared for the process:
c
c       Q(-p1) + Qbar(-p2) --> q(p3) + qbar(p4) + g(p5)
C
c       summed over initial and final spins and colors including
c       masses for the incoming quark-antiquark pair.
c       The final state quark-antiquark pair are massless.
c       No averaging is performed for initial spins or colors.
c
c
c       1.      p1...p5 are the outgoing four momenta of the
c               partons, satisfying
c                       p1 + p2 + p3 + p4 + p5 = 0
c
c       2.      xm2 is the mass^2 of the heavy Quark and AntiQuark lines
c               with
c                       p1**2 = p2**2 = xm2
c               where
c                       pi**2 = pi(4)**2-pi(1)**2-pi(2)**2-pi(3)**2
c
c Abbiamo:
c             p12 = p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
c             etc...
c
c Per il processo
c       Q(-p1) + Qbar(-p2) --> q(p3) + qbar(p4) + g(p5)
c chiamare
c       q4G1M(XM2,P12,P13,P14,P15,P23,P24,P25,P34,P35,P45)*g**6/4/N**2
c
c Per
c       q(-p1) + qbar(-p2) --> Q(p3) + Qbar(p4) + g(p5)
c che equivarrebbe a
c       Q(-p4)+Qbar(-p3) --> q(p2)+qbar(p1)+g(p5)
c chiamare
c       q4G1M(XM2,P34,P24,P14,P45,P23,P13,P35,P12,P25,P15)*g**6/4/N**2
c
c

      FUNCTION  q4G1M(XM2,P12,P13,P14,P15,P23,P24,P25,P34,P35,P45)
      double precision  dL1, dL2, dL3, dL4
      double precision  P12, P13, P14, P15, P23, P24, P25, P34, P35, P45
      double precision  q4G1M, XM2, RES, S, XN, XV
      PARAMETER         (XN = 3.d0)

      XV  = XN**2 - 1
      S = 2*(P12+XM2)

      dL1 = P13/P25-2*P35/S
      dL2 = P14/P25-2*P45/S
      dL3 = P23/P15-2*P35/S
      dL4 = P24/P15-2*P45/S


      RES =
     & +(P13**2+P23**2+P14**2+P24**2+XM2*(P12+P34+XM2))/2/S/P34
     & *(4*XV**2/XN*(P13/P15/P35+P24/P25/P45)
     & +4*XV/XN*(2*P14/P15/P45+2*P23/P25/P35
     & -P13/P15/P35-P24/P25/P45-P12/P15/P25-P34/P35/P45) )

      RES = RES
     & -XV*(XN**2-4)/XN
     & *2*XM2/S/P34*((P13-P14)/P25-(P23-P24)/P15)
     & +4*XV**2/XN*XM2*( (P35**2+P45**2)/P35/P45/S**2
     & -0.5*(1/P15+1/P25+1/P35+1/P45)/S
     & -0.25*(1/P15+1/P25+XM2/P15**2+XM2/P25**2+4/S)/P34
     & -(dL1**2+dL2**2+dL3**2+dL4**2)/4/P34**2)
     & -2*XV/XN
     & *XM2/S/P34*(1+2*P34/S+XM2/P15+XM2/P25+(P35**2+P45**2)/P15/P25
     & +(P13-P14)*(dL1-dL2)/P34+(P23-P24)*(dL3-dL4)/P34 )

      q4G1M = RES


      RETURN
      END



      FUNCTION FBB1(XM2,P12,P13,P14,P15,P23,P24,P25,P34,P35,P45)
      implicit double precision (A-Z)
      n = 3
      n2 = n**2
      v = n2-1
      s = 2*xm2+2*p12
      ans = 0
      tmp0 = v**3*xm2*((xm2**2+4*p24*xm2-2*p24*p25)/(p13**2*p24**2)+(p34
     1   *s+2*p15*p25-p14*p24)/(p13*p14*p23*p24))/n2/2+v*xm2*(xm2**2/(
     2   p13*p14*p23*p25)-(5*p13*s+6*p13*p25+8*p13*p23-4*p13**2)/(p13*p1
     3   4*p23*p24)/4)/n2
      tmp0 = v**2*xm2*(-4*xm2**2/(p13*p14*p25**2)+2*s*xm2/(p13*p14*p23*p
     1   25)+(-(p23+p13)*p45+p23**2+p13**2)/(p13*p14*p23*p25))/n2/4+((
     2   -p13*s+p15*p25+p23**2+p15**2+p13**2)/(p13*p15*p24*p25*p34)+(p24
     3   +p15)/(p13*p14*p23*p25))*v*xm2**2-v*((2*(p45+p35)*xm2+p25**2+p1
     4   5**2)/(p14*p23*p34)/2+p13*(2*(p35+p34)*xm2+p23**2+p13**2)/(p1
     5   4*p15*p25*p34))+(n2+1)*p12*v*(2*(p45+p35)*xm2+p25**2+p15**2)/(n
     6   2*p13*p14*p23*p24)/4+v**2*xm2*((p45+p34)*xm2/(p13*p25*p34*p45
     7   )+2*xm2/(p15**2*p34)+(p13/p25-p23/p15)**2/p34**2)+tmp0
      ans = -2*n2*(s**2/4+p45**2+p34**2)*v*xm2**2/(p13*p25*p34*p45*s)+
     1   n2*p14*p24*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(p13*p25*p34*p45*s
     2   )+2*n2*p23*v*(2*(p35+p34)*xm2+p23**2+p13**2)/(p25*p34*p45*s)-(n
     3   2+1)*v*xm2*(4*p12*xm2+s**2+2*p23*p24+2*p13*p14)/(n2*p13*p14*p23
     4   *p24)/8+n2*v*((p34-2*xm2)/(p13*p24*s)-s/(p15*p25*p34)/4)*xm
     5   2+v*((p12-3*xm2)/(p15*p25*p34)+(5*p13*p25+3*p13*p23-2*p13**2)/(
     6   p13*p14*p23*p24))*xm2+2*n2*(-s*(s+2*p34)/8+4*p15*p25+p14*p24+
     7   p13*p23)*v*xm2/(p15*p25*p34*s)+4*n2*(2*p23**2/(p15*p34**2*s)+(p
     8   45**2+p35**2+p34**2)/(p34**2*s**2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p23*xm2-2*p23*p25)/(p14**2*p23**2)+(p34
     1   *s+2*p15*p25-p13*p23)/(p13*p14*p23*p24))/n2/2
      tmp0 = v*xm2*(xm2**2/(p13*p14*p24*p25)-(5*p14*s+6*p14*p25+8*p14*p2
     1   4-4*p14**2)/(p13*p14*p23*p24)/4)/n2+v**2*xm2*(-4*xm2**2/(p13*
     2   p14*p25**2)+2*s*xm2/(p13*p14*p24*p25)+(-(p24+p14)*p35+p24**2+p1
     3   4**2)/(p13*p14*p24*p25))/n2/4+((-p14*s+p15*p25+p24**2+p15**2+
     4   p14**2)/(p14*p15*p23*p25*p34)+(p23+p15)/(p13*p14*p24*p25))*v*xm
     5   2**2-v*((2*(p45+p35)*xm2+p25**2+p15**2)/(p13*p24*p34)/2+p14*(
     6   2*(p45+p34)*xm2+p24**2+p14**2)/(p13*p15*p25*p34))+(n2+1)*p12*v*
     7   (2*(p45+p35)*xm2+p25**2+p15**2)/(n2*p13*p14*p23*p24)/4+2*n2*p
     8   24*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(p25*p34*p35*s)+tmp0
      ans = -2*n2*(s**2/4+p35**2+p34**2)*v*xm2**2/(p14*p25*p34*p35*s)+
     1   v**2*xm2*((p35+p34)*xm2/(p14*p25*p34*p35)+2*xm2/(p15**2*p34)+(p
     2   14/p25-p24/p15)**2/p34**2)+n2*p13*p23*v*(2*(p35+p34)*xm2+p23**2
     3   +p13**2)/(p14*p25*p34*p35*s)-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p23
     4   *p24+2*p13*p14)/(n2*p13*p14*p23*p24)/8+n2*v*((p34-2*xm2)/(p14
     5   *p23*s)-s/(p15*p25*p34)/4)*xm2+v*((p12-3*xm2)/(p15*p25*p34)+(
     6   5*p14*p25+3*p14*p24-2*p14**2)/(p13*p14*p23*p24))*xm2+2*n2*(-s*(
     7   s+2*p34)/8+4*p15*p25+p14*p24+p13*p23)*v*xm2/(p15*p25*p34*s)+4
     8   *n2*(2*p24**2/(p15*p34**2*s)+(p45**2+p35**2+p34**2)/(p34**2*s**
     9   2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p25*xm2-2*p23*p25)/(p14**2*p25**2)+(p45
     1   *s-p15*p25+2*p13*p23)/(p14*p15*p24*p25))/n2/2
      tmp0 = v*xm2*(xm2**2/(p14*p15*p23*p24)-(5*p14*s+8*p14*p24+6*p14*p2
     1   3-4*p14**2)/(p14*p15*p24*p25)/4)/n2+v**2*xm2*(-4*xm2**2/(p14*
     2   p15*p23**2)+2*s*xm2/(p14*p15*p23*p24)+(-(p24+p14)*p35+p24**2+p1
     3   4**2)/(p14*p15*p23*p24))/n2/4-2*n2*(s**2/4+p45**2+p35**2)*v
     4   *xm2**2/(p14*p23*p35*p45*s)+((-p14*s+p24**2+p13*p23+p14**2+p13*
     5   *2)/(p13*p14*p23*p25*p45)+(p25+p13)/(p14*p15*p23*p24))*v*xm2**2
     6   -v*(p14*(2*(p45+p34)*xm2+p24**2+p14**2)/(p13*p15*p23*p45)+(2*(p
     7   35+p34)*xm2+p23**2+p13**2)/(p15*p24*p45)/2)+v**2*xm2*((p45+p3
     8   5)*xm2/(p14*p23*p35*p45)+2*xm2/(p13**2*p45)+(p14/p23-p24/p13)**
     9   2/p45**2)+n2*p15*p25*v*(2*(p45+p35)*xm2+p25**2+p15**2)/(p14*p23
     :   *p35*p45*s)+tmp0
      ans = 2*n2*p24*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(p23*p35*p45*s)+(
     1   n2+1)*p12*v*(2*(p35+p34)*xm2+p23**2+p13**2)/(n2*p14*p15*p24*p25
     2   )/4-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p24*p25+2*p14*p15)/(n2*p14
     3   *p15*p24*p25)/8+n2*v*((p45-2*xm2)/(p14*p25*s)-s/(p13*p23*p45)
     4   /4)*xm2+v*((p12-3*xm2)/(p13*p23*p45)+(3*p14*p24+5*p14*p23-2*p
     5   14**2)/(p14*p15*p24*p25))*xm2+2*n2*(-s*(s+2*p45)/8+p15*p25+p1
     6   4*p24+4*p13*p23)*v*xm2/(p13*p23*p45*s)+4*n2*(2*p24**2/(p13*p45*
     7   *2*s)+(p45**2+p35**2+p34**2)/(p45**2*s**2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p24*xm2-2*p23*p24)/(p15**2*p24**2)+(p45
     1   *s-p14*p24+2*p13*p23)/(p14*p15*p24*p25))/n2/2
      tmp0 = v*xm2*(xm2**2/(p14*p15*p23*p25)-(5*p15*s+8*p15*p25+6*p15*p2
     1   3-4*p15**2)/(p14*p15*p24*p25)/4)/n2+v**2*xm2*(-4*xm2**2/(p14*
     2   p15*p23**2)+2*s*xm2/(p14*p15*p23*p25)+(-(p25+p15)*p34+p25**2+p1
     3   5**2)/(p14*p15*p23*p25))/n2/4-2*n2*(s**2/4+p45**2+p34**2)*v
     4   *xm2**2/(p15*p23*p34*p45*s)+((-p15*s+p25**2+p13*p23+p15**2+p13*
     5   *2)/(p13*p15*p23*p24*p45)+(p24+p13)/(p14*p15*p23*p25))*v*xm2**2
     6   -v*(p15*(2*(p45+p35)*xm2+p25**2+p15**2)/(p13*p14*p23*p45)+(2*(p
     7   35+p34)*xm2+p23**2+p13**2)/(p14*p25*p45)/2)+2*n2*p25*v*(2*(p4
     8   5+p35)*xm2+p25**2+p15**2)/(p23*p34*p45*s)+v**2*xm2*((p45+p34)*x
     9   m2/(p15*p23*p34*p45)+2*xm2/(p13**2*p45)+(p15/p23-p25/p13)**2/p4
     :   5**2)+tmp0
      ans = n2*p14*p24*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(p15*p23*p34*p4
     1   5*s)+(n2+1)*p12*v*(2*(p35+p34)*xm2+p23**2+p13**2)/(n2*p14*p15*p
     2   24*p25)/4-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p24*p25+2*p14*p15)/(
     3   n2*p14*p15*p24*p25)/8+n2*v*((p45-2*xm2)/(p15*p24*s)-s/(p13*p2
     4   3*p45)/4)*xm2+v*((p12-3*xm2)/(p13*p23*p45)+(3*p15*p25+5*p15*p
     5   23-2*p15**2)/(p14*p15*p24*p25))*xm2+2*n2*(-s*(s+2*p45)/8+p15*
     6   p25+p14*p24+4*p13*p23)*v*xm2/(p13*p23*p45*s)+4*n2*(2*p25**2/(p1
     7   3*p45**2*s)+(p45**2+p35**2+p34**2)/(p45**2*s**2))*v*xm2+tmp0+an
     8   s
      tmp0 = v**3*xm2*((xm2**2+4*p23*xm2-2*p23*p24)/(p15**2*p23**2)+(p35
     1   *s+2*p14*p24-p13*p23)/(p13*p15*p23*p25))/n2/2
      tmp0 = v*xm2*(xm2**2/(p13*p15*p24*p25)-(5*p15*s+8*p15*p25+6*p15*p2
     1   4-4*p15**2)/(p13*p15*p23*p25)/4)/n2+v**2*xm2*(-4*xm2**2/(p13*
     2   p15*p24**2)+2*s*xm2/(p13*p15*p24*p25)+(-(p25+p15)*p34+p25**2+p1
     3   5**2)/(p13*p15*p24*p25))/n2/4+((-p15*s+p25**2+p14*p24+p15**2+
     4   p14**2)/(p14*p15*p23*p24*p35)+(p23+p14)/(p13*p15*p24*p25))*v*xm
     5   2**2-v*(p15*(2*(p45+p35)*xm2+p25**2+p15**2)/(p13*p14*p24*p35)+(
     6   2*(p45+p34)*xm2+p24**2+p14**2)/(p13*p25*p35)/2)+2*n2*p25*v*(2
     7   *(p45+p35)*xm2+p25**2+p15**2)/(p24*p34*p35*s)+(n2+1)*p12*v*(2*(
     8   p45+p34)*xm2+p24**2+p14**2)/(n2*p13*p15*p23*p25)/4+tmp0
      ans = -2*n2*(s**2/4+p35**2+p34**2)*v*xm2**2/(p15*p24*p34*p35*s)+
     1   v**2*xm2*((p35+p34)*xm2/(p15*p24*p34*p35)+2*xm2/(p14**2*p35)+(p
     2   15/p24-p25/p14)**2/p35**2)+n2*p13*p23*v*(2*(p35+p34)*xm2+p23**2
     3   +p13**2)/(p15*p24*p34*p35*s)-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p23
     4   *p25+2*p13*p15)/(n2*p13*p15*p23*p25)/8+n2*v*((p35-2*xm2)/(p15
     5   *p23*s)-s/(p14*p24*p35)/4)*xm2+v*((p12-3*xm2)/(p14*p24*p35)+(
     6   3*p15*p25+5*p15*p24-2*p15**2)/(p13*p15*p23*p25))*xm2+2*n2*(-s*(
     7   s+2*p35)/8+p15*p25+4*p14*p24+p13*p23)*v*xm2/(p14*p24*p35*s)+4
     8   *n2*(2*p25**2/(p14*p35**2*s)+(p45**2+p35**2+p34**2)/(p35**2*s**
     9   2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p25*xm2-2*p24*p25)/(p13**2*p25**2)+(p35
     1   *s-p15*p25+2*p14*p24)/(p13*p15*p23*p25))/n2/2
      tmp0 = v*xm2*(xm2**2/(p13*p15*p23*p24)-(5*p13*s+6*p13*p24+8*p13*p2
     1   3-4*p13**2)/(p13*p15*p23*p25)/4)/n2+v**2*xm2*(-4*xm2**2/(p13*
     2   p15*p24**2)+2*s*xm2/(p13*p15*p23*p24)+(-(p23+p13)*p45+p23**2+p1
     3   3**2)/(p13*p15*p23*p24))/n2/4-2*n2*(s**2/4+p45**2+p35**2)*v
     4   *xm2**2/(p13*p24*p35*p45*s)+((-p13*s+p14*p24+p23**2+p14**2+p13*
     5   *2)/(p13*p14*p24*p25*p35)+(p25+p14)/(p13*p15*p23*p24))*v*xm2**2
     6   -v*((2*(p45+p34)*xm2+p24**2+p14**2)/(p15*p23*p35)/2+p13*(2*(p
     7   35+p34)*xm2+p23**2+p13**2)/(p14*p15*p24*p35))+v**2*xm2*((p45+p3
     8   5)*xm2/(p13*p24*p35*p45)+2*xm2/(p14**2*p35)+(p13/p24-p23/p14)**
     9   2/p35**2)+n2*p15*p25*v*(2*(p45+p35)*xm2+p25**2+p15**2)/(p13*p24
     :   *p35*p45*s)+tmp0
      ans = (n2+1)*p12*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(n2*p13*p15*p23
     1   *p25)/4+2*n2*p23*v*(2*(p35+p34)*xm2+p23**2+p13**2)/(p24*p35*p
     2   45*s)-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p23*p25+2*p13*p15)/(n2*p13
     3   *p15*p23*p25)/8+n2*v*((p35-2*xm2)/(p13*p25*s)-s/(p14*p24*p35)
     4   /4)*xm2+v*((p12-3*xm2)/(p14*p24*p35)+(5*p13*p24+3*p13*p23-2*p
     5   13**2)/(p13*p15*p23*p25))*xm2+2*n2*(-s*(s+2*p35)/8+p15*p25+4*
     6   p14*p24+p13*p23)*v*xm2/(p14*p24*p35*s)+4*n2*(2*p23**2/(p14*p35*
     7   *2*s)+(p45**2+p35**2+p34**2)/(p35**2*s**2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p14*xm2-2*p14*p15)/(p14**2*p23**2)+(p34
     1   *s+2*p15*p25-p14*p24)/(p13*p14*p23*p24))/n2/2+v*xm2*(xm2**2/(
     2   p13*p15*p23*p24)-(5*p23*s-4*p23**2+6*p15*p23+8*p13*p23)/(p13*p1
     3   4*p23*p24)/4)/n2
      tmp0 = v**2*xm2*(-4*xm2**2/(p15**2*p23*p24)+2*s*xm2/(p13*p15*p23*p
     1   24)+(-(p23+p13)*p45+p23**2+p13**2)/(p13*p15*p23*p24))/n2/4+((
     2   -p23*s+p25**2+p15*p25+p23**2+p13**2)/(p14*p15*p23*p25*p34)+(p25
     3   +p14)/(p13*p15*p23*p24))*v*xm2**2-v*((2*(p45+p35)*xm2+p25**2+p1
     4   5**2)/(p13*p24*p34)/2+p23*(2*(p35+p34)*xm2+p23**2+p13**2)/(p1
     5   5*p24*p25*p34))+(n2+1)*p12*v*(2*(p45+p35)*xm2+p25**2+p15**2)/(n
     6   2*p13*p14*p23*p24)/4+v**2*xm2*((p45+p34)*xm2/(p15*p23*p34*p45
     7   )+2*xm2/(p25**2*p34)+(p23/p15-p13/p25)**2/p34**2)+tmp0
      ans = -2*n2*(s**2/4+p45**2+p34**2)*v*xm2**2/(p15*p23*p34*p45*s)+
     1   n2*p14*p24*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(p15*p23*p34*p45*s
     2   )+2*n2*p13*v*(2*(p35+p34)*xm2+p23**2+p13**2)/(p15*p34*p45*s)-(n
     3   2+1)*v*xm2*(4*p12*xm2+s**2+2*p23*p24+2*p13*p14)/(n2*p13*p14*p23
     4   *p24)/8+n2*v*((p34-2*xm2)/(p14*p23*s)-s/(p15*p25*p34)/4)*xm
     5   2+v*((p12-3*xm2)/(p15*p25*p34)+(-2*p23**2+5*p15*p23+3*p13*p23)/
     6   (p13*p14*p23*p24))*xm2+2*n2*(-s*(s+2*p34)/8+4*p15*p25+p14*p24
     7   +p13*p23)*v*xm2/(p15*p25*p34*s)+4*n2*(2*p13**2/(p25*p34**2*s)+(
     8   p45**2+p35**2+p34**2)/(p34**2*s**2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p13*xm2-2*p13*p15)/(p13**2*p24**2)+(p34
     1   *s+2*p15*p25-p13*p23)/(p13*p14*p23*p24))/n2/2
      tmp0 = v*xm2*(xm2**2/(p14*p15*p23*p24)-(5*p24*s-4*p24**2+6*p15*p24
     1   +8*p14*p24)/(p13*p14*p23*p24)/4)/n2+v**2*xm2*(-4*xm2**2/(p15*
     2   *2*p23*p24)+2*s*xm2/(p14*p15*p23*p24)+(-(p24+p14)*p35+p24**2+p1
     3   4**2)/(p14*p15*p23*p24))/n2/4+((-p24*s+p25**2+p15*p25+p24**2+
     4   p14**2)/(p13*p15*p24*p25*p34)+(p25+p13)/(p14*p15*p23*p24))*v*xm
     5   2**2-v*((2*(p45+p35)*xm2+p25**2+p15**2)/(p14*p23*p34)/2+p24*(
     6   2*(p45+p34)*xm2+p24**2+p14**2)/(p15*p23*p25*p34))+(n2+1)*p12*v*
     7   (2*(p45+p35)*xm2+p25**2+p15**2)/(n2*p13*p14*p23*p24)/4+2*n2*p
     8   14*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(p15*p34*p35*s)+tmp0
      ans = -2*n2*(s**2/4+p35**2+p34**2)*v*xm2**2/(p15*p24*p34*p35*s)+
     1   v**2*xm2*((p35+p34)*xm2/(p15*p24*p34*p35)+2*xm2/(p25**2*p34)+(p
     2   24/p15-p14/p25)**2/p34**2)+n2*p13*p23*v*(2*(p35+p34)*xm2+p23**2
     3   +p13**2)/(p15*p24*p34*p35*s)-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p23
     4   *p24+2*p13*p14)/(n2*p13*p14*p23*p24)/8+n2*v*((p34-2*xm2)/(p13
     5   *p24*s)-s/(p15*p25*p34)/4)*xm2+v*((p12-3*xm2)/(p15*p25*p34)+(
     6   -2*p24**2+5*p15*p24+3*p14*p24)/(p13*p14*p23*p24))*xm2+2*n2*(-s*
     7   (s+2*p34)/8+4*p15*p25+p14*p24+p13*p23)*v*xm2/(p15*p25*p34*s)+
     8   4*n2*(2*p14**2/(p25*p34**2*s)+(p45**2+p35**2+p34**2)/(p34**2*s*
     9   *2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p15*xm2-2*p13*p15)/(p15**2*p24**2)+(p45
     1   *s-p15*p25+2*p13*p23)/(p14*p15*p24*p25))/n2/2
      tmp0 = v*xm2*(xm2**2/(p13*p14*p24*p25)-(5*p24*s-4*p24**2+8*p14*p24
     1   +6*p13*p24)/(p14*p15*p24*p25)/4)/n2+v**2*xm2*(-4*xm2**2/(p13*
     2   *2*p24*p25)+2*s*xm2/(p13*p14*p24*p25)+(-(p24+p14)*p35+p24**2+p1
     3   4**2)/(p13*p14*p24*p25))/n2/4-2*n2*(s**2/4+p45**2+p35**2)*v
     4   *xm2**2/(p13*p24*p35*p45*s)+((-p24*s+p24**2+p23**2+p13*p23+p14*
     5   *2)/(p13*p15*p23*p24*p45)+(p23+p15)/(p13*p14*p24*p25))*v*xm2**2
     6   -v*(p24*(2*(p45+p34)*xm2+p24**2+p14**2)/(p13*p23*p25*p45)+(2*(p
     7   35+p34)*xm2+p23**2+p13**2)/(p14*p25*p45)/2)+v**2*xm2*((p45+p3
     8   5)*xm2/(p13*p24*p35*p45)+2*xm2/(p23**2*p45)+(p24/p13-p14/p23)**
     9   2/p45**2)+n2*p15*p25*v*(2*(p45+p35)*xm2+p25**2+p15**2)/(p13*p24
     :   *p35*p45*s)+tmp0
      ans = 2*n2*p14*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(p13*p35*p45*s)+(
     1   n2+1)*p12*v*(2*(p35+p34)*xm2+p23**2+p13**2)/(n2*p14*p15*p24*p25
     2   )/4-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p24*p25+2*p14*p15)/(n2*p14
     3   *p15*p24*p25)/8+n2*v*((p45-2*xm2)/(p15*p24*s)-s/(p13*p23*p45)
     4   /4)*xm2+v*((p12-3*xm2)/(p13*p23*p45)+(-2*p24**2+3*p14*p24+5*p
     5   13*p24)/(p14*p15*p24*p25))*xm2+2*n2*(-s*(s+2*p45)/8+p15*p25+p
     6   14*p24+4*p13*p23)*v*xm2/(p13*p23*p45*s)+4*n2*(2*p14**2/(p23*p45
     7   **2*s)+(p45**2+p35**2+p34**2)/(p45**2*s**2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p14*xm2-2*p13*p14)/(p14**2*p25**2)+(p45
     1   *s-p14*p24+2*p13*p23)/(p14*p15*p24*p25))/n2/2
      tmp0 = v*xm2*(xm2**2/(p13*p15*p24*p25)-(5*p25*s-4*p25**2+8*p15*p25
     1   +6*p13*p25)/(p14*p15*p24*p25)/4)/n2+v**2*xm2*(-4*xm2**2/(p13*
     2   *2*p24*p25)+2*s*xm2/(p13*p15*p24*p25)+(-(p25+p15)*p34+p25**2+p1
     3   5**2)/(p13*p15*p24*p25))/n2/4-2*n2*(s**2/4+p45**2+p34**2)*v
     4   *xm2**2/(p13*p25*p34*p45*s)+((-p25*s+p25**2+p23**2+p13*p23+p15*
     5   *2)/(p13*p14*p23*p25*p45)+(p23+p14)/(p13*p15*p24*p25))*v*xm2**2
     6   -v*(p25*(2*(p45+p35)*xm2+p25**2+p15**2)/(p13*p23*p24*p45)+(2*(p
     7   35+p34)*xm2+p23**2+p13**2)/(p15*p24*p45)/2)+2*n2*p15*v*(2*(p4
     8   5+p35)*xm2+p25**2+p15**2)/(p13*p34*p45*s)+v**2*xm2*((p45+p34)*x
     9   m2/(p13*p25*p34*p45)+2*xm2/(p23**2*p45)+(p25/p13-p15/p23)**2/p4
     :   5**2)+tmp0
      ans = n2*p14*p24*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(p13*p25*p34*p4
     1   5*s)+(n2+1)*p12*v*(2*(p35+p34)*xm2+p23**2+p13**2)/(n2*p14*p15*p
     2   24*p25)/4-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p24*p25+2*p14*p15)/(
     3   n2*p14*p15*p24*p25)/8+n2*v*((p45-2*xm2)/(p14*p25*s)-s/(p13*p2
     4   3*p45)/4)*xm2+v*((p12-3*xm2)/(p13*p23*p45)+(-2*p25**2+3*p15*p
     5   25+5*p13*p25)/(p14*p15*p24*p25))*xm2+2*n2*(-s*(s+2*p45)/8+p15
     6   *p25+p14*p24+4*p13*p23)*v*xm2/(p13*p23*p45*s)+4*n2*(2*p15**2/(p
     7   23*p45**2*s)+(p45**2+p35**2+p34**2)/(p45**2*s**2))*v*xm2+tmp0+a
     8   ns
      tmp0 = v**3*xm2*((xm2**2+4*p13*xm2-2*p13*p14)/(p13**2*p25**2)+(p35
     1   *s+2*p14*p24-p13*p23)/(p13*p15*p23*p25))/n2/2
      tmp0 = v*xm2*(xm2**2/(p14*p15*p23*p25)-(5*p25*s-4*p25**2+8*p15*p25
     1   +6*p14*p25)/(p13*p15*p23*p25)/4)/n2+v**2*xm2*(-4*xm2**2/(p14*
     2   *2*p23*p25)+2*s*xm2/(p14*p15*p23*p25)+(-(p25+p15)*p34+p25**2+p1
     3   5**2)/(p14*p15*p23*p25))/n2/4+((-p25*s+p25**2+p24**2+p14*p24+
     4   p15**2)/(p13*p14*p24*p25*p35)+(p24+p13)/(p14*p15*p23*p25))*v*xm
     5   2**2-v*(p25*(2*(p45+p35)*xm2+p25**2+p15**2)/(p14*p23*p24*p35)+(
     6   2*(p45+p34)*xm2+p24**2+p14**2)/(p15*p23*p35)/2)+2*n2*p15*v*(2
     7   *(p45+p35)*xm2+p25**2+p15**2)/(p14*p34*p35*s)+(n2+1)*p12*v*(2*(
     8   p45+p34)*xm2+p24**2+p14**2)/(n2*p13*p15*p23*p25)/4+tmp0
      ans = -2*n2*(s**2/4+p35**2+p34**2)*v*xm2**2/(p14*p25*p34*p35*s)+
     1   v**2*xm2*((p35+p34)*xm2/(p14*p25*p34*p35)+2*xm2/(p24**2*p35)+(p
     2   25/p14-p15/p24)**2/p35**2)+n2*p13*p23*v*(2*(p35+p34)*xm2+p23**2
     3   +p13**2)/(p14*p25*p34*p35*s)-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p23
     4   *p25+2*p13*p15)/(n2*p13*p15*p23*p25)/8+n2*v*((p35-2*xm2)/(p13
     5   *p25*s)-s/(p14*p24*p35)/4)*xm2+v*((p12-3*xm2)/(p14*p24*p35)+(
     6   -2*p25**2+3*p15*p25+5*p14*p25)/(p13*p15*p23*p25))*xm2+2*n2*(-s*
     7   (s+2*p35)/8+p15*p25+4*p14*p24+p13*p23)*v*xm2/(p14*p24*p35*s)+
     8   4*n2*(2*p15**2/(p24*p35**2*s)+(p45**2+p35**2+p34**2)/(p35**2*s*
     9   *2))*v*xm2+tmp0+ans
      tmp0 = v**3*xm2*((xm2**2+4*p15*xm2-2*p14*p15)/(p15**2*p23**2)+(p35
     1   *s-p15*p25+2*p14*p24)/(p13*p15*p23*p25))/n2/2
      tmp0 = v*xm2*(xm2**2/(p13*p14*p23*p25)-(5*p23*s-4*p23**2+6*p14*p23
     1   +8*p13*p23)/(p13*p15*p23*p25)/4)/n2+v**2*xm2*(-4*xm2**2/(p14*
     2   *2*p23*p25)+2*s*xm2/(p13*p14*p23*p25)+(-(p23+p13)*p45+p23**2+p1
     3   3**2)/(p13*p14*p23*p25))/n2/4-2*n2*(s**2/4+p45**2+p35**2)*v
     4   *xm2**2/(p14*p23*p35*p45*s)+((-p23*s+p24**2+p14*p24+p23**2+p13*
     5   *2)/(p14*p15*p23*p24*p35)+(p24+p15)/(p13*p14*p23*p25))*v*xm2**2
     6   -v*((2*(p45+p34)*xm2+p24**2+p14**2)/(p13*p25*p35)/2+p23*(2*(p
     7   35+p34)*xm2+p23**2+p13**2)/(p14*p24*p25*p35))+v**2*xm2*((p45+p3
     8   5)*xm2/(p14*p23*p35*p45)+2*xm2/(p24**2*p35)+(p23/p14-p13/p24)**
     9   2/p35**2)+n2*p15*p25*v*(2*(p45+p35)*xm2+p25**2+p15**2)/(p14*p23
     :   *p35*p45*s)+tmp0
      ans = (n2+1)*p12*v*(2*(p45+p34)*xm2+p24**2+p14**2)/(n2*p13*p15*p23
     1   *p25)/4+2*n2*p13*v*(2*(p35+p34)*xm2+p23**2+p13**2)/(p14*p35*p
     2   45*s)-(n2+1)*v*xm2*(4*p12*xm2+s**2+2*p23*p25+2*p13*p15)/(n2*p13
     3   *p15*p23*p25)/8+n2*v*((p35-2*xm2)/(p15*p23*s)-s/(p14*p24*p35)
     4   /4)*xm2+v*((p12-3*xm2)/(p14*p24*p35)+(-2*p23**2+5*p14*p23+3*p
     5   13*p23)/(p13*p15*p23*p25))*xm2+2*n2*(-s*(s+2*p35)/8+p15*p25+4
     6   *p14*p24+p13*p23)*v*xm2/(p14*p24*p35*s)+4*n2*(2*p13**2/(p24*p35
     7   **2*s)+(p45**2+p35**2+p34**2)/(p35**2*s**2))*v*xm2+tmp0+ans
      fbb1 = ans
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC    ALTERNATIVE ROUTINES (courtesy of V. Yundin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function virtuality(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3),virtuality
      integer mu
      virtuality=(p1(0)+p2(0))**2
      do mu=1,3
         virtuality=virtuality-(p1(mu)+p2(mu))**2
      enddo
      end
      

      subroutine alt_ttqqg(p,flav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      double precision p(0:3,nlegborn),pp(0:3,nlegborn),amp2
      integer flav(nlegborn),mu
      double precision s12,s15,s23,s34,s45,mqt,ave_fac
      double precision virtuality,dotp
      external virtuality,dotp
      
c     crossing
      if ((flav(1)*flav(2).ne.0) .and. (flav(1)+flav(2).eq.0)) then
         if (flav(1).gt.0) then
c     q qb -> t tb g
         do mu=0,3
            pp(mu,1) = p(mu,3)
            pp(mu,2) = p(mu,4)
            pp(mu,3) = -p(mu,2)
            pp(mu,4) = -p(mu,1)
            pp(mu,5) = p(mu,5)
         enddo
         else
c     qb q -> t tb g
         do mu=0,3
            pp(mu,1) = p(mu,3)
            pp(mu,2) = p(mu,4)
            pp(mu,3) = -p(mu,1)
            pp(mu,4) = -p(mu,2)
            pp(mu,5) = p(mu,5)
         enddo
         endif
         ave_fac=36d0
      elseif(flav(1).gt.0) then
c     q g  -> t tb q
         do mu=0,3
            pp(mu,1) = p(mu,3)
            pp(mu,2) = p(mu,4)
            pp(mu,3) = p(mu,5)
            pp(mu,4) = -p(mu,1)
            pp(mu,5) = -p(mu,2)
         enddo
c The minus sign come from the crossing of an odd number of fermions
       ave_fac=-96d0
 
      elseif(flav(1).lt.0) then
c     qb g  -> t tb qb
         do mu=0,3
            pp(mu,1) = p(mu,3)
            pp(mu,2) = p(mu,4)
            pp(mu,3) = -p(mu,1)
            pp(mu,4) = p(mu,5)
            pp(mu,5) = -p(mu,2)
         enddo
c The minus sign come from the crossing of an odd number of fermions
         ave_fac=-96d0    
      elseif(flav(2).gt.0) then
c     g q  -> t tb q
         do mu=0,3
            pp(mu,1) = p(mu,3)
            pp(mu,2) = p(mu,4)
            pp(mu,3) = p(mu,5)
            pp(mu,4) = -p(mu,2)
            pp(mu,5) = -p(mu,1)
         enddo
c The minus sign come from the crossing of an odd number of fermions
         ave_fac=-96d0 
      elseif(flav(2).lt.0) then
c     g qb  -> t tb qb
         do mu=0,3
            pp(mu,1) = p(mu,3)
            pp(mu,2) = p(mu,4)
            pp(mu,3) = -p(mu,2)
            pp(mu,4) = p(mu,5)
            pp(mu,5) = -p(mu,1)
         enddo
c The minus sign come from the crossing of an odd number of fermions
         ave_fac=-96d0 

      else
         print*, " not yet done"
      endif
         

      
c     compute  invariants
      s12=virtuality(pp(0,1),pp(0,2))
      s15=virtuality(pp(0,1),pp(0,5))
      s23=virtuality(pp(0,2),pp(0,3))
      s34=virtuality(pp(0,3),pp(0,4))
      s45=virtuality(pp(0,4),pp(0,5))
      
            
c     mass
      mqt=kn_masses(3)

      amp2=(16*(-2*mqt**12*s12*s34*(s34 + 9*s45) +
     $     2*mqt**10*(s12**2*s34*(9*s15 - 9*s23 + 2*s45) +
     $     2*(s34**2 + 2*s34*s45 + 2*s45**2)*
     $     (4*s34**2 + 9*s34*s45 + 9*s45**2) +
     $     s12*((-5*s15 + 11*s23 - 9*s34)*s34**2 +
     $     6*(3*s15 + 6*s23 - 7*s34)*s34*s45 - 36*s34*s45**2 -
     $     36*s45**3)) +
     $     mqt**8*(s12**3*s34*(-32*s15 + 32*s23 + s34 - 47*s45) -
     $     8*(2*s15 - s34)*(s34**2 + 2*s34*s45 + 2*s45**2)*
     $     (4*s34**2 + 9*s34*s45 + 9*s45**2) +
     $     2*s12**2*(s34*
     $     (-27*s15**2 + 9*s15*s23 + 18*s23**2 + 41*s15*s34 -
     $     42*s23*s34 + 14*s34**2) +
     $     s34*(12*s15 - 49*s23 + 104*s34)*s45 +
     $     2*(-36*s15 + 36*s23 + 73*s34)*s45**2 + 36*s45**3) +
     $     s12*(s34**2*(40*s15**2 - 40*s23**2 + 48*s23*s34 -
     $     63*s34**2 + 6*s15*(-5*s23 + 4*s34)) -
     $     s34*(18*(s15**2 + 8*s15*s23 + 6*s23**2) -
     $     2*(142*s15 + 43*s23)*s34 + 305*s34**2)*s45 +
     $     4*(108*s15 - 27*s23 - 133*s34)*s34*s45**2 +
     $     4*(108*s15 - 36*s23 - 113*s34)*s45**3 - 72*s45**4)) +
     $     mqt**6*(s12**4*s34*(23*s15 - 23*s23 + 56*s45) +
     $     4*(6*s15**2 - 6*s15*s34 + s34**2)*
     $     (s34**2 + 2*s34*s45 + 2*s45**2)*
     $     (4*s34**2 + 9*s34*s45 + 9*s45**2) +
     $     s12**3*(s34*(64*s15**2 + 18*s15*s23 - 82*s23**2 -
     $     83*s15*s34 + 79*s23*s34 - 14*s34**2) +
     $     2*(-36*(s15 - s23)**2 + (77*s15 + s23)*s34 -
     $     135*s34**2)*s45 +
     $     (144*(s15 - s23) - 251*s34)*s45**2 - 32*s45**3) +
     $     s12**2*(s34*(72*s15**3 - 36*s23**3 - 200*s15**2*s34 +
     $     128*s23**2*s34 - 107*s23*s34**2 + 44*s34**3 +
     $     s15*(-36*s23**2 + 80*s23*s34 + 27*s34**2)) +
     $     2*s34*(44*s15**2 - 82*s15*s23 + 166*s23**2 -
     $     342*s15*s34 + 16*s23*s34 + 221*s34**2)*s45 +
     $     2*(252*s15**2 + 36*s23**2 + 143*s23*s34 +
     $     341*s34**2 - 24*s15*(12*s23 + 25*s34))*s45**2 +
     $     36*(-10*s15 + 4*s23 + 9*s34)*s45**3 + 32*s45**4) +
     $     s12*(s34**2*(-56*s15**3 + 36*s23**3 - 46*s23**2*s34 +
     $     51*s23*s34**2 - 46*s34**3 +
     $     4*s15**2*(2*s23 + 5*s34) +
     $     s15*(52*s23**2 - 82*s23*s34 + 137*s34**2)) +
     $     2*s34*(36*s23*(s15**2 + 3*s15*s23 + s23**2) -
     $     2*(108*s15**2 + 9*s15*s23 + 43*s23**2)*s34 +
     $     (467*s15 + 9*s23)*s34**2 - 144*s34**3)*s45 -
     $     s34*(900*s15**2 + 108*s23**2 + 154*s23*s34 +
     $     563*s34**2 - 4*s15*(108*s23 + 455*s34))*s45**2 -
     $     16*(54*s15**2 + s34*(9*s23 + 28*s34) -
     $     3*s15*(9*s23 + 32*s34))*s45**3 +
     $     18*(12*s15 - 5*s34)*s45**4)) -
     $     mqt**4*(-(s12**5*s34*(-16*s15 + 16*s23 + s34 - 35*s45)) +
     $     8*s15*(2*s15**2 - 3*s15*s34 + s34**2)*
     $     (s34**2 + 2*s34*s45 + 2*s45**2)*
     $     (4*s34**2 + 9*s34*s45 + 9*s45**2) +
     $     s12**4*(s34*(41*s15**2 + 19*s15*s23 - 60*s23**2 -
     $     46*s15*s34 + 44*s23*s34 + 2*s34**2) +
     $     3*(-24*(s15 - s23)**2 + (33*s15 + 5*s23)*s34 -
     $     53*s34**2)*s45 + (64*s15 - 64*s23 - 119*s34)*s45**2
     $     ) + s12**3*
     $     (s34*(64*s15**3 + s15**2*(54*s23 - 163*s34) +
     $     s15*(-54*s23**2 + 45*s23*s34 + 48*s34**2) -
     $     2*(32*s23**3 - 56*s23**2*s34 + 36*s23*s34**2 +
     $     s34**3)) +
     $     (-144*s15*(s15 - s23)**2 +
     $     2*(186*s15**2 - 244*s15*s23 + 187*s23**2)*s34 +
     $     9*(-63*s15 + 5*s23)*s34**2 + 279*s34**3)*s45 +
     $     (360*s15**2 - 432*s15*s23 + 72*s23**2 - 649*s15*s34 +
     $     82*s23*s34 + 340*s34**2)*s45**2 +
     $     8*(-16*s15 + 8*s23 + 15*s34)*s45**3) +
     $     s12**2*(s34*(36*s15*(s15 - s23)*(2*s15**2 + 3*s23**2) +
     $     4*(-59*s15**3 + 30*s15**2*s23 + 32*s23**3)*s34 +
     $     (173*s15**2 - 201*s15*s23 - 44*s23**2)*s34**2 +
     $     2*(5*s15 + 38*s23)*s34**3 + 2*s34**4) +
     $     s34*(268*s15**3 + 108*s23**3 - 404*s23**2*s34 -
     $     61*s23*s34**2 - 211*s34**3 -
     $     68*s15**2*(9*s23 + 14*s34) +
     $     s15*(600*s23**2 + 576*s23*s34 + 961*s34**2))*s45 +
     $     (144*s15*(4*s15**2 - 5*s15*s23 + s23**2) -
     $     2*(804*s15**2 - 397*s15*s23 + 124*s23**2)*s34 +
     $     4*(374*s15 - 9*s23)*s34**2 - 331*s34**3)*s45**2 -
     $     2*(252*s15**2 - 144*s15*s23 - 377*s15*s34 +
     $     7*s23*s34 + 102*s34**2)*s45**3 +
     $     2*(32*s15 - 17*s34)*s45**4) +
     $     s12*(s34**2*(-58*s15**4 + 8*s15**3*(8*s23 + 11*s34) -
     $     3*s15**2*(28*s23**2 + 50*s23*s34 - 19*s34**2) +
     $     s15*(108*s23**3 + 54*s23**2*s34 + 137*s23*s34**2 -
     $     60*s34**3) -
     $     s34*(64*s23**3 + 8*s23**2*s34 + 32*s23*s34**2 +
     $     s34**3)) +
     $     s34*(-18*s15*(s15**3 - 4*s15**2*s23 - 12*s23**3) -
     $     12*(27*s15**3 - 7*s15**2*s23 + 16*s15*s23**2 +
     $     9*s23**3)*s34 +
     $     2*(500*s15**2 - 38*s15*s23 + 51*s23**2)*s34**2 +
     $     (-605*s15 + s23)*s34**3 + 56*s34**4)*s45 +
     $     s34*(-828*s15**3 + 8*s15**2*(81*s23 + 263*s34) +
     $     2*s34*(90*s23**2 + 11*s23*s34 + 56*s34**2) -
     $     s15*(324*s23**2 + 478*s23*s34 + 1203*s34**2))*
     $     s45**2 - 2*
     $     (72*s15**2*(5*s15 - 3*s23) +
     $     12*s15*(-70*s15 + 9*s23)*s34 +
     $     (455*s15 + 27*s23)*s34**2 - 41*s34**3)*s45**3 +
     $     18*(3*s15 - 2*s34)*(4*s15 - s34)*s45**4)) -
     $     s12*s15*(s12 + s15 - s34)*s34*
     $     (s12**2 + 2*s15**2 + 4*s23**2 + 2*s23*s34 + s34**2 +
     $     2*s12*(s23 - s45) - 2*s15*(2*s23 + s34 - s45) -
     $     4*s23*s45 + 2*s45**2)*
     $     (s12**2*(7*s15 - 7*s23 + 8*s45) -
     $     s34*(s23 - s45)*(7*s34 + 9*s45) -
     $     s15**2*(8*s34 + 9*s45) +
     $     s15*(s34*(9*s23 + 8*s34) + 2*(9*s23 + s34)*s45 -
     $     9*s45**2) +
     $     s12*(9*s15**2 + 14*s23*s34 + 9*s23*s45 - 15*s34*s45 -
     $     8*s45**2 + s15*(-9*s23 - 15*s34 + 2*s45))) +
     $     mqt**2*(s12**6*s34*(7*s15 - 7*s23 + 8*s45) +
     $     4*s15**2*(s15 - s34)**2*(s34**2 + 2*s34*s45 + 2*s45**2)*
     $     (4*s34**2 + 9*s34*s45 + 9*s45**2) +
     $     s12**5*(s34*(32*s15**2 + 7*s23*(-2*s23 + 3*s34) -
     $     s15*(18*s23 + 23*s34)) +
     $     (-32*(s15 - s23)**2 + (31*s15 + 7*s23)*s34 -
     $     39*s34**2)*s45 - 24*s34*s45**2) +
     $     s12**4*(s34*(41*s15**3 + 23*s15**2*(s23 - 4*s34) -
     $     28*s23*(s23**2 - s23*s34 + s34**2) +
     $     4*s15*(-9*s23**2 + 15*s23*s34 + 8*s34**2)) +
     $     2*(-36*s15*(s15 - s23)**2 +
     $     (65*s15**2 - 86*s15*s23 + 87*s23**2)*s34 +
     $     (-73*s15 + 3*s23)*s34**2 + 39*s34**3)*s45 +
     $     (32*(3*s15**2 - 4*s15*s23 + s23**2) -
     $     6*(21*s15 + 8*s23)*s34 + 79*s34**2)*s45**2 +
     $     32*s34*s45**3) +
     $     s12**3*(s34*(64*s15**4 - s15**3*(74*s23 + 165*s34) +
     $     28*s23*s34*(3*s23**2 + s34**2) +
     $     3*s15**2*(46*s23**2 + 35*s23*s34 + 46*s34**2) -
     $     4*s15*(32*s23**3 + 7*s23**2*s34 + 36*s23*s34**2 +
     $     8*s34**3)) -
     $     2*(36*s15**2*(s15 - s23)**2 -
     $     (145*s15**3 - 293*s15**2*s23 + 248*s15*s23**2 +
     $     18*s23**3)*s34 +
     $     (214*s15**2 - 170*s15*s23 + 145*s23**2)*s34**2 +
     $     (-147*s15 + 19*s23)*s34**3 + 39*s34**4)*s45 +
     $     (72*s15*(3*s15**2 - 4*s15*s23 + s23**2) +
     $     (-421*s15**2 + 100*s15*s23 - 132*s23**2)*s34 +
     $     2*(157*s15 + 62*s23)*s34**2 - 93*s34**3)*s45**2 +
     $     2*(-48*s15**2 + 5*(5*s23 - 8*s34)*s34 +
     $     s15*(32*s23 + 79*s34))*s45**3 - 16*s34*s45**4) +
     $     s12**2*(s34*(18*s15**2*(s15 - s23)*
     $     (3*s15**2 - 4*s15*s23 + 6*s23**2) -
     $     8*s15*(25*s15**3 - 42*s15**2*s23 + 48*s15*s23**2 -
     $     32*s23**3)*s34 +
     $     (237*s15**3 - 333*s15**2*s23 + 164*s15*s23**2 -
     $     84*s23**3)*s34**2 -
     $     4*(28*s15**2 - 45*s15*s23 + 7*s23**2)*s34**3 +
     $     (25*s15 - 21*s23)*s34**4) +
     $     2*s34*(2*s15*
     $     (39*s15**3 - 85*s15**2*s23 + 51*s15*s23**2 +
     $     54*s23**3) -
     $     2*(135*s15**3 - 174*s15**2*s23 + 130*s15*s23**2 +
     $     18*s23**3)*s34 +
     $     (316*s15**2 - 136*s15*s23 + 93*s23**2)*s34**2 +
     $     (-127*s15 + 15*s23)*s34**3 + 19*s34**4)*s45 +
     $     (72*s15**2*(3*s15**2 - 4*s15*s23 + s23**2) -
     $     2*s15*(356*s15**2 - 293*s15*s23 + 212*s23**2)*
     $     s34 + 2*(399*s15**2 - 12*s15*s23 + 86*s23**2)*
     $     s34**2 - 2*(151*s15 + 50*s23)*s34**3 + 47*s34**4)*
     $     s45**2 + 2*
     $     (-108*s15**3 + 2*s15*(27*s23 - 59*s34)*s34 +
     $     8*s15**2*(9*s23 + 25*s34) +
     $     s34**2*(-52*s23 + 31*s34))*s45**3 +
     $     2*(16*s15**2 - 34*s15*s34 + 17*s34**2)*s45**4) +
     $     s12*(-2*s15**5*s34*(23*s34 + 18*s45) +
     $     s34**3*(s23 - s45)*(7*s34 + 9*s45)*
     $     (4*s23**2 + s34**2 + 2*s23*(s34 - 2*s45) + 2*s45**2)
     $     + s15**3*(-(s34**2*
     $     (164*s23**2 + 246*s23*s34 + 73*s34**2)) +
     $     2*s34*(-108*s23**2 - 34*s23*s34 + 191*s34**2)*
     $     s45 + 12*s34*(36*s23 + 73*s34)*s45**2 +
     $     16*(9*s23 + 35*s34)*s45**3 + 72*s45**4) +
     $     2*s15**4*(55*s34**3 - 34*s34**2*s45 -
     $     162*s34*s45**2 - 108*s45**3 +
     $     3*s23*s34*(19*s34 + 24*s45)) +
     $     s15**2*s34*
     $     (18*s34**4 - 346*s34**3*s45 - 661*s34**2*s45**2 -
     $     404*s34*s45**3 - 126*s45**4 +
     $     108*s23**3*(s34 + 2*s45) +
     $     s23*s34*(205*s34**2 - 98*s34*s45 - 494*s45**2) +
     $     6*s23**2*(41*s34**2 + 22*s34*s45 - 54*s45**2)) +
     $     s15*s34**2*
     $     (-9*s34**4 + 75*s34**3*s45 + 118*s34**2*s45**2 +
     $     74*s34*s45**3 + 72*s45**4 -
     $     8*s23**3*(16*s34 + 27*s45) +
     $     4*s23**2*(-25*s34**2 + 24*s34*s45 + 90*s45**2) +
     $     s23*(-78*s34**3 + 40*s34**2*s45 + 60*s34*s45**2 -
     $     180*s45**3))))))/
     $     (3.*s12**2*(mqt**2 - s15)**2*s34**2*
     $     (mqt**2 - s12 - s15 + s34)**2*(s12 - s34 - s45)*s45)


c     add couplings and averaging factors
      amp2=amp2* (4*pi*st_alpha)**3/ave_fac 
      end

