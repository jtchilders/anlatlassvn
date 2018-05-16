      subroutine setborn(p,bflav,born,bornjk,bornmunu)
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_math.h'
      include '../include/pwhg_flst.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      integer bflav(nlegs)
      double precision p(0:3,nlegs),born,bornjk(nlegs,nlegs),
     $     bornmunu(0:3,0:3,nlegs)
      double precision bborn,bmunu(0:3,0:3,nlegs),bjk(nlegs,nlegs)
      integer j,k,mu,nu

cccccccccccccccccccc
c$$$c     to check PS
c$$$      born=1d0 !:
c$$$      return
cccccccccccccccccccc

      call compborn(p,bflav,bborn,bmunu,bjk)
      do j=1,nlegs
         if(abs(bflav(j)).le.6) then
c     Spin correlated
            if(bflav(j).eq.0) then
               do mu=0,3
                  do nu=0,3
                     bornmunu(mu,nu,j)=bmunu(mu,nu,j)
                  enddo
               enddo
            endif
c     Colour linked
            do k=j+1,nlegs
               if(abs(bflav(k)).le.6) then
                  bornjk(j,k)=bjk(j,k)
                  bornjk(k,j)=bornjk(j,k)
               endif
            enddo
         endif
      enddo
      born=bborn
      end


      subroutine compborn(p,bflav,born,bmunu,bjk)
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_math.h'
      include '../include/pwhg_flst.h'
      include '../include/pwhg_st.h'
      integer bflav(nlegborn)
      double precision p(0:3,nlegborn),born,bmunu(0:3,0:3,nlegborn),
     $     bjk(nlegborn,nlegborn)
      double precision amp2,bbmunu(0:3,0:3,nlegborn),
     $     bbjk(nlegborn,nlegborn)
      integer ileg,ileg2,mu,nu
      if (abs(bflav(3)).le.6.or.abs(bflav(4)).le.6) then
         write(*,*) 'born_ampsq: ERROR in flavor assignement'
         stop
      endif
c     if present, the gluon is assumed to be the last particle
      call born_ampsq_g_last(p,bflav,amp2,bbmunu,bbjk)
cccccccccccccccccccccccccc
      if(amp2.lt.0.) then
         write(*,*) 'WARNING: in compborn, amp2=0'
         write(*,*) 'wrong permutation ?'
         write(*,*) 'PROGRAM STOPS'
         call exit(1) 
      endif
cccccccccccccccccccccccccc
c     Assign born
      born=amp2
c     Assign bmunu
      do ileg=1,nlegborn
         if(bflav(ileg).eq.0) then
            do mu=0,3
               do nu=0,3
                  bmunu(mu,nu,ileg)=bbmunu(mu,nu,ileg)
               enddo
            enddo
         else
            do mu=0,3
               do nu=0,3
                  bmunu(mu,nu,ileg)=0.
               enddo
            enddo
         endif
c     Assign bjk (here all the matrix is filled)
         do ileg2=1,nlegborn
            bjk(ileg,ileg2)=bbjk(ileg,ileg2)
         enddo
      enddo
      end


c     Select and call the proper Born subroutine.
c     NB: If present, a final-state gluon must be the last particle
c     in the list.
      subroutine born_ampsq_g_last(p,bornflav,amp2,bmunu,bjk)
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_flst.h'
      integer bornflav(nlegborn)
      double precision p(0:3,nlegborn)
      double precision amp2,bmunu(0:3,0:3,nlegborn),
     $     bjk(nlegborn,nlegborn)
      integer ferm_type(nlegborn)
      double precision ferm_charge(nlegborn),charge(-5:5)
      integer i,j,k,l,mu,nu
cccccccccccccccccccc
      double precision bmunu_1(0:3,0:3),bmunu_2(0:3,0:3)
      double precision b12,b13,b14,b23,b24,b34
      charge(-5)=     1d0/3.
      charge(-4)=    -2d0/3.  
      charge(-3)=     1d0/3. 
      charge(-2)=    -2d0/3. 
      charge(-1)=     1d0/3.  
      charge(0)=      0d0    
      charge(1)=     -1d0/3. 
      charge(2)=      2d0/3.  
      charge(3)=     -1d0/3. 
      charge(4)=      2d0/3.  
      charge(5)=     -1d0/3.
      do j=1,6
         do k=1,6
            bjk(j,k)=0d0
         enddo
      enddo
cccccccccccccccccccc

c     lepton-antilepton from Z decay
      ferm_type(3) = +1
      ferm_type(4) = -1
      ferm_charge(3) = -1d0
      ferm_charge(4) = +1d0
      
c     i: flavour index of first incoming parton
c     j: flavour index of second incoming parton
c     k,l: flavours of outgoing partons
      i = bornflav(1)
      j = bornflav(2)
      k = bornflav(5)
      l = bornflav(6)      
      ferm_charge(1) = charge(i)
      ferm_charge(2) = charge(j)
      ferm_charge(5) = charge(k)
      ferm_charge(6) = charge(l)

c     assign ferm_type for QCD particles
      if (i.eq.0) then
         ferm_type(1) = 0
      else 
         ferm_type(1) = i/abs(i)
      endif 
      if (j.eq.0) then
         ferm_type(2) = 0
      else 
         ferm_type(2) = j/abs(j)
      endif   
      if (k.eq.0) then
         ferm_type(5) = 0
      else 
         ferm_type(5) = k/abs(k)
      endif   
      if (l.eq.0) then
         ferm_type(6) = 0
      else 
         ferm_type(6) = l/abs(l)
      endif   

c     choose subprocess, assign bjk and bmunu (when needed)
      if ((i.eq.(-j)).and.(k.eq.0).and.(l.eq.0).and.(i.ne.0)) then
c     q q~ -> g g
         call q_aq_to_l_al_g_g(p,ferm_type,ferm_charge,amp2,bmunu_1
     $        ,bmunu_2,b12,b13,b14,b23,b24,b34)
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,5)=bmunu_1(mu,nu)
               bmunu(mu,nu,6)=bmunu_2(mu,nu)
            enddo
         enddo
      elseif ((i.eq.k).and.(j.eq.l).and.(abs(i).ne.abs(j))
     $        .and. ((i*j).ne.0)) then
c     q qp -> q qp (DIFF FLAVOURS)
         call q_qp_to_l_al_q_qp(p,ferm_type,ferm_charge,amp2,
     $        b12,b13,b14,b23,b24,b34)
      elseif ((i.eq.l).and.(j.eq.k).and.(abs(i).ne.abs(j))
     $        .and. ((i*j).ne.0)) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ADDED IN LH, to check.
c     WITH MY FLAVOUR ORDERING, IT SHOULD NOT BE CALLED.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     q qp -> qp q (DIFF FLAVOURS)
         call q_qp_to_l_al_qp_q(p,ferm_type,ferm_charge,amp2,
     $        b12,b13,b14,b23,b24,b34)
         write(*,*) 'Called the LH flavour config'
      elseif ((i.eq.j).and.(k.eq.l).and.(i.eq.k)
     $        .and. ((i*k).ne.0)) then
c     q q -> q q (SAME FLAVOURS)
         call q_q_to_l_al_q_q(p,ferm_type,ferm_charge,amp2,
     $        b12,b13,b14,b23,b24,b34)
      elseif ((i.eq.(-j)).and.(k.eq.(-l)).and.
     $        (abs(i).ne.abs(k)).and.((i*k).ne.0)
     $        .and.(k.gt.0))      then
c     q q~ -> qp qp~ (DIFF FLAVOURS)
         call q_aq_to_l_al_qp_aqp(p,ferm_type,ferm_charge,amp2,
     $        b12,b13,b14,b23,b24,b34)
      elseif ((i.eq.(-j)).and.(k.eq.(-l)).and.(j.eq.k) 
     $        .and.((i*k).ne.0).and.(i.lt.0) )      then
c     q~ q -> q q~ (SAME FLAVOURS)
         call aq_q_to_l_al_q_aq(p,ferm_type,ferm_charge,amp2,
     $        b12,b13,b14,b23,b24,b34)
      elseif ((i.eq.(-j)).and.(k.eq.(-l)).and.(i.eq.k) 
     $        .and.((i*k).ne.0).and.(i.gt.0) )      then
c     q q~ -> q q~ (SAME FLAVOURS)
         call q_aq_to_l_al_q_aq(p,ferm_type,ferm_charge,amp2,
     $        b12,b13,b14,b23,b24,b34)
      elseif ((i.eq.k).and.(j.eq.0).and.(l.eq.0) 
     $        .and.(i.ne.0)) then
c     q g -> q g  
         call q_g_to_l_al_q_g(p,ferm_type,ferm_charge,amp2,bmunu_1
     $        ,bmunu_2,b12,b13,b14,b23,b24,b34)
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,2)=bmunu_1(mu,nu)
               bmunu(mu,nu,6)=bmunu_2(mu,nu)
            enddo
         enddo
      elseif ((j.eq.k).and.(i.eq.0).and.(l.eq.0).and.(j.ne.0)) then
c     g q -> q g  
         call g_q_to_l_al_q_g(p,ferm_type,ferm_charge,amp2,bmunu_1
     $        ,bmunu_2,b12,b13,b14,b23,b24,b34)
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,1)=bmunu_1(mu,nu)
               bmunu(mu,nu,6)=bmunu_2(mu,nu)
            enddo
         enddo
      elseif ((i.eq.0).and.(j.eq.0).and.(k.eq.(-l))
     #        .and.((k*l).ne.0).and.(k.gt.0))        then
c     g g  -> q q~  
         call g_g_to_l_al_q_aq(p,ferm_type,ferm_charge,amp2,bmunu_1
     $        ,bmunu_2,b12,b13,b14,b23,b24,b34)
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,1)=bmunu_1(mu,nu)
               bmunu(mu,nu,2)=bmunu_2(mu,nu)
            enddo
         enddo
      else
         write(*,*) 'WARNING: in born_ampsq_g_last'
         write(*,*) 'did not match flavour string'
         amp2 = -1.d0
      endif

      bjk(1,2)=b12
      bjk(1,5)=b13
      bjk(1,6)=b14
      bjk(2,5)=b23
      bjk(2,6)=b24
      bjk(5,6)=b34
c     
      bjk(2,1)=b12
      bjk(5,1)=b13
      bjk(6,1)=b14
      bjk(5,2)=b23
      bjk(6,2)=b24
      bjk(6,5)=b34
      end


c     Compute the tree-level squared amplitude for the process
c     q(p1) aq(p2) -> Z(p3+p4) g(p5) g(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     It also returns the spin-correlated and the color-linked
c     squared amplitudes.
c
c     q   --->---ggggggggg  g
c                \  
c                 \ 
c                  ZZZZZZZ  Z  +  u-channel + 3g vertex
c                 /
c                /  
c     aq  ---<---ggggggggg  g
c     
c     fermion_type = +1 fermion
c     fermion_type = -1 antifermion
c     fermion_charge = +2/3, -1/3, -2/3, +1/3
      subroutine q_aq_to_l_al_g_g(pphy,fermion_type,fermion_charge,amp2
     $     ,bmunu_1,bmunu_2,born12,born13,born14,born23,born24,born34)
      implicit none
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg),itmp
      double precision fermion_charge(nleg),ferm_charge(nleg),rtmp
      double precision pphy(0:3,nleg),p(0:3,nleg)
      double precision amp2,amp2check
ccccccccccccc
      double precision bmunu_1(0:3,0:3),bmunu_2(0:3,0:3)
      double precision born12,born13,born14,born23,born24,born34
      double precision born12tmp,born13tmp,born14tmp,
     $     born23tmp,born24tmp,born34tmp
      logical cross12
ccccccccccccc
      include '../include/pwhg_st.h'
      include '../include/pwhg_math.h'
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'planar.h'
      integer ic,ctmp
      double precision afact,bfact
      parameter (afact=CF*CF*nc,bfact=-CF/2d0) !c. factors for tree level
      double complex unit_I
      parameter (unit_I=(0,1))
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      double precision pp1(0:3),pp2(0:3),pp5(0:3),pp6(0:3)
      double complex psi1(2,-1:1),psi2(2,-1:1),psi3(2,-1:1),psi4(2,-1:1)
      double precision epsg5(0:3,1:2),epsg6(0:3,1:2)
      double precision p34,dotp
      double complex ccdotp
      double complex jlep(0:3,-1:1),jqua(0:3,-1:1)
      double complex amp_Ta_Tb(-1:1,-1:1),amp_Tb_Ta(-1:1,-1:1),
     $     amp_3gv(-1:1,-1:1)
      integer mu,i,pol5,pol6,hel_lep,hel_qua
      integer utype_q,utype_l
      double precision q_q,v_q,a_q,L_q,R_q
      double precision q_l,v_l,a_l,L_l,R_l
      double precision Zcoup_q(-1:1),Zcoup_l(-1:1)
      double complex prop34V, prop34gamma
      double precision tiny
      parameter (tiny=1.d-5)
ccccccccccccc
      double complex amp_Ta_Tb_vec56(2,2,-1:1,-1:1)
      double complex amp_Tb_Ta_vec56(2,2,-1:1,-1:1)
      double complex amp_3gv_vec56(2,2,-1:1,-1:1)
      integer nu,ipol1,ipol2
      double precision a12,a13,a14,a23,a24,a34
      double precision b12,b13,b14,b23,b24,b34
      parameter (
     $     a12= + (CA-2*CF)*CF/4.    ,
     $     a13= + CA**2 * CF**2 /2.  ,
     $     a14= - CA*CF/4.           ,
     $     a23= - CA*CF/4.           ,
     $     a24= + CA**2 * CF**2 /2.  ,
     $     a34= + CA**3 * CF/4.
     $     )
      parameter (
     $     b12= + (CA-2*CF)*CF/4. *(CA**2 +1),
     $     b13= - CA*CF/4.                   ,
     $     b14= - CA*CF/4.                   ,
     $     b23= - CA*CF/4.                   ,
     $     b24= - CA*CF/4.                   ,
     $     b34= 0.
     $     )
      born12=0.
      born13=0.
      born14=0.
      born23=0.
      born24=0.
      born34=0.
cccccccccccc
c     planar
      do mu=1,ncstructmax
         sigmaofncstruct(mu)=0d0
         do i=1,6
            clineofpart(mu,i)=0
         enddo
      enddo
c     first colour order: 5 before 6 (corresponds to ta_tb)
c     struct 1, 3 clines, (25),(16) (and 56, which is a gluon-gluon)
c     sequence is: psi2 g5 g6 psi1

      ncstruct=1
      nclines=3
      clineofpart(ncstruct,1)=2
      clineofpart(ncstruct,2)=1
      clineofpart(ncstruct,5)=1
      clineofpart(ncstruct,6)=2

c     second colour order: 6 before 5 (corresponds to tb_ta)
c     struct 2, 3 clines, (26),(15) (and 56, which is a gluon-gluon)
c     sequence is: psi2 g6 g5 psi1
      ncstruct=2
      nclines=3
      clineofpart(ncstruct,1)=2
      clineofpart(ncstruct,2)=1
      clineofpart(ncstruct,5)=2
      clineofpart(ncstruct,6)=1
cccccccccccc


c     check Z decay products
      if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with Z decay'
         stop
      endif
     
c     local copy of variables
      do i=1,nleg
         do mu=0,3
            p(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     now only p, ferm_charge and ferm_type should be modified, if
c     needed

c     exchance particle 1 and 2 if needed
      cross12=.false.
      if (ferm_type(1).eq.-1) then
         if (ferm_type(2).eq.1) then
            cross12=.true.
            call exchange_momenta(p(0,1),p(0,2))
            rtmp = ferm_charge(1)
            ferm_charge(1)=-ferm_charge(2)
            ferm_charge(2)=-rtmp
            itmp = ferm_type(1)
            ferm_type(1)=ferm_type(2)
            ferm_type(2)=itmp
         else
            write(*,*) 'Error in quark type 1-2 (q_aq_to_l_al_g_g)'
            stop
         endif
      endif

c     utype = +1 if up-type quark (u,c,ubar,cbar)
c     utype = -1 otherwise
      if (abs(abs(ferm_charge(1))-2d0/3).lt.tiny) then
         utype_q = +1
         q_q = 2d0/3
      elseif (abs(abs(ferm_charge(1))-1d0/3).lt.tiny) then
         utype_q = -1
         q_q = -1d0/3
      else
         write(*,*) 'Wrong charge in q_aq_to_l_al_g_g ',ferm_charge(1)
         stop
      endif

c     as before, for lepton current
      if (abs(abs(ferm_charge(3))-1d0).lt.tiny) then
         utype_l = -1
         q_l = -1d0
      elseif (abs(abs(ferm_charge(3))-0d0).lt.tiny) then
         utype_l = +1
         q_l = 0d0
      else
         write(*,*) 'Wrong charge in q_aq_to_l_al_g_g ',ferm_charge(3)
         stop
      endif
      
      v_q = utype_q*1.d0/2 - 2d0*q_q*ph_sthw2
      a_q = utype_q*1.d0/2d0
      v_l = utype_l*1.d0/2 - 2d0*q_l*ph_sthw2
      a_l = utype_l*1.d0/2d0
      L_q = v_q + a_q
      R_q = v_q - a_q
      L_l = v_l + a_l
      R_l = v_l - a_l
      
c     copy of z couplings useful for do loops
      Zcoup_q(-1)=L_q
      Zcoup_q(1)=R_q
      Zcoup_l(-1)=L_l
      Zcoup_l(1)=R_l

c     define momenta according to fermionic lines
c     gluon momenta always outgoing
      do mu=0,3
         p1(mu) = ferm_type(1)*p(mu,1)
         p2(mu) = ferm_type(2)*p(mu,2)
         p3(mu) = ferm_type(3)*p(mu,3)
         p4(mu) = ferm_type(4)*p(mu,4)
         p5(mu) = p(mu,5)
         p6(mu) = p(mu,6)
      enddo
      p34=dotp(p3,p4)
c     build wave functions from px (bra/ket need positive energies)
c     q
      if (p(0,1).lt.0d0) then
         do mu=0,3
            pp1(mu) = -p(mu,1)
         enddo         
         call ket(pp1,ferm_type(1),psi1)
      else
         call ket(p(0,1),ferm_type(1),psi1)
      endif
c     aq
      if (p(0,2).lt.0d0) then
         do mu=0,3
            pp2(mu) = -p(mu,2)
         enddo         
         call bra(pp2,ferm_type(2),psi2)
      else
         call bra(p(0,2),ferm_type(2),psi2)
      endif   
c     em
      call bra(p(0,3),ferm_type(3),psi3)
c     ep
      call ket(p(0,4),ferm_type(4),psi4)
c     build leptonic current
      do i=-1,1,2
         call bra_gamma_ket(psi3,psi4,i,jlep(0,i))
      enddo
      
      amp2=0d0
      amp2check=0d0
      do pol5=1,2
         do pol6=1,2
c     build gluon wave functions
            if (p(0,5).lt.0d0) then
               do mu=0,3
                  pp5(mu) = -p(mu,5)
               enddo         
               call polvec(pp5,pol5,epsg5(0,pol5))
            else   
               call polvec(p(0,5),pol5,epsg5(0,pol5))
            endif
            if (p(0,6).lt.0d0) then
               do mu=0,3
                  pp6(mu) = -p(mu,6)
               enddo         
               call polvec(pp6,pol6,epsg6(0,pol6))
            else
               call polvec(p(0,6),pol6,epsg6(0,pol6))
            endif
            do hel_lep=-1,1,2            
               do hel_qua=-1,1,2

cccccccccccccccccccccccccccccccccccc
c     first colour order: 5 before 6
cccccccccccccccccccccccccccccccccccc             
                  call bra_slash_ket_g1_g2(psi2,psi1,hel_qua,p2,p1,
     $ p5,epsg5(0,pol5),p6,epsg6(0,pol6),jlep(0,hel_lep),
     $ amp_Ta_Tb(hel_lep,hel_qua))

ccccccccccccccccccccccccccccccccccccc
c     second colour order: 6 before 5
ccccccccccccccccccccccccccccccccccccc         
                  call bra_slash_ket_g1_g2(psi2,psi1,hel_qua,p2,p1,
     $ p6,epsg6(0,pol6),p5,epsg5(0,pol5),jlep(0,hel_lep),
     $ amp_Tb_Ta(hel_lep,hel_qua))

cccccccccccccccccccccccccccccccccccc
c     third configuration: 3g vertex
cccccccccccccccccccccccccccccccccccc
                  call bra_gamma_ket_3gv(psi2,psi1,hel_qua,p2,p1,p5,
     $ epsg5(0,pol5),p6,epsg6(0,pol6),jqua(0,hel_qua))
                  
                  amp_3gv(hel_lep,hel_qua)=
     $ ccdotp(jlep(0,hel_lep),jqua(0,hel_qua))
                  
cccccccccccccccccccccccccc
c     Z/gamma interference
cccccccccccccccccccccccccc
                  prop34V = 1d0/dcmplx(-2*p34-ph_Zmass2,ph_ZmZw) 
                  prop34gamma = 1d0/(-2*p34)
                  
                  amp_Ta_Tb(hel_lep,hel_qua) = ((q_q*q_l*prop34gamma) + 
     $ (1/(2*ph_sthw*ph_cthw)**2 * 
     $ Zcoup_q(hel_qua)*Zcoup_l(hel_lep)*prop34V)) * 
     $ amp_Ta_Tb(hel_lep,hel_qua)
                  
                  amp_Tb_Ta(hel_lep,hel_qua) = ((q_q*q_l*prop34gamma) + 
     $ (1/(2*ph_sthw*ph_cthw)**2 *
     $ Zcoup_q(hel_qua)*Zcoup_l(hel_lep)*prop34V)) *
     $ amp_Tb_Ta(hel_lep,hel_qua)  
                  
                  amp_3gv(hel_lep,hel_qua) = ((q_q*q_l*prop34gamma) + 
     $ (1/(2*ph_sthw*ph_cthw)**2 * 
     $ Zcoup_q(hel_qua)*Zcoup_l(hel_lep)*prop34V)) * 
     $ amp_3gv(hel_lep,hel_qua)              

cccccccccccccccccccccccc
c     final coherent sum
cccccccccccccccccccccccc
c     the new one, which is more compact.
c     (use that 3gv is |ab>-|ba>)
c     it has been checked, it gives the same result as the old one
                  amp2=amp2
     #+ afact *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
     #+ afact *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ bfact *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ bfact *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
c$$$c     the original one, that works fine
c$$$                  amp2check=amp2check + CF*CF*nc*
c$$$     $           ( amp_Ta_Tb(hel_lep,hel_qua)*
c$$$     $           DCONJG(amp_Ta_Tb(hel_lep,hel_qua))
c$$$     $           + amp_Tb_Ta(hel_lep,hel_qua)*
c$$$     $           DCONJG(amp_Tb_Ta(hel_lep,hel_qua)))
c$$$     $           - CF/2*
c$$$     $           ( amp_Ta_Tb(hel_lep,hel_qua)*
c$$$     $           DCONJG(amp_Tb_Ta(hel_lep,hel_qua))
c$$$     $           + amp_Tb_Ta(hel_lep,hel_qua)*
c$$$     $           DCONJG(amp_Ta_Tb(hel_lep,hel_qua)))
c$$$     $           + CF*CA*nc*
c$$$     $           ( amp_3gv(hel_lep,hel_qua)*
c$$$     $           DCONJG(amp_3gv(hel_lep,hel_qua)))            
c$$$     $           + CF*CA/2*nc*
c$$$     $           ( amp_3gv(hel_lep,hel_qua)*
c$$$     $           DCONJG(amp_Ta_Tb(hel_lep,hel_qua) - 
c$$$     $           amp_Tb_Ta(hel_lep,hel_qua))
c$$$     $           + (amp_Ta_Tb(hel_lep,hel_qua) - 
c$$$     $           amp_Tb_Ta(hel_lep,hel_qua)) * 
c$$$     $           DCONJG(amp_3gv(hel_lep,hel_qua)))
c$$$                  print *, amp2/amp2check

cccccccccccc
c     planar
cccccccccccc
c     first colour order: 5 before 6 (corresponds to ta_tb)
c     struct 1, 3 clines, (25),(16) (and 56, which is a gluon-gluon)
c     sequence is: psi2 g5 g6 psi1
         sigmaofncstruct(1)=sigmaofncstruct(1)+
     $   (amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     $   DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))

c     second colour order: 6 before 5 (corresponds to tb_ta)
c     struct 2, 3 clines, (26),(15) (and 56, which is a gluon-gluon)
c     sequence is: psi2 g6 g5 psi1
         sigmaofncstruct(2)=sigmaofncstruct(2)+
     $   (amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     $   DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))

ccccccccccccccccccccccccccccccccccccc
c     needed to build spin projection
ccccccccccccccccccccccccccccccccccccc
                  amp_Ta_Tb_vec56(pol5,pol6,hel_lep,hel_qua)=0.
                  amp_Tb_Ta_vec56(pol5,pol6,hel_lep,hel_qua)=0.
                  amp_3gv_vec56(pol5,pol6,hel_lep,hel_qua)=0.
                  
                  amp_Ta_Tb_vec56(pol5,pol6,hel_lep,hel_qua)=
     $                 amp_Ta_Tb_vec56(pol5,pol6,hel_lep,hel_qua)+
     $                 amp_Ta_Tb(hel_lep,hel_qua)
                  
                  amp_Tb_Ta_vec56(pol5,pol6,hel_lep,hel_qua)=
     $                 amp_Tb_Ta_vec56(pol5,pol6,hel_lep,hel_qua)+
     $                 amp_Tb_Ta(hel_lep,hel_qua)
                  
                  amp_3gv_vec56(pol5,pol6,hel_lep,hel_qua)=
     $                 amp_3gv_vec56(pol5,pol6,hel_lep,hel_qua)+
     $                 amp_3gv(hel_lep,hel_qua)

cccccccccccccccccccccccccccccccccc
c     needed to build color-linked
cccccccccccccccccccccccccccccccccc

cccccccccccccc
c     Using |3g>= |ab> - |ba>
c$$$            amp2=amp2 
c$$$     #+ CF*CF*nc*
c$$$     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
c$$$     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
c$$$     #+ CF*CF*nc*
c$$$     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
c$$$     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
c$$$     #- CF/2*
c$$$     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
c$$$     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
c$$$     #- CF/2*
c$$$     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
c$$$     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))


                  born12=born12 
     #+ a12 *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
     #+ a12 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b12 *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b12 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))

                  born13=born13
     #+ a14 *  
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
     #+ a13 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b13 *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b13 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))

                  born14=born14
     #+ a13 * 
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
     #+ a14 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b14 *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b14 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))

                  born23=born23
     #+ a24 * 
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
     #+ a23 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b23 *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b23 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
                  
                  born24=born24
     #+ a23 * 
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
     #+ a24 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b24 *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b24 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))

                  born34=born34
     #+ a34 * 
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))
     #+ a34 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b34 *
     #(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))
     #+ b34 *
     #(amp_Tb_Ta(hel_lep,hel_qua)-amp_3gv(hel_lep,hel_qua))*
     #DCONJG(amp_Ta_Tb(hel_lep,hel_qua)+amp_3gv(hel_lep,hel_qua))


               enddo      
            enddo
         enddo
      enddo ! close all the helicity loops

ccccccccccccccccccccc
c     spin projection
ccccccccccccccccccccc
c     bmunu_1 (epsg5)
      do mu=0,3
      do nu=0,3
         bmunu_1(mu,nu) = 0d0
         do ipol1=1,2
         do ipol2=1,2
         do hel_lep=-1,1,2             
         do hel_qua=-1,1,2   
         do pol6=1,2
            bmunu_1(mu,nu) = bmunu_1(mu,nu) +
     $      (epsg5(mu,ipol1) * epsg5(nu,ipol2))*
     $       (
     $           CF*CF*nc*
     $           ( amp_Ta_Tb_vec56(ipol1,pol6,hel_lep,hel_qua)*
     $      DCONJG(amp_Ta_Tb_vec56(ipol2,pol6,hel_lep,hel_qua))
     $           + amp_Tb_Ta_vec56(ipol1,pol6,hel_lep,hel_qua)*
     $      DCONJG(amp_Tb_Ta_vec56(ipol2,pol6,hel_lep,hel_qua)))
     $           - CF/2*
     $           ( amp_Ta_Tb_vec56(ipol1,pol6,hel_lep,hel_qua)*
     $      DCONJG(amp_Tb_Ta_vec56(ipol2,pol6,hel_lep,hel_qua))
     $           + amp_Tb_Ta_vec56(ipol1,pol6,hel_lep,hel_qua)*
     $      DCONJG(amp_Ta_Tb_vec56(ipol2,pol6,hel_lep,hel_qua)))
     $           + CF*CA*nc*
     $           ( amp_3gv_vec56(ipol1,pol6,hel_lep,hel_qua)*
     $      DCONJG(amp_3gv_vec56(ipol2,pol6,hel_lep,hel_qua)))            
     $           + CF*CA/2*nc*
     $           ( amp_3gv_vec56(ipol1,pol6,hel_lep,hel_qua)*
     $    DCONJG(amp_Ta_Tb_vec56(ipol2,pol6,hel_lep,hel_qua) - 
     $           amp_Tb_Ta_vec56(ipol2,pol6,hel_lep,hel_qua))
     $           + (amp_Ta_Tb_vec56(ipol1,pol6,hel_lep,hel_qua) - 
     $              amp_Tb_Ta_vec56(ipol1,pol6,hel_lep,hel_qua)) * 
     $         DCONJG(amp_3gv_vec56(ipol2,pol6,hel_lep,hel_qua)))
     $       )    
         enddo
         enddo
         enddo
         enddo
         enddo
      enddo
      enddo

c     bmunu_2 (epsg6)
      do mu=0,3
      do nu=0,3
         bmunu_2(mu,nu) = 0d0
         do ipol1=1,2
         do ipol2=1,2
         do hel_lep=-1,1,2             
         do hel_qua=-1,1,2   
         do pol5=1,2
            bmunu_2(mu,nu) = bmunu_2(mu,nu) +
     $      (epsg6(mu,ipol1) * epsg6(nu,ipol2))*
     $       (
     $           CF*CF*nc*
     $           ( amp_Ta_Tb_vec56(pol5,ipol1,hel_lep,hel_qua)*
     $      DCONJG(amp_Ta_Tb_vec56(pol5,ipol2,hel_lep,hel_qua))
     $           + amp_Tb_Ta_vec56(pol5,ipol1,hel_lep,hel_qua)*
     $      DCONJG(amp_Tb_Ta_vec56(pol5,ipol2,hel_lep,hel_qua)))
     $           - CF/2*
     $           ( amp_Ta_Tb_vec56(pol5,ipol1,hel_lep,hel_qua)*
     $      DCONJG(amp_Tb_Ta_vec56(pol5,ipol2,hel_lep,hel_qua))
     $           + amp_Tb_Ta_vec56(pol5,ipol1,hel_lep,hel_qua)*
     $      DCONJG(amp_Ta_Tb_vec56(pol5,ipol2,hel_lep,hel_qua)))
     $           + CF*CA*nc*
     $           ( amp_3gv_vec56(pol5,ipol1,hel_lep,hel_qua)*
     $      DCONJG(amp_3gv_vec56(pol5,ipol2,hel_lep,hel_qua)))            
     $           + CF*CA/2*nc*
     $           ( amp_3gv_vec56(pol5,ipol1,hel_lep,hel_qua)*
     $    DCONJG(amp_Ta_Tb_vec56(pol5,ipol2,hel_lep,hel_qua) - 
     $           amp_Tb_Ta_vec56(pol5,ipol2,hel_lep,hel_qua))
     $           + (amp_Ta_Tb_vec56(pol5,ipol1,hel_lep,hel_qua) - 
     $              amp_Tb_Ta_vec56(pol5,ipol1,hel_lep,hel_qua)) * 
     $         DCONJG(amp_3gv_vec56(pol5,ipol2,hel_lep,hel_qua)))
     $       )
         enddo
         enddo
         enddo
         enddo
         enddo
      enddo
      enddo

ccccccccccccccccccccccccccccccccccccccccccccc
c     coupling costants and averaging factors
ccccccccccccccccccccccccccccccccccccccccccccc 
      amp2 = amp2*ph_unit_e**4 * (4*pi*st_alpha)**2 
      amp2=  amp2/4/nc/nc/2

      born12 = born12*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born12tmp = born12/4/nc/nc/2

      born13 = born13*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born13tmp = born13/4/nc/nc/2

      born14 = born14*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born14tmp = born14/4/nc/nc/2

      born23 = born23*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born23tmp = born23/4/nc/nc/2

      born24 = born24*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born24tmp = born24/4/nc/nc/2

      born34 = born34*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born34tmp = born34/4/nc/nc/2

cccccccccccccccccccccccccccccccccccccccccc
      if(cross12) then
c     exchance 1 and 2 color linked if 1 and 2 were
c     exchanged at the beginning
c     (1q  2aq 3g 4g )
c            |
c            |
c            v
c     (2aq 1q  3g 4g )
         born12=born12tmp
         born23=born13tmp
         born24=born14tmp
         born13=born23tmp
         born14=born24tmp
         born34=born34tmp

c     planar: exchange 1-2
         if(ncstruct.ne.2) then
            write(*,*) 'q_aq_to_l_al_g_g, planar'
            call exit(1)
         endif
         do ic=1,ncstruct
            ctmp=clineofpart(ic,1)
            clineofpart(ic,1)=clineofpart(ic,2)
            clineofpart(ic,2)=ctmp
         enddo

      else
         born12 = born12tmp
         born13 = born13tmp
         born14 = born14tmp
         born23 = born23tmp
         born24 = born24tmp
         born34 = born34tmp
      endif
ccccccccccccccccccccccccccccccccccccccccc

      do mu=0,3
         do nu=0,3
            bmunu_1(mu,nu) = bmunu_1(mu,nu) * 
     $           ph_unit_e**4 * (4*pi*st_alpha)**2 /4/nc/nc/2
            bmunu_2(mu,nu) = bmunu_2(mu,nu) * 
     $           ph_unit_e**4 * (4*pi*st_alpha)**2 /4/nc/nc/2
         enddo
      enddo

c$$$c      check
c$$$      if((ferm_type(5).eq.0).and.(ferm_type(6).eq.0)) then
c$$$         print*, (born12+born13+born14)/amp2
c$$$      endif
      end
      

c     Compute the tree-level squared amplitude for the process
c     q(p1) g(p2) -> Z(p3+p4) q(p5) g(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     It also returns the spin-correlated and the color-linked
c     squared amplitudes.
c     Perform crossing and call q_aq_to_l_al_g_g.
      subroutine q_g_to_l_al_q_g(pphy,fermion_type,fermion_charge,amp2
     $     ,bmunu_1,bmunu_2,b12,b13,b14,b23,b24,b34)      
      implicit none
      include 'nlegborn.h'
      include 'planar.h'
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg)
      double precision fermion_charge(nleg),ferm_charge(nleg)
      double precision pphy(0:3,nleg),pp(0:3,nleg)
      double precision amp2
      integer mu,nu,i
cccccccccccccccccccccccccccc
      double precision bmunu_1(0:3,0:3),bmunu_2(0:3,0:3)
      double precision b12,b13,b14,b23,b24,b34
      double precision b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp
c$$$      double complex amp_Ta_Tb_vec56(2,2,-1:1,-1:1)
c$$$      double complex amp_Tb_Ta_vec56(2,2,-1:1,-1:1)
c$$$      double complex amp_3gv_vec56(2,2,-1:1,-1:1)
c$$$      integer ipol1,ipol2
cccccccccccccccccccccccccccc
      integer ic,ctmp

c     copy of local variables     
      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo
      
c     exchange initial gluon <-> final quark
      do mu=0,3
         pp(mu,5) = -pphy(mu,2)
         pp(mu,2) = -pphy(mu,5)
      enddo
      
c     assign type and charge informations to do the crossing
c     NB: no useful information is in ferm_type(1) and ferm_charge(1),
c     since they represented gluons.
c     NOTE the MINUS sign!!!
      ferm_type(2) = -ferm_type(5)
      ferm_charge(2) = -ferm_charge(5)
      call q_aq_to_l_al_g_g(pp,ferm_type,ferm_charge,amp2,bmunu_1
     $     ,bmunu_2,b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp)
      
c     correct for color average
      amp2 = amp2 *2d0 *3d0/8d0
cccccccccccccccc
      do mu=0,3
         do nu=0,3
            bmunu_1(mu,nu) = bmunu_1(mu,nu) * 
     $           2d0 *3d0/8d0
            bmunu_2(mu,nu) = bmunu_2(mu,nu) * 
     $           2d0 *3d0/8d0
         enddo
      enddo
c     (1q 2aq 3g 4g )
c            |
c            |
c            v
c     (1q 3g  2q 4g)
      b13=b12tmp *2d0 *3d0/8d0
      b12=b13tmp *2d0 *3d0/8d0
      b14=b14tmp *2d0 *3d0/8d0
      b23=b23tmp *2d0 *3d0/8d0
      b34=b24tmp *2d0 *3d0/8d0
      b24=b34tmp *2d0 *3d0/8d0
cccccccccccccccc

c$$$c     check
c$$$      if((fermion_type(1).ne.0).and.(fermion_type(2).eq.0)) then
c$$$         print*, (b12+b13+b14)/amp2
c$$$      endif

c     planar: exchange 2-5
      if(ncstruct.ne.2) then
         write(*,*) 'problem in q_g_to_l_al_q_g'
         call exit(1)
      endif
      do ic=1,ncstruct
         ctmp=clineofpart(ic,2)
         clineofpart(ic,2)=clineofpart(ic,5)
         clineofpart(ic,5)=ctmp
      enddo

      end
      
      
c     Compute the tree-level squared amplitude for the process
c     g(p1) q(p2) -> Z(p3+p4) q(p5) g(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     It also returns the spin-correlated and the color-linked
c     squared amplitudes.
c     Perform crossing and call q_aq_to_l_al_g_g.      
      subroutine g_q_to_l_al_q_g(pphy,fermion_type,fermion_charge,amp2
     $     ,bmunu_1,bmunu_2,b12,b13,b14,b23,b24,b34)
      implicit none
      include 'nlegborn.h'
      include 'planar.h'
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg)
      double precision fermion_charge(nleg),ferm_charge(nleg)
      double precision pphy(0:3,nleg),pp(0:3,nleg)
      double precision amp2
      integer mu,nu,i
cccccccccccccccccccccccccccc
      double precision bmunu_1(0:3,0:3),bmunu_2(0:3,0:3)
      double precision b12,b13,b14,b23,b24,b34
      double precision b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp
c$$$      double complex amp_Ta_Tb_vec56(2,2,-1:1,-1:1)
c$$$      double complex amp_Tb_Ta_vec56(2,2,-1:1,-1:1)
c$$$      double complex amp_3gv_vec56(2,2,-1:1,-1:1)
c$$$      integer ipol1,ipol2
cccccccccccccccccccccccccccc
      integer ic,ctmp

c     copy of local variables     
      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     exchange initial gluon <-> final quark
      do mu=0,3
         pp(mu,5) = -pphy(mu,1)
         pp(mu,1) = -pphy(mu,5)
      enddo

c     assign type and charge informations to do the crossing
c     NB: no useful information is in ferm_type(1) and ferm_charge(1),
c     since they represented gluons.
c     NOTE the MINUS sign!!!
      ferm_type(1) = -ferm_type(5)
      ferm_charge(1) = -ferm_charge(5)      
      call q_aq_to_l_al_g_g(pp,ferm_type,ferm_charge,amp2,bmunu_1
     $     ,bmunu_2,b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp)
c     correct for color average
      amp2 = amp2 *2d0* 3d0/8d0

cccccccccccccccc
      do mu=0,3
         do nu=0,3
            bmunu_1(mu,nu) = bmunu_1(mu,nu) * 
     $           2d0 *3d0/8d0
            bmunu_2(mu,nu) = bmunu_2(mu,nu) * 
     $           2d0 *3d0/8d0
         enddo
      enddo
c     (1q 2aq 3g 4g )
c            |
c            |
c            v
c     (3g 2q  1q 4g)
      b23=b12tmp *2d0 *3d0/8d0
      b13=b13tmp *2d0 *3d0/8d0
      b34=b14tmp *2d0 *3d0/8d0
      b12=b23tmp *2d0 *3d0/8d0
      b24=b24tmp *2d0 *3d0/8d0
      b14=b34tmp *2d0 *3d0/8d0
cccccccccccccccc

c$$$c     check
c$$$      if((fermion_type(1).eq.0).and.(fermion_type(2).ne.0)) then
c$$$         print*, (b12+b13+b14)/amp2
c$$$      endif

c     planar: exchange 1-5
      if(ncstruct.ne.2) then
         write(*,*) 'problem in g_q_to_l_al_q_g'
         call exit(1)
      endif
      do ic=1,ncstruct
         ctmp=clineofpart(ic,1)
         clineofpart(ic,1)=clineofpart(ic,5)
         clineofpart(ic,5)=ctmp
      enddo

      end
      

c     Compute the tree-level squared amplitude for the process
c     g(p1) g(p2) -> Z(p3+p4) q(p5) aq(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     It also returns the spin-correlated and the color-linked
c     squared amplitudes.
c     Perform crossing and call q_aq_to_l_al_g_g.
      subroutine g_g_to_l_al_q_aq(pphy,fermion_type,fermion_charge,amp2
     $     ,bmunu_1,bmunu_2,b12,b13,b14,b23,b24,b34)      
      implicit none
      include 'nlegborn.h'
      include 'planar.h'
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg)
      double precision fermion_charge(nleg),ferm_charge(nleg)
      double precision pphy(0:3,nleg),pp(0:3,nleg)
      double precision amp2
      integer mu,nu,i
cccccccccccccccccccccccccccc
      double precision bmunu_1(0:3,0:3),bmunu_2(0:3,0:3)
      double precision b12,b13,b14,b23,b24,b34
      double precision b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp
c$$$      double complex amp_Ta_Tb_vec56(2,2,-1:1,-1:1)
c$$$      double complex amp_Tb_Ta_vec56(2,2,-1:1,-1:1)
c$$$      double complex amp_3gv_vec56(2,2,-1:1,-1:1)
c$$$      integer nu,ipol1,ipol2
cccccccccccccccccccccccccccc
      integer ic,ileg
      double precision tmpclineofpart(nleg)

c     copy of local variables     
      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     !: original
c$$$c     exchange initial gluons <-> final quarks
c$$$      do mu=0,3
c$$$         pp(mu,5) = -pphy(mu,1)
c$$$         pp(mu,1) = -pphy(mu,5)
c$$$         pp(mu,6) = -pphy(mu,2)
c$$$         pp(mu,2) = -pphy(mu,6)         
c$$$      enddo
c$$$
c$$$c     assign type and charge informations to do the crossing
c$$$c     NB: no useful information is in ferm_type(1) and ferm_charge(1),
c$$$c     since they represented gluons.
c$$$c     NOTE the MINUS sign!!!
c$$$      ferm_type(1) = -ferm_type(5)
c$$$      ferm_type(2) = -ferm_type(6)
c$$$      ferm_charge(1) = -ferm_charge(5)
c$$$      ferm_charge(2) = -ferm_charge(6)
c$$$
c$$$      ferm_type(5) = 0
c$$$      ferm_type(6) = 0
c$$$      ferm_charge(5) = 0
c$$$      ferm_charge(6) = 0


c     !:
c     exchange initial gluons <-> final quarks
      do mu=0,3
         pp(mu,5) = -pphy(mu,1)
         pp(mu,2) = -pphy(mu,5)
         pp(mu,6) = -pphy(mu,2)
         pp(mu,1) = -pphy(mu,6)         
      enddo

c     assign type and charge informations to do the crossing
c     NB: no useful information is in ferm_type(1) and ferm_charge(1),
c     since they represented gluons.
c     NOTE the MINUS sign!!!
      ferm_type(2) = -ferm_type(5)
      ferm_type(1) = -ferm_type(6)
      ferm_charge(2) = -ferm_charge(5)
      ferm_charge(1) = -ferm_charge(6)

      ferm_type(5) = 0
      ferm_type(6) = 0
      ferm_charge(5) = 0
      ferm_charge(6) = 0



      call q_aq_to_l_al_g_g(pp,ferm_type,ferm_charge,amp2,bmunu_1
     $     ,bmunu_2,b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp)

c     correct for color average
      amp2 = amp2 *2d0 *3d0/8 *3d0/8
ccccccccccccccccccccccccc     
      do mu=0,3
         do nu=0,3
            bmunu_1(mu,nu) = bmunu_1(mu,nu) * 
     $           2d0 *3d0/8 *3d0/8
            bmunu_2(mu,nu) = bmunu_2(mu,nu) * 
     $           2d0 *3d0/8 *3d0/8
         enddo
      enddo
c     (1q 2aq 3g 4g )
c            |
c            |
c            v
c     (3g 4g  2q 1aq)
      b34=b12tmp *2d0 *3d0/8d0 *3d0/8d0 
      b14=b13tmp *2d0 *3d0/8d0 *3d0/8d0
      b24=b14tmp *2d0 *3d0/8d0 *3d0/8d0
      b13=b23tmp *2d0 *3d0/8d0 *3d0/8d0
      b23=b24tmp *2d0 *3d0/8d0 *3d0/8d0
      b12=b34tmp *2d0 *3d0/8d0 *3d0/8d0
ccccccccccccccccccccccccc

c$$$c     check
c$$$      if((fermion_type(1).eq.0).and.(fermion_type(2).eq.0)) then
c$$$         print*, (b12+b13+b14)/amp2
c$$$      endif

c     planar: exchange as follows
c     1->5, 5->2, 2->6, 6->1

      if(ncstruct.ne.2) then
         write(*,*) 'problem in g_g_to_l_al_q_aq'
         call exit(1)
      endif
      do ic=1,ncstruct
c     do a copy of uncrossed planar links
         do ileg=1,nlegborn
            tmpclineofpart(ileg)=clineofpart(ic,ileg)
         enddo
c     do the crossings
         clineofpart(ic,5)=tmpclineofpart(1)
         clineofpart(ic,6)=tmpclineofpart(2)
         clineofpart(ic,2)=tmpclineofpart(5)
         clineofpart(ic,1)=tmpclineofpart(6)
      enddo

      end
      

      
c     Compute the tree-level squared amplitude for the process
c     q(p1) q(p2) -> Z(p3+p4) q(p5) q(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     NB: flavour of q(p1) MUST be equal from that of q(p2).
c     If it's not the case, use q_qp_to_l_al_q_qp.
c     It also returns the color-linked squared amplitudes.
c
c     q   --->---g--->--- q
c                g
c                g                +  u-channel
c                g
c     q   --->---g--->--- q
c     
c     Z/gamma current inserted on all quark legs
c
c     fermion_type = +1 fermion
c     fermion_type = -1 antifermion
c     fermion_charge = +2/3, -1/3, -2/3, +1/3
      subroutine q_q_to_l_al_q_q(pphy,fermion_type,fermion_charge,amp2,
     $     born12,born13,born14,born23,born24,born34)   
      implicit none
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg)
      double precision fermion_charge(nleg),ferm_charge(nleg)
      double precision pphy(0:3,nleg),p(0:3,nleg)
      double precision amp2
ccccccccccccc
      double precision born12,born13,born14,born23,born24,born34
ccccccccccccc
      include '../include/pwhg_st.h'
      include '../include/pwhg_math.h'
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'planar.h'
      double complex unit_I
      parameter (unit_I=(0,1))
      double precision px1(0:3,nleg),px2(0:3,nleg)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      double precision pp1(0:3),pp2(0:3),pp5(0:3),pp6(0:3)
      double complex psi1(2,-1:1),psi2(2,-1:1),psi3(2,-1:1),
     $     psi4(2,-1:1),psi5(2,-1:1),psi6(2,-1:1)      
      double precision p34,dotp
      double complex ccdotp
      double complex z1,z2
      double complex jlep(0:3,-1:1),jqua15(0:3,-1:1),jqua26(0:3,-1:1),
     $     jqua16(0:3,-1:1),jqua25(0:3,-1:1),jtemp(0:3,-1:1)
      double complex amp_ljj(4,-1:1,-1:1,-1:1)
      integer mu,i,hel_1,hel_2,hel_lep
      double precision pcurr15(0:3),pcurr26(0:3),
     $     pcurr16(0:3),pcurr25(0:3)
      integer utype_q,utype_l
      double precision q_q,v_q,a_q,L_q,R_q
      double precision q_l,v_l,a_l,L_l,R_l
      double precision Zcoup_q(-1:1),Zcoup_l(-1:1)
      double complex prop34V, prop34gamma
      double precision tiny
      parameter (tiny=1.d-5)
ccccccccccccc
      double precision f12,f13,f14,f23,f24,f34
      double precision g12,g13,g14,g23,g24,g34
      parameter (
     $     f12= + CF/2.            ,
     $     f13= - CF/4.            ,
     $     f14= - CF/4. *(2.-CA**2),
     $     f23= - CF/4. *(2.-CA**2),
     $     f24= - CF/4.            ,
     $     f34= + CF/2.
     $     )
      parameter (
     $     g12= - (CA-2*CF)*CF/4. *(CA**2 +1),
     $     g13= + (CA-2*CF)*CF/4.            ,
     $     g14= + (CA-2*CF)*CF/4.            ,
     $     g23= + (CA-2*CF)*CF/4.            ,
     $     g24= + (CA-2*CF)*CF/4.            ,
     $     g34= - (CA-2*CF)*CF/4. *(CA**2 +1)
     $     )

      born12=0.
      born13=0.
      born14=0.
      born23=0.
      born24=0.
      born34=0.
cccccccccccc
c     planar
      do mu=1,ncstructmax
         sigmaofncstruct(mu)=0d0
         do i=1,6
            clineofpart(mu,i)=0
         enddo
      enddo
c$$$c     struct 1, 2 clines, (15),(26)
c$$$      ncstruct=1
c$$$      nclines=2
c$$$      clineofpart(ncstruct,1)=1
c$$$      clineofpart(ncstruct,2)=2
c$$$      clineofpart(ncstruct,5)=1
c$$$      clineofpart(ncstruct,6)=2
c$$$
c$$$c     struct 2, 2 clines, (16),(25)
c$$$      ncstruct=2
c$$$      nclines=2
c$$$      clineofpart(ncstruct,1)=1
c$$$      clineofpart(ncstruct,2)=2
c$$$      clineofpart(ncstruct,5)=2
c$$$      clineofpart(ncstruct,6)=1


c     struct 1, 2 clines, 
c     currents are (15),(26) [t-ch]
c     colour flow as follows, i.e. (16) (25)
c     1-->---     ---->----5
c             | |
c             v ^
c             | |
c     2-->---     ---->---6

      ncstruct=1
      nclines=2
      clineofpart(ncstruct,1)=1
      clineofpart(ncstruct,2)=2
      clineofpart(ncstruct,5)=2
      clineofpart(ncstruct,6)=1

c     struct 2, 2 clines, 
c     currents are (16),(25) [u-ch]
c     colour flow as follows, i.e. (15) (26)
      ncstruct=2
      nclines=2
      clineofpart(ncstruct,1)=1
      clineofpart(ncstruct,2)=2
      clineofpart(ncstruct,5)=1
      clineofpart(ncstruct,6)=2


cccccccccccc

c     check Z decay products
      if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with Z decay'
         stop
      endif
      
c     local copy of variables
      do i=1,6
         do mu=0,3
            p(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     now only p, ferm_charge and ferm_type should be modified, if
c     needed
      
c     utype = +1 if up-type quark (u,c,ubar,cbar)
c     utype = -1 otherwise
      if (abs(abs(ferm_charge(1))-2d0/3).lt.tiny) then
         utype_q = +1
         q_q = 2d0/3
      elseif (abs(abs(ferm_charge(1))-1d0/3).lt.tiny) then
         utype_q = -1
         q_q = -1d0/3
      else
         write(*,*) 'Wrong charge in q_q_to_l_al_q_q ',ferm_charge(1)
         stop
      endif

c     as before, for lepton current
      if (abs(abs(ferm_charge(3))-1d0).lt.tiny) then
         utype_l = -1
         q_l = -1d0
      elseif (abs(abs(ferm_charge(3))-0d0).lt.tiny) then
         utype_l = +1
         q_l = 0d0
      else
         write(*,*) 'Wrong charge in q_q_to_l_al_q_q ',ferm_charge(3)
         stop
      endif      
      
      v_q = utype_q*1.d0/2d0 - 2d0*q_q*ph_sthw2
      a_q = utype_q*1.d0/2d0
      v_l = utype_l*1.d0/2d0 - 2d0*q_l*ph_sthw2
      a_l = utype_l*1.d0/2d0
      L_q = v_q + a_q
      R_q = v_q - a_q
      L_l = v_l + a_l
      R_l = v_l - a_l
      
c     copy of z couplings useful for do loops      
      Zcoup_q(-1)=L_q
      Zcoup_q(1)=R_q
      Zcoup_l(-1)=L_l
      Zcoup_l(1)=R_l

c     px1: momenta for t-ch topology
c     do the crossings, if needed
c     exchance particle 1 and 5
      if (ferm_type(1).eq.-1) then
         if (ferm_type(5).eq.-1) then
            call exchange_mom(p,1,5,6,px1)
         else
            write(*,*) 'Error in the type of the quark 1-5'
            stop
         endif
      else  
         call exchange_mom(p,1,1,6,px1)    
      endif
c     exchance particle 2 and 6
      if (ferm_type(2).eq.-1) then
         if (ferm_type(6).eq.-1) then
            call exchange_mom(px1,2,6,6,px1)
         else
            write(*,*) 'Error in the type of the quark 2-6'
            stop
         endif
      endif

c     px2: momenta for u-ch topology
c     do the crossings, if needed
c     exchance particle 1 and 6
      if (ferm_type(1).eq.-1) then
         if (ferm_type(6).eq.-1) then
            call exchange_mom(p,1,6,6,px2)
         else
            write(*,*) 'Error in the type of the quark 1-6'
            stop
         endif
      else  
         call exchange_mom(p,1,1,6,px2)   
      endif
c     exchance particle 2 and 5
      if (ferm_type(2).eq.-1) then
         if (ferm_type(5).eq.-1) then
            call exchange_mom(px2,2,5,6,px2)
         else
            write(*,*) 'Error in the type of the quark 2-5'
            stop
         endif
      endif
      
      amp2=0d0      
      do hel_lep=-1,1,2         
      do hel_1=-1,1,2            
      do hel_2=-1,1,2
ccccccccccccccccccc
c     t-ch topology
ccccccccccccccccccc
c     define momenta according to fermionic lines
c     current momenta are built from physical vectors always outgoing
         do mu=0,3
            p1(mu) = ferm_type(1)*px1(mu,1)
            p2(mu) = ferm_type(2)*px1(mu,2)
            p3(mu) = ferm_type(3)*px1(mu,3)
            p4(mu) = ferm_type(4)*px1(mu,4)
            p5(mu) = ferm_type(5)*px1(mu,5)
            p6(mu) = ferm_type(6)*px1(mu,6)
            pcurr15(mu) = p5(mu)-p1(mu)
            pcurr26(mu) = p6(mu)-p2(mu)
         enddo               
         p34=dotp(p3,p4)
c     build wave functions from px (bra/ket need positive energies)
c     q
         if (px1(0,1).lt.0d0) then
            do mu=0,3
               pp1(mu) = -px1(mu,1)
            enddo         
            call ket(pp1,ferm_type(1),psi1)
         else
            call ket(px1(0,1),ferm_type(1),psi1)
         endif                
c     q
         if (px1(0,5).lt.0d0) then
            do mu=0,3
               pp5(mu) = -px1(mu,5)
            enddo         
            call bra(pp5,ferm_type(5),psi5)
         else
            call bra(px1(0,5),ferm_type(5),psi5)
         endif
c     em
         call bra(px1(0,3),ferm_type(3),psi3)
c     ep
         call ket(px1(0,4),ferm_type(4),psi4)
c     q
         if (px1(0,6).lt.0d0) then
            do mu=0,3
               pp6(mu) = -px1(mu,6)
            enddo         
            call bra(pp6,ferm_type(6),psi6)
         else
            call bra(px1(0,6),ferm_type(6),psi6)
         endif
c     q
         if (px1(0,2).lt.0d0) then
            do mu=0,3
               pp2(mu) = -px2(mu,2)
            enddo         
            call ket(pp2,ferm_type(2),psi2)
         else
            call ket(px1(0,2),ferm_type(2),psi2)
         endif
c     build currents (15,34,26)
         do i=-1,1,2
            call bra_gamma_ket(psi3,psi4,i,jlep(0,i))
            call bra_gamma_ket(psi6,psi2,i,jqua26(0,i))
            call bra_gamma_ket(psi5,psi1,i,jqua15(0,i))
         enddo
         
         call bra_gamma_ket_curr(psi6,psi2,hel_2,p6,p2,pcurr15,
     $        jqua15(0,hel_1),jtemp(0,hel_2))
         
         amp_ljj(1,hel_lep,hel_1,hel_2) = 
     $        ccdotp(jlep(0,hel_lep),jtemp(0,hel_2))
         
         call bra_gamma_ket_curr(psi5,psi1,hel_1,p5,p1,pcurr26,
     $        jqua26(0,hel_2),jtemp(0,hel_1))
         
         amp_ljj(2,hel_lep,hel_1,hel_2) = 
     $        ccdotp(jlep(0,hel_lep),jtemp(0,hel_1))               

ccccccccccccccccccc
c     u-ch topology
ccccccccccccccccccc
c     define momenta according to fermionic lines
c     current momenta are built from physical vectors always outgoing
         do mu=0,3
            p1(mu) = ferm_type(1)*px2(mu,1)
            p2(mu) = ferm_type(2)*px2(mu,2)
            p3(mu) = ferm_type(3)*px2(mu,3)
            p4(mu) = ferm_type(4)*px2(mu,4)
            p5(mu) = ferm_type(5)*px2(mu,5)
            p6(mu) = ferm_type(6)*px2(mu,6)
            pcurr16(mu) = p6(mu)-p1(mu)
            pcurr25(mu) = p5(mu)-p2(mu)         
         enddo               
         p34=dotp(p3,p4)
c     build wave functions from px (bra/ket need positive energies)              
c     q
         if (px2(0,1).lt.0d0) then
            do mu=0,3
               pp1(mu) = -px2(mu,1)
            enddo         
            call ket(pp1,ferm_type(1),psi1)
         else
            call ket(px2(0,1),ferm_type(1),psi1)
         endif               
c     q        
         if (px2(0,5).lt.0d0) then
            do mu=0,3
               pp5(mu) = -px2(mu,5)
            enddo         
            call bra(pp5,ferm_type(5),psi5)
         else
            call bra(px2(0,5),ferm_type(5),psi5)
         endif
c     em
         call bra(px2(0,3),ferm_type(3),psi3)
c     ep
         call ket(px2(0,4),ferm_type(4),psi4)
c     q
         if (px2(0,6).lt.0d0) then
            do mu=0,3
               pp6(mu) = -px2(mu,6)
            enddo         
            call bra(pp6,ferm_type(6),psi6)
         else
            call bra(px2(0,6),ferm_type(6),psi6)
         endif
c     q 
         if (px1(0,2).lt.0d0) then
            do mu=0,3
               pp2(mu) = -px2(mu,2)
            enddo         
            call ket(pp2,ferm_type(2),psi2)
         else
            call ket(px2(0,2),ferm_type(2),psi2)
         endif
c     build currents (16,34,25)
         do i=-1,1,2
            call bra_gamma_ket(psi3,psi4,i,jlep(0,i))
            call bra_gamma_ket(psi6,psi1,i,jqua16(0,i))
            call bra_gamma_ket(psi5,psi2,i,jqua25(0,i))
         enddo
         
         call bra_gamma_ket_curr(psi6,psi1,hel_1,p6,p1,pcurr25,
     $        jqua25(0,hel_2),jtemp(0,hel_1))
         
         amp_ljj(3,hel_lep,hel_1,hel_2) = 
     $        ccdotp(jlep(0,hel_lep),jtemp(0,hel_1))
         
         call bra_gamma_ket_curr(psi5,psi2,hel_2,p5,p2,pcurr16,
     $        jqua16(0,hel_1),jtemp(0,hel_2))
         
         amp_ljj(4,hel_lep,hel_1,hel_2) = 
     $        ccdotp(jlep(0,hel_lep),jtemp(0,hel_2))
         
cccccccccccccccccccccccccc
c     Z/gamma interference
cccccccccccccccccccccccccc
         prop34V = 1d0/dcmplx(-2*p34-ph_Zmass2,ph_ZmZw) 
         prop34gamma = 1d0/(-2*p34)                    

         amp_ljj(1,hel_lep,hel_1,hel_2) = 
     $        (((q_q*q_l)*prop34gamma) + 
     $        (1/(2*ph_sthw*ph_cthw)**2 * 
     $        Zcoup_q(hel_2)*Zcoup_l(hel_lep)*prop34V))*
     $        amp_ljj(1,hel_lep,hel_1,hel_2)
         
         amp_ljj(2,hel_lep,hel_1,hel_2) = 
     $        (((q_q*q_l)*prop34gamma) + 
     $        (1/(2*ph_sthw*ph_cthw)**2 * 
     $        Zcoup_q(hel_1)*Zcoup_l(hel_lep)*prop34V))*
     $        amp_ljj(2,hel_lep,hel_1,hel_2)
               
         amp_ljj(3,hel_lep,hel_1,hel_2) = 
     $        (((q_q*q_l)*prop34gamma) + 
     $        (1/(2*ph_sthw*ph_cthw)**2 * 
     $        Zcoup_q(hel_1)*Zcoup_l(hel_lep)*prop34V))*
     $        amp_ljj(3,hel_lep,hel_1,hel_2)
         
         amp_ljj(4,hel_lep,hel_1,hel_2) = 
     $        (((q_q*q_l)*prop34gamma) + 
     $        (1/(2*ph_sthw*ph_cthw)**2 * 
     $        Zcoup_q(hel_2)*Zcoup_l(hel_lep)*prop34V))*
     $        amp_ljj(4,hel_lep,hel_1,hel_2)
         
cccccccccccccccccccccccc
c     final coherent sum
cccccccccccccccccccccccc
c     minus sign for fermion statistic
         z1=amp_ljj(1,hel_lep,hel_1,hel_2)+
     $        amp_ljj(2,hel_lep,hel_1,hel_2)
         
         z2=-amp_ljj(3,hel_lep,hel_1,hel_2)-
     $        amp_ljj(4,hel_lep,hel_1,hel_2)
         
         amp2 = amp2 + 
     $        CF*nc/2*(z1*DCONJG(z1)+z2*DCONJG(z2)) 
         
         if (hel_1.eq.hel_2) then
            amp2 = amp2
     $           -CF/2*(z1*DCONJG(z2)+z2*DCONJG(z1))               
         endif   

cccccccccccc
c     planar
cccccccccccc
c     struct 1,
c     z1*z1 corresponds to (15),(26) as currents [tch]
c     2 clines, (16),(25)
         sigmaofncstruct(1)=sigmaofncstruct(1)+z1*DCONJG(z1)

c     struct 2,
c     z2*z2 corresponds to (16),(25) as currents [uch]
c     2 clines, (15),(26)
         sigmaofncstruct(2)=sigmaofncstruct(2)+z2*DCONJG(z2)

cccccccccccccccccccccccccccccccccc
c     needed to build color-linked
cccccccccccccccccccccccccccccccccc
                  born12=born12 
     #+ f12 *(z1*DCONJG(z1)+z2*DCONJG(z2))
                  if (hel_1.eq.hel_2) born12=born12
     #+ g12 *(z1*DCONJG(z2)+z2*DCONJG(z1))

                  born13=born13 
     #+ f13 * z1*DCONJG(z1)
     #+ f14 * z2*DCONJG(z2)
                  if (hel_1.eq.hel_2) born13=born13
     #+ g13 *(z1*DCONJG(z2)+z2*DCONJG(z1))

                  born14=born14 
     #+ f14 * z1*DCONJG(z1)
     #+ f13 * z2*DCONJG(z2)
                  if (hel_1.eq.hel_2) born14=born14
     #+ g14 *(z1*DCONJG(z2)+z2*DCONJG(z1))

                  born23=born23 
     #+ f23 * z1*DCONJG(z1)
     #+ f24 * z2*DCONJG(z2)
                  if (hel_1.eq.hel_2) born23=born23
     #+ g23 *(z1*DCONJG(z2)+z2*DCONJG(z1))

                  born24=born24 
     #+ f24 * z1*DCONJG(z1)
     #+ f23 * z2*DCONJG(z2)
                  if (hel_1.eq.hel_2) born24=born24
     #+ g24 *(z1*DCONJG(z2)+z2*DCONJG(z1))

                  born34=born34 
     #+ f34 *(z1*DCONJG(z1)+z2*DCONJG(z2))
                  if (hel_1.eq.hel_2) born34=born34
     #+ g34 *(z1*DCONJG(z2)+z2*DCONJG(z1))
            
      enddo      
      enddo         
      enddo ! close all the helicity loops

ccccccccccccccccccccccccccccccccccccccccccccc
c     coupling costants and averaging factors
ccccccccccccccccccccccccccccccccccccccccccccc 
      amp2 = amp2*ph_unit_e**4 * (4*pi*st_alpha)**2 
      amp2=  amp2/nc/nc/4/2

      born12 = born12*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born12=  born12/nc/nc/4/2

      born13 = born13*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born13=  born13/nc/nc/4/2

      born14 = born14*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born14=  born14/nc/nc/4/2

      born23 = born23*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born23=  born23/nc/nc/4/2

      born24 = born24*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born24=  born24/nc/nc/4/2

      born34 = born34*ph_unit_e**4 * (4*pi*st_alpha)**2 
      born34=  born34/nc/nc/4/2

c$$$c     check
c$$$      if((ferm_type(1).eq.1).and.(ferm_type(2).eq.1).and.
c$$$     $     (ferm_type(5).eq.1)) then
c$$$         print*, (born12+born13+born14)/amp2
c$$$      endif
      end
      

c     Compute the tree-level squared amplitude for the process
c     aq(p1) q(p2) -> Z(p3+p4) q(p5) aq(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     It also returns the color-linked squared amplitudes.
c     Perform crossing and call q_q_to_l_al_q_q.      
      subroutine aq_q_to_l_al_q_aq(pphy,fermion_type,fermion_charge,
     $     amp2,b12,b13,b14,b23,b24,b34)      
      implicit none
      include 'nlegborn.h'
      include 'planar.h'
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg)
      double precision fermion_charge(nleg),ferm_charge(nleg)
      double precision pphy(0:3,nleg),pp(0:3,nleg)
      double precision amp2
      integer mu,i
cccccccccccccccccccccccccccc
      double precision b12,b13,b14,b23,b24,b34
      double precision b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp
cccccccccccccccccccccccccccc
      integer ic,ctmp

c     copy of local variables          
      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     exchange initial antiquark <-> final antiquark
      do mu=0,3
         pp(mu,6) = -pphy(mu,1)
         pp(mu,1) = -pphy(mu,6)         
      enddo
      
c     assign type and charge informations to do the crossing
c     NOTE the MINUS sign!!!
      ferm_type(1) = -ferm_type(6)
      ferm_type(6) = -fermion_type(1)
      ferm_charge(1) = -ferm_charge(6)
      ferm_charge(6) = -fermion_charge(1)
      call q_q_to_l_al_q_q(pp,ferm_type,ferm_charge,amp2,
     $     b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp)

c     correct for different symmetry factor
      amp2= amp2*2d0
cccccccccccccccccccccccccccccc
c     (1q 2q 3q 4q )
c            |
c            |
c            v
c     (4q 2q  3q 1q)
      b24=b12tmp *2d0 
      b34=b13tmp *2d0
      b14=b14tmp *2d0
      b23=b23tmp *2d0
      b12=b24tmp *2d0
      b13=b34tmp *2d0
cccccccccccccccccccccccccccccc

c$$$c     check
c$$$      if((fermion_type(1).lt.0).and.(fermion_type(2).gt.0)) then
c$$$         print*, '1',(b12+b13+b14)/amp2
c$$$         print*, '2',(b12+b23+b24)/amp2
c$$$         print*, '3',(b13+b23+b34)/amp2
c$$$         print*, '4',(b14+b24+b34)/amp2
c$$$      endif

c     planar: exchange 1-6
      if(ncstruct.ne.2) then
         write(*,*) 'problem in aq_q_to_l_al_q_aq'
         call exit(1)
      endif
      do ic=1,ncstruct
         ctmp=clineofpart(ic,1)
         clineofpart(ic,1)=clineofpart(ic,6)
         clineofpart(ic,6)=ctmp
      enddo

      end


c     Compute the tree-level squared amplitude for the process
c     q(p1) aq(p2) -> Z(p3+p4) q(p5) aq(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     It also returns the color-linked squared amplitudes.
c     Perform crossing and call q_q_to_l_al_q_q.
      subroutine q_aq_to_l_al_q_aq(pphy,fermion_type,fermion_charge,
     #     amp2,b12,b13,b14,b23,b24,b34)      
      implicit none
      include 'nlegborn.h'
      include 'planar.h'
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg)
      double precision fermion_charge(nleg),ferm_charge(nleg)
      double precision pphy(0:3,nleg),pp(0:3,nleg)
      double precision amp2
      integer mu,i
cccccccccccccccccccccccccccc
      double precision b12,b13,b14,b23,b24,b34
      double precision b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp
cccccccccccccccccccccccccccc
      integer ic,ctmp

c     copy of local variables           
      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     exchange initial antiquark <-> final antiquark
      do mu=0,3
         pp(mu,6) = -pphy(mu,2)
         pp(mu,2) = -pphy(mu,6)         
      enddo

c     assign type and charge informations to do the crossing
c     NOTE the MINUS sign!!!
      ferm_type(2) = -ferm_type(6)
      ferm_type(6) = -fermion_type(2)
      ferm_charge(2) = -ferm_charge(6)
      ferm_charge(6) = -fermion_charge(2)
      call q_q_to_l_al_q_q(pp,ferm_type,ferm_charge,amp2,
     $     b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp)

c     correct for different symmetry factor
      amp2= amp2*2d0
cccccccccccccccccccccccccccccc
c     (1q 2q 3q 4q )
c            |
c            |
c            v
c     (1q 4q  3q 2q)
      b14=b12tmp *2d0 
      b13=b13tmp *2d0
      b12=b14tmp *2d0
      b34=b23tmp *2d0
      b24=b24tmp *2d0
      b23=b34tmp *2d0
cccccccccccccccccccccccccccccc

c$$$c     check
c$$$      if((fermion_type(1).gt.0).and.(fermion_type(2).lt.0)) then
c$$$         print*, (b12+b13+b14)/amp2
c$$$      endif

c     planar: exchange 2-6
      if(ncstruct.ne.2) then
         write(*,*) 'q_aq_to_l_al_q_aq'
         call exit(1)
      endif
      do ic=1,ncstruct
         ctmp=clineofpart(ic,2)
         clineofpart(ic,2)=clineofpart(ic,6)
         clineofpart(ic,6)=ctmp
      enddo

      end
      
c     Compute the tree-level squared amplitude for the process
c     q(p1) qp(p2) -> Z(p3+p4) q(p5) qp(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     NB: flavour of q(p1) MUST be DIFFFERENT to that of q(p2).
c     If it's not the case, use q_q_to_l_al_q_q.
c     It also returns the color-linked squared amplitudes.
c
c     q   --->---g--->--- q
c                g
c                g          
c                g
c     qp  --->---g--->--- qp
c     
c     Z/gamma current inserted on all quark legs
c
c     fermion_type = +1 fermion
c     fermion_type = -1 antifermion
c     fermion_charge = +2/3, -1/3, -2/3, +1/3      
      subroutine q_qp_to_l_al_q_qp(pphy,fermion_type,fermion_charge,
     $     amp2,born12,born13,born14,born23,born24,born34)
      implicit none
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg),itmp
      double precision fermion_charge(nleg),ferm_charge(nleg),rtmp
      double precision pphy(0:3,nleg),p(0:3,nleg)
      double precision amp2
ccccccccccccc
      double precision born12,born13,born14,born23,born24,born34
      double precision born12tmp,born13tmp,born14tmp,
     $     born23tmp,born24tmp,born34tmp
      logical cross15,cross26
ccccccccccccc
      include '../include/pwhg_st.h'
      include '../include/pwhg_math.h'
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'planar.h'
      integer ic,ctmp
      double complex unit_I
      parameter (unit_I=(0,1))
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      double complex psi1(2,-1:1),psi2(2,-1:1),psi3(2,-1:1),
     $     psi4(2,-1:1),psi5(2,-1:1),psi6(2,-1:1)      
      double precision p34,dotp
      double complex ccdotp
      double complex jlep(0:3,-1:1),jqua15(0:3,-1:1),jqua26(0:3,-1:1),
     $     jtemp(0:3,-1:1)
      double complex amp_ljj(2,-1:1,-1:1,-1:1)
      integer mu,i,hel_15,hel_26,hel_lep
      double precision pcurr15(0:3),pcurr26(0:3)
      integer utype_q,utype_qp,utype_l
      double precision q_q,v_q,a_q,L_q,R_q
      double precision q_qp,v_qp,a_qp,L_qp,R_qp
      double precision q_l,v_l,a_l,L_l,R_l
      double precision Zcoup_q(-1:1),Zcoup_qp(-1:1),Zcoup_l(-1:1)
      double complex prop34V, prop34gamma
      double precision tiny
      parameter (tiny=1.d-5)
ccccccccccccc
      double precision f12,f13,f14,f23,f24,f34
      parameter (
     $     f12= + CF/2.            ,
     $     f13= - CF/4.            ,
     $     f14= - CF/4. *(2.-CA**2),
     $     f23= - CF/4. *(2.-CA**2),
     $     f24= - CF/4.            ,
     $     f34= + CF/2.
     $     )
      born12=0.
      born13=0.
      born14=0.
      born23=0.
      born24=0.
      born34=0.
ccccccccccccc
c     planar
      do mu=1,ncstructmax
         sigmaofncstruct(mu)=0d0
         do i=1,6
            clineofpart(mu,i)=0
         enddo
      enddo
c$$$c     struct 1, 2 clines, (15),(26)
c$$$      ncstruct=1
c$$$      nclines=2
c$$$      clineofpart(ncstruct,1)=1
c$$$      clineofpart(ncstruct,2)=2
c$$$      clineofpart(ncstruct,5)=1
c$$$      clineofpart(ncstruct,6)=2
c     gluon in t-channel!!!
c     color-connected are 16 and 25
      ncstruct=1
      nclines=2
      clineofpart(ncstruct,1)=1
      clineofpart(ncstruct,2)=2
      clineofpart(ncstruct,5)=2
      clineofpart(ncstruct,6)=1
cccccccccccc


c     check Z decay products
      if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with Z decay'
         stop
      endif

c     local copy of variables
      do i=1,nleg
         do mu=0,3
            p(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     now only p, ferm_charge and ferm_type should be modified, if
c     needed

c     exchance particle 1 and 5 if needed
      cross15=.false.
      if (ferm_type(1).eq.-1) then
        if (ferm_type(5).eq.-1) then
           cross15=.true.
            call exchange_momenta(p(0,1),p(0,5))
            rtmp = ferm_charge(1)
            ferm_charge(1)=-ferm_charge(5)
            ferm_charge(5)=-rtmp
            itmp = ferm_type(1)
            ferm_type(1)=ferm_type(5)
            ferm_type(5)=itmp
         else
            write(*,*) 'Error in quark type 1-5 (q_qp_to_l_al_q_qp)'
            stop
         endif
      endif
      
c     exchance particle 2 and 6 if needed
      cross26=.false.
      if (ferm_type(2).eq.-1) then
         if (ferm_type(6).eq.-1) then
            cross26=.true.
            call exchange_momenta(p(0,2),p(0,6))
            rtmp = ferm_charge(2)
            ferm_charge(2)=-ferm_charge(6)
            ferm_charge(6)=-rtmp
            itmp = ferm_type(2)
            ferm_type(2)=ferm_type(6)
            ferm_type(6)=itmp
         else
            write(*,*) 'Error in quark type 2-6 (q_qp_to_l_al_q_qp)'
            stop
         endif
      endif

c     utype = +1 if up-type quark (u,c,ubar,cbar)
c     utype = -1 otherwise
      if (abs(abs(ferm_charge(1))-2d0/3).lt.tiny) then
         utype_q = +1
         q_q = 2d0/3
      elseif (abs(abs(ferm_charge(1))-1d0/3).lt.tiny) then
         utype_q = -1
         q_q = -1d0/3
      else
         write(*,*) 'Wrong charge in q_qp_to_l_al_q_qp ', ferm_charge(1)
         stop
      endif

c     as before, for qp      
      if (abs(abs(ferm_charge(2))-2d0/3).lt.tiny) then
         utype_qp = +1
         q_qp = 2d0/3
      elseif (abs(abs(ferm_charge(2))-1d0/3).lt.tiny) then
         utype_qp = -1
         q_qp = -1d0/3
      else
         write(*,*) 'Wrong charge in q_qp_to_l_al_q_qp ', ferm_charge(2)
         stop
      endif

c     as before, for lepton current
      if (abs(abs(ferm_charge(3))-1d0).lt.tiny) then
         utype_l = -1
         q_l = -1d0
      elseif (abs(abs(ferm_charge(3))-0d0).lt.tiny) then
         utype_l = +1
         q_l = 0d0
      else
         write(*,*) 'Wrong charge in q_qp_to_l_al_q_qp ',ferm_charge(3)
         stop
      endif

      v_q = utype_q*1.d0/2 - 2*q_q*ph_sthw**2 
      a_q = utype_q*1.d0/2
      v_qp = utype_qp*1.d0/2 - 2*q_qp*ph_sthw**2 
      a_qp = utype_qp*1.d0/2
      v_l = utype_l*1.d0/2 - 2*q_l*ph_sthw**2 
      a_l = utype_l*1.d0/2
      L_q = v_q + a_q
      R_q = v_q - a_q
      L_qp = v_qp + a_qp
      R_qp = v_qp - a_qp
      L_l = v_l + a_l
      R_l = v_l - a_l
      
c     copy of z couplings useful for do loops
      Zcoup_q(-1)=L_q
      Zcoup_q(1)=R_q
      Zcoup_qp(-1)=L_qp
      Zcoup_qp(1)=R_qp
      Zcoup_l(-1)=L_l
      Zcoup_l(1)=R_l

c     define momenta according to fermionic lines
c     current momenta are built from physical vectors always outgoing
      do mu=0,3
         p1(mu) = ferm_type(1)*p(mu,1)
         p2(mu) = ferm_type(2)*p(mu,2)
         p3(mu) = ferm_type(3)*p(mu,3)
         p4(mu) = ferm_type(4)*p(mu,4)
         p5(mu) = ferm_type(5)*p(mu,5)
         p6(mu) = ferm_type(6)*p(mu,6)
         pcurr15(mu) = p5(mu)-p1(mu)
         pcurr26(mu) = p6(mu)-p2(mu)
      enddo
      p34=dotp(p3,p4)

c     build wave functions from p (p here are always physical - positive
c     energies)
c     q
      call ket(p(0,1),ferm_type(1),psi1)
c     q
      call bra(p(0,5),ferm_type(5),psi5)
c     em
      call bra(p(0,3),ferm_type(3),psi3)
c     ep
      call ket(p(0,4),ferm_type(4),psi4)
c     qp
      call bra(p(0,6),ferm_type(6),psi6)
c     qp
      call ket(p(0,2),ferm_type(2),psi2)
c     build currents (15,34,26)
      do i=-1,1,2
         call bra_gamma_ket(psi3,psi4,i,jlep(0,i))
         call bra_gamma_ket(psi6,psi2,i,jqua26(0,i))
         call bra_gamma_ket(psi5,psi1,i,jqua15(0,i))
      enddo
      prop34V = 1d0/dcmplx(-2*p34-ph_Zmass2,ph_ZmZw) 
      prop34gamma = 1d0/(-2*p34)
      
      amp2=0d0
      do hel_lep=-1,1,2         
      do hel_15=-1,1,2         
      do hel_26=-1,1,2
         call bra_gamma_ket_curr(psi6,psi2,hel_26,p6,p2,pcurr15,
     $        jqua15(0,hel_15),jtemp(0,hel_26))
         amp_ljj(1,hel_lep,hel_15,hel_26) = 
     $        ccdotp(jlep(0,hel_lep),jtemp(0,hel_26))

         call bra_gamma_ket_curr(psi5,psi1,hel_15,p5,p1,pcurr26,
     $        jqua26(0,hel_26),jtemp(0,hel_15))
         
         amp_ljj(2,hel_lep,hel_15,hel_26) = 
     $        ccdotp(jlep(0,hel_lep),jtemp(0,hel_15))
         
cccccccccccccccccccccccccc
c     Z/gamma interference
cccccccccccccccccccccccccc   
         amp_ljj(1,hel_lep,hel_15,hel_26) = 
     $        ((q_qp*q_l*prop34gamma) + 
     $        (1/(2*ph_sthw*ph_cthw)**2 * 
     $        Zcoup_qp(hel_26)*Zcoup_l(hel_lep)*prop34V))*
     $        amp_ljj(1,hel_lep,hel_15,hel_26)
         
         amp_ljj(2,hel_lep,hel_15,hel_26) = 
     $        ((q_q*q_l*prop34gamma) + 
     $        (1/(2*ph_sthw*ph_cthw)**2 *
     $        Zcoup_q(hel_15)*Zcoup_l(hel_lep)*prop34V))*
     $        amp_ljj(2,hel_lep,hel_15,hel_26)

cccccccccccccccccccccccc
c     final coherent sum
cccccccccccccccccccccccc
         amp2 = amp2 + CF*nc/2*
     $        (amp_ljj(1,hel_lep,hel_15,hel_26)+
     $        amp_ljj(2,hel_lep,hel_15,hel_26))
     $        * DCONJG(amp_ljj(1,hel_lep,hel_15,hel_26)+
     $          amp_ljj(2,hel_lep,hel_15,hel_26))  


cccccccccccc
c     planar
cccccccccccc
c     struct 1, 2 clines, (15),(26)
         sigmaofncstruct(1)=sigmaofncstruct(1)+amp2


cccccccccccccccccccccccccccccccccc
c     needed to build color-linked
cccccccccccccccccccccccccccccccccc         
         born12 = born12 + f12 *
     $        (amp_ljj(1,hel_lep,hel_15,hel_26)+
     $        amp_ljj(2,hel_lep,hel_15,hel_26))
     $        * DCONJG(amp_ljj(1,hel_lep,hel_15,hel_26)+
     $          amp_ljj(2,hel_lep,hel_15,hel_26))  

         born13 = born13 + f13 *
     $        (amp_ljj(1,hel_lep,hel_15,hel_26)+
     $        amp_ljj(2,hel_lep,hel_15,hel_26))
     $        * DCONJG(amp_ljj(1,hel_lep,hel_15,hel_26)+
     $          amp_ljj(2,hel_lep,hel_15,hel_26))  

         born14 = born14 + f14 *
     $        (amp_ljj(1,hel_lep,hel_15,hel_26)+
     $        amp_ljj(2,hel_lep,hel_15,hel_26))
     $        * DCONJG(amp_ljj(1,hel_lep,hel_15,hel_26)+
     $          amp_ljj(2,hel_lep,hel_15,hel_26))  

         born23 = born23 + f23 *
     $        (amp_ljj(1,hel_lep,hel_15,hel_26)+
     $        amp_ljj(2,hel_lep,hel_15,hel_26))
     $        * DCONJG(amp_ljj(1,hel_lep,hel_15,hel_26)+
     $          amp_ljj(2,hel_lep,hel_15,hel_26))  

         born24 = born24 + f24 *
     $        (amp_ljj(1,hel_lep,hel_15,hel_26)+
     $        amp_ljj(2,hel_lep,hel_15,hel_26))
     $        * DCONJG(amp_ljj(1,hel_lep,hel_15,hel_26)+
     $          amp_ljj(2,hel_lep,hel_15,hel_26))  

         born34 = born34 + f34 *
     $        (amp_ljj(1,hel_lep,hel_15,hel_26)+
     $        amp_ljj(2,hel_lep,hel_15,hel_26))
     $        * DCONJG(amp_ljj(1,hel_lep,hel_15,hel_26)+
     $          amp_ljj(2,hel_lep,hel_15,hel_26))  

      enddo      
      enddo         
      enddo ! close all the helicity loops

ccccccccccccccccccccccccccccccccccccccccccccc
c     coupling costants and averaging factors
ccccccccccccccccccccccccccccccccccccccccccccc 
      amp2 = amp2 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      amp2=  amp2 /nc/nc/4

      born12 = born12 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born12tmp = born12 /nc/nc/4

      born13 = born13 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born13tmp = born13 /nc/nc/4

      born14 = born14 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born14tmp = born14 /nc/nc/4

      born23 = born23 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born23tmp = born23 /nc/nc/4

      born24 = born24 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born24tmp = born24 /nc/nc/4

      born34 = born34 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born34tmp = born34 /nc/nc/4


cccccccccccccccccccccccccccccccccccccccccc
      if(cross15.and.(.not.cross26)) then
c     exchance 1 and 5 color linked if 1 and 5 were
c     exchanged at the beginning
c     (1q  2qp  3q  4qp )
c            |
c            |
c            v
c     (3aq 2qp  1aq 4qp)
         born23=born12tmp
         born13=born13tmp
         born34=born14tmp
         born12=born23tmp
         born24=born24tmp
         born14=born34tmp

c     planar: exchange 1-5
         if(ncstruct.ne.1) then
            write(*,*) 'q_qp_to_l_al_q_qp, 15'
            call exit(1)
         endif
         do ic=1,ncstruct
            ctmp=clineofpart(ic,1)
            clineofpart(ic,1)=clineofpart(ic,5)
            clineofpart(ic,5)=ctmp
         enddo


      elseif(cross26.and.(.not.cross15)) then
c     exchance 2 and 6 color linked if 2 and 6 were
c     exchanged at the beginning
c     (1q  2qp  3q  4qp )
c            |
c            |
c            v
c     (1q  4aqp 3q 2aqp)
         born14=born12tmp
         born13=born13tmp
         born12=born14tmp
         born34=born23tmp
         born24=born24tmp
         born23=born34tmp

c     planar: exchange 2-6
         if(ncstruct.ne.1) then
            write(*,*) 'q_qp_to_l_al_q_qp, 26'
            call exit(1)
         endif
         do ic=1,ncstruct
            ctmp=clineofpart(ic,2)
            clineofpart(ic,2)=clineofpart(ic,6)
            clineofpart(ic,6)=ctmp
         enddo

      else
         born12 = born12tmp
         born13 = born13tmp
         born14 = born14tmp
         born23 = born23tmp
         born24 = born24tmp
         born34 = born34tmp
      endif
ccccccccccccccccccccccccccccccccccccccccc


c$$$c     check
c$$$      print*, (born12+born13+born14)/amp2
      end
      
c     Compute the tree-level squared amplitude for the process
c     q(p1) aq(p2) -> Z(p3+p4) qp(p5) aqp(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     NB: flavour of q(p1) MUST be DIFFFERENT to that of qp(p5).
c     If it's not the case, use q_q_to_l_al_q_q.
c     It also returns the color-linked squared amplitudes.
c
c     q   --->---\        --->---  qp
c                 \      /
c                  ggggg/ 
c                 /     \
c     aq  ---<---/       \---<---  aqp
c     
c     Z/gamma current inserted on all quark legs
c
c     fermion_type = +1 fermion
c     fermion_type = -1 antifermion
c     fermion_charge = +2/3, -1/3, -2/3, +1/3      
      subroutine q_aq_to_l_al_qp_aqp(pphy,fermion_type,fermion_charge,
     $     amp2,born12,born13,born14,born23,born24,born34)   
      implicit none
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg),itmp
      double precision fermion_charge(nleg),ferm_charge(nleg),rtmp
      double precision pphy(0:3,nleg),p(0:3,nleg)
      double precision amp2
ccccccccccccc
      double precision born12,born13,born14,born23,born24,born34
      double precision born12tmp,born13tmp,born14tmp,
     $     born23tmp,born24tmp,born34tmp
      logical cross12,cross56
ccccccccccccc  
      include '../include/pwhg_st.h'
      include '../include/pwhg_math.h'
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'planar.h'
      integer ic,ctmp
      double complex unit_I
      parameter (unit_I=(0,1))
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      double complex psi1(2,-1:1),psi2(2,-1:1),psi3(2,-1:1),
     $     psi4(2,-1:1),psi5(2,-1:1),psi6(2,-1:1)    
      double precision p34,dotp
      double complex ccdotp
      double complex jlep(0:3,-1:1),jqua56(0:3,-1:1),jqua12(0:3,-1:1),
     $     jtemp(0:3,-1:1)
      double complex amp_ljj(2,-1:1,-1:1,-1:1)
      integer mu,i,hel_12,hel_56,hel_lep
      double precision pcurr56(0:3),pcurr12(0:3)
      integer utype_q,utype_qp,utype_l
      double precision q_q,v_q,a_q,L_q,R_q
      double precision q_qp,v_qp,a_qp,L_qp,R_qp
      double precision q_l,v_l,a_l,L_l,R_l
      double precision Zcoup_q(-1:1),Zcoup_qp(-1:1),Zcoup_l(-1:1)
      double complex prop34V, prop34gamma
      double precision tiny
      parameter (tiny=1.d-5)
ccccccccccccc
      double precision h12,h13,h14,h23,h24,h34
      parameter (
     $     h12= - CF/4.            ,
     $     h13= - CF/4. *(2.-CA**2),
     $     h14= + CF/2.            ,
     $     h23= + CF/2.            ,
     $     h24= - CF/4. *(2.-CA**2),
     $     h34= - CF/4.
     $     )
      born12=0.
      born13=0.
      born14=0.
      born23=0.
      born24=0.
      born34=0.
ccccccccccccc
c     planar
      do mu=1,ncstructmax
         sigmaofncstruct(mu)=0d0
         do i=1,6
            clineofpart(mu,i)=0
         enddo
      enddo
c$$$c     struct 1, 2 clines, (12),(56)
c$$$      ncstruct=1
c$$$      nclines=2
c$$$      clineofpart(ncstruct,1)=1
c$$$      clineofpart(ncstruct,2)=1
c$$$      clineofpart(ncstruct,5)=2
c$$$      clineofpart(ncstruct,6)=2
      ncstruct=1
      nclines=2
      clineofpart(ncstruct,1)=1
      clineofpart(ncstruct,2)=2
      clineofpart(ncstruct,5)=1
      clineofpart(ncstruct,6)=2
cccccccccccc


c     check Z decay products
      if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with Z decay'
         stop
      endif

c     local copy of variables
      do i=1,nleg
         do mu=0,3
            p(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     now only p, ferm_charge and ferm_type should be modified, if
c     needed
      
c     exchance particle 1 and 2 if needed
      cross12=.false.
      if (ferm_type(1).eq.-1) then
         if (ferm_type(2).eq.1) then
            cross12=.true.
            call exchange_momenta(p(0,1),p(0,2))
            rtmp = ferm_charge(1)
            ferm_charge(1)=-ferm_charge(2)
            ferm_charge(2)=-rtmp
            itmp = ferm_type(1)
            ferm_type(1)=ferm_type(2)
            ferm_type(2)=itmp
         else
            write(*,*) 'Error in quark type 1-2 (q_aq_to_l_al_qp_aqp)'
            stop
         endif
      endif
      
c     exchance particle 5 and 6 if needed
      cross56=.false.
      if (ferm_type(5).eq.-1) then
         if (ferm_type(6).eq.1) then
            cross56=.true.
            call exchange_momenta(p(0,5),p(0,6))
            rtmp = ferm_charge(5)
            ferm_charge(5)=-ferm_charge(6)
            ferm_charge(6)=-rtmp
            itmp = ferm_type(5)
            ferm_type(5)=ferm_type(6)
            ferm_type(6)=itmp
         else
            write(*,*) 'Error in quark type 5-6 (q_aq_to_l_al_qp_aqp)'
            stop
         endif
      endif
      
c     utype = +1 if up-type quark (u,c,ubar,cbar)
c     utype = -1 otherwise
      if (abs(abs(ferm_charge(1))-2d0/3).lt.tiny) then
         utype_q = +1
         q_q = 2d0/3
      elseif (abs(abs(ferm_charge(1))-1d0/3).lt.tiny) then
         utype_q = -1
         q_q = -1d0/3
      else
         write(*,*) 'Wrong charge in q_aq_to_l_al_qp_aqp ',
     $        ferm_charge(1)
         stop
      endif

c     as before, for qp            
      if (abs(abs(ferm_charge(5))-2d0/3).lt.tiny) then
         utype_qp = +1
         q_qp = 2d0/3
      elseif (abs(abs(ferm_charge(5))-1d0/3).lt.tiny) then
         utype_qp = -1
         q_qp = -1d0/3
      else
         write(*,*) 'Wrong charge in q_aq_to_l_al_qp_aqp ',
     $        ferm_charge(5)
         stop
      endif

c     as before, for lepton current      
      if (abs(abs(ferm_charge(3))-1d0).lt.tiny) then
         utype_l = -1
         q_l = -1d0
      elseif (abs(abs(ferm_charge(3))-0d0).lt.tiny) then
         utype_l = +1
         q_l = 0d0
      else
         write(*,*) 'Wrong charge in q_aq_to_l_al_qp_aqp ',
     $        ferm_charge(3)
         stop
      endif                 
      
      v_q = utype_q*1.d0/2 - 2*q_q*ph_sthw2
      a_q = utype_q*1.d0/2
      v_qp = utype_qp*1.d0/2 - 2*q_qp*ph_sthw2
      a_qp = utype_qp*1.d0/2
      v_l = utype_l*1.d0/2 - 2*q_l*ph_sthw2
      a_l = utype_l*1.d0/2
      L_q = v_q + a_q
      R_q = v_q - a_q
      L_qp = v_qp + a_qp
      R_qp = v_qp - a_qp
      L_l = v_l + a_l
      R_l = v_l - a_l
      
c     copy of z couplings useful for do loops      
      Zcoup_q(-1)=L_q
      Zcoup_q(1)=R_q
      Zcoup_qp(-1)=L_qp
      Zcoup_qp(1)=R_qp
      Zcoup_l(-1)=L_l
      Zcoup_l(1)=R_l

c     define momenta according to fermionic lines
c     current momenta are built from physical vectors always outgoing
      do mu=0,3
         p1(mu) = ferm_type(1)*p(mu,1)
         p2(mu) = ferm_type(2)*p(mu,2)
         p3(mu) = ferm_type(3)*p(mu,3)
         p4(mu) = ferm_type(4)*p(mu,4)
         p5(mu) = ferm_type(5)*p(mu,5)
         p6(mu) = ferm_type(6)*p(mu,6)
         pcurr56(mu) = p5(mu)-p6(mu)
         pcurr12(mu) = p2(mu)-p1(mu)
      enddo
      p34=dotp(p3,p4)

c     build wave functions from p (p here are always physical - positive
c     energies)
c     q
      call ket(p(0,1),ferm_type(1),psi1)
c     aq
      call bra(p(0,2),ferm_type(2),psi2)
c     em
      call bra(p(0,3),ferm_type(3),psi3)
c     ep
      call ket(p(0,4),ferm_type(4),psi4)
c     qp
      call bra(p(0,5),ferm_type(5),psi5)
c     aqp
      call ket(p(0,6),ferm_type(6),psi6)
c     build currents (12,34,56)
      do i=-1,1,2
         call bra_gamma_ket(psi3,psi4,i,jlep(0,i))
         call bra_gamma_ket(psi5,psi6,i,jqua56(0,i))
         call bra_gamma_ket(psi2,psi1,i,jqua12(0,i))
      enddo
      prop34V = 1d0/dcmplx(-2*p34-ph_Zmass2,ph_ZmZw) 
      prop34gamma = 1d0/(-2*p34)
      
      amp2=0d0
      do hel_lep=-1,1,2         
         do hel_12=-1,1,2         
            do hel_56=-1,1,2                                             
               call bra_gamma_ket_curr(psi2,psi1,hel_12,p2,p1,pcurr56,
     $              jqua56(0,hel_56),jtemp(0,hel_12))
               
               amp_ljj(1,hel_lep,hel_12,hel_56) = 
     $              ccdotp(jlep(0,hel_lep),jtemp(0,hel_12))
               
               call bra_gamma_ket_curr(psi5,psi6,hel_56,p5,p6,pcurr12,
     $              jqua12(0,hel_12),jtemp(0,hel_56))
               
               amp_ljj(2,hel_lep,hel_12,hel_56) = 
     $              ccdotp(jlep(0,hel_lep),jtemp(0,hel_56))
               
cccccccccccccccccccccccccc
c     Z/gamma interference
cccccccccccccccccccccccccc       
               amp_ljj(1,hel_lep,hel_12,hel_56) = 
     $              ((q_q*q_l*prop34gamma) + 
     $              (1/(2*ph_sthw*ph_cthw)**2 * 
     $              Zcoup_q(hel_12)*Zcoup_l(hel_lep)*prop34V))*
     $              amp_ljj(1,hel_lep,hel_12,hel_56)
               
               amp_ljj(2,hel_lep,hel_12,hel_56) = 
     $              ((q_qp*q_l*prop34gamma) + 
     $              (1/(2*ph_sthw*ph_cthw)**2 *
     $              Zcoup_qp(hel_56)*Zcoup_l(hel_lep)*prop34V))*
     $              amp_ljj(2,hel_lep,hel_12,hel_56)

cccccccccccccccccccccccc
c     final coherent sum
cccccccccccccccccccccccc
               amp2 = amp2 + CF*nc/2*
     $              (amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))
     $              * DCONJG(amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))           

cccccccccccc
c     planar
cccccccccccc
c     struct 1, 2 clines, (12),(56)
         sigmaofncstruct(1)=sigmaofncstruct(1)+amp2

cccccccccccccccccccccccccccccccccc
c     needed to build color-linked
cccccccccccccccccccccccccccccccccc   
               born12 = born12 + h12 *
     $              (amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))
     $              * DCONJG(amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))           

               born13 = born13 + h13 *
     $              (amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))
     $              * DCONJG(amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))           

               born14 = born14 + h14 *
     $              (amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))
     $              * DCONJG(amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))           

               born23 = born23 + h23 *
     $              (amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))
     $              * DCONJG(amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))           

               born24 = born24 + h24 *
     $              (amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))
     $              * DCONJG(amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))           

               born34 = born34 + h34 *
     $              (amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))
     $              * DCONJG(amp_ljj(1,hel_lep,hel_12,hel_56)+
     $              amp_ljj(2,hel_lep,hel_12,hel_56))           

            enddo      
         enddo         
      enddo ! close all the helicity loops

ccccccccccccccccccccccccccccccccccccccccccccc
c     coupling costants and averaging factors
ccccccccccccccccccccccccccccccccccccccccccccc 
      amp2 = amp2 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      amp2=  amp2 /nc/nc/4

      born12 = born12 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born12tmp = born12 /nc/nc/4

      born13 = born13 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born13tmp = born13 /nc/nc/4

      born14 = born14 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born14tmp = born14 /nc/nc/4

      born23 = born23 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born23tmp = born23 /nc/nc/4

      born24 = born24 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born24tmp = born24 /nc/nc/4

      born34 = born34 *ph_unit_e**4 * (4*pi*st_alpha)**2 
      born34tmp = born34 /nc/nc/4

cccccccccccccccccccccccccccccccccccccccccc
      if(cross12) then
c     exchance 1 and 2 color linked if 1 and 2 were
c     exchanged at the beginning
c     (1q  2aq 3qp 4aqp )
c            |
c            |
c            v
c     (2aq 1q  3qp 4aqp)
         born12=born12tmp
         born23=born13tmp
         born24=born14tmp
         born13=born23tmp
         born14=born24tmp
         born34=born34tmp

c     planar: exchange 1-2
         if(ncstruct.ne.1) then
            write(*,*) 'q_aq_to_l_al_qp_aqp, planar'
            call exit(1)
         endif
         do ic=1,ncstruct
            ctmp=clineofpart(ic,1)
            clineofpart(ic,1)=clineofpart(ic,2)
            clineofpart(ic,2)=ctmp
         enddo

      elseif(cross56) then
         write(*,*) 'c. linked problem in q_aq_to_l_al_qp_aqp '
         stop
      else
         born12 = born12tmp
         born13 = born13tmp
         born14 = born14tmp
         born23 = born23tmp
         born24 = born24tmp
         born34 = born34tmp
      endif
ccccccccccccccccccccccccccccccccccccccccc

c$$$c     check
c$$$      print*, (born12+born13+born14)/amp2
      end
      


c     Compute the tree-level squared amplitude for the process
c     q(p1) qp(p2) -> Z(p3+p4) qp(p5) q(p6), with Z -> l-(p3) l+(p4)
c     It uses the bra/ket formalism (HagZep).
c     It also returns the color-linked squared amplitudes.
c     Perform crossing and call q_qp_to_l_al_q_qp.
      subroutine q_qp_to_l_al_qp_q(pphy,fermion_type,fermion_charge,
     #     amp2,b12,b13,b14,b23,b24,b34)      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ADDED IN LH, to check.
c     WITH MY FLAVOUR ORDERING, IT SHOULD NOT BE CALLED.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer nleg
      parameter (nleg=6)
      integer fermion_type(nleg),ferm_type(nleg)
      double precision fermion_charge(nleg),ferm_charge(nleg)
      double precision pphy(0:3,nleg),pp(0:3,nleg)
      double precision amp2
      integer mu,i
cccccccccccccccccccccccccccc
      double precision b12,b13,b14,b23,b24,b34
      double precision b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp
cccccccccccccccccccccccccccc

c     copy of local variables           
      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     exchange final state particles
      do mu=0,3
         pp(mu,5) = pphy(mu,6)
         pp(mu,6) = pphy(mu,5)         
      enddo

c     assign type and charge informations to do the crossing
      ferm_type(5) = ferm_type(6)
      ferm_type(6) = fermion_type(5)
      ferm_charge(5) = ferm_charge(6)
      ferm_charge(6) = fermion_charge(5)
      call q_qp_to_l_al_q_qp(pp,ferm_type,ferm_charge,amp2,
     $     b12tmp,b13tmp,b14tmp,b23tmp,b24tmp,b34tmp)

cccccccccccccccccccccccccccccc
c     (1q 2qp 3q 4qp )
c            |
c            |
c            v
c     (1q 2qp 4qp 3q)
      b12=b12tmp 
      b14=b13tmp
      b13=b14tmp
      b24=b23tmp
      b23=b24tmp
      b34=b34tmp
cccccccccccccccccccccccccccccc
      write(*,*) 'WARNING: called q_qp_to_l_al_qp_q'
      end


      subroutine borncolour_lh
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structure, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'planar.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      character *4 proc
      integer col(3)
      data col/501,502,503/
      save col
      integer ubflav(nlegborn),ileg,ic,ctag,ig1,ig2,iq,ia,iq1,iq2,mu
      double precision ubmom(0:3,nlegborn)
      double precision dummyamp2,dummymunu(0:3,0:3,nlegborn),
     $     dummyjk(nlegborn,nlegborn)
      double precision r,random,sigtot,siguptoic(ncstructmax)
      external random

      logical debug
      parameter (debug=.false.)

c     blank all colours (not only neutral particles)
      do ileg=1,nup
         icolup(1,ileg)=0
         icolup(2,ileg)=0
      enddo

      if(debug) then
         write(*,*) '***************************'
         write(*,*) 'TESTING borncolour_lh'
         if(nup.ne.6) then
            write(*,*) 'nup is not 6 in borncolour_lh ',nup
         endif
      endif

      do ileg=1,nup
         ubflav(ileg)=idup(ileg)
         if(ubflav(ileg).eq.21) ubflav(ileg)=0
c     use kn_cmpborn to compute planars, as in the dijet case
         do mu=0,3
            ubmom(mu,ileg)=kn_cmpborn(mu,ileg)
         enddo
      enddo

      if(debug) then
         write(*,*) 'UBFLAV: ',ubflav
      endif

      call born_ampsq_g_last(ubmom,ubflav,dummyamp2,dummymunu,dummyjk)

c     now needed informations needed in common/planar/
      sigtot=0d0
      do ic=1,ncstruct
         sigtot=sigtot+sigmaofncstruct(ic)
         siguptoic(ic)=sigtot
         if(debug) write(*,*) 'sigma_ic= ',sigmaofncstruct(ic)
      enddo

      r=random()*sigtot
      do ic=1,ncstruct
         if(r.lt.siguptoic(ic)) then
            ctag=ic
            goto 123
         endif
      enddo

 123  continue

      if(debug) write(*,*) 'chosen ic = ',ctag



      if((idup(1).eq.21).or.(idup(2).eq.21).or.
     1   (idup(5).eq.21).or.(idup(6).eq.21)) then
         proc='qqgg'
      elseif((abs(idup(1)).eq.abs(idup(2))).and.
     $       (abs(idup(1)).eq.abs(idup(3))).and.
     $       (abs(idup(1)).eq.abs(idup(4)))) then
         proc='qqqq'
      else
         proc='qqQQ'
      endif

      if(debug) write(*,*) proc


c q qb g g or permutations-crossing
      if(proc.eq.'qqgg') then
c     find the quarks and gluons
         ig1=-1
         iq1=-1
         do ileg=1,nlegborn
            if(ileg.ne.3.and.ileg.ne.4) then ! exclude leptons
               if(idup(ileg).eq.21) then
                  if(ig1.lt.0) then
                     ig1=ileg
                  else
                     ig2=ileg
                  endif
               else
                  if(iq1.lt.0) then
                     iq1=ileg
                  else
                     iq2=ileg
                  endif
               endif
            endif
         enddo


         if((clineofpart(ctag,iq1).eq.clineofpart(ctag,ig2)).and.
     $      (clineofpart(ctag,iq2).eq.clineofpart(ctag,ig1))) then
            if(debug) print*, 'CC are ',iq1,ig2
            if(idup(iq1).gt.0) then
c     q ---->---------- g
c             |---<----  
               icolup(1,iq1)=col(clineofpart(ctag,iq1))
               icolup(1,ig2)=col(clineofpart(ctag,ig2))
               icolup(2,ig2)=666
            else
c     qbar ----<------- g
c             |--- >---- 
               icolup(2,iq1)=col(clineofpart(ctag,iq1))
               icolup(2,ig2)=col(clineofpart(ctag,ig2))
               icolup(1,ig2)=666
            endif
            if(idup(iq2).gt.0) then
c     q ---->---------- g
c             |---<----  
               icolup(1,iq2)=col(clineofpart(ctag,iq2))
               icolup(1,ig1)=col(clineofpart(ctag,ig1))
               icolup(2,ig1)=666
            else
c     qbar ----<------- g
c             |--- >---- 
               icolup(2,iq2)=col(clineofpart(ctag,iq2))
               icolup(2,ig1)=col(clineofpart(ctag,ig1))
               icolup(1,ig1)=666
            endif

c     Until here we assumed that the Q-G color linked system was IF or FI.
c     Therefore the color (or anticolor) of the quark is the
c     color (or anticolor) of the gluon (my clineofpart vectors don't
c     know about color/anticolor, but only about who is linked to who!).
c     If the Q-G are II or FF, I have to exchange the color and anticolor of the
c     gluons, i.e. do this change (here II, but the same for FF):
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     WHAT I HAVE      | WHAT I MUST HAVE EVENTUALLY
c     Q c1 --->---     |       Q c1 --->---
c                      |                  
c                      |                  
c                      |                  
c     G c1 --->---     |       G c1 ---<---
c       c2 ---<---     |         c2 --->---
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if((istup(iq1)*istup(ig2).gt.0).and.
     $           (istup(iq2)*istup(ig1).gt.0)) then
               call colour_conj(icolup(1,ig1))
               call colour_conj(icolup(1,ig2))
               if(debug) write(*,*) 'color conj trick 1'
            endif

         elseif((clineofpart(ctag,iq1).eq.clineofpart(ctag,ig1)).and.
     $      (clineofpart(ctag,iq2).eq.clineofpart(ctag,ig2))) then
            if(debug) print*, 'CC are ',iq1,ig1
            if(idup(iq1).gt.0) then
c     q ---->---------- g
c             |---<----  
               icolup(1,iq1)=col(clineofpart(ctag,iq1))
               icolup(1,ig1)=col(clineofpart(ctag,ig1))
               icolup(2,ig1)=666
            else
c     qbar ----<------- g
c             |--- >---- 
               icolup(2,iq1)=col(clineofpart(ctag,iq1))
               icolup(2,ig1)=col(clineofpart(ctag,ig1))
               icolup(1,ig1)=666
            endif
            if(idup(iq2).gt.0) then
c     q ---->---------- g
c             |---<----  
               icolup(1,iq2)=col(clineofpart(ctag,iq2))
               icolup(1,ig2)=col(clineofpart(ctag,ig2))
               icolup(2,ig2)=666
            else
c     qbar ----<------- g
c             |--- >---- 
               icolup(2,iq2)=col(clineofpart(ctag,iq2))
               icolup(2,ig2)=col(clineofpart(ctag,ig2))
               icolup(1,ig2)=666
            endif

            if((istup(iq1)*istup(ig1).gt.0).and.
     $           (istup(iq2)*istup(ig2).gt.0)) then
               call colour_conj(icolup(1,ig1))
               call colour_conj(icolup(1,ig2))
               if(debug) write(*,*) 'color conj trick 2'
            endif
            
         else
            print*, 'problem in gg planar connections'
            stop
         endif

      elseif(proc.eq.'qqqq') then
c q q q q or permutations-crossing
         do ileg=1,nlegborn
            if(ileg.ne.3.and.ileg.ne.4) then ! exclude leptons         
               icolup(1,ileg)=col(clineofpart(ctag,ileg))
               icolup(2,ileg)=0
               if(idup(ileg).lt.0) call colour_conj(icolup(1,ileg))
            endif
         enddo

      elseif(proc.eq.'qqQQ') then
c q qp q qp or permutations-crossing
         do ileg=1,nlegborn
            if(ileg.ne.3.and.ileg.ne.4) then ! exclude leptons         
               icolup(1,ileg)=col(clineofpart(ctag,ileg))
               icolup(2,ileg)=0
               if(idup(ileg).lt.0) call colour_conj(icolup(1,ileg))
            endif
         enddo
      else
         print*, 'not enter in any flav struct while doing planars'
         stop
      endif

c     Probably this is not needed 
c     (these occurences were already taken care of in the above code).
      do ileg=1,nlegborn
         if(ileg.ne.3.and.ileg.ne.4) then ! exclude leptons         
            if((idup(ileg).lt.0).and.(icolup(2,ileg).eq.0)) then
               call colour_conj(icolup(1,ileg))
            endif
         endif
      enddo

      end



      subroutine colourjoin2g(icol1,icol2,icol3,icol4)
c                             q     qbar  g     g
c perform a planar colour connection on the planar sequence
c q qbar g g: it does the following
c     c1 --->---||--->---||--->--- c2 (antiquark)
c               ||       ||     
c               ||       ||     
c               v^       v^     
c               ||       ||     
c               ||       ||
c               c4       c3

      implicit none
      integer icol1(2),icol2(2),icol3(2),icol4(2)
      integer newcolor
      icol1(2)=0
      icol2(1)=0
      call getnewcolor(newcolor)
      icol1(1)=newcolor
      icol4(2)=newcolor
      call getnewcolor(newcolor)
      icol4(1)=newcolor
      icol3(2)=newcolor
      call getnewcolor(newcolor)
      icol3(1)=newcolor
      icol2(2)=newcolor
      end



      subroutine finalize_lh
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface.
c     Resonance Z -> e-(3) e+(4)
      call add_resonance(23,3,4)

c     The general reshuffling procedure.
      call lhefinitemasses

      end




c$$$c     if I have 2 gluons, then there are 3 color lines,
c$$$c     and I tagged only the lines that correspond to the qq
c$$$c     pair (c1,c2 in the following)
c$$$c     c1 --->---||--->---||--->--- c2
c$$$c               ||       ||     
c$$$c               ||       ||     
c$$$c               v^       v^     
c$$$c               ||       ||     
c$$$c               ||       ||
c$$$c             c1  c3   c3  c2
c$$$
c$$$c     check colourjoin2g in the jj code.
c$$$c     as it is now it connects the 1st argument with the 4th
c$$$c     (i.e. the 1st quark with the 2nd gluon)
c$$$         if(clineofpart(ctag,iq1).eq.clineofpart(ctag,ig2)) then
c$$$            print*, 'cc ',iq1,ig2
c$$$            call colourjoin2g(icolup(1,iq1),icolup(1,iq2),
c$$$     $           icolup(1,ig1),icolup(1,ig2))
c$$$         elseif(clineofpart(ctag,iq1).eq.clineofpart(ctag,ig1)) then
c$$$            print*, 'cc ',iq1,ig1
c$$$            call colourjoin2g(icolup(1,iq1),icolup(1,iq2),
c$$$     $           icolup(1,ig2),icolup(1,ig1))
