      subroutine setborn(p,bflav,born,bornjk,bmunu)
      use consts_MCFM; use dpinitialization 
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      include 'cvecbos.h'
      integer nlegs,nf
      parameter (nlegs=nlegborn)
      parameter (nf=5)
      real * 8 p(0:3,nlegs),p0(0:3,nlegs),bornjk(nlegs,nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),born
      integer bflav0(nlegs),bflav(nlegs)
!      integer i,j
      if(idvecbos.eq.24) then
         bflav0=bflav
         p0=p
      else
c Apply CP to the kinematics
         bflav0=-bflav
         p0=p
         p0(1,:)=-p(1,:)
      endif
      call compborn(p0,bflav0,born,bmunu,bornjk)
C No bmunu to worry about here
C      if(idvecbos.ne.24) then
Cc     Apply P also to bmunu tensor
C         bmunu(1,:,:) = -bmunu(1,:,:) 
C         bmunu(:,1,:) = -bmunu(:,1,:) 
C      endif

!      call qqb_wpwp_qqb_colborn(bflav,p,bornjk2)

c     Colour factors for colour-correlated Born amplitudes;
C     Rule from 2.98 in FNO2007, leads to \sum_i B_ij=Cj * B, where i/=j
C      write(*,*) 'born', born  
C      do j=1,nlegborn
C         write(*,*) 'sum (eq. 2.98)', j, sum(bornjk(:,j))/born, 
C     .        sum(bornjk2(:,j))/born
C      enddo

      end


      subroutine compborn(pin,bflav,born,bmunu,bornjk)
      use consts_MCFM; use dpinitialization; use define_ampl  
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h' ! for alphas 
      include 'PhysPars.h'
      integer nlegs,nf
      parameter (nlegs=nlegborn)
      parameter (nf=5)
      real * 8 pin(0:3,nlegs)
      real * 8 p(12,1:4)  ! 12 = mxpart in MCFM 
      real * 8 p1(12,1:4) ! 12 = mxpart in MCFM 
      real * 8 bmunu(0:3,0:3),born, bornjk(nlegs,nlegs)
      real * 8 msqB_str(-nf:nf,-nf:nf,3)
      real * 8 bornjk_cp(nlegs,nlegs)
      double precision msq(-5:5,-5:5)
      double complex virt_q_bq(3), virt_bq_q(3) 
      double complex virt_q_q(3), virt_bq_bq(3) 
      double complex xl(nlegborn,nlegborn)
      integer i,j,k,ud, cs, icol1, icol2   
      integer bflav(nlegs)
      logical :: identical, iswap  
      logical, save :: firsttime = .true. 
      character(len=3) :: chn 
      include 'cvecbos.h'

      if (firsttime) then 
         call dpinitialize(6) 

c-----setting MCFMconstants
	 MCFMwmass  = ph_wmass
	 MCFMwwidth = ph_wwidth
	 MCFMzmass  = ph_zmass
	 MCFMzwidth = ph_zwidth
	 MCFMgw     = ph_unit_e/ph_sthw
         firsttime = .false. 
      endif

      Npoint = 6 
      MCFMgsq    = st_alpha*4d0*pi 

      do i=1,nlegborn
         p(i,4)   = pin(0,i) 
         p(i,1:3) = pin(1:3,i) 
      enddo
      p(1,:) = -p(1,:) 
      p(2,:) = -p(2,:) 

C     get incoming channel 
      if (bflav(1) > 0 .and. bflav(2) < 0) then 
         chn = 'qqb' 
      elseif (bflav(1) < 0 .and. bflav(2) > 0) then 
         chn = 'qbq' 
      elseif (bflav(1) > 0 .and. bflav(2) > 0) then 
         chn = 'qqq' 
      elseif (bflav(1) < 0 .and. bflav(2) < 0) then 
         chn = 'qbb' 
      else
         write(*,*) 'bflav', bflav
         stop 'setborn: incoming flavour not OK' 
      endif

C     count number of families involved (whether quark
C     pairs are identical) 
      ud = 0;  cs = 0 
      do i = 1,8 
         if (i < 3) then ! skip leptons in counting 
            if (bflav(i) == -1 .or. bflav(i) == 2) then 
               ud = ud +1
            elseif (bflav(i) == -3 .or. bflav(i) == 4) then 
               cs =cs +1
            endif
         elseif (i > 6) then 
            if (bflav(i) == 1 .or. bflav(i) == -2) then 
               ud = ud +1
            elseif (bflav(i) == 3 .or. bflav(i) == -4) then 
               cs =cs +1
            endif
         endif
      enddo
      if (ud == 2 .and. cs == 2) then 
         identical = .false. 
      elseif (ud == 4 .or. cs == 4) then 
         identical = .true. 
      else
         write(*,*) 'ud, cs', ud, cs 
         stop 'setborn: number of quark pairs not ok'
      endif
      p1 = p 

C     now make sure that routines are called with momenta
C     in the right order 
      iswap = .false. 
      if (chn == 'qbb') then 
         if (bflav(1) - 1 .ne. bflav(7)) then 
            p1(8,:) = p(7,:) 
            p1(7,:) = p(8,:) 
            iswap = .true. 
         endif
      elseif (chn == 'qqq') then 
         if (bflav(1) - 1 .ne. bflav(7)) then 
            p1(8,:) = p(7,:) 
            p1(7,:) = p(8,:) 
            iswap = .true. 
         endif
      endif

      call qqb_wpwp_qqb(p1,msq,chn,identical) 
      call qqb_wpwp_qqb_str(p1,msqB_str,chn,identical) 
      msq = msq * vsymfact 
      msqB_str = msqB_str * vsymfact 
      born = msq(bflav(1),bflav(2))

C     -- no gluons, so no spin correlated Born  
      do i=0,3
         do j=0,3
            bmunu(i,j)=0d0
         enddo
      enddo

      bornjk = 0d0 
C     -- compute bornjk 
      do icol1=1,nlegborn
         do icol2=icol1+1,nlegborn
            xl = 0d0
C     factor 1/2 has the following origin: (2.92) of 0709.2092
C     sums over all i,j (unordered), the pole routine had only ordered pairs
C     so setting xl = 1d0 one would get e.g. B'_12 (for B'_21 = 0), we want the
C     symmetric B_12 = B_21 = B'12/2d0, so need 1/2 here 
C     now e.g. (2.98) is satisfied 
            xl(icol1,icol2) = 1d0/2d0
            
            j = bflav(1) 
            k = bflav(2) 
C     -- chn qqb 
            if (chn == 'qqb') then 
               virt_q_bq(1) = (
     .              ((-4d0*Cf+2d0*Nc)*(xl(1,7)+xl(2,8))
     .              +        ( 4d0*Cf-1d0*Nc)*(xl(1,2)+xl(7,8)) 
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,8)+xl(2,7))
     .              ))
               
               virt_q_bq(2) = (
     .              ((-4d0*Cf+2d0*Nc)*(xl(1,7)+xl(2,8))
     .              +        ( 4d0*Cf-1d0*Nc)*(xl(1,8)+xl(2,7))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,2)+xl(7,8)) 
     .              ))
               
               
               virt_q_bq(3) = (
     .              ((-2d0*Cf+2d0*Nc)*(xl(1,7)+xl(2,8))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,8)+xl(2,7))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,2)+xl(7,8)) 
     .              ))
               
               if ((j == 2 .and. k == -3) .or.
     .              (j == 4 .and. k == -1)) then 
                  bornjk(icol1,icol2) = msqB_str(j,k,1)*virt_q_bq(1) ! s
               else
                  bornjk(icol1,icol2) = 
     .                 + msqB_str(j,k,2)*virt_q_bq(1) ! t 
     .                 + msqB_str(j,k,1)*virt_q_bq(2) ! s
     .                 + msqB_str(j,k,3)*virt_q_bq(3) ! s*t 
               endif
               
               
               
C     -- chn qbq
            elseif (chn == 'qbq') then 
               virt_bq_q(1) = (
     .              ((-4d0*Cf+2d0*Nc)*(xl(2,7)+xl(1,8))
     .              +        ( 4d0*Cf-1d0*Nc)*(xl(1,2)+xl(7,8))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(2,8)+xl(1,7))
     .              ))
               
               virt_bq_q(2) = (
     .              ((-4d0*Cf+2d0*Nc)*(xl(2,7)+xl(1,8))
     .              +        ( 4d0*Cf-1d0*Nc)*(xl(2,8)+xl(1,7))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,2)+xl(7,8)) 
     .              ))
               
               
               virt_bq_q(3) = (
     .              ((-2d0*Cf+2d0*Nc)*(xl(2,7)+xl(1,8))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(2,8)+xl(1,7))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,2)+xl(7,8)) 
     .              ))
               
               if ((j == -3 .and. k == 2) .or.
     .              (j == -1 .and. k == 4)) then  
                  bornjk(icol1,icol2) = msqB_str(j,k,1)*virt_bq_q(1) ! s
               else
                  bornjk(icol1,icol2) = 
     .                 + msqB_str(j,k,2)*virt_bq_q(1)
     .                 + msqB_str(j,k,1)*virt_bq_q(2)
     .                 + msqB_str(j,k,3)*virt_bq_q(3)
               endif

               
C     -- chn qqq 
            elseif (chn == 'qqq') then 
               virt_q_q(1) = ( 
     .              ((-4d0*Cf+2d0*Nc)*(xl(1,2)+xl(7,8))
     .              +        ( 4d0*Cf-1d0*Nc)*(xl(1,8)+xl(2,7))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,7)+xl(2,8))
     .              ))
               
               virt_q_q(2) = ( 
     .              ((-4d0*Cf+2d0*Nc)*(xl(1,2)+xl(7,8))  
     .              +        ( 4d0*Cf-1d0*Nc)*(xl(1,7)+xl(2,8))  
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(2,7)+xl(1,8))  
     .              ))
               
               virt_q_q(3) = ( 
     .              ((-2d0*Cf+2d0*Nc)*(xl(1,2)+xl(7,8))  
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,8)+xl(2,7))  
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,7)+xl(2,8))  
     .              ))
               
               if (j == k) then 
                  bornjk(icol1,icol2) = 
     .                 + msqB_str(j,k,1)*virt_q_q(1)
     .                 + msqB_str(j,k,2)*virt_q_q(2)
     .                 + msqB_str(j,k,3)*virt_q_q(3)
               else
                  bornjk(icol1,icol2) = msqB_str(j,k,1)*virt_q_q(1)
               endif
               
               
            
C     -- chn qbb 
            elseif (chn == 'qbb') then 
               virt_bq_bq(1) = (
     .              ((-4d0*Cf+2d0*Nc)*(xl(7,8)+xl(1,2))
     .              +        ( 4d0*Cf-1d0*Nc)*(xl(2,7)+xl(1,8))
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,7)+xl(2,8))
     .              ))
               
               virt_bq_bq(2) = ( 
     .              ((-4d0*Cf+2d0*Nc)*(xl(1,2)+xl(7,8))  
     .              +        ( 4d0*Cf-1d0*Nc)*(xl(1,7)+xl(2,8))  
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(2,7)+xl(1,8))  
     .              ))
               
               virt_bq_bq(3) = ( 
     .              ((-2d0*Cf+2d0*Nc)*(xl(1,2)+xl(7,8))  
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,8)+xl(2,7))  
     .              +        ( 2d0*Cf-1d0*Nc)*(xl(1,7)+xl(2,8))  
     .              ))
               
               if (j == k) then 
                  bornjk(icol1,icol2) = 
     .                 + msqB_str(j,k,1)*virt_bq_bq(1)
     .                 + msqB_str(j,k,2)*virt_bq_bq(2)
     .                 + msqB_str(j,k,3)*virt_bq_bq(3)
               else
                  bornjk(icol1,icol2) = msqB_str(j,k,1)*virt_bq_bq(1)
               endif
                              
            else
               write(*,*) 'chn', chn
               stop 'qqb_wpw_qqb_colborn: undefined channel' 
            endif

            bornjk(icol2,icol1) = bornjk(icol1,icol2)
         enddo
         bornjk(icol1,icol1) = 0d0 
      enddo

      if(iswap) then
         bornjk_cp = bornjk
         bornjk(1,7) = bornjk_cp(1,8)
         bornjk(1,8) = bornjk_cp(1,7)
         bornjk(2,7) = bornjk_cp(2,8)
         bornjk(2,8) = bornjk_cp(2,7)
         bornjk(7,1) = bornjk(1,7)
         bornjk(8,1) = bornjk(1,8)
         bornjk(7,2) = bornjk(2,7)
         bornjk(8,2) = bornjk(2,8)
      endif

      end



      
      subroutine borncolour_lh
c     Sets up the colour for the given flavour configuration
c     already filled in the Les Houches interface.
c     In case there are several colour structure, one
c     should pick one with a probability proportional to
c     the value of the corresponding cross section, for the
c     kinematics defined in the Les Houches interface
      implicit none 
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      
      integer iq1,iq2,iq3,iq4,i
      integer idpar(8)
      character(len=3) :: chn 
      real *8 w1,w2,r,p(12,4), random  

C     -- neutral particles
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0
      icolup(1,5)=0
      icolup(2,5)=0
      icolup(1,6)=0
      icolup(2,6)=0

c     -- colored particles
      icolup(1,1)=0
      icolup(2,1)=0
      icolup(1,2)=0
      icolup(2,2)=0
      icolup(1,7)=0
      icolup(2,7)=0
      icolup(1,8)=0
      icolup(2,8)=0


C     id of all outgoing particles 
      idpar(1)   = -idup(1)  
      idpar(2)   = -idup(2)  
      idpar(3:8) =  idup(3:8)  
 
      
C     -- qq -> qq 
      if (idup(1) > 0 .and. idup(2) > 0) then
         iq1 = 1; iq2 = 2; iq3 = 7; iq4 = 8;
         chn = 'qqq' 
C     -- qbqb -> qbqb 
      elseif (idup(1) < 0 .and. idup(2) < 0) then
C         iq1 = 2; iq2 = 1; iq3 = 7; iq4 = 8;
         iq1 = 2; iq2 = 1; iq3 = 8; iq4 = 7;
         chn = 'qbb' 
C     -- qqb -> qqb 
      elseif (idup(1) > 0 .and. idup(2) < 0 .and. idup(7) > 0) then
         iq1 = 1; iq2 = 8; iq3 = 2; iq4 = 7;
         chn = 'qqb' 
C     -- qqb -> qbq 
      elseif (idup(1) > 0 .and. idup(2) < 0 .and. idup(8) > 0) then
         iq1 = 1; iq2 = 7; iq3 = 2; iq4 = 8;
         chn = 'qqb' 
C     -- qbq -> qbq 
      elseif (idup(1) < 0 .and. idup(2) > 0 .and. idup(8) > 0) then
         iq1 = 1; iq2 = 8; iq3 = 2; iq4 = 7;
         chn = 'qbq' 
C     -- qbq -> qqb 
      elseif (idup(1) < 0 .and. idup(2) > 0 .and. idup(8) < 0) then
         iq1 = 1; iq2 = 7; iq3 = 2; iq4 = 8;
         chn = 'qbq' 
      else
         write(*,*) 'IDUP', idup(1:8) 
         stop 'borncolor_lc: idup out of range'
      endif

C     distinct families, 1 diag, LC color flow fixed  
         if (idpar(iq1)+idpar(iq3)==-1 .and.
     .     idpar(iq1)/=idpar(iq2)) then
            if (idup(iq1) > 0) then 
               icolup(1,iq1) = 501
            else
               icolup(2,iq1) = 501
            endif               
            if (idup(iq4) > 0) then 
               icolup(1,iq4) = 501 
            else
               icolup(2,iq4) = 501 
            endif               
            if (idup(iq2) > 0) then 
               icolup(1,iq2) = 502 
            else
               icolup(2,iq2) = 502 
            endif
            if (idup(iq3) > 0) then 
               icolup(1,iq3) = 502 
            else
               icolup(2,iq3) = 502 
            endif
C     distinct families, 1 diag, LC color flow fixed  
         elseif (idpar(iq1)+idpar(iq4)==-1 .and.
     .           idpar(iq1)/=idpar(iq2)) then
            if (idup(iq1) > 0) then 
               icolup(1,iq1) = 501
            else
               icolup(2,iq1) = 501
            endif               
            if (idup(iq3) > 0) then 
               icolup(1,iq3) = 501 
            else
               icolup(2,iq3) = 501 
            endif               
            if (idup(iq2) > 0) then 
               icolup(1,iq2) = 502 
            else
               icolup(2,iq2) = 502 
            endif
            if (idup(iq4) > 0) then 
               icolup(1,iq4) = 502 
            else
               icolup(2,iq4) = 502 
            endif

C     identical case, two LC possible color flows  
         elseif (idpar(iq1)+idpar(iq3)==-1 .and.
     .           idpar(iq1)==idpar(iq2)) then

C     need to pick one of the two above possibilities based on |M|^2 at LC 
            p = 0d0 
            do i=1,nlegborn
               p(i,4)   = kn_cmpborn(0,i) 
               p(i,1:3) = kn_cmpborn(1:3,i) 
            enddo
            p(1,:) = -p(1,:) 
            p(2,:) = -p(2,:) 
            
            call qqb_wpwp_qqb_w1w2(p,chn,w1,w2)
            r = random() 
            if (r < w1/(w1+w2) ) then
            if (idup(iq1) > 0) then 
               icolup(1,iq1) = 501
            else
               icolup(2,iq1) = 501
            endif               
            if (idup(iq4) > 0) then 
               icolup(1,iq4) = 501 
            else
               icolup(2,iq4) = 501 
            endif               
            if (idup(iq2) > 0) then 
               icolup(1,iq2) = 502 
            else
               icolup(2,iq2) = 502 
            endif
            if (idup(iq3) > 0) then 
               icolup(1,iq3) = 502 
            else
               icolup(2,iq3) = 502 
            endif
            else
               
               if (idup(iq1) > 0) then 
               icolup(1,iq1) = 501
            else
               icolup(2,iq1) = 501
            endif               
            if (idup(iq3) > 0) then 
               icolup(1,iq3) = 501 
            else
               icolup(2,iq3) = 501 
            endif               
            if (idup(iq2) > 0) then 
               icolup(1,iq2) = 502 
            else
               icolup(2,iq2) = 502 
            endif
            if (idup(iq4) > 0) then 
               icolup(1,iq4) = 502 
            else
               icolup(2,iq4) = 502 
            endif
            endif
            
         else
            write(*,*) 'IQ', iq1,iq2,iq3,iq4 
            stop 'borncolour_lc: idup not recognised' 
         endif

      end

      subroutine finalize_lh
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface.
c     
c     vector boson id and decay
      include 'cvecbos.h'
c     lepton masses
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass

      call add_resonance(idvecbos,3,4)
C     need to shift (56) to (67) since previous res adds a label 
      call add_resonance(idvecbos,6,7)

C     c     The following routine also performs the reshuffling of momenta if
Cc     a massive decay is chosen
C      call momenta_reshuffle(3,4,5,decmass)
      end



c     i1<i2
      subroutine momenta_reshuffle(ires,i1,i2,decmass)
      implicit none
      include 'LesHouches.h'
      integer ires,i1,i2,j
      real * 8 ptemp(0:3),ptemp1(0:3),beta(3),betainv(3),modbeta,decmass
      if (i1.ge.i2) then
         write(*,*) 'wrong sequence in momenta_reshuffle'
         stop
      endif
cccccccccccccccccccccccccccccc
c construct boosts from/to vector boson rest frame 
      do j=1,3
         beta(j)=-pup(j,ires)/pup(4,ires)
      enddo
      modbeta=sqrt(beta(1)**2+beta(2)**2+beta(3)**2)
      do j=1,3
         beta(j)=beta(j)/modbeta
         betainv(j)=-beta(j)
      enddo
cccccccccccccccccccccccccccccccccccccccc
c first decay product (massive)
      ptemp(0)=pup(4,i1)
      do j=1,3
         ptemp(j)=pup(j,i1)
      enddo
      call mboost(1,beta,modbeta,ptemp,ptemp)
      ptemp1(0)=0.5d0*(pup(5,ires)+(decmass**2)/pup(5,ires))
      do j=1,3
         ptemp1(j)=ptemp(j)/ptemp(0)*sqrt(ptemp1(0)**2 -decmass**2)
      enddo
      call mboost(1,betainv,modbeta,ptemp1,ptemp)
      do j=1,3
         pup(j,i1)=ptemp(j)
      enddo
      pup(4,i1)=ptemp(0)
      pup(5,i1)=sqrt(pup(4,i1)**2-pup(1,i1)**2
     $     -pup(2,i1)**2-pup(3,i1)**2)
      
c second decay product (massless)

      ptemp(0)=pup(4,i2)
      do j=1,3
         ptemp(j)=pup(j,i2)
      enddo
      call mboost(1,beta,modbeta,ptemp,ptemp)
      ptemp1(0)=0.5d0*(pup(5,ires)-(decmass**2)/pup(5,ires))
      do j=1,3
         ptemp1(j)=ptemp(j)/ptemp(0)*ptemp1(0)
      enddo
      call mboost(1,betainv,modbeta,ptemp1,ptemp)
      do j=1,3
         pup(j,i2)=ptemp(j)
      enddo
      pup(4,i2)=ptemp(0)
c abs to avoid tiny negative values
      pup(5,i2)=sqrt(abs(pup(4,i2)**2-pup(1,i2)**2
     $     -pup(2,i2)**2-pup(3,i2)**2))
cccccccccccccccccccccccccccccccccccccccc
      end

 
