ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute leptonic tensors for VBF pp->W+W-jj 
c
      subroutine compute_tensors_wpm(pin)

c     fully leptonic decays
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)
c
c vbfnlo stuff:
      include 'global.inc'
      include 'coupl.inc'
      include 'tensor_born.inc'
      include 'higgs_graphs.h'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 q12(0:4,3),q34(0:4,3),qww(0:3),qq

      complex*16 zero
      parameter (zero=(0d0,0d0))

      double complex dotrc
      external dotrc
c
c declare local variables
      integer i,j,mu
      integer ttype
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ttype = 1 ! tensors are Born-type

      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,4
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,7)
         p(mu,4) = pin(mu,8) 
         p(mu,5) = 0d0   
      enddo ! mu

c compute invariants:
      call calc_invariants(p,v,q12,q34,1)
c
c reset leptonic tensors:
      call vtoww_born_reset
c compute leptonic tensors:
      CALL IXXXXX(v(0,2),ZERO ,1,-1,wep_born) !W(1,5))          !e+       
      CALL OXXXXX(v(0,1),ZERO ,-1,1,wve_born) !W(1,6))          !ve 
      CALL OXXXXX(v(0,4),ZERO ,-1,1,wmu_born) !W(1,7))          !mu-      
      CALL IXXXXX(v(0,3),ZERO ,1,-1,wvm_born) !W(1,8))          !vm~
      CALL JIOXXX(wep_born,wve_born,GWF ,WMASS,WWIDTH,wp_born) !W(1,10))!W+
      CALL JIOXXX(wvm_born,wmu_born,GWF ,WMASS,WWIDTH,wm_born) !W(1,11))!W-

      do mu = 0,3
         qp_born(mu) = v(mu,1)+v(mu,2)
         qm_born(mu) = v(mu,3)+v(mu,4)
         qww(mu) = qp_born(mu) + qm_born(mu)
      enddo
      qp_born(4)=qp_born(0)**2-qp_born(1)**2-qp_born(2)**2-qp_born(3)**2
      qm_born(4)=qm_born(0)**2-qm_born(1)**2-qm_born(2)**2-qm_born(3)**2

c determine "reduced" W+- polarization vectors 
c wpp(mu) = wp(mu) - xp*qp(mu) and wmp(mu) = wm(mu) - xm*qm(mu)
            qq = qp_born(0)*qww(0)-qp_born(1)*qww(1)-
     &           qp_born(2)*qww(2)-qp_born(3)*qww(3)
            xp_born = dotrc(qww,wp_born)/qq
            qq = qm_born(0)*qww(0)-qm_born(1)*qww(1)-
     &           qm_born(2)*qww(2)-qm_born(3)*qww(3)
            xm_born = dotrc(qww,wm_born)/qq
            do mu = 0,3
               wpp_born(mu+1) = wp_born(mu+1) - xp_born*qp_born(mu)
               wmp_born(mu+1) = wm_born(mu+1) - xm_born*qm_born(mu)
            enddo
            do mu = 5,6
               wpp_born(mu) = wp_born(mu)
               wmp_born(mu) = wm_born(mu)
            enddo

c leptonic tensors
C for W+W-
      if (.not.higgs_only) then
         call atoww(v,aww_born,ttype)
         call ztoww(v,zww_born,ttype)
      endif

      j = 1

      if (.not.higgs_only) then
         call aatoww(q12(0,j),q34(0,j),v,aaww_born(0,0,j),ttype)
         call aztoww(q12(0,j),q34(0,j),v,azww_born(0,0,j),ttype)
         call aztoww(q34(0,j),q12(0,j),v,zaww_born(0,0,j),ttype)
      endif
      call zztoww(q12(0,j),q34(0,j),v,zzww_born(0,0,j),ttype)
      call wwtoww(q12(0,j),q34(0,j),v,wwww6_born(0,0,j),ttype)  ! q12 = W-
      call wwtoww(q34(0,j),q12(0,j),v,wwww5_born(0,0,j),ttype)  ! q12 = W+

      if (.not.higgs_only) then
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
      call WVtoWP(2,q34(0,j),v,NCwpa_born(0,0,1,j),
     &        NCwpz_born(0,0,1,j),ttype) !emit W- on upper
      call WVtoWP(2,q12(0,j),v,NCwpa_born(0,0,2,j),
     &        NCwpz_born(0,0,2,j),ttype)!emit W- on lower
      call WVtoWP(1,q34(0,j),v,CCwpa_born(0,0,1,j),
     &        CCwpz_born(0,0,1,j),ttype)!emit W- on upper
      call WVtoWP(1,q12(0,j),v,CCwpa_born(0,0,2,j),
     &        CCwpz_born(0,0,2,j),ttype)!emit W- on lower
C for WV --> mu- nu_mu
      call WVtoWM(2,q34(0,j),v,NCwma_born(0,0,1,j),
     &        NCwmz_born(0,0,1,j),ttype)!emit W+ on upper
      call WVtoWM(2,q12(0,j),v,NCwma_born(0,0,2,j),
     &        NCwmz_born(0,0,2,j),ttype)!emit W+ on lower
      call WVtoWM(1,q34(0,j),v,CCwma_born(0,0,1,j),
     &        CCwmz_born(0,0,1,j),ttype)!emit W+ on upper
      call WVtoWM(1,q12(0,j),v,CCwma_born(0,0,2,j),
     &        CCwmz_born(0,0,2,j),ttype)!emit W+ on lower
         
      endif ! H->WW only

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute tensors for VBF pp->W+W-jj with semi-leptonic decay
c
      subroutine compute_tensors_wpm_slp(pin)
c
c     "slp": (e+ve) is replaced with (ud~) in final state
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)
c
c vbfnlo stuff:
      include 'global.inc'
      include 'coupl.inc'
      include 'tensor_born.inc'
      include 'higgs_graphs.h'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 q12(0:4,3),q34(0:4,3),qww(0:3),qq

      complex*16 zero
      parameter (zero=(0d0,0d0))

      double complex dotrc
      external dotrc
c
c declare local variables
      integer i,j,mu
      integer ttype
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ttype = 1 ! tensors are Born-type
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,4
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,7)
         p(mu,4) = pin(mu,8) 
         p(mu,5) = 0d0   
      enddo ! mu

c compute invariants:
      call calc_invariants(p,v,q12,q34,1)
c
c reset (semi-)leptonic tensors:
      call vtoww_born_reset
c compute (semi-)leptonic tensors:
      CALL IXXXXX(v(0,2),ZERO ,1,-1,wve_born)  !b~               
      CALL OXXXXX(v(0,1),ZERO ,-1,1,wep_born)  !t        
      CALL OXXXXX(v(0,4),ZERO ,-1,1,wmu_born)  !mu-      
      CALL IXXXXX(v(0,3),ZERO ,1,-1,wvm_born)  !vm~

      CALL JIOXXX(wve_born,wep_born,GWF ,WMASS,WWIDTH,wp_born) !W+
      CALL JIOXXX(wvm_born,wmu_born,GWF ,WMASS,WWIDTH,wm_born) !W-

      do mu = 0,3
         qp_born(mu) = v(mu,1)+v(mu,2)
         qm_born(mu) = v(mu,3)+v(mu,4)
         qww(mu) = qp_born(mu) + qm_born(mu)
      enddo
      qp_born(4)=qp_born(0)**2-qp_born(1)**2-
     &           qp_born(2)**2-qp_born(3)**2
      qm_born(4)=qm_born(0)**2-qm_born(1)**2-
     &           qm_born(2)**2-qm_born(3)**2

c determine "reduced" W+- polarization vectors 
c wpp(mu) = wp(mu) - xp*qp(mu) and wmp(mu) = wm(mu) - xm*qm(mu)
      qq = qp_born(0)*qww(0)-qp_born(1)*qww(1)-
     &     qp_born(2)*qww(2)-qp_born(3)*qww(3)
      xp_born = dotrc(qww,wp_born)/qq
      qq = qm_born(0)*qww(0)-qm_born(1)*qww(1)-
     &     qm_born(2)*qww(2)-qm_born(3)*qww(3)
      xm_born = dotrc(qww,wm_born)/qq
      do mu = 0,3
         wpp_born(mu+1) = wp_born(mu+1) - xp_born*qp_born(mu)
         wmp_born(mu+1) = wm_born(mu+1) - xm_born*qm_born(mu)
      enddo
      do mu = 5,6
         wpp_born(mu) = wp_born(mu)
         wmp_born(mu) = wm_born(mu)
      enddo

c (semi-)leptonic tensors
C for W+W-
        if (.not.higgs_only) then
           call atoww_slp(v,aww_born,ttype)
           call ztoww_slp(v,zww_born,ttype)
        endif
        j = 1

        if (.not.higgs_only) then
           call aatoww_slp(q12(0,j),q34(0,j),v,aaww_born(0,0,j),ttype)
           call aztoww_slp(q12(0,j),q34(0,j),v,azww_born(0,0,j),ttype)
           call aztoww_slp(q34(0,j),q12(0,j),v,zaww_born(0,0,j),ttype)
        endif
        call zztoww_slp(q12(0,j),q34(0,j),v,zzww_born(0,0,j),ttype)
        call wwtoww_slp(q12(0,j),q34(0,j),v,wwww6_born(0,0,j),ttype)! q12=W-
        call wwtoww_slp(q34(0,j),q12(0,j),v,wwww5_born(0,0,j),ttype)! q12=W+

        if (.not.higgs_only) then
cC for WV --> tb
cc NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
        call WVtoWP_slp(2,q34(0,j),v,NCwpa_born(0,0,1,j),
     &                  NCwpz_born(0,0,1,j),ttype) !emit W- on upper
        call WVtoWP_slp(2,q12(0,j),v,NCwpa_born(0,0,2,j),
     &                  NCwpz_born(0,0,2,j),ttype)     !emit W- on lower
        call WVtoWP_slp(1,q34(0,j),v,CCwpa_born(0,0,1,j),
     &                  CCwpz_born(0,0,1,j),ttype)     !emit W- on upper
        call WVtoWP_slp(1,q12(0,j),v,CCwpa_born(0,0,2,j),
     &                  CCwpz_born(0,0,2,j),ttype)     !emit W- on lower
C for WV --> mu- nu_mu
        call WVtoWM(2,q34(0,j),v,NCwma_born(0,0,1,j),
     &        NCwmz_born(0,0,1,j),ttype)     !emit W+ on upper
        call WVtoWM(2,q12(0,j),v,NCwma_born(0,0,2,j),
     &        NCwmz_born(0,0,2,j),ttype)     !emit W+ on lower
        call WVtoWM(1,q34(0,j),v,CCwma_born(0,0,1,j),
     &        CCwmz_born(0,0,1,j),ttype)     !emit W+ on upper
        call WVtoWM(1,q12(0,j),v,CCwma_born(0,0,2,j),
     &        CCwmz_born(0,0,2,j),ttype)     !emit W+ on lower
         
        endif ! H->WW only

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute tensors for VBF pp->W+W-jj with semi-leptonic decay
c
      subroutine compute_tensors_wpm_slm(pin)
c
c     "slm": (vm~ mu) is replaced with (u~ d) in final state
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)
c
c vbfnlo stuff:
      include 'global.inc'
      include 'coupl.inc'
      include 'tensor_born.inc'
      include 'higgs_graphs.h'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 q12(0:4,3),q34(0:4,3),qww(0:3),qq

      complex*16 zero
      parameter (zero=(0d0,0d0))

      double complex dotrc
      external dotrc
c
c declare local variables
      integer i,j,mu
      integer ttype
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ttype = 1 ! tensors are Born-type
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,4
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,7)
         p(mu,4) = pin(mu,8) 
         p(mu,5) = 0d0   
      enddo ! mu

c compute invariants:
      call calc_invariants(p,v,q12,q34,1)
c
c reset (semi-)leptonic tensors:
      call vtoww_born_reset
c compute (semi-)leptonic tensors:
      CALL IXXXXX(v(0,2),ZERO ,1,-1,wep_born)        !e+       
      CALL OXXXXX(v(0,1),ZERO ,-1,1,wve_born)        !ve 
      CALL OXXXXX(v(0,4),ZERO ,-1,1,wmu_born)        !b    
      CALL IXXXXX(v(0,3),ZERO ,1,-1,wvm_born)        !t~
      CALL JIOXXX(wep_born,wve_born,GWF ,WMASS,WWIDTH,wp_born) !W+
      CALL JIOXXX(wvm_born,wmu_born,GWF ,WMASS,WWIDTH,wm_born) !W-

      do mu = 0,3
         qp_born(mu) = v(mu,1)+v(mu,2)
         qm_born(mu) = v(mu,3)+v(mu,4)
         qww(mu) = qp_born(mu) + qm_born(mu)
      enddo
      qp_born(4)=qp_born(0)**2-qp_born(1)**2-
     &           qp_born(2)**2-qp_born(3)**2
      qm_born(4)=qm_born(0)**2-qm_born(1)**2-
     &           qm_born(2)**2-qm_born(3)**2

c determine "reduced" W+- polarization vectors 
c wpp(mu) = wp(mu) - xp*qp(mu) and wmp(mu) = wm(mu) - xm*qm(mu)
      qq = qp_born(0)*qww(0)-qp_born(1)*qww(1)-
     &     qp_born(2)*qww(2)-qp_born(3)*qww(3)
      xp_born = dotrc(qww,wp_born)/qq
      qq = qm_born(0)*qww(0)-qm_born(1)*qww(1)-
     &     qm_born(2)*qww(2)-qm_born(3)*qww(3)
      xm_born = dotrc(qww,wm_born)/qq
      do mu = 0,3
         wpp_born(mu+1) = wp_born(mu+1) - xp_born*qp_born(mu)
         wmp_born(mu+1) = wm_born(mu+1) - xm_born*qm_born(mu)
      enddo
      do mu = 5,6
         wpp_born(mu) = wp_born(mu)
         wmp_born(mu) = wm_born(mu)
      enddo

c (semi-)leptonic tensors
C for W+W-
         if (.not.higgs_only) then
            call atoww_slm(v,aww_born,ttype)
            call ztoww_slm(v,zww_born,ttype)
         endif

         j = 1
         if (.not.higgs_only) then
            call aatoww_slm(q12(0,j),q34(0,j),v,aaww_born(0,0,j),ttype)
            call aztoww_slm(q12(0,j),q34(0,j),v,azww_born(0,0,j),ttype)
            call aztoww_slm(q34(0,j),q12(0,j),v,zaww_born(0,0,j),ttype)
         endif
         call zztoww_slm(q12(0,j),q34(0,j),v,zzww_born(0,0,j),ttype)
         call wwtoww_slm(q12(0,j),q34(0,j),v,wwww6_born(0,0,j),ttype) !q12=W-
         call wwtoww_slm(q34(0,j),q12(0,j),v,wwww5_born(0,0,j),ttype) !q12=W+

         if (.not.higgs_only) then
cC for WV --> tb
cc NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
         call WVtoWP(2,q34(0,j),v,NCwpa_born(0,0,1,j),
     &               NCwpz_born(0,0,1,j),ttype)     !emit W- on upper
         call WVtoWP(2,q12(0,j),v,NCwpa_born(0,0,2,j),
     &               NCwpz_born(0,0,2,j),ttype)     !emit W- on lower
         call WVtoWP(1,q34(0,j),v,CCwpa_born(0,0,1,j),
     &               CCwpz_born(0,0,1,j),ttype)     !emit W- on upper
         call WVtoWP(1,q12(0,j),v,CCwpa_born(0,0,2,j),
     &               CCwpz_born(0,0,2,j),ttype)     !emit W- on lower
C for WV --> mu- nu_mu
         call WVtoWM_slm(2,q34(0,j),v,NCwma_born(0,0,1,j),
     &               NCwmz_born(0,0,1,j),ttype)     !emit W+ on upper
         call WVtoWM_slm(2,q12(0,j),v,NCwma_born(0,0,2,j),
     &               NCwmz_born(0,0,2,j),ttype)     !emit W+ on lower
         call WVtoWM_slm(1,q34(0,j),v,CCwma_born(0,0,1,j),
     &               CCwmz_born(0,0,1,j),ttype)     !emit W+ on upper
         call WVtoWM_slm(1,q12(0,j),v,CCwma_born(0,0,2,j),
     &               CCwmz_born(0,0,2,j),ttype)     !emit W+ on lower
         
         endif ! H->WW only

      return
      end
