      subroutine qqb_z2jet_v_pwhg(pin,vflav,res)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     July 2012.                                                       *
*                                                                      *
*     Calculate the virtual matrix element squared and                 *
*     subtraction terms for the process                                *
*                                                                      *
*     q(-p1) + qbar(-p2) --> Z + j(p5) + j(p6)                         *
*                            |                                         *
*                            --> e^-(p3) + e^+(p4)                     *
*                                                                      *
*     where the partons are either q(p5) and qbar(p6) [Qflag = .true.] *
*                               or g(p5) and g(p6)    [Gflag = .true.] *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'scale.f'
      include 'flags.f'
      include 'msq_cs.f'
      double precision res,res0
      double precision, save :: msqv(-nf:nf,-nf:nf),
     . mmsq_qqb(2,2),mmsq_qbq(2,2),mmsq_gq(2,2),mmsq_qg(2,2),
     . mmsq_qbg(2,2),mmsq_gqb(2,2),mmsq_gg(2,2)
      double precision
     . p(mxpart,4),pswap(mxpart,4),fac,pin(mxpart,4)
      double complex, save ::  mmsq_qqb_ax(2,2),mmsq_qbq_ax(2,2),
     . mmsq_gq_ax(2,2),mmsq_qg_ax(2,2),mmsq_qbg_ax(2,2),
     . mmsq_gqb_ax(2,2),mmsq_gg_ax(2,2),
     . mmsq_qqb_vec(2,2),mmsq_qbq_vec(2,2),
     . mmsq_gq_vec(2,2),mmsq_qg_vec(2,2),mmsq_qbg_vec(2,2),
     . mmsq_gqb_vec(2,2),mmsq_gg_vec(2,2)
      double complex atreez,a61z,a62z,a63z,prop,vcouple(2),
     . tamp,lamp,tampup,lampup,tampdo,lampdo,tamps,lamps,lampx,lampsx
      integer nu,j,k,cs,polq,polb,polz
      integer jj(-nf:nf),nup,ndo
      double precision subuv
      double precision faclo,v2(2),vQ(nf,2),idfac
      double precision mqq(0:2,fn:nf,fn:nf)
      double complex, save ::
     .               atreez_526143(2,2,2),atreez_251643(2,2,2),
     .               atreez_625143(2,2,2),atreez_261543(2,2,2),
     .               atreez_162543(2,2,2),atreez_615243(2,2,2),
     .               atreez_152643(2,2,2),atreez_516243(2,2,2),
     .               atreez_562143(2,2,2),atreez_651243(2,2,2),
     .               atreez_265143(2,2,2),atreez_621543(2,2,2),
     .               atreez_561243(2,2,2),atreez_652143(2,2,2),
     .               atreez_165243(2,2,2),atreez_612543(2,2,2),
     .               a61z_526143(2,2,2),a61z_251643(2,2,2),
     .               a61z_625143(2,2,2),a61z_261543(2,2,2),
     .               a61z_162543(2,2,2),a61z_615243(2,2,2),
     .               a61z_152643(2,2,2),a61z_516243(2,2,2),
     .               a61z_562143(2,2,2),a61z_651243(2,2,2),
     .               a61z_265143(2,2,2),a61z_621543(2,2,2),
     .               a61z_561243(2,2,2),a61z_652143(2,2,2),
     .               a61z_165243(2,2,2),a61z_612543(2,2,2),
     .               a62z_625143(2,2,2),a62z_261543(2,2,2),
     .               a62z_526143(2,2,2),a62z_251643(2,2,2),
     .               a62z_152643(2,2,2),a62z_516243(2,2,2),
     .               a62z_162543(2,2,2),a62z_615243(2,2,2),
     .               a62z_265143(2,2,2),a62z_621543(2,2,2),
     .               a62z_562143(2,2,2),a62z_651243(2,2,2),
     .               a62z_165243(2,2,2),a62z_612543(2,2,2),
     .               a62z_561243(2,2,2),a62z_652143(2,2,2),
     .               a63z_526143(2,2,2),a63z_625143(2,2,2),
     .               a63z_162543(2,2,2),a63z_152643(2,2,2),
     .               a63z_562143(2,2,2),a63z_265143(2,2,2),
     .               a63z_561243(2,2,2),a63z_165243(2,2,2)

      logical compare,checkvector,checkaxial
      common/mqq/mqq

      parameter (compare=.false.)
c--- These parameters allow for a point-by-point comparison of the
c--- vector and axial pieces with MadLoop. They should normally
c--- both be set to .false.
      parameter(checkvector=.false.,checkaxial=.false.)
!      data scheme/'dred'/
      data jj/-1,-2,-1,-2,-1,0,1,2,1,2,1/
!      save first
      integer vflav(6),nq,i
      logical same 
      real *8 born
C     variable needed to avoid recalculating same stuff 
      logical ::  recalc_g, recalc_q 
      real * 8, save :: opin_g(mxpart,4),scale_g 
      real * 8, save :: opin_q(mxpart,4),scale_q 

!     Set up momenta 
      p = pin 

C-----calculate born
      call qqb_z2jet_pwhg(p,vflav,res0)

      ! now understand if it's a two-quark or four-quark process 
      nq = 0 
      do i=1,6
         if (i < 3 .or. i > 4) then 
            if (abs(vflav(i)) .gt. 0) then 
               nq = nq+1
            endif
         endif
      enddo
      if (nq == 2) then 
         Gflag = .true. 
         Qflag = .false. 
      elseif (nq == 4) then 
         Qflag = .true. 
         Gflag = .false. 
      else 
         write(*,*) 'nq',nq 
         write(*,*) 'vflav',vflav
         stop 'nq out of range' 
      endif

!     decide if recalc
      recalc_q = .false. 
      recalc_g = .false. 
      

      if (Gflag) then
         if (scale_g .ne. scale) then
            recalc_g = .true.
         else
            do i=1,6
               do nu=1,4
                  if(opin_g(i,nu).ne.p(i,nu)) then
                     recalc_g = .true.
                     goto 10
                  endif
               enddo
            enddo
         endif
      else 
         if (scale_q .ne. scale) then
            recalc_q = .true.
         else
            do i=1,6
               do nu=1,4
                  if(opin_q(i,nu).ne.p(i,nu)) then
                     recalc_q = .true.
                     goto 10
                  endif
               enddo
            enddo
         endif
      endif

 10    continue 
       
       if (recalc_g) then
          opin_g = p
          scale_g = scale
       endif
       if (recalc_q) then
          opin_q = p
          scale_q = scale
       endif  
      


      same = .false. 
C     Gflag case 
      if (vflav(1) == 0 .and. vflav(2) == 0) then 
         if (abs(vflav(5)) == 2 .or. abs(vflav(5)) == 4) then 
            nup = 1
            ndo = 0 
         else
            nup = 0
            ndo = 1
         endif

C     Qflav 
C     qqb -> qqb or qbq -> qbq  
      elseif ((vflav(1).eq.-vflav(2)).and.(vflav(5).eq.-vflav(6))) then 

         if (abs(vflav(1)) .eq. abs(vflav(5))) then 
            same = .true. 
         else
            same = .false. 
         endif

         if ((vflav(5)/2)*2 == vflav(5)) then 
            nup = 1
            ndo = 0 
         else
            nup = 0
            ndo = 1
         endif
      endif


      if (Qflag .and. Gflag) then
        write(6,*) 'Both Qflag and Gflag cannot be true'
        write(6,*) 'They are set in file options.DAT'
        write(6,*) 'Failed in qqb_z2jet_v.f'
        stop
      endif
      
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

      res=0d0

      s(3,4) = (p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .     -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
      prop=s(3,4)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      
      v2(1)=l1
      v2(2)=r1

      do j=1,nf
        vQ(j,1)=L(j)
        vQ(j,2)=R(j)
      enddo


************************************************************************
*     Endpoint contributions from QQGG matrix elements                 *
*     Loop matrix elements are also initialized here                   *
************************************************************************
      j=vflav(1)
      k=vflav(2)

      if (Gflag) then
c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme. 
      subuv=2d0*xn*(epinv*(11d0-2d0*dble(nf)/xn)-1d0)/6d0

c--- Now calculate the relevant lowest-order matrix elements
c--- for each possible initial state from the QQGG contribution

      if (recalc_g) then 

c---  calculate the qqb terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p2)+g(p6)+g(p5)+qb(p1)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq_qqb,mmsq_qqb_vec,mmsq_qqb_ax)
     
c--- obtain qbq from qqb by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_qbq(polq,polz)=mmsq_qqb(3-polq,polz)
        mmsq_qbq_vec(polq,polz)=mmsq_qqb_vec(3-polq,polz)
        mmsq_qbq_ax(polq,polz)=-mmsq_qqb_ax(3-polq,polz)
      enddo
      enddo     

c---  calculate the gq terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p5)+g(p6)+g(p1)+qb(p2)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq_gq,mmsq_gq_vec,mmsq_gq_ax)

c--- obtain gqb from gq by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_gqb(polq,polz)=mmsq_gq(3-polq,polz)
        mmsq_gqb_vec(polq,polz)=mmsq_gq_vec(3-polq,polz)
        mmsq_gqb_ax(polq,polz)=-mmsq_gq_ax(3-polq,polz)
      enddo
      enddo
c---  calculate the qg terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p5)+g(p6)+g(p2)+qb(p1)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq_qg,mmsq_qg_vec,mmsq_qg_ax)

c--- obtain qbg from qg by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_qbg(polq,polz)=mmsq_qg(3-polq,polz)
        mmsq_qbg_vec(polq,polz)=mmsq_qg_vec(3-polq,polz)
        mmsq_qbg_ax(polq,polz)=-mmsq_qg_ax(3-polq,polz)
      enddo
      enddo
c--- calculate the gg terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p6)+g(p1)+g(p2)+qb(p5)+lbar(p4)+l(p3)
      do nu=1,4
       pswap(1,nu)=p(6,nu)
       pswap(2,nu)=p(1,nu)
       pswap(3,nu)=p(2,nu)
       pswap(4,nu)=p(5,nu)
       pswap(5,nu)=p(4,nu)
       pswap(6,nu)=p(3,nu)
 
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq_gg,mmsq_gg_vec,mmsq_gg_ax)
!      write(*,*) 'calc mmsq_gg', mmsq_gg
C     endif 

      endif ! recalc_g 
      endif ! Gflag 
************************************************************************
*     Endpoint contributions from QQQQ matrix elements                 *
************************************************************************            
      if (Qflag) then
c--- UV counter-term is already included in a6routine.f
      subuv=0d0
      
      
c--- transfer lowest order matrix elements
c--- NB: this breaks the routine if Qflag = Gflag = .true.

      do cs=0,2
        do j=-nf,nf
        do k=-nf,nf
        msq_cs(cs,j,k)=mqq(cs,j,k)
        enddo
        enddo
      enddo

      endif

c--- Add VIRTUAL terms

************************************************************************
*     Include loop contributions from QQGG matrix elements             *
************************************************************************
      if (Gflag) then

c--- compute correct vector-like coupling for diagrams with Z coupled to a loop
      vcouple(1)=czip
      vcouple(2)=czip
      do j=1,nf
      do polz=1,2
      vcouple(polz)=vcouple(polz)
     & +Q(j)*q1+0.5d0*(vQ(j,1)+vQ(j,2))*v2(polz)*prop
      enddo
      enddo

      do polq=1,2
      do polz=1,2

      j=vflav(1) 
      k=vflav(2) 

c--- quark-antiquark
      if ((j .gt. 0) .and. (k .lt. 0)) then
        if (j .eq. -k) 
     .  msqv(j,k)=msqv(j,k)+half*(aveqq/avegg)*(mmsq_qqb(polq,polz)*(
     .          cdabs(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)**2)
     .                     +dble(mmsq_qqb_vec(polq,polz)
     .        *dconjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     .  		   +dble(mmsq_qqb_ax(polq,polz)
     .        *dconjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .  	  *(v2(polz)*prop)/sin2w))
c--- antiquark-quark
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        if (j .eq. -k)
     .  msqv(j,k)=msqv(j,k)+half*(aveqq/avegg)*(mmsq_qbq(polq,polz)*(
     .          cdabs(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)**2)
     .                     +dble(mmsq_qbq_vec(polq,polz)
     .        *dconjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*vcouple(polz))
     .  		   +dble(mmsq_qbq_ax(polq,polz)
     .        *dconjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .  	  *(v2(polz)*prop)/sin2w))


c--- quark-gluon
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
      
************************** BEGIN MADLOOP CHECKING CODE ************************
        if (checkvector) then
c--- MadLoop check: vector couplings only - no Z coupling to leptons and quarks
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(
     .                   +dble(mmsq_qg_vec(polq,polz)
     .          *dconjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     .                   )
        endif

        if (checkaxial) then
c--- MadLoop check: axial couplings only - photon contribution removed
c------- full Z in Born (no photon)
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(
     .  		+dble(mmsq_qg_ax(polq,polz)
     .        *dconjg(Q(j)*q1*0+vQ(j,polq)*v2(polz)*prop)
     .         *(v2(polz)*prop)/sin2w))

c------- axial Z everywehere (Born and virt, quarks and leptons)
c        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(
c     .  		+dble(mmsq_qg_ax(polq,polz)
c     .         *dconjg(Q(j)*q1*0+
c     .           (vQ(j,polq)-vQ(j,3-polq))*(v2(polz)-v2(3-polz))*prop)
c     .         *((v2(polz)-v2(3-polz))*prop)/sin2w))/8d0
	endif
*************************** END MADLOOP CHECKING CODE *************************
	
	if ((checkvector.eqv..false.).and.(checkaxial.eqv..false.)) then
c--- normal case
	msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_qg(polq,polz)*(
     .  	cdabs(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)**2)
     .                     +dble(mmsq_qg_vec(polq,polz)
     .        *dconjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     .  		   +dble(mmsq_qg_ax(polq,polz)
     .        *dconjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .  	  *(v2(polz)*prop)/sin2w))


        endif
	
c--- antiquark-gluon
      elseif ((j .lt. 0) .and. (k .eq. 0)) then

        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_qbg(polq,polz)*(
     .          cdabs(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)**2)
     .                     +dble(mmsq_qbg_vec(polq,polz)
     .        *dconjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)*vcouple(polz))
     .  		   +dble(mmsq_qbg_ax(polq,polz)
     .        *dconjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     .  	  *(v2(polz)*prop)/sin2w))

c--- gluon-quark
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_gq(polq,polz)*(
     .          cdabs(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)**2)
     .                     +dble(mmsq_gq_vec(polq,polz)
     .        *dconjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*vcouple(polz))
     .  		   +dble(mmsq_gq_ax(polq,polz)
     .        *dconjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .  	  *(v2(polz)*prop)/sin2w))

c--- gluon-antiquark
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_gqb(polq,polz)*(
     .          cdabs(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)**2)
     .                     +dble(mmsq_gqb_vec(polq,polz)
     .        *dconjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)*vcouple(polz))
     .  		   +dble(mmsq_gqb_ax(polq,polz)
     .        *dconjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)
     .  	  *(v2(polz)*prop)/sin2w))

c--- gluon-gluon
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
        msqv(j,k)=msqv(j,k)+dfloat(ndo)*(mmsq_gg(polq,polz)*(
     .             cdabs(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)**2)
     .                     +dble(mmsq_gg_vec(polq,polz)
     .        *dconjg(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)*vcouple(polz))
     .  		   +dble(mmsq_gg_ax(polq,polz)
     .        *dconjg(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)
     .  	  *(v2(polz)*prop)/sin2w))
        msqv(j,k)=msqv(j,k)+dfloat(nup)*(mmsq_gg(polq,polz)*(
     .          cdabs(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)**2)
     .                     +dble(mmsq_gg_vec(polq,polz)
     .        *dconjg(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)*vcouple(polz))
     .  		   +dble(mmsq_gg_ax(polq,polz)
     .        *dconjg(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)
     .  	  *(v2(polz)*prop)/sin2w))
!
      endif
      enddo
      enddo
!      if (j.eq. 0 .and. k.eq. 0) stop  
      endif
      
************************************************************************
*     Include loop contributions from QQQQ matrix elements             *
************************************************************************

      if (Qflag) then
      
      j=vflav(1) 
      k=vflav(2) 

      faclo=4d0*V*aveqq*esq**2*gsq**2
      fac=faclo*xn*0.5d0*ason2pi

      if (recalc_q) then 
      call spinoru(6,p,za,zb)
c--- Set-up the desired amplitudes, in order to minimize number of calls
      do polq=1,2
      do polz=1,2
      do polb=1,2
c--- atreez
      atreez_526143(polq,polb,polz)=
     .       atreez(polq,polb,polz,5,2,6,1,4,3,za,zb)
      atreez_251643(polq,polb,polz)=
     .       atreez(polq,polb,polz,2,5,1,6,4,3,za,zb)
      atreez_625143(polq,polb,polz)=
     .       atreez(polq,polb,polz,6,2,5,1,4,3,za,zb)
      atreez_261543(polq,polb,polz)=
     .       atreez(polq,polb,polz,2,6,1,5,4,3,za,zb)
      atreez_162543(polq,polb,polz)=
     .       atreez(polq,polb,polz,1,6,2,5,4,3,za,zb)
      atreez_615243(polq,polb,polz)=
     .       atreez(polq,polb,polz,6,1,5,2,4,3,za,zb)
      atreez_152643(polq,polb,polz)=
     .       atreez(polq,polb,polz,1,5,2,6,4,3,za,zb)
      atreez_516243(polq,polb,polz)=
     .       atreez(polq,polb,polz,5,1,6,2,4,3,za,zb)
      atreez_562143(polq,polb,polz)=
     .       atreez(polq,polb,polz,5,6,2,1,4,3,za,zb)
      atreez_651243(polq,polb,polz)=
     .       atreez(polq,polb,polz,6,5,1,2,4,3,za,zb)
      atreez_265143(polq,polb,polz)=
     .       atreez(polq,polb,polz,2,6,5,1,4,3,za,zb)
      atreez_621543(polq,polb,polz)=
     .       atreez(polq,polb,polz,6,2,1,5,4,3,za,zb)
      atreez_561243(polq,polb,polz)=
     .       atreez(polq,polb,polz,5,6,1,2,4,3,za,zb)
      atreez_652143(polq,polb,polz)=
     .       atreez(polq,polb,polz,6,5,2,1,4,3,za,zb)
      atreez_165243(polq,polb,polz)=
     .       atreez(polq,polb,polz,1,6,5,2,4,3,za,zb)
      atreez_612543(polq,polb,polz)=
     .       atreez(polq,polb,polz,6,1,2,5,4,3,za,zb)  
c--- a61z
      a61z_526143(polq,polb,polz)=
     .     a61z(polq,polb,polz,5,2,6,1,4,3,za,zb)
      a61z_251643(polq,polb,polz)=
     .     a61z(polq,polb,polz,2,5,1,6,4,3,za,zb)
      a61z_625143(polq,polb,polz)=
     .     a61z(polq,polb,polz,6,2,5,1,4,3,za,zb)
      a61z_261543(polq,polb,polz)=
     .     a61z(polq,polb,polz,2,6,1,5,4,3,za,zb)
      a61z_162543(polq,polb,polz)=
     .     a61z(polq,polb,polz,1,6,2,5,4,3,za,zb)
      a61z_615243(polq,polb,polz)=
     .     a61z(polq,polb,polz,6,1,5,2,4,3,za,zb)
      a61z_152643(polq,polb,polz)=
     .     a61z(polq,polb,polz,1,5,2,6,4,3,za,zb)
      a61z_516243(polq,polb,polz)=
     .     a61z(polq,polb,polz,5,1,6,2,4,3,za,zb)
      a61z_562143(polq,polb,polz)=
     .     a61z(polq,polb,polz,5,6,2,1,4,3,za,zb)
      a61z_651243(polq,polb,polz)=
     .     a61z(polq,polb,polz,6,5,1,2,4,3,za,zb)
      a61z_265143(polq,polb,polz)=
     .     a61z(polq,polb,polz,2,6,5,1,4,3,za,zb)
      a61z_621543(polq,polb,polz)=
     .     a61z(polq,polb,polz,6,2,1,5,4,3,za,zb)
      a61z_561243(polq,polb,polz)=
     .     a61z(polq,polb,polz,5,6,1,2,4,3,za,zb)
      a61z_652143(polq,polb,polz)=
     .     a61z(polq,polb,polz,6,5,2,1,4,3,za,zb)
      a61z_165243(polq,polb,polz)=
     .     a61z(polq,polb,polz,1,6,5,2,4,3,za,zb)
      a61z_612543(polq,polb,polz)=
     .     a61z(polq,polb,polz,6,1,2,5,4,3,za,zb)
c---  a62z
      a62z_625143(polq,polb,polz)=
     .     a62z(polq,polb,polz,6,2,5,1,4,3,za,zb)
      a62z_261543(polq,polb,polz)=
     .     a62z(polq,polb,polz,2,6,1,5,4,3,za,zb)
      a62z_526143(polq,polb,polz)=
     .     a62z(polq,polb,polz,5,2,6,1,4,3,za,zb)
      a62z_251643(polq,polb,polz)=
     .     a62z(polq,polb,polz,2,5,1,6,4,3,za,zb)
      a62z_152643(polq,polb,polz)=
     .     a62z(polq,polb,polz,1,5,2,6,4,3,za,zb)
      a62z_516243(polq,polb,polz)=
     .     a62z(polq,polb,polz,5,1,6,2,4,3,za,zb)
      a62z_162543(polq,polb,polz)=
     .     a62z(polq,polb,polz,1,6,2,5,4,3,za,zb)
      a62z_615243(polq,polb,polz)=
     .     a62z(polq,polb,polz,6,1,5,2,4,3,za,zb)
      a62z_265143(polq,polb,polz)=
     .     a62z(polq,polb,polz,2,6,5,1,4,3,za,zb)
      a62z_621543(polq,polb,polz)=
     .     a62z(polq,polb,polz,6,2,1,5,4,3,za,zb)
      a62z_562143(polq,polb,polz)=
     .     a62z(polq,polb,polz,5,6,2,1,4,3,za,zb)
      a62z_651243(polq,polb,polz)=
     .     a62z(polq,polb,polz,6,5,1,2,4,3,za,zb)
      a62z_165243(polq,polb,polz)=
     .     a62z(polq,polb,polz,1,6,5,2,4,3,za,zb)
      a62z_612543(polq,polb,polz)=
     .     a62z(polq,polb,polz,6,1,2,5,4,3,za,zb)
      a62z_561243(polq,polb,polz)=
     .     a62z(polq,polb,polz,5,6,1,2,4,3,za,zb)
      a62z_652143(polq,polb,polz)=
     .     a62z(polq,polb,polz,6,5,2,1,4,3,za,zb)
c---  a63z
      a63z_526143(polq,polb,polz)=
     .     a63z(polq,polb,polz,5,2,6,1,4,3,za,zb)
      a63z_625143(polq,polb,polz)=
     .     a63z(polq,polb,polz,6,2,5,1,4,3,za,zb)
      a63z_162543(polq,polb,polz)=
     .     a63z(polq,polb,polz,1,6,2,5,4,3,za,zb)
      a63z_152643(polq,polb,polz)=
     .     a63z(polq,polb,polz,1,5,2,6,4,3,za,zb)
      a63z_562143(polq,polb,polz)=
     .     a63z(polq,polb,polz,5,6,2,1,4,3,za,zb)
      a63z_265143(polq,polb,polz)=
     .     a63z(polq,polb,polz,2,6,5,1,4,3,za,zb)
      a63z_561243(polq,polb,polz)=
     .     a63z(polq,polb,polz,5,6,1,2,4,3,za,zb)
      a63z_165243(polq,polb,polz)=
     .     a63z(polq,polb,polz,1,6,5,2,4,3,za,zb)
      enddo
      enddo
      enddo
      endif 

c---Desired formula =(Att*(A61+A61o+(A63-A63s/n+A62s+A62os)/n) 
c                     +Atts*(A61s+A61os+(A63s-A63/n +A62+A62o)/n))
c                   =tamp*(lamp+lampx)
c                   +tamps*(lamps+lampsx)
c   where lamp,lampsx are the pieces without the "s" (5<->6 swap)
c     and lamps,lampx are the pieces that have the "s"
      do polq=1,2
      do polz=1,2
      do polb=1,2
        tamp=czip
        lamp=czip
        lampx=czip
        tamps=czip
        lamps=czip
        lampsx=czip
c--- Q-Q
        if ((j .gt. 0) .and. (k .gt. 0)) then
           idfac=1d0
           tamp=atreez_526143(polq,polb,polz)
     .          *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .          -atreez_251643(3-polb,3-polq,polz)
     .          *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
           lamp=a61z_526143(polq,polb,polz)
     .          *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .          -a61z_251643(3-polb,3-polq,polz)
     .          *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
     .          +a63z_526143(polq,polb,polz)/xn
     .          *v2(polz)/sin2w*prop
           if (j .eq. k) then
c---  statistical factor
              idfac=0.5d0
              tamps=-(atreez_625143(polq,polb,polz)
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -atreez_261543(3-polb,3-polq,polz)
     .             *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop))
              lampx=-(
     .             -a63z_625143(polq,polb,polz)/xn**2
     .             *v2(polz)/sin2w*prop
     .             +a62z_625143(polq,polb,polz)/xn
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -a62z_261543(3-polb,3-polq,polz)/xn
     .             *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop))
              lamps=-(a61z_625143(polq,polb,polz)
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -a61z_261543(3-polb,3-polq,polz)
     .             *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
     .             +a63z_625143(polq,polb,polz)/xn
     .             *v2(polz)/sin2w*prop)
              lampsx=
     .             -a63z_526143(polq,polb,polz)/xn**2
     .             *v2(polz)/sin2w*prop
     .             +a62z_526143(polq,polb,polz)/xn
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -a62z_251643(3-polb,3-polq,polz)/xn
     .             *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
           endif
           born = born 
     .          +idfac*faclo*dble(tamp*dconjg(tamp))
     .          +idfac*faclo*dble(tamps*dconjg(tamps))
           if (polq .eq. polb) born=born
     .          +idfac*faclo*dble(tamp*dconjg(tamps))*(-2d0/xn)

           if (compare) then
              msqv(j,k) = born 
           else
              msqv(j,k)=msqv(j,k)
     .             +idfac*fac*2d0*dble(tamp*dconjg(lamp))
     .             +idfac*fac*2d0*dble(tamps*dconjg(lamps))
              if (polq .eq. polb) msqv(j,k)=msqv(j,k)
     .             +idfac*fac*2d0*dble(tamp*dconjg(lampx))
     .             +idfac*fac*2d0*dble(tamps*dconjg(lampsx))
           endif
c---  Qbar-Qbar
        elseif ((j .lt. 0) .and. (k .lt. 0)) then
           idfac=1d0
           tamp=atreez_162543(polq,polb,polz)
     .          *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     .          -atreez_615243(3-polb,3-polq,polz)
     .          *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
           lamp=a61z_162543(polq,polb,polz)
     .          *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     .          -a61z_615243(3-polb,3-polq,polz)
     .          *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     .          +a63z_162543(polq,polb,polz)/xn
     .          *v2(polz)/sin2w*prop
           if (j .eq. k) then
c---  statistical factor
              idfac=0.5d0
              tamps=-(atreez_152643(polq,polb,polz)
     .             *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     .             -atreez_516243(3-polb,3-polq,polz)
     .             *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop))
              lampx=-(
     .             -a63z_152643(polq,polb,polz)/xn**2
     .             *v2(polz)/sin2w*prop
     .             +a62z_152643(polq,polb,polz)/xn
     .             *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     .             -a62z_516243(3-polb,3-polq,polz)/xn
     .             *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop))
              lamps=-(a61z_152643(polq,polb,polz)
     .             *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     .             -a61z_516243(3-polb,3-polq,polz)
     .             *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     .             +a63z_152643(polq,polb,polz)/xn
     .             *v2(polz)/sin2w*prop)
              lampsx=
     .             -a63z_162543(polq,polb,polz)/xn**2
     .             *v2(polz)/sin2w*prop
     .             +a62z_162543(polq,polb,polz)/xn
     .             *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     .             -a62z_615243(3-polb,3-polq,polz)/xn
     .             *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
           endif
           born=born
     .          +idfac*faclo*dble(tamp*dconjg(tamp))
     .          +idfac*faclo*dble(tamps*dconjg(tamps))
           if (polq .eq. polb) born=born
     .          +idfac*faclo*dble(tamp*dconjg(tamps))*(-2d0/xn)
           if (compare) then
              msqv(j,k) = born 
           else
              msqv(j,k)=msqv(j,k)
     .             +idfac*fac*2d0*dble(tamp*dconjg(lamp))
     .             +idfac*fac*2d0*dble(tamps*dconjg(lamps))
              if (polq .eq. polb) msqv(j,k)=msqv(j,k)
     .             +idfac*fac*2d0*dble(tamp*dconjg(lampx))
     .             +idfac*fac*2d0*dble(tamps*dconjg(lampsx))
           endif
c---  Q-Qbar
        elseif ((j .gt. 0) .and. (k .lt. 0)) then
           tamp=atreez_562143(polq,polb,polz)
     .          *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .          -atreez_651243(3-polb,3-polq,polz)
     .          *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
           lamp=a61z_562143(polq,polb,polz)
     .          *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .          -a61z_651243(3-polb,3-polq,polz)
     .          *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     .          +a63z_562143(polq,polb,polz)/xn
     .          *v2(polz)/sin2w*prop
           if (j .eq. -k) then
              tamps=-(atreez_265143(polq,polb,polz)
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -atreez_621543(3-polb,3-polq,polz)
     .             *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop))
              lampx=-(
     .             -a63z_265143(polq,polb,polz)/xn**2
     .             *v2(polz)/sin2w*prop
     .             +a62z_265143(polq,polb,polz)/xn
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -a62z_621543(3-polb,3-polq,polz)/xn
     .             *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop))
              lamps=-(a61z_265143(polq,polb,polz)
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -a61z_621543(3-polb,3-polq,polz)
     .             *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     .             +a63z_265143(polq,polb,polz)/xn
     .             *v2(polz)/sin2w*prop)
              lampsx=
     .             -a63z_562143(polq,polb,polz)/xn**2
     .             *v2(polz)/sin2w*prop
     .             +a62z_562143(polq,polb,polz)/xn
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -a62z_651243(3-polb,3-polq,polz)/xn
     .             *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
              tampdo=atreez_265143(polq,polb,polz)
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -atreez_621543(3-polb,3-polq,polz)
     .             *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
              lampdo=a61z_265143(polq,polb,polz)
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -a61z_621543(3-polb,3-polq,polz)
     .             *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
     .             +a63z_265143(polq,polb,polz)/xn
     .             *v2(polz)/sin2w*prop
              tampup=atreez_265143(polq,polb,polz)
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -atreez_621543(3-polb,3-polq,polz)
     .             *(Q(2)*q1+vQ(2,polb)*v2(polz)*prop)
              lampup=a61z_265143(polq,polb,polz)
     .             *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .             -a61z_621543(3-polb,3-polq,polz)
     .             *(Q(2)*q1+vQ(2,polb)*v2(polz)*prop)
     .             +a63z_265143(polq,polb,polz)/xn
     .             *v2(polz)/sin2w*prop
           endif

           if (same .or. j.ne. -k) then 
           born=born
     .          +faclo*dble(tamp*dconjg(tamp))
           if (same) then 
              born= born
     .             +faclo*dble(tamps*dconjg(tamps))
              if (polq .eq. polb) born=born
     .             +faclo*dble(tamp*dconjg(tamps))*(-2d0/xn)
           endif

           if (compare) then
              msqv(j,k) = born 
           else
              msqv(j,k)=msqv(j,k)+fac*2d0*(
     .             +dble(tamp*dconjg(lamp)))
              if (same) then 
                 msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                +dble(tamps*dconjg(lamps)))
                 if (polq .eq. polb) msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                +dble(tamp*dconjg(lampx))
     .                +dble(tamps*dconjg(lampsx)))
              endif
           endif
        endif ! endif same 
c---  add additional annihilation diagrams if necessary
           if (j .eq. -k .and. .not.same) then
              if (jj(j) .eq. 1) then
                 born=born+faclo*(
     .                +dble(tampup*dconjg(tampup))*dfloat(nup)             
     .                +dble(tampdo*dconjg(tampdo))*dfloat(ndo))
                 
                 if (compare) then
                    msqv(j,k) = born 
                 else
                    msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                   +dble(tampup*dconjg(lampup))*dfloat(nup)             
     .                   +dble(tampdo*dconjg(lampdo))*dfloat(ndo))
                 endif
              elseif (jj(j) .eq. 2) then             
                 born=born+faclo*(
     .                +dble(tampup*dconjg(tampup))*dfloat(nup)             
     .                +dble(tampdo*dconjg(tampdo))*dfloat(ndo))
                 if (compare) then
                    msqv(j,k) = born 
                 else
                    msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                   +dble(tampup*dconjg(lampup))*dfloat(nup)             
     .                   +dble(tampdo*dconjg(lampdo))*dfloat(ndo))
                 endif
              endif
           endif
c---  Qbar-Q
        elseif ((j .lt. 0) .and. (k .gt. 0)) then
           tamp=atreez_561243(polq,polb,polz)
     .          *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .          -atreez_652143(3-polb,3-polq,polz)
     .          *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
           lamp=a61z_561243(polq,polb,polz)
     .          *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .          -a61z_652143(3-polb,3-polq,polz)
     .          *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
     .          +a63z_561243(polq,polb,polz)/xn
     .          *v2(polz)/sin2w*prop
           if (j .eq. -k) then
              tamps=-(atreez_165243(polq,polb,polz)
     .             *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .             -atreez_612543(3-polb,3-polq,polz)
     .             *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop))
              lampx=-(
     .             -a63z_165243(polq,polb,polz)/xn**2
     .             *v2(polz)/sin2w*prop
     .             +a62z_165243(polq,polb,polz)/xn
     .             *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .             -a62z_612543(3-polb,3-polq,polz)/xn
     .             *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop))
              lamps=-(a61z_165243(polq,polb,polz)
     .             *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .             -a61z_612543(3-polb,3-polq,polz)
     .             *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
     .             +a63z_165243(polq,polb,polz)/xn
     .             *v2(polz)/sin2w*prop)
              lampsx=
     .             -a63z_561243(polq,polb,polz)/xn**2
     .             *v2(polz)/sin2w*prop
     .             +a62z_561243(polq,polb,polz)/xn
     .             *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .             -a62z_652143(3-polb,3-polq,polz)/xn
     .             *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
              tampdo=atreez_165243(polq,polb,polz)
     .             *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .             -atreez_612543(3-polb,3-polq,polz)
     .             *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
              lampdo=a61z_165243(polq,polb,polz)
     .             *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .             -a61z_612543(3-polb,3-polq,polz)
     .             *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
     .             +a63z_165243(polq,polb,polz)/xn
     .             *v2(polz)/sin2w*prop
              tampup=atreez_165243(polq,polb,polz)
     .             *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .             -atreez_612543(3-polb,3-polq,polz)
     .             *(Q(2)*q1+vQ(2,polb)*v2(polz)*prop)
              lampup=a61z_165243(polq,polb,polz)
     .             *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .             -a61z_612543(3-polb,3-polq,polz)
     .             *(Q(2)*q1+vQ(2,polb)*v2(polz)*prop)
     .             +a63z_165243(polq,polb,polz)/xn
     .             *v2(polz)/sin2w*prop
           endif
           if (same .or. j.ne. -k) then 
              born=born
     .             +faclo*dble(tamp*dconjg(tamp))
              if (same) then 
                 born= born
     .                +faclo*dble(tamps*dconjg(tamps))
                 if (polq .eq. polb) born=born
     .                +faclo*dble(tamp*dconjg(tamps))*(-2d0/xn)
              endif
              if (compare) then
                 msqv(j,k) = born 
              else
                 msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                +dble(tamp*dconjg(lamp)))
                 if (same) then 
                    msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                   +dble(tamps*dconjg(lamps)))
                    if (polq .eq. polb) msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                   +dble(tamp*dconjg(lampx))
     .                   +dble(tamps*dconjg(lampsx)))
                 endif
              endif
           endif ! same or ... 
c---  add additional annihilation diagrams if necessary
           if (j .eq. -k .and. .not. same ) then
             if (jj(j) .eq. -1) then
                born=born+faclo*(
     .               +dble(tampup*dconjg(tampup))*dfloat(nup)             
     .               +dble(tampdo*dconjg(tampdo))*dfloat(ndo))
                if (compare) then
                   msqv(j,k) = born 
                else
                   msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                  +dble(tampup*dconjg(lampup))*dfloat(nup)             
     .                  +dble(tampdo*dconjg(lampdo))*dfloat(ndo))
                endif
             elseif (jj(j) .eq. -2) then             
                born=born+faclo*(
     .               +dble(tampup*dconjg(tampup))*dfloat(nup)             
     .               +dble(tampdo*dconjg(tampdo))*dfloat(ndo))
                if (compare) then
                   msqv(j,k) = born 
                else
                   msqv(j,k)=msqv(j,k)+fac*2d0*(
     .                  +dble(tampup*dconjg(lampup))*dfloat(nup)             
     .                  +dble(tampdo*dconjg(lampdo))*dfloat(ndo))
                endif
             endif
          endif
       endif
      enddo
      enddo
      enddo
      
      endif

c--- write out ug virtual amplitude when checking
      if (checkvector .or. checkaxial) then      
        write(6,*) 'Madloop check: ug Virt',msqv(2,0)
c	pause
      endif
      
************************************************************************
*     UV contributions are included here                               *
************************************************************************

C     This is the correction to put the answer in UV renormalized
C     dred scheme with msbar coupling

      msqv(j,k)=msqv(j,k)-ason2pi*subuv*res0

************************************************************************
C    Ths correction to MSbar scheme
      if (scheme .eq. 'msbr') then
      if (Qflag) then 
         msqv(j,k) = msqv(j,k) + 
     &        ason2pi*res0*(-4d0*(cf/2d0)- 0d0*Nc/6d0)
      else
         msqv(j,k) = msqv(j,k) + 
     &        ason2pi*res0*(-2d0*(cf/2d0)- 2d0*Nc/6d0)
      endif
      endif


C-----divide by ason2pi as wanted by powheg
      res = msqv(j,k)/ason2pi 
      return
      end
     
     
 
