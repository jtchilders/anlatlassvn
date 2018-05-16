module dimen
!
  integer, parameter :: d_offsh=1,d_offshcol=2,d_ampfld=3,d_ampcol=4, &
                      d_psnew=5,d_ordtmp=6,d_impul=7
!
end module dimen
!
module incnst
!
  integer*4, parameter :: nmax=10     !maximum number of particles
!  integer*4, parameter :: nxop=10    !maximum number of particles, optimized version
  integer*4, parameter :: nlor=6      !maximum number of lorentz degrees of freedom
  integer*4, parameter :: npmax=43    !maximum number of different particles "flavour"
  integer*4, parameter :: nmrgmx=5001 !
!
end module incnst
!
module extsrc
!
  use incnst
  implicit none
!
  complex*16 src(nlor)
  logical :: exflg=.false.
!
end module extsrc
!
module alptyp
!
  use incnst
  type proc
    integer*4, dimension(nmax) :: flv, hel   !FLV(j), HEL(j) flavour and helicity of j-th particle
    integer*4, dimension (2*nmax) :: col     !(COL(2*j-1),COL(2*j)) coulor of j-th particle (gluon == (qb,q))
    integer*4 :: nfrm                        !NFRM number of fermions for the selected process
    real*8, dimension (4*nmax) :: mom        !momenta of external particle: mom(4*j-3) enrgy mom(4*j-l),
!                                             l=2,1,0 px,py,pz of j-th particle
  end type proc
!
  type particle
!
    integer :: chcg          !chcg(j) is the charge conjugate of J
    integer :: cnschg        !cnschg(j) discriminate wether two particles shares the same
! conserve charges i.e. only Z,gluon and photon are equal (if 0 cabibbo mixing); for 
! optimization pourposes only 
!
    integer :: lortyp        !lorentz representation
    integer :: nlordof       !number of lorentz degrees of freedom
    integer :: colrep        !coulor representation
    integer :: zcoup         !couplings with the Z boson (massive fermions only)
!
  end type particle
!
end module alptyp
!
module naming
!
  use incnst
  use alptyp
! particle naming scheme internal to alpha
  integer*4, parameter :: gluon=1,wm=4,wp=3,z=2,photon=5,qdbr=6,qbtbr=7,qupbrl=8,qdbrl=9,     &
                          ebrl=10,nuebr=11,qupbrr=12,qdbrr=13,qcbr=14,qcbrl=15,qcbrr=16,      &
                          ebrr=17,qsbrl=18,qsbrr=19,qtopbr=20,qupbr=21,qup=22,qd=23,qbt=24,   &
                          qupl=25,qdl=26,qupr=27,qdr=28,el=29,nue=30,qc=31,qcr=32,er=33,      &
                          qsl=34,qsr=35,qtop=36,qcl=37,xpp=38,ypp=39,zzww=40,higgs=41,        &
                          xglu=42,yglu=43
  type(particle), dimension (npmax) :: prtcl
  integer*4, dimension (npmax) :: chcg,cnschg                    !chcg(j) is the charged conjugate of j
  integer, parameter :: glu=1,frm=3,frmbr=4,msgbs=2,frmhl=5,frmbrhl=6,frmhr=7,frmbrhr=8,axww=9,axzzww=10,pht=11,scal=12,xgl=13     !to disentangle lorentz representation
  integer, parameter :: up = 1, dn = 2, lpt = 3
  integer, dimension(npmax) :: conv
!
end module naming
!
!
!********************************************************************************
!
! process inizialization section
!
module prcint
!
  use naming
!
  integer*4, parameter :: nprcmx=1, nintmax=50                ! NPRCMX is the number of processes which needs
!                      to be computed simultaneously and NINTMAX the mazimal number of interaction terms
  integer*4 :: tabint(3*nintmax+1,nprcmx)                    ! interactions used for the computation
  integer*4 :: tabintmrg(9*nintmax+1,nprcmx)                 ! interactions used for the computation 
  integer*4 :: opam(nintmax,nprcmx),opfs(3*nintmax,nprcmx),nvasfs(3*nintmax+1,nprcmx),nvasam(nintmax+1,nprcmx)
  real*8    :: cpam(nintmax,nprcmx),cpfs(3*nintmax,nprcmx)
  integer*4 :: nprc
!
end module prcint
!
module couplings
!
  use incnst
  implicit none
!
  real*8, dimension (npmax) :: masses = -99.d0, width = -99.d0         ! MASSES(j) mass of j-th particle
!
  integer*4, dimension (npmax,npmax,npmax) :: oper     ! OPER(j,k,l) labels the operator to be used for 
!                                                        particles j,k and l interaction
  real*8, dimension (npmax,npmax,npmax) :: coup        ! COUP(j,k,l) coupling constant for OPER(j,k,l)
!
  real*8 :: gstrong,g2weak,gbar,ctw,stw,glup,gldn,gllp,grup,grdn,grlp,emch,trihcp,yuk   ! strong coupling constant 
  real*8 :: qrtc(3,4),qrtca(3,4),qrtc2(2,4),qrtc3(2,4)  !for quartic gauge couplings: 1-3==wp  4-6==wm  7-9==z 
  real*8 :: qrtchh(3),qrtchg(2) !for quartic higgs and higgs-gauge couplings
!
  real*8 :: ggh                 ! effective higgs-glu-glu coupling
!
end module couplings
!
subroutine initprc
!
  use prcint
  use couplings
  implicit none
!
  integer*4 :: j1,j2,pos,j3,j4,tmp,nit,nitm,x1,x2,x3              !auxiliary variables
  integer*4, dimension (3) :: tmpa                         !auxiliary variables
!
! reordering fields in ascendig order, according to alpha naming scheme, TABINT is used for the amplitude
!
  do j1=1,nprcmx
    nit=tabint(1,j1)
    do j2=1,nit
      pos=3*(j2-1)+2
      do j3=pos,pos+2
       do j4=j3+1,pos+2
        if (tabint(j3,j1).gt.tabint(j4,j1)) then
          tmp=tabint(j3,j1)
          tabint(j3,j1)=tabint(j4,j1)
          tabint(j4,j1)=tmp
        endif
       enddo
      enddo      
    enddo
  enddo
!
!  generating TABINTMRG which is used for the iterative steps before amplitude evaluation: to each 
!  three-particle vertex do correspond three equation of motion (< three if equal particles are involved
!  again TABINTMRG is ordered according to ALPHA NAMING scheme
!
  do j1=1,nprcmx
    nitm=0
    nit=tabint(1,j1)
    do j2=1,nit
      tmpa=tabint(3*(j2-1)+2:3*(j2-1)+4,j1)
      tabintmrg(3*nitm+2:3*nitm+4,j1)=tmpa
      nitm=nitm+1
      if (tmpa(2).ne.tmpa(1)) then
        tabintmrg(3*nitm+2:3*nitm+4,j1)=(/tmpa(2),tmpa(1),tmpa(3)/)
        nitm=nitm+1
      endif
      if (tmpa(3).ne.tmpa(1).and.tmpa(3).ne.tmpa(2)) then
        tabintmrg(3*nitm+2:3*nitm+4,j1)=(/tmpa(3),tmpa(1),tmpa(2)/)
        nitm=nitm+1
      endif
    enddo
    tabintmrg(1,j1)=nitm
  enddo
!
  do j1=1,nprcmx
    nitm=tabintmrg(1,j1)
    do j2=1,nitm
      do j3=j2+1,nitm
        if(tabintmrg(3*(j2-1)+2,j1).gt.tabintmrg(3*(j3-1)+2,j1)) then
          tmpa=tabintmrg(3*(j2-1)+2:3*(j2-1)+4,j1)
          tabintmrg(3*(j2-1)+2:3*(j2-1)+4,j1)=tabintmrg(3*(j3-1)+2:3*(j3-1)+4,j1)
          tabintmrg(3*(j3-1)+2:3*(j3-1)+4,j1)=tmpa
        endif
      enddo
    enddo
    nvasfs(1,j1)=0
    do j2=1,nitm
      x1=tabintmrg(3*(j2-1)+2,j1); x2=tabintmrg(3*(j2-1)+3,j1); x3=tabintmrg(3*(j2-1)+4,j1)
      opfs(j2,j1)=oper(x1,x2,x3); cpfs(j2,j1)=coup(x1,x2,x3); 
      if(x1.eq.gluon.or.x2.eq.gluon.or.x3.eq.gluon) then
        if(.not.(x1.eq.higgs.or.x2.eq.higgs.or.x3.eq.higgs))  then
          nvasfs(1,j1)= nvasfs(1,j1)+1; nvasfs(nvasfs(1,j1)+1,j1)=j2
        endif
      endif
    enddo
  enddo
!
  do j1=1,nprcmx
    nitm=tabint(1,j1)
    do j2=1,nitm
      do j3=j2+1,nitm
        if(tabint(3*(j2-1)+2,j1).gt.tabint(3*(j3-1)+2,j1)) then
          tmpa=tabint(3*(j2-1)+2:3*(j2-1)+4,j1)
          tabint(3*(j2-1)+2:3*(j2-1)+4,j1)=tabint(3*(j3-1)+2:3*(j3-1)+4,j1)
          tabint(3*(j3-1)+2:3*(j3-1)+4,j1)=tmpa
        endif
      enddo
    enddo
    nvasam(1,j1)=0
    do j2=1,nitm
      x1=tabint(3*(j2-1)+2,j1); x2=tabint(3*(j2-1)+3,j1); x3=tabint(3*(j2-1)+4,j1)
      opam(j2,j1)=oper(x1,x2,x3); cpam(j2,j1)=coup(x1,x2,x3);
      if(x1.eq.gluon.or.x2.eq.gluon.or.x3.eq.gluon) then
        if(.not.(x1.eq.higgs.or.x2.eq.higgs.or.x3.eq.higgs)) then
          nvasam(1,j1)= nvasam(1,j1)+1; nvasam(nvasam(1,j1)+1,j1)=j2
        endif
      endif
    enddo
  enddo
!
end subroutine initprc
!
! end process inizialization section
!
!*******************************************************************************
!
module strfld
!
! P (j,n) provides the binary representation of N : 
! (i.e.  P(J,N)=(1,3,7,0,0,...,0), N=2^(1-1)+2^(3-1)+2^(7-1) )
! DEC (n) reports how many J=1 are present in P(J,N) corresponding to the number of momenta contributing 
! to the given one, EVN(N)=0,1 if N is odd, even,  for optimization pourposes
  use incnst
!
  integer*4, dimension (2**(nmax+1)-1) :: momd      !number of momenta contributing 
!
! LORTYP(j) lables the lorentz representation of j; NLORDOF(j) returns the number of lorentz d.o.f. of j
! COLREP(j) lables the coulor representation of j
  integer*4, dimension(npmax) :: lortyp=-9999,nlordof=-9999,colrep=-9999
! FLDPTR is used to select the field amplitudes in the array OFFSH. FLDPTR(j1,j2,j3) returns :
! if j1=1 the number of off-shell fields of type J2 constructed out of J3 on shell external particles,
! if j1=2 the position in the array OFFSHCOL where relevant information for the fields are stored
! are stored (momenta, color ...)
  integer*4, dimension(2,npmax,nmax/2) :: fldptr=-9999
! NEXC,NEXCSM arrays for the iteration
  integer*4, parameter :: nexcmax=13
  integer*4, dimension(nexcmax,2:nmax/2) :: nexc,nexcsm
!
end module strfld
!
module incol
!
! the following array are logically two or three index arrays, are rearranged as one dimensional 
! for optimization pourposes
  integer*4, dimension(1058) :: tbmom   ! TBMOM(I,J) gives the coulor of the combination of particles i and j
  integer*4, dimension(529) :: tbnew    ! TBNEW(I,J) gives the number of different coulors from i and j 
!                                         combination
  real*8,  dimension(1058) :: tbcoeff   ! TBCOEFF(I,J) SU3 clebsh coefficient               
!
end module incol
!
!
subroutine initmom
!
  use strfld
  implicit none
! auxiliary variables
  integer*4 :: j1,i1,j2,aux,zr,j3,nf,nf0
  integer*4, dimension(nmax) :: s
  real*8 pm
  integer*4 popcnt0
!
  momd=0
  do j1=1,2**nmax
    momd(j1)=popcnt0(j1)
  enddo
!
  do j1=2,nmax/2
    aux=0                   ! number of possible J2+J3=J1 non simmetric case
    do j2=1,nmax/2
      do j3=1,nmax/2
        if (j3+j2.ne.j1) cycle
        aux=aux+1
        nexc(2*aux:2*aux+1,j1)=(/j2,j3/)
      enddo
    enddo
    nexc(1,j1)=aux
  enddo
!
  do j1=2,nmax/2
    aux=0                   ! number of possible J2+J3=J1  simmetric case
    do j2=1,nmax/2
      do j3=j2,nmax/2
        if (j3+j2.ne.j1) cycle
        aux=aux+1
        nexcsm(2*aux:2*aux+1,j1)=(/j2,j3/)
      enddo
    enddo
    nexcsm(1,j1)=aux
  enddo
!
end subroutine initmom
!
subroutine fillpr(i1,i2,i3,i4,i5,i6,tmp)
!
  use alptyp
  implicit none
  integer, intent(in) :: i1,i2,i3,i4,i5,i6
  type (particle), intent (out) :: tmp
!
  tmp%chcg=i1; tmp%cnschg=i2; tmp%lortyp=i3; tmp%nlordof=i4; tmp%colrep=i5; tmp%zcoup=i6;
!
end subroutine fillpr
!
subroutine inittab
!
  use naming
  use strfld, only : lortyp,nlordof,colrep
!
  call fillpr(gluon,0,glu,6,3,-99,prtcl(gluon)); 
  call fillpr(photon,0,pht,4,1,-99,prtcl(photon)); 
  call fillpr(qupbr,1,frm,4,2,up,prtcl(qup));  
  call fillpr(qcbr,2,frm,4,2,up,prtcl(qc));  
  call fillpr(qdbr,3,frm,4,2,dn,prtcl(qd));
  call fillpr(qbtbr,4,frm,4,2,dn,prtcl(qbt));
  call fillpr(qup,5,frmbr,4,2,up,prtcl(qupbr));  
  call fillpr(qc,6,frmbr,4,2,up,prtcl(qcbr));  
  call fillpr(qd,7,frmbr,4,2,dn,prtcl(qdbr));
  call fillpr(qbt,8,frmbr,4,2,dn,prtcl(qbtbr));
  call fillpr(qupbrl,9,frmhl,2,2,-99,prtcl(qupl));  
  call fillpr(qcbrl,10,frmhl,2,2,-99,prtcl(qcl));  
  call fillpr(qdbrl,11,frmhl,2,2,-99,prtcl(qdl));
  call fillpr(qupbrr,12,frmhr,2,2,-99,prtcl(qupr));  
  call fillpr(qcbrr,13,frmhr,2,2,-99,prtcl(qcr));  
  call fillpr(qdbrr,14,frmhr,2,2,-99,prtcl(qdr));
  call fillpr(qupl,15,frmbrhl,2,2,-99,prtcl(qupbrl));  
  call fillpr(qcl,16,frmbrhl,2,2,-99,prtcl(qcbrl));  
  call fillpr(qdl,17,frmbrhl,2,2,-99,prtcl(qdbrl));
  call fillpr(qupr,18,frmbrhr,2,2,up,prtcl(qupbrr));  
  call fillpr(qcr,19,frmbrhr,2,2,up,prtcl(qcbrr));  
  call fillpr(qdr,20,frmbrhr,2,2,dn,prtcl(qdbrr));
  call fillpr(wm,21,msgbs,5,1,-99,prtcl(wp));
  call fillpr(wp,22,msgbs,5,1,-99,prtcl(wm));
  call fillpr(z,0,msgbs,5,1,-99,prtcl(z));
  call fillpr(el,23,frmbrhl,2,1,-99,prtcl(ebrl));
  call fillpr(ebrl,24,frmhl,2,1,-99,prtcl(el));
  call fillpr(er,25,frmbrhr,2,1,-99,prtcl(ebrr));
  call fillpr(ebrr,26,frmhr,2,1,-99,prtcl(er));
  call fillpr(nue,27,frmbrhl,2,1,-99,prtcl(nuebr));
  call fillpr(nuebr,28,frmhl,2,1,-99,prtcl(nue));
  call fillpr(qsbrl,29,frmhl,2,2,-99,prtcl(qsl));
  call fillpr(qsbrr,30,frmhr,2,2,-99,prtcl(qsr));
  call fillpr(qsl,31,frmbrhl,2,2,-99,prtcl(qsbrl));
  call fillpr(qsr,32,frmbrhr,2,2,-99,prtcl(qsbrr));
  call fillpr(qtopbr,33,frm,4,2,up,prtcl(qtop));
  call fillpr(qtop,34,frmbr,4,2,up,prtcl(qtopbr));
  call fillpr(ypp,35,axww,1,1,-99,prtcl(xpp));
  call fillpr(xpp,36,axww,1,1,-99,prtcl(ypp));
  call fillpr(zzww,0,axzzww,2,1,-99,prtcl(zzww));
  call fillpr(higgs,0,scal,3,1,-99,prtcl(higgs));
  call fillpr(xglu,0,xgl,3,3,-99,prtcl(yglu));
  call fillpr(yglu,0,xgl,3,3,-99,prtcl(xglu));
!
  chcg=prtcl%chcg; cnschg=prtcl%cnschg; lortyp=prtcl%lortyp; nlordof=prtcl%nlordof
  colrep=prtcl%colrep; conv=prtcl%zcoup 
!
end subroutine inittab
!
subroutine initcoup
!
  use couplings
  use naming
  implicit none
!
  real*8 :: stw2,sqd2,ch(npmax)
  integer, dimension(npmax) :: tmp1,tmp2
  integer :: j1,nt
!
  call inital
!
  sqd2=sqrt(2.d0)
  oper(gluon,gluon,gluon)=1;  coup(gluon,gluon,gluon)=gstrong/sqd2
!
  oper(xglu,higgs,xglu)=47; coup(xglu,higgs,xglu)=ggh
  oper(higgs,xglu,xglu)=48; coup(higgs,xglu,xglu)=ggh
!
  oper(gluon,gluon,yglu)=49; coup(gluon,gluon,yglu)=gstrong/sqd2
  oper(yglu,gluon,gluon)=50; coup(yglu,gluon,gluon)=gstrong/sqd2
!
  oper(gluon,higgs,xglu)=51; coup(gluon,higgs,xglu)=ggh
  oper(higgs,gluon,xglu)=52; coup(higgs,gluon,xglu)=ggh
  oper(xglu,gluon,higgs)=53; coup(xglu,gluon,higgs)=ggh
!
  oper(gluon,gluon,higgs)=45; coup(gluon,gluon,higgs)=ggh
  oper(higgs,gluon,gluon)=46; coup(higgs,gluon,gluon)=ggh
!
  nt=5; tmp1(1:nt)=(/qdbr,qbtbr,qcbr,qtopbr,qupbr/); tmp2(1:nt)=(/qd,qbt,qc,qtop,qup/)
  do j1=1,nt
    oper(higgs,tmp1(j1),tmp2(j1))=40; oper(tmp1(j1),tmp2(j1),higgs) =41; oper(tmp2(j1),tmp1(j1),higgs) =41; 
    yuk=-masses(tmp1(j1))*gbar/(sqd2*2.d0*masses(z))
    coup(higgs,tmp1(j1),tmp2(j1))=yuk; coup(tmp1(j1),tmp2(j1),higgs) =yuk; coup(tmp2(j1),tmp1(j1),higgs) =yuk; 
  enddo
!
  trihcp=-masses(higgs)**2/masses(z)*0.25*gbar/sqd2
  qrtchh(1)=6.d0; qrtchh(2)=gbar/(8.d0*masses(z)*sqd2)*2.d0; qrtchh(3)=-2.d0/trihcp
  oper(higgs,higgs,higgs)=39; coup(higgs,higgs,higgs)=trihcp
!
  qrtchg(1)=g2weak/(4.d0*masses(wp)*sqd2); qrtchg(2)=gbar/(4.d0*masses(z)*sqd2);
  oper(wp,wm,higgs)=38; coup(wp,wm,higgs)=masses(wp)*g2weak/sqd2; oper(wm,wp,higgs)=38; 
  coup(wm,wp,higgs)=masses(wp)*g2weak/sqd2; oper(higgs,wp,wm)=37; coup(higgs,wp,wm)=masses(wp)*g2weak/sqd2; 
!
  oper(z,z,higgs)=38; coup(z,z,higgs)=masses(z)*gbar/sqd2;
  oper(higgs,z,z)=37; coup(higgs,z,z)=masses(z)*gbar/sqd2;
!
  qrtc=0.d0; qrtc2=0.d0; qrtc3 =0.d0
  qrtc2(1,1)=0.d0; qrtc2(2,1)=2.d0; qrtc2(1,2)=1.d0; qrtc2(2,2)=0.d0; qrtc2(1,3)=1.d0; 
  qrtc2(2,3)=0.d0; qrtc2(1,4)=0.d0; qrtc2(2,4)=2.d0; qrtc3(1,1)=2.d0; qrtc3(2,1)=0.d0;
  qrtc3(1,2)=0.d0; qrtc3(2,2)=1.d0;  qrtc3(1,4)=2.d0; qrtc3(2,4)=0.d0;

  oper(wp,wm,zzww)=28; oper(wm,wp,zzww)=28; oper(zzww,wp,wm)=29;
  coup(wp,wm,zzww)=-g2weak/sqd2; coup(wm,wp,zzww)=-g2weak/sqd2; coup(zzww,wp,wm)=-g2weak/sqd2;
  oper(z,z,zzww)=28; oper(zzww,z,z)=29; coup(z,z,zzww)=-g2weak/sqd2*ctw*ctw; coup(zzww,z,z)=-g2weak/sqd2*ctw*ctw; 
!
  oper(z,photon,zzww)=33; oper(photon,z,zzww)=32; oper(zzww,z,photon)=35; 
  coup(z,photon,zzww)=-g2weak/sqd2*ctw*stw; coup(zzww,z,photon)=-g2weak/sqd2*ctw*stw; 
  coup(photon,z,zzww)=-g2weak/sqd2*ctw*stw; 
!
  oper(photon,photon,zzww)=34;  oper(zzww,photon,photon)=36; 
  coup(photon,photon,zzww)=-g2weak/sqd2*stw*stw;  coup(zzww,photon,photon)=-g2weak/sqd2*stw*stw; 
!
  oper(wp,wp,ypp)=27; oper(ypp,wp,wp)=26; coup(wp,wp,ypp)=-g2weak*sqd2; coup(ypp,wp,wp)=-g2weak*sqd2
  oper(wm,wm,xpp)=27; oper(xpp,wm,wm)=26; coup(wm,wm,xpp)=g2weak/sqd2; coup(xpp,wm,wm)=g2weak/sqd2
!
  oper(z,wp,wm)=25;  oper(wp,z,wm)=25; oper(wm,z,wp)=25; 
  coup(z,wp,wm)=-g2weak*ctw/sqd2; coup(wp,z,wm)=g2weak*ctw/sqd2; coup(wm,z,wp)=-g2weak*ctw/sqd2; 
  oper(photon,wp,wm)=30;  oper(wp,wm,photon)=31; oper(wm,wp,photon)=31; 
  coup(photon,wp,wm)=-g2weak*stw/sqd2; coup(wp,wm,photon)=-g2weak*stw/sqd2; coup(wm,wp,photon)=g2weak*stw/sqd2; 
  qrtc=0.d0;  qrtc(1,1)=-1.d0/ctw; qrtc(2,2)=-1.d0/ctw; qrtc(2,3)=1.d0/ctw;  !Z5 W+ W-
              qrtc(2,1)=1.d0;      qrtc(1,2)=1.d0;      qrtc(3,3)=1.d0;     !W+5 Z W- 
              qrtc(3,1)=-1.d0;     qrtc(3,2)=1.d0;      qrtc(1,3)=1.d0;     !W-5 Z W+ 
  qrtca=0.d0; ! qrtca(1,4)=-1.d0/stw; qrtca(3,2)=1.d0/stw; qrtca(3,3)=-1.d0/stw;  !Z5 W+ W-
              qrtca(2,4)=1.d0;     qrtca(1,2)=-1.d0;      qrtca(2,3)=-1.d0;     !W+5 Z W- 
              qrtca(3,4)=-1.d0;     qrtca(2,2)=-1.d0;      qrtca(1,3)=-1.d0;     !W-5 Z W+ 
!
  nt=5; tmp1(1:nt)=(/qupbr,qdbr,qcbr,qbtbr,qtopbr/); tmp2(1:nt)=(/qup,qd,qc,qbt,qtop/); 
  do j1=1,nt  
    oper(gluon,tmp1(j1),tmp2(j1))=4; oper(tmp2(j1),gluon,tmp1(j1))=5;
    oper(tmp1(j1),gluon,tmp2(j1))=6; coup(gluon,tmp1(j1),tmp2(j1))=gstrong/sqd2; 
    coup(tmp2(j1),gluon,tmp1(j1))=gstrong/sqd2; coup(tmp1(j1),gluon,tmp2(j1))=gstrong/sqd2; 
  enddo
!
  nt=4; tmp1(1:nt)=(/qupbrl,qdbrl,qcbrl,qsbrl/); tmp2(1:nt)=(/qupl,qdl,qcl,qsl/); 
  do j1=1,nt  
    oper(gluon,tmp1(j1),tmp2(j1))=10; oper(tmp2(j1),gluon,tmp1(j1))=11
    oper(tmp1(j1),gluon,tmp2(j1))=12; coup(gluon,tmp1(j1),tmp2(j1))=gstrong/sqd2; 
    coup(tmp2(j1),gluon,tmp1(j1))=gstrong/sqd2; coup(tmp1(j1),gluon,tmp2(j1))=gstrong/sqd2; 
  enddo
!
  nt=4; tmp1(1:nt)=(/qupbrr,qdbrr,qcbrr,qsbrr/); tmp2(1:nt)=(/qupr,qdr,qcr,qsr/); 
  do j1=1,nt  
    oper(gluon,tmp1(j1),tmp2(j1))=16; oper(tmp2(j1),gluon,tmp1(j1))=17;
    oper(tmp1(j1),gluon,tmp2(j1))=18; coup(gluon,tmp1(j1),tmp2(j1))=gstrong/sqd2; 
    coup(tmp2(j1),gluon,tmp1(j1))=gstrong/sqd2; coup(tmp1(j1),gluon,tmp2(j1))=gstrong/sqd2; 
  enddo
!
  nt=2; tmp1(1:nt)=(/qdbr,qbtbr/); tmp2(1:nt)=(/qup,qtop/); 
  do j1=1,nt  
    oper(wm,tmp1(j1),tmp2(j1))=7; oper(tmp2(j1),wm,tmp1(j1))=8;
    oper(tmp1(j1),wm,tmp2(j1))=9; coup(wm,tmp1(j1),tmp2(j1))=g2weak/2.d0; 
    coup(tmp2(j1),wm,tmp1(j1))=g2weak/2.d0; coup(tmp1(j1),wm,tmp2(j1))=g2weak/2.d0; 
  enddo
!
  nt=2; tmp1(1:nt)=(/qupbr,qtopbr/); tmp2(1:nt)=(/qd,qbt/); 
  do j1=1,nt  
    oper(wp,tmp1(j1),tmp2(j1))=7; oper(tmp2(j1),wp,tmp1(j1))=8;
    oper(tmp1(j1),wp,tmp2(j1))=9; coup(wp,tmp1(j1),tmp2(j1))=g2weak/2.d0; 
    coup(tmp2(j1),wp,tmp1(j1))=g2weak/2.d0; coup(tmp1(j1),wp,tmp2(j1))=g2weak/2.d0; 
  enddo
!
  nt=3; tmp1(1:nt)=(/qdbrl,ebrl,qsbrl/); tmp2(1:nt)=(/qupl,nue,qcl/); 
  do j1=1,nt  
    oper(wm,tmp1(j1),tmp2(j1))=13; oper(tmp2(j1),wm,tmp1(j1))=14;
    oper(tmp1(j1),wm,tmp2(j1))=15; coup(wm,tmp1(j1),tmp2(j1))=g2weak/2.d0; 
    coup(tmp2(j1),wm,tmp1(j1))=g2weak/2.d0; coup(tmp1(j1),wm,tmp2(j1))=g2weak/2.d0; 
  enddo
!
  nt=3; tmp1(1:nt)=(/qupbrl,nuebr,qcbrl/); tmp2(1:nt)=(/qdl,el,qsl/); 
  do j1=1,nt  
    oper(wp,tmp1(j1),tmp2(j1))=13; oper(tmp2(j1),wp,tmp1(j1))=14;
    oper(tmp1(j1),wp,tmp2(j1))=15; coup(wp,tmp1(j1),tmp2(j1))=g2weak/2.d0; 
    coup(tmp2(j1),wp,tmp1(j1))=g2weak/2.d0; coup(tmp1(j1),wp,tmp2(j1))=g2weak/2.d0; 
  enddo
!
  stw2=stw**2
  glup=-gbar/sqd2*(1.d0-4.d0/3.d0*stw2);   gldn=gbar/sqd2*(1.d0-2.d0/3.d0*stw2);
  gllp=gbar/sqd2*(1.d0-2.d0*stw2); grup=gbar*2.d0*sqd2/3.d0*stw2; 
  grdn=-gbar*sqd2/3.d0*stw2;  grlp=-gbar*stw2*sqd2;  emch= -g2weak*stw/sqd2
!
  nt=6; tmp1(1:nt)=(/qupbrl,qcbrl,qdbrl,ebrl,qsbrl,nuebr/); tmp2(1:nt)=(/qupl,qcl,qdl,el,qsl,nue/); 
  ch(1:nt)=(/glup/2.d0,glup/2.d0,gldn/2.d0,gllp/2.d0,gldn/2.d0,-gbar/sqd2/2.d0/)
  do j1=1,nt  
    oper(z,tmp1(j1),tmp2(j1))=13; oper(tmp2(j1),z,tmp1(j1))=14;
    oper(tmp1(j1),z,tmp2(j1))=15; coup(z,tmp1(j1),tmp2(j1))=ch(j1); 
    coup(tmp2(j1),z,tmp1(j1))=ch(j1); coup(tmp1(j1),z,tmp2(j1))=ch(j1); 
  enddo
!
  nt=5; tmp1(1:nt)=(/qupbrr,qcbrr,qdbrr,ebrr,qsbrr/); tmp2(1:nt)=(/qupr,qcr,qdr,er,qsr/); 
  ch(1:nt)=(/grup/2.d0,grup/2.d0,grdn/2.d0,grlp/2.d0,grdn/2.d0/)
  do j1=1,nt  
    oper(z,tmp1(j1),tmp2(j1))=19; oper(tmp2(j1),z,tmp1(j1))=20;
    oper(tmp1(j1),z,tmp2(j1))=21; coup(z,tmp1(j1),tmp2(j1))=ch(j1); 
    coup(tmp2(j1),z,tmp1(j1))=ch(j1); coup(tmp1(j1),z,tmp2(j1))=ch(j1); 
  enddo
!
  nt=4; tmp1(1:nt)=(/qupbr,qdbr,qbtbr,qtopbr/); tmp2(1:nt)=(/qup,qd,qbt,qtop/); 
  do j1=1,nt 
    oper(z,tmp1(j1),tmp2(j1))=22; oper(tmp2(j1),z,tmp1(j1))=23;
    oper(tmp1(j1),z,tmp2(j1))=24; coup(z,tmp1(j1),tmp2(j1))=1.d0/2.d0; 
    coup(tmp2(j1),z,tmp1(j1))=1.d0/2.d0; coup(tmp1(j1),z,tmp2(j1))=1.d0/2.d0; 
  enddo
!
  nt=5; tmp1(1:nt)=(/qupbr,qdbr,qcbr,qbtbr,qtopbr/); tmp2(1:nt)=(/qup,qd,qc,qbt,qtop/); 
  ch(1:nt)=-(/2.d0/3.d0*emch,-1.d0/3.d0*emch,2.d0/3.d0*emch,-1.d0/3.d0*emch,2.d0/3.d0*emch/)
  do j1=1,nt  
    oper(photon,tmp1(j1),tmp2(j1))=4; oper(tmp2(j1),photon,tmp1(j1))=5;
    oper(tmp1(j1),photon,tmp2(j1))=6; coup(photon,tmp1(j1),tmp2(j1))=ch(j1); 
    coup(tmp2(j1),photon,tmp1(j1))=ch(j1); coup(tmp1(j1),photon,tmp2(j1))=ch(j1);
  enddo
!
  nt=5; tmp1(1:nt)=(/qupbrl,qdbrl,qcbrl,ebrl,qsbrl/); tmp2(1:nt)=(/qupl,qdl,qcl,el,qsl/); 
  ch(1:nt)=-(/2.d0/3.d0*emch,-1.d0/3.d0*emch,2.d0/3.d0*emch,-emch,-1.d0/3.d0*emch/)
  do j1=1,nt  
    oper(photon,tmp1(j1),tmp2(j1))=10; oper(tmp2(j1),photon,tmp1(j1))=11;
    oper(tmp1(j1),photon,tmp2(j1))=12; coup(photon,tmp1(j1),tmp2(j1))=ch(j1); 
    coup(tmp2(j1),photon,tmp1(j1))=ch(j1); coup(tmp1(j1),photon,tmp2(j1))=ch(j1);
  enddo
!
  nt=5; tmp1(1:nt)=(/qupbrr,qdbrr,qcbrr,ebrr,qsbrr/); tmp2(1:nt)=(/qupr,qdr,qcr,er,qsr/); 
  ch(1:nt)=-(/2.d0/3.d0*emch,-1.d0/3.d0*emch,2.d0/3.d0*emch,-emch,-1.d0/3.0*emch/)
  do j1=1,nt  
    oper(photon,tmp1(j1),tmp2(j1))=16; oper(tmp2(j1),photon,tmp1(j1))=17;
    oper(tmp1(j1),photon,tmp2(j1))=18; coup(photon,tmp1(j1),tmp2(j1))=ch(j1); 
    coup(tmp2(j1),photon,tmp1(j1))=ch(j1); coup(tmp1(j1),photon,tmp2(j1))=ch(j1);
  enddo
!
end subroutine initcoup
!
module utilities
!
  implicit none
  integer :: fl1,fl2,fl3
  logical :: diagonal
!
end module utilities
!
module running
!
  real*8 rnas         !alpha strong, to be changed on event by event basis
!
end module running
!
subroutine rncpl
!
  use running
  implicit none
  real*8 :: apar(100)
  common/alphapar/apar
!
  rnas=apar(56)       !alpha strong, to be changed on event by event basis   
!
end subroutine rncpl
!
subroutine amp (prcss,dim,mtel)
!
  use alptyp
  use strfld
  use prcint
  use couplings, only : oper, coup 
  use incol
  use dimen
  use utilities
  use running
  implicit none
!
  type (proc), intent(in) :: prcss
  integer, intent(in) :: dim(7)
  complex*16, intent(out) :: mtel             !amplitude
!
  integer*4, parameter :: nca=3
  integer*4 :: nprt,nfld,flvold,padd,nconf,j5,padd0
  integer*4 :: j1,j2,iter,steptm,stepex,ex2,ex3,nf2,nf3,psfl2,psfl3
  integer*4 :: pscfl2,pscfl3,nl1,nl2,nl3,j3,j4,beg3,nc,psfl3i,pscfl3i
  integer*4 :: p1,p2,p3,posfld,poscol,nlold
  integer*4 :: nf1,psfl1,pscfl1,pref,beg
  integer*4, dimension (2) :: cl,clnw
  integer*4 :: psnew(dim(d_psnew))
  integer*4, save :: fststp=0
  integer*4, dimension (9*nintmax+1) :: tabmrgpr
  integer*4, dimension(nexcmax) :: exc
  integer*4 :: offshcol(dim(d_offshcol)),ampcol(dim(d_ampcol)),ordtmp(dim(d_ordtmp))
  complex*16  ::  offsh(dim(d_offsh)),ampfld(dim(d_ampfld))
  complex*16, dimension (nlor) :: fusion,polar
  complex*16 :: amplitude,tmp
  real*8, dimension (2) :: colcff
  real*8 :: impul(4,dim(d_impul))
  real*8, dimension (4) :: imp
  real*8 :: perm
  integer*4, save, dimension(nmax) :: colconv(0:3,0:3)
  integer, save, dimension(-11:11,-11:11) :: cmbcol
!
  integer*4 :: nphf,even,nmom,ot,npmit,tmpi,ci,cf,j5b,ptr,nfr2
  integer*4 :: operot,jop,op(3*nintmax),sim
  logical :: choose,choose2,choose3,choose4,cnsold(0:npmax)
  real*8, dimension(4) :: imp3,imp2
  real*8 :: coupot,cpot1,cp(3*nintmax)
  integer,save :: psfit=0,pscit=0,psfa=0,psca=0
  logical, parameter :: dbg=.false.
!
! once for all inizialization
  if (fststp.eq.0) then
    fststp=1
! inizialization of arrays
    call initcoup
    call initmom
    call inittab
    call colprdsu3
!
! inizialization of merging tables
    call initprc
!
!   converting from (qb,q) coulor represenation to internal coulor representation
!   1,8 gluon; (-) 9,11 (anti) quarks; 0 coulorless
    colconv=1000; colconv(1,1)=1; colconv(2,2)=2; colconv(1,2)=3; colconv(1,3)=4
    colconv(2,3)=5; colconv(2,1)=6; colconv(3,1)=7; colconv(3,2)=8; colconv(0,1)=9
    colconv(0,2)=10; colconv(0,3)=11; colconv(1,0)=-9; colconv(2,0)=-10
    colconv(3,0)=-11; colconv(0,0)=0
!
! cmbcol(i,j) returns 1 if the particle i and j are coulor charged conjugate 
    cmbcol=0; cmbcol(1,1)=1; cmbcol(2,2)=1; cmbcol(3,6)=1; cmbcol(6,3)=1
    cmbcol(4,7)=1; cmbcol(7,4)=1; cmbcol(5,8)=1; cmbcol(8,5)=1; cmbcol(9,-9)=1
    cmbcol(10,-10)=1; cmbcol(11,-11)=1; cmbcol(-9,9)=1; cmbcol(-10,10)=1
    cmbcol(-11,11)=1; cmbcol(0,0)=1
!
  endif
!
  fldptr=-9999
  nprt=0                                   ! number of external particle
  do j1=1,nmax
    if (prcss%flv(j1).eq.1001) exit
    nprt=nprt+1
  enddo
!
  nphf=nprt/2              ! for optmization
  even=mod(nprt,2)         ! for optmization
  nmom=2**nprt-1           ! maximum number of different combinations of nprt 
!                            momenta, each momenta enters at most once
! inizialization step
!
  call rncpl         !loading running couplings (alpha-strong)
!
  posfld=1           ! each off shell amplitude will be collected into a single 
  poscol=1           ! array OFFSH, POSFLD and POSCOL allows to find the position
  flvold=-99          
  do j1=1,nprt
    if(prcss%flv(j1).ne.flvold) then
      if(flvold.gt.0) fldptr(1,flvold,1)=nfld 
      nfld=1
      flvold=prcss%flv(j1)
      fldptr(2,flvold,1)=poscol
    else
      nfld=nfld+1
    endif
    p1=2**(j1-1)
! for convenience if a diagonal gluon (1,1) or (2,2) is assigned as external
! particle we assume that also the other one is present with 0 polarization
    if (colconv(prcss%col(2*j1-1),prcss%col(2*j1)).eq.2) then
      offshcol(poscol:poscol+nca-1)=(/p1,posfld,1/)
      offsh(posfld:posfld+nlordof(prcss%flv(j1))-1)=0.d0
      posfld=posfld+nlordof(prcss%flv(j1))
      poscol=poscol+nca
      nfld=nfld+1
    endif
    call initpol(prcss,j1,offsh(posfld)) !OFFSH(posfld,....) contains the 
!                                        polarization of the jth particle
    offshcol(poscol:poscol+nca-1)=(/p1,posfld,colconv(prcss%col(2*j1-1),prcss%col(2*j1))/)
    impul(1:4,offshcol(poscol))=prcss%mom(4*(j1-1)+1:4*j1)
    posfld=posfld+nlordof(prcss%flv(j1))
    poscol=poscol+nca
    if (colconv(prcss%col(2*j1-1),prcss%col(2*j1)).eq.1) then
      offshcol(poscol:poscol+nca-1)=(/p1,posfld,2/)
      offsh(posfld:posfld+nlordof(prcss%flv(j1))-1)=0
      posfld=posfld+nlordof(prcss%flv(j1))
      poscol=poscol+nca
      nfld=nfld+1
    endif
  enddo
  fldptr(1,flvold,1)=nfld
!
! iteration 
! we merge on/off shell field2 and field3 into field1 via equation of motion up to the step required by the algorithm
  tabmrgpr=tabintmrg(1:,nprc)
  cp(1:tabintmrg(1,nprc))=cpfs(1:tabintmrg(1,nprc),nprc)
  op(1:tabintmrg(1,nprc))=opfs(1:tabintmrg(1,nprc),nprc)
  do j1=1,nvasfs(1,nprc)
    cp(nvasfs(j1+1,nprc))=cp(nvasfs(j1+1,nprc))*rnas
  enddo
  psnew(1:nmom)=0        
  do iter=2,nphf
    sim=0
    cnsold=.false.
    steptm=2                  !for optimization
    flvold=-99   
    do j1=1,tabmrgpr(1)
      fl1=tabmrgpr(steptm)    !labels the flavour of the daughter offshell particle
      fl2=tabmrgpr(steptm+1)  !fl2 and fl3 label the flavour of the parents offshell
      fl3=tabmrgpr(steptm+2)  !particles to be merged
!
      coupot=cp(j1)                                    !coup(fl1,fl2,fl3)
      operot=op(j1)                                    !oper(fl1,fl2,fl3)
!
      nl1=nlordof(fl1)        !nl1,nl2,nl3 are the corresponding numbers of lorentz
      nl2=nlordof(fl2)        !degrees of freedom
      nl3=nlordof(fl3)        
      if(fl1.ne.flvold) then
        if(flvold.gt.0) fldptr(1,chcg(flvold),iter)=nconf
        flvold=fl1
        if (cnsold(cnschg(fl1))) then
          psnew(1:nmom)=0                      !for optimization
          cnsold=.false.
        endif
        cnsold(cnschg(fl1))=.true.
        nconf=0
        fldptr(2,chcg(fl1),iter)=poscol
      endif
      steptm=steptm+3
      if(fl2.ne.fl3) then
        if(sim.ne.1) exc=nexc(1:,iter); sim=1
      else
        if (sim.ne.2) exc=nexcsm(1:,iter); sim=2
      endif
      stepex=2
      do j2=1,exc(1)
        ex2=exc(stepex)
        ex3=exc(stepex+1)
        stepex=stepex+2
        nf2=fldptr(1,fl2,ex2)                               !number of excitation of field2
        if(nf2.le.0) cycle
!
        nf3=fldptr(1,fl3,ex3)                               !number of excitation of field3
!
        if(nf3.le.0) cycle
!
        pscfl2=fldptr(2,fl2,ex2)                            !starting position of field2 in array OFFSHCOL
        pscfl3=fldptr(2,fl3,ex3)                            !starting position of field2 in array OFFSHCOL
        psfl2=offshcol(pscfl2+1)                             !starting position of field2 in array OFFSH
        psfl3=offshcol(pscfl3+1)                             !starting position of field3 in array OFFSH
        do j3=1,nf2
          if (ex2.eq.ex3.and.fl2.eq.fl3) then               !if the flavour of field 2 and 3 are the same don't 
            beg3=j3+1                                       !repeat twice the same computation 
            psfl3i=psfl3+nl3*(beg3-2)
            pscfl3i=pscfl3+nca*(beg3-2)
          else
           beg3=1
           psfl3i=psfl3-nl3
           pscfl3i=pscfl3-nca
          endif
!
          p2=offshcol(pscfl2)                               !momentum of field2 off-shell amplitude 
          cl(1)=offshcol(pscfl2+2); 
          imp2=impul(1:4,p2)
!
          choose=iter.eq.nphf.and.even.eq.0          !for optimization: true for even number of external particle and
!                                                    !the resulting off shell field is the merging of half this number
!
          if(prcss%nfrm.gt.2) then
            nfr2=ibits(p2,0,prcss%nfrm)
          else
            nfr2=0
          endif
!
          do j4=beg3,nf3
!
            pscfl3i=pscfl3i+nca
            psfl3i=offshcol(pscfl3i+1)                             !starting position of field3 in array OFFSH

!
            p3=offshcol(pscfl3i)                         !momentum of field3
            p1=p2+p3                                     !momentum of field1: resulting out of the merger of field2,3
!            if(p1.gt.nmom) cycle
            if(momd(p1).ne.iter) cycle               !if this is not true some momenta enter more than once in p1
            if(choose.and.btest(p1+1,0)) cycle       !for even number of particle we accept merger of npart/2 
!                                                        !momenta only if momenta 1 contribute 
            cl(2)=offshcol(pscfl3i+2)  ! cl=(/offshcol(pscfl2+2),offshcol(pscfl3i+2)/)  !contains the coulors of field2 and 3
            colprdit : select case (colrep(fl1))         
              case(2,3)                                  !quarks and gluons
                ptr=(cl(1)+12)+23*(cl(2)+11)   !pointer to deal with a one dimensional array as two dimensional, opt
                nc=tbnew(ptr)                  !number of new coulored object: 1,2 2 only for neutral gluons
                if(nc.eq.0) cycle
                ptr=2*ptr-1
                colcff(1)=tbcoeff(ptr); colcff(2)=tbcoeff(ptr+1); ! colcff=(/tbcoeff(ptr),tbcoeff(ptr+1)/)  !coulor clebsch
                clnw(1)=tbmom(ptr); clnw(2)=tbmom(ptr+1); !      clnw=tbmom(ptr:ptr+1)                   !coulor of the new object
              case(1)                                    !coulorless object
                if(cmbcol(cl(1),cl(2)).eq.0) cycle
                nc=1
                colcff(1)=1.d0; colcff(2)=0.d0;      ! colcff=(/1.,0./)
                clnw(1)=0; clnw(2)=0;              ! clnw=0
              case default
                write(*,*)'something wrong in color assignment'
            end select colprdit
            imp3=impul(1:4,p3)
            imp=imp2+imp3             !momentum of the new fields
            impul (1:4,p1)=imp                          !array to store momenta combinations
            call fuse (-imp,imp2,offsh(psfl2),imp3,  &     !merging the fields 2 and 3 into fusion
                         offsh(psfl3i),operot,fusion)         !via the proper trilinear interaction (oper)
            call prpgt (chcg(fl1),imp,fusion,polar)                      !multiplying by the propagator -> polar
            if(nfr2.gt.0) then
              call permtn(prcss%nfrm,p2,p3,perm)                            !returning proper fermi statistics in perm
            else
              perm=1.d0
            endif
            padd=psnew(p1)                    !returning the position, in OFFSHCOL, of field1
            if (padd.eq.0) then             !new configuration
              psnew(p1)=poscol
              if(diagonal) then;
                if(clnw(1).eq.1) then; nc=2; colcff(2)=0.d0; clnw(2)=2; endif
              endif
              do j5=1,nc
                cpot1=coupot*colcff(j5)*perm
                offsh(posfld:posfld+nl1-1)=cpot1*polar(1:nl1)            !storing field1
                offshcol(poscol:poscol+nca-1)=(/p1,posfld,clnw(j5)/)     !storing pointers to field1
                nconf=nconf+1                      !number of new excitation of the same flavour and number of momenta
                poscol=poscol+nca
                posfld=posfld+nl1
              enddo
            else                                                         !new contribution to an old configuration
              j5b=0
              if(diagonal) then;
                if(clnw(1).eq.2) j5b=nl1
              endif
              do j5=1,nc
                cpot1=coupot*colcff(j5)*perm
                padd0=offshcol(padd+1)+j5b                                                !position of field1 in OFFSH
                offsh(padd0:padd0+nl1-1)=offsh(padd0:padd0+nl1-1)+cpot1*polar(1:nl1)      !adding the new contribution
                j5b=j5b+nl1
              enddo
            endif
!
          enddo
          pscfl2=pscfl2+nca
          psfl2=offshcol(pscfl2+1)
        enddo
      enddo
!
    enddo
!
    fldptr(1,chcg(fl1),iter)=nconf
  enddo
!
  if (fststp.eq.1.and.dbg) then
    fststp=2
    if (poscol.gt.pscit) then
      pscit=poscol
      write(*,*)'iteration,poscol,posfld',poscol,posfld
    endif
    if (posfld.gt.psfit) then
      psfit=posfld
      write(*,*)'iteration,poscol,posfld',poscol,posfld
    endif
  endif
!
  if(poscol.gt.dim(d_offshcol).or.posfld.gt.dim(d_offsh)) then
    write(*,*)'ill dimensioned array OFFSH and/or OFFSHCOL in AMP'
    write(*,*)'iteration,poscol,posfld',poscol,posfld,dim(d_offshcol),dim(d_offsh)
    stop
  endif
!
!
! amplitude 
!
  amplitude=(0.,0.)
!
  tabmrgpr(1:3*nintmax+1)=tabint(1:3*nintmax+1,nprc)      !most of the steps are the same as in iteration
  cp(1:tabint(1,nprc))=cpam(1:tabint(1,nprc),nprc)
  op(1:tabint(1,nprc))=opam(1:tabint(1,nprc),nprc)
  do j1=1,nvasam(1,nprc)
    cp(nvasam(j1+1,nprc))=cp(nvasam(j1+1,nprc))*rnas
  enddo    
!                                            we first merge field2 and field3 and finally contract with field1
!                                            only the list of acceptable fusions is different
  cnsold=.false.
  psnew(1:nmom)=0             !none of these mergin has been computed before
  ordtmp(1:nmom)=0            !reporting the order of field1 momenta to avoid repeating the same computation for
!                                  interaction quadratic or cubic in the same field
  steptm=2 
  flvold=-99   
  do j1=1,tabmrgpr(1)+1
!
    if(j1.eq.tabmrgpr(1)+1) then
      flvold=9999                                !last step compute the final contribution to the amplitude and exit
    else
      fl1=tabmrgpr(steptm)
      fl2=tabmrgpr(steptm+1)
      fl3=tabmrgpr(steptm+2)
      steptm=steptm+3    
      nl1=nlordof(fl1)
      nl2=nlordof(fl2)
      nl3=nlordof(fl3)
!
      coupot=cp(j1)          !coup(fl1,fl2,fl3)
      operot=op(j1)          !oper(fl1,fl2,fl3)
!
    endif
!
    if(fl1.ne.flvold) then             !we have computed all 2,3 merging compatible with contraction with old 1
!
  if (fststp.eq.2.and.j1.ne.1) then
!
    if(poscol.gt.dim(d_ampcol).or.posfld.gt.dim(d_ampfld)) then
      write(*,*)'dim',dim
      write(*,*)'ill dimensioned array AMPCOL and/or AMPFLD in AMP'
      write(*,*)'amplitude,poscol,posfld',poscol,posfld,dim(d_ampcol),dim(d_ampfld)
      stop
    endif
!
    fststp=2
    if(dbg) then
      if (poscol.gt.psca) then
        psca=poscol
        write(*,*)'amplitude,poscol,posfld',poscol,posfld
      endif
      if (posfld.gt.psfa) then
        psfa=posfld
        write(*,*)'amplitude,poscol,posfld',poscol,posfld
      endif
    endif
  endif
!
      if (flvold.gt.0) then            !compute the contribution to the amplitude (flvold < 0 iteration just finished)
        if(flvold.eq.9999)flvold=fl1   !flvold=9999 last contribution
        cnsold(cnschg(flvold))=.true.
        nlold=nlordof(flvold)          !number of lorentz degrees of freedom          
        do iter=1,nphf
          nf1=fldptr(1,flvold,iter)       !same as in iteration
          if(nf1.le.0) cycle
          pscfl1=fldptr(2,flvold,iter)
          psfl1=offshcol(pscfl1+1)
          pref=0                          !for optimization
          do j2=1,nf1
            psfl1=offshcol(pscfl1+1)
            p1=offshcol(pscfl1)
            p1=nmom-p1
            padd=psnew(p1)                !locating position of field1
            if (padd.eq.0) then           !no 2,3 merging compatible with field1
              pscfl1=pscfl1+nca
              cycle
            endif
            if (p1.ne.pref) then
              pref=p1
            else
              padd=padd+nca               !there are two neutral gluons
            endif
            cl=(/offshcol(pscfl1+2),ampcol(padd+2)/)
            if(cmbcol(cl(1),cl(2)).ne.1) cycle
            padd=ampcol(padd+1)
            tmp=sum(offsh(psfl1:psfl1+nlold-1)*ampfld(padd:padd+nlold-1))
            if(prcss%nfrm.gt.2) then
              call permtn(prcss%nfrm,nmom-p1,p1,perm)                            !returning proper fermi statistics in perm
            else
              perm=1.d0
            endif
            amplitude=amplitude+tmp*perm
            pscfl1=pscfl1+nca
          enddo
        enddo
        if(j1.eq.tabmrgpr(1)+1) cycle                             !end of computation
        if(cnsold(cnschg(fl1))) then
          psnew(1:nmom)=0             !none of these mergin has been computed before
          ordtmp(1:nmom)=0            !reporting the order of field1 momenta to avoid repeating the same computation for
!                                      interaction quadratic or cubic in the same field

        endif
!
      endif
! 
      flvold=fl1                  !a new class of fusion 2,3 beginning, compatible to match with new flavour fl1
      posfld=1                    !pointer to offsh
      poscol=1                    !pointer to offshcol    
!
    endif
!
    do iter=1,nphf
!
      npmit=nprt-iter
      choose2=iter.eq.nphf.and.even.eq.0
      nf1=fldptr(1,fl1,iter)
      if (nf1.le.0) cycle
      pscfl1=fldptr(2,fl1,iter)
!
      if(fl1.eq.fl2) then                     !again to deal with equal flavour
        beg=iter
      else
        beg=1
      endif
!
      do ex2=beg,nphf
!
        choose=fl1.eq.fl2.and.iter.eq.ex2     !again to deal with equal flavour  
!
        ex3=nprt-iter-ex2                     !the total number of momenta must match the number of external particle
        if(ex3.lt.1.or.ex3.gt.nphf) cycle
        if (fl2.eq.fl3.and.ex3.lt.ex2) cycle  !again to deal with equal flavour  
        nf2=fldptr(1,fl2,ex2)
        if(nf2.le.0) cycle
        nf3=fldptr(1,fl3,ex3)
        if(nf3.le.0) cycle
        pscfl2=fldptr(2,fl2,ex2)
        pscfl3=fldptr(2,fl3,ex3)
!
        choose3=ex2.eq.ex3.and.fl2.eq.fl3     !again to deal with equal flavour for optimization
!
          j5=pscfl1
          do j4=1,nf1
            ordtmp(offshcol(j5))=j4           !now ordtmp register the ordering of momenta of field1
            j5=j5+nca
          enddo
!
        do j3=1,nf2
!
          cl(1)=offshcol(pscfl2+2)  !ci=offshcol(pscfl2+2)
          psfl2=offshcol(pscfl2+1)
!
          if (choose3) then
            beg3=j3+1
            pscfl3i=pscfl3+nca*(beg3-2)
          else
            beg3=1
            pscfl3i=pscfl3-nca
          endif
!
          p2=offshcol(pscfl2)
          imp2=impul(1:4,p2)
          if(prcss%nfrm.gt.2) then
            nfr2=ibits(p2,0,prcss%nfrm)
          else
            nfr2=0
          endif
!
          do j4=beg3,nf3
            pscfl3i=pscfl3i+nca
!
            p3=offshcol(pscfl3i) 
            p1=p2+p3
!            if(p1.gt.nmom) cycle
            if(momd(p1).ne.npmit) cycle          !momenta appearing more than once
            ot=ordtmp(nmom-p1)
            if(ot.eq.0) cycle              !field1 does not contain the frequency to match the proposed 2,3 fusion
            if(choose.and.ordtmp(p2).le.ot)cycle     !to avoid repeating the same computation twice
!            if(choose2.and.momb%evn(p1).eq.1) cycle
            if(choose2.and.btest(p1,0)) cycle
            cl(2)=offshcol(pscfl3i+2) ! cl=(/ci,offshcol(pscfl3i+2)/)
            colprdam : select case (colrep(fl1))
              case(2,3)
                ptr=(cl(1)+12)+23*(cl(2)+11)
                nc=tbnew(ptr)
                if(nc.eq.0) cycle
                ptr=2*ptr-1
                colcff=tbcoeff(ptr:ptr+1)
                clnw=tbmom(ptr:ptr+1)
              case(1)
                if(cmbcol(cl(1),cl(2)).eq.0) cycle
                nc=1
                colcff(1)=1.d0; colcff(2)=0.d0;        !  colcff=(/1.,0./)
                clnw(1)=0; clnw(2)=0                   !  clnw=0
              case default
                write(*,*)'something wrong in color assignment'
              end select colprdam
!
              imp3=impul(1:4,p3)
              imp=imp2+imp3
              psfl3i=offshcol(pscfl3i+1)
              call fuse (-imp,imp2,offsh(psfl2),imp3,offsh(psfl3i),operot,fusion)
!              call permt(prcss%nfrm,p2,p3,perm)
              if(nfr2.gt.0) then
!                if(ibits(p3,0,prcss%nfrm).ne.0) then
                  call permtn(prcss%nfrm,p2,p3,perm)
!                else
!                  perm=1.d0
!                endif
              else
                perm=1.d0
              endif
!
              padd=psnew(p1)
              if (padd.eq.0) then             !new configuration
                psnew(p1)=poscol
                if(diagonal) then;
                  if(clnw(1).eq.1) then; nc=2; colcff(2)=0.d0; clnw(2)=2; endif
                endif
                do j5=1,nc
                  cpot1=coupot*colcff(j5)*perm
                  ampfld(posfld:posfld+nl1-1)=cpot1*fusion(1:nl1)
                  ampcol(poscol:poscol+nca-1)=(/p1,posfld,clnw(j5)/)
                  nconf=nconf+1
                  posfld=posfld+nl1
                  poscol=poscol+nca
                enddo
              else                            !new contribution to an old configuration
                j5b=0
                if(diagonal) then;
                  if(clnw(1).eq.2) j5b=nl1
                endif
                do j5=1,nc
                  padd0=ampcol(padd+1)+j5b
                  cpot1=coupot*colcff(j5)*perm
                  ampfld(padd0:padd0+nl1-1)=ampfld(padd0:padd0+nl1-1)+cpot1*fusion(1:nl1)
                  j5b=j5b+nl1
                enddo
              endif
!
          enddo
          pscfl2=pscfl2+nca
        enddo
      enddo
!
    enddo
!
  enddo
!
  mtel=amplitude
  fststp=1
!
end subroutine amp
!
subroutine initpol(prcss,prt,fld)              
!
  use strfld
  use alptyp
  use couplings, only : masses
  use naming
  use extsrc
  implicit none
!
  type (proc), intent(in) :: prcss                       !input   process
  integer, intent(in) :: prt                             !prt-th to be initialized
  complex*16, intent(out), dimension(nlor) :: fld        !polarization vector of particle number prt
!
  integer tmp
  real*8 :: im(4)
  complex*16 :: tmpc(4)
!
  integer, parameter :: maxpar=10
  integer :: vcnt
  character*1 idecay
  common/tdecflag/idecay
  complex*16, dimension(4) :: vtp,utpb
  common/tspinors/vtp,utpb
  complex*16 eps(4,maxpar)
  common/vdec/eps
  save vcnt
!
  if(prt.eq.1) vcnt=0
!
  tmp=prcss%hel(prt)
  im=prcss%mom(4*(prt-1)+1:4*prt)
  if(im(1).lt.0) im=-im
  polarization: select case (lortyp(prcss%flv(prt)))
    case (glu,pht)                                     !gluon source
!      if (prt.eq.1) then; fld(1:3)=im(2:4); fld(4:6)=0.d0; return; endif
      call sourcemassboson(tmp,tmpc,im,masses(prcss%flv(prt)))
      fld(1:3)=tmpc(2:4)
      fld(4:6)=0.d0
    case (frm)                                         !fermion source (incoming fermion outcoming anti fermion)
      if(prcss%flv(prt).eq.qtop.and.idecay.eq.'y') then
        fld(1:4)=vtp
      else
        call fermionsources(tmp,fld,im,masses(prcss%flv(prt)))
      endif
    case (frmbr)                                       !fermion bar source
      if(prcss%flv(prt).eq.qtopbr.and.idecay.eq.'y') then
        fld(1:4)=utpb
      else
        call fermionbarsources(tmp,fld,im,masses(prcss%flv(prt)))
      endif
    case (frmhl)                                       !fermion source (incoming fermion outcoming anti fermion)
      if(tmp.lt.0) then
        write(*,*)'wrong helicity assignment',tmp
        stop
      endif
      call fermionsourceshl(tmp,fld,im)
    case (frmbrhl)                                     !fermion bar source
      if(tmp.lt.0) then
        write(*,*)'wrong helicity assignment',tmp
        stop
      endif
      call fermionbarsourceshl(tmp,fld,im)
    case (frmhr)                                       !fermion source (incoming fermion outcoming anti fermion)
      if(tmp.gt.0) then
        write(*,*)'wrong helicity assignment',tmp
        stop
      endif
      call fermionsourceshl(tmp,fld,im)
    case (frmbrhr)                                     !fermion bar source
      if(tmp.gt.0) then
        write(*,*)'wrong helicity assignment',tmp
        stop
      endif
      call fermionbarsourceshl(tmp,fld,im)
    case (msgbs)
      if (exflg) then                                  !w,z source
        fld(1:4)=src(1:4)
        fld(5)=0.d0
      else      
        if(idecay.eq.'y') then
          vcnt=vcnt+1
          tmpc=eps(1:4,vcnt)
        else
          call sourcemassboson(tmp,tmpc,im,masses(prcss%flv(prt)))
        endif
        fld(1:4)=tmpc
        fld(5)=0.d0
      endif
    case (scal)
        fld(1)=1.d0
        fld(2:3)=0.d0
    case default
       write(6,*)'wrong lorenz type in INITPOL',lortyp(prcss%flv(prt)),prcss%flv(prt)
       stop
  end select polarization
!
end  subroutine initpol
!
subroutine fuse (p1,p2,amp2,p3,amp3,oper,nwam)
!
  use incnst
  use naming, only : z,wp,wm
  use couplings, only : qrtc,qrtc2,qrtc3,qrtca,qrtchg
  use utilities, only : fl1,fl2,fl3,diagonal
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3           !on/off-shell amplitudes to be merged
  real*8, intent(in), dimension (4) :: p1,p2,p3                   !particle momenta
  integer*4, intent(in) :: oper                               !labelling the interaction type and the flavour of the 3rd field
  complex*16, intent(out), dimension (nlor) :: nwam               !result of the merging
!
  real*8 cp
!
  diagonal=.false.
!
  operator: select case (oper)
    case(1)
      call triglu (p1(2),p2(2),amp2,p3(2),amp3,nwam)    !trilinear gluon interaction
    case (4)
      call fbftog (amp2,amp3,nwam)             !fusion of fermions into gluon
    case (5)
      call gfbtofb (amp2,amp3,nwam)            !fusion of f-bar  gluon into f-bar
    case (6)
      call gftof (amp2,amp3,nwam)              !fusion of f-bar  gluon into f-bar
    case (7)
      call fbftow (amp2,amp3,nwam)             !fusion of fermions into gluon
    case (8)
      call wfbtofb (amp2,amp3,nwam)            !fusion of f-bar  gluon into f-bar
    case (9)
      call wftof (amp2,amp3,nwam)              !fusion of f-bar  gluon into f-bar
    case (10)
      call fbftogl (amp2,amp3,nwam)             !fusion of fermions into gluon
    case (11)
      call gfbtofbl (amp2,amp3,nwam)            !fusion of f-bar  gluon into f-bar
    case (12)
      call gftofl (amp2,amp3,nwam)              !fusion of f-bar  gluon into f-bar
    case (13)
      call fbftowl (amp2,amp3,nwam)             !fusion of fermions into gluon
    case (14)
      call wfbtofbl (amp2,amp3,nwam)            !fusion of f-bar  gluon into f-bar
    case (15)
      call wftofl (amp2,amp3,nwam)              !fusion of f-bar  gluon into f-bar
    case (16)
      call fbftogr (amp2,amp3,nwam)             !fusion of fermions into gluon
    case (17)
      call gfbtofbr (amp2,amp3,nwam)            !fusion of f-bar  gluon into f-bar
    case (18)
      call gftofr (amp2,amp3,nwam)              !fusion of f-bar  gluon into f-bar
    case (19)
      call fbftowr (amp2,amp3,nwam)             !fusion of fermions into gluon
    case (20)
      call wfbtofbr (amp2,amp3,nwam)            !fusion of f-bar  gluon into f-bar
    case (21)
      call wftofr (amp2,amp3,nwam)              !fusion of f-bar  gluon into f-bar
    case (22)
      call fbftoz (fl3,amp2,amp3,nwam)             !fusion of fermions into gluon
    case (23)
      call zfbtofb (fl3,amp2,amp3,nwam)            !fusion of f-bar  gluon into f-bar
    case (24)
      call zftof (fl3,amp2,amp3,nwam)              !fusion of f-bar  gluon into f-bar
    case(25)
      call triw (qrtc(1,fl1-1),p1,p2,amp2,p3,amp3,nwam)    !trilinear gluon interaction 
    case (26)
      call wwtoax (amp2,amp3,nwam)              !fusion W+ W+ into X++
    case (27)
      call waxtow (amp2,amp3,nwam)              !fusion W+ Y++ into W+
    case (28)
      call wax2tow (qrtc2(1,fl1-1),amp2,amp3,nwam)              !fusion W+ Y++ into W+
    case (29)
      call wwtoax2 (qrtc3(1,fl2-1),amp2,amp3,nwam)              !fusion W+ Y++ into W+
    case(30)
      call triphw (qrtca(1,fl1-1),p1,p2,amp2,p3,amp3,nwam)     ! fusion W+ W- into photon
    case(31)
      call triwph (qrtca(1,fl1-1),p1,p2,amp2,p3,amp3,nwam)     ! fusion W photot into W
    case (32)
      call wax2top (qrtc2(1,fl1-1),amp2,amp3,nwam)              !fusion W+ Y++ into W+
    case (33)
      call pax2tow (qrtc2(1,fl1-1),amp2,amp3,nwam)              !fusion W+ Y++ into W+
    case (34)
      call pax2top (qrtc2(1,fl1-1),amp2,amp3,nwam)              !fusion W+ Y++ into W+
    case (35)
      call zptoax2 (qrtc3(1,fl2-1),amp2,amp3,nwam)              !fusion W+ Y++ into W+
    case (36)
      call pptoax2 (qrtc3(1,fl2-1),amp2,amp3,nwam)              !fusion W+ Y++ into W+
    case (37)
      gb0 : select case (fl2)
        case (wp,wm)
          cp=qrtchg(1)
        case(z)
          cp=qrtchg(2)
        case default
          write(*,*)'wrong flavour in fuse'; stop
      end select gb0
      call wwtoh (cp,amp2,amp3,nwam)              !fusion W+ W+ into X++
    case (38)
      gb1 : select case (fl2)
        case (wp,wm)
          cp=qrtchg(1)
        case(z)
          cp=qrtchg(2)
        case default
          write(*,*)'wrong flavour in fuse'; stop
      end select gb1
      call whtow (cp,amp2,amp3,nwam)              !fusion W+ W+ into X++
    case (39)
      call trih(amp2,amp3,nwam)
    case (40)
      call fbftoh(amp2,amp3,nwam)
    case (41)
      call fhtof(amp2,amp3,nwam)
    case (45)
      call gluhtoglu(p1,p2,amp2,amp3,nwam)
      diagonal=.true.
    case (46)
      call gluglutoh(p2,amp2,p3,amp3,nwam)
    case (47)
      call hxglutoxglu(amp2,amp3,nwam)
      diagonal=.true.
    case (48)
      call xgluxglutoh(amp2,amp3,nwam)
    case (49)
      call gluyglutoglu(amp2,amp3,nwam)
    case (50)
      call gluglutoyglu(amp2,amp3,nwam)
    case (51)
      call hxglutoglu(p1,amp2,amp3,nwam)
      diagonal=.true.
    case (52)
      call gluxglutoh(p2,amp2,amp3,nwam)
    case (53)
      call gluhtoxglu(p2,amp2,amp3,nwam)
      diagonal=.true.
    case default
      write (*,*) 'WRONG lorentz operator selection',oper
      stop
  end select operator
!
end subroutine fuse
!
subroutine gluyglutoglu(amp2,amp3,nwam)
!
  use incnst, only : nlor
  implicit none
!
  complex*16, intent(in) :: amp2(nlor),amp3(nlor)
  complex*16, intent(out) :: nwam(nlor)
!
  nwam(1)=-(amp2(2)*amp3(3)-amp2(3)*amp3(2))
  nwam(2)=-(amp2(3)*amp3(1)-amp2(1)*amp3(3))
  nwam(3)=-(amp2(1)*amp3(2)-amp2(2)*amp3(1))
  nwam(4:6)=0.d0
!
end subroutine gluyglutoglu
!
subroutine gluglutoyglu(amp2,amp3,nwam)
!
  use incnst, only : nlor
  implicit none
!
  complex*16, intent(in) :: amp2(nlor),amp3(nlor)
  complex*16, intent(out) :: nwam(nlor)
!
  nwam(1)=(amp3(2)*amp2(3)-amp3(3)*amp2(2))
  nwam(2)=(amp3(3)*amp2(1)-amp3(1)*amp2(3))
  nwam(3)=(amp3(1)*amp2(2)-amp3(2)*amp2(1))
!
end subroutine gluglutoyglu
!
subroutine hxglutoglu(p1,amp2,amp3,nwam)
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in) :: p1(4)
  complex*16, intent(in) :: amp2(nlor),amp3(nlor)
  complex*16, intent(out) :: nwam(nlor)
!
  nwam(4:6)=0.d0
  nwam(1)=-(amp3(2)*p1(4)-amp3(3)*p1(3))
  nwam(2)=-(amp3(3)*p1(2)-amp3(1)*p1(4))
  nwam(3)=-(amp3(1)*p1(3)-amp3(2)*p1(2))
  nwam(1:3)=nwam(1:3)*amp2(1)
!
end subroutine hxglutoglu
!
subroutine gluxglutoh(p2,amp2,amp3,nwam)
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in) :: p2(4)
  complex*16, intent(in) :: amp2(nlor),amp3(nlor)
  complex*16, intent(out) :: nwam(nlor)
!
  nwam(1)=       -amp3(1)*(p2(3)*amp2(3)-p2(4)*amp2(2))
  nwam(1)=nwam(1)-amp3(2)*(p2(4)*amp2(1)-p2(2)*amp2(3))
  nwam(1)=nwam(1)-amp3(3)*(p2(2)*amp2(2)-p2(3)*amp2(1))
  nwam(2:3)=(0.d0,0.d0)
!
end subroutine gluxglutoh
!
subroutine gluhtoxglu(p2,amp2,amp3,nwam)
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in) :: p2(4)
  complex*16, intent(in) :: amp2(nlor),amp3(nlor)
  complex*16, intent(out) :: nwam(nlor)
!
  nwam(1)=-p2(3)*amp2(3)+p2(4)*amp2(2)
  nwam(2)=-p2(4)*amp2(1)+p2(2)*amp2(3)
  nwam(3)=-p2(2)*amp2(2)+p2(3)*amp2(1)
  nwam(1:3)=nwam(1:3)*amp3(1)
!
end subroutine gluhtoxglu
!
subroutine xgluxglutoh(amp2,amp3,nwam)
!
  use incnst, only : nlor
  implicit none
!
  complex*16, intent(in) :: amp2(nlor),amp3(nlor)
  complex*16, intent(out) :: nwam(nlor)
!
  nwam(1)=sum(amp2(1:3)*amp3(1:3))
  nwam(2:3)=(0.d0,0.d0)
!
end subroutine xgluxglutoh
!
subroutine hxglutoxglu(amp2,amp3,nwam)
!
  use incnst, only : nlor
  implicit none
!
  complex*16, intent(in) :: amp2(nlor),amp3(nlor)
  complex*16, intent(out) :: nwam(nlor)
!
  nwam(1:3)=amp3(1:3); nwam(1:3)=amp2(1)*nwam(1:3)
!
end subroutine hxglutoxglu
!
subroutine gluhtoglu (p1,p2,amp2,amp3,nwam)
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in), dimension (4) :: p1,p2
  complex*16, intent(out), dimension (nlor) :: nwam
!
  complex*16 tmp
!
  tmp=p1(1)*p2(1)-sum(p1(2:4)*p2(2:4))
  nwam(1:3)=-tmp*amp2(1:3)
  tmp=-sum(p1(2:4)*amp2(1:3))
  nwam(1:3)=nwam(1:3)+tmp*p2(2:4)
  nwam(1:3)=nwam(1:3)*amp3(1)
  nwam(4:6)=0.d0
!
end subroutine gluhtoglu
!
subroutine gluglutoh (p2,amp2,p3,amp3,nwam) 
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in), dimension (4) :: p3,p2
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=(p2(1)*p3(1)-sum(p2(2:4)*p3(2:4)))*(-sum(amp2(1:3)*amp3(1:3)))
  nwam(1)=nwam(1)-sum(p2(2:4)*amp3(1:3))*sum(amp2(1:3)*p3(2:4))
  nwam(2:3)=(0.d0,0.d0)
!
end subroutine gluglutoh
!
subroutine triglu (p1,p2,amp2,p3,amp3,nwam)     !trilinear gluon interaction   !
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in), dimension (3) :: p1,p2,p3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  complex*16 :: a2a3,a2p,a3p
!
  integer j1
!
  a2a3=sum(amp2(1:3)*amp3(1:3))
  a2p=sum(amp2(1:3)*(p3-p1))
  a3p=sum(amp3(1:3)*(p1-p2))
  nwam(1:3)=a2a3*(p2-p3)+a2p*amp3(1:3)+a3p*amp2(1:3)
  nwam(1)=nwam(1)+amp3(5)*amp2(3)-amp3(6)*amp2(2)-amp2(5)*amp3(3)+amp2(6)*amp3(2)
  nwam(2)=nwam(2)+amp3(6)*amp2(1)-amp3(4)*amp2(3)-amp2(6)*amp3(1)+amp2(4)*amp3(3)
  nwam(3)=nwam(3)+amp3(4)*amp2(2)-amp3(5)*amp2(1) -amp2(4)*amp3(2)+amp2(5)*amp3(1)
  nwam(4)=-amp2(2)*amp3(3)+amp2(3)*amp3(2)
  nwam(5)=-amp2(3)*amp3(1)+amp2(1)*amp3(3)
  nwam(6)=-amp2(1)*amp3(2)+amp2(2)*amp3(1)
!
end subroutine triglu 
!
subroutine triphw (qrt,p1,p2,amp2,p3,amp3,nwam)   !trilinear WWZ interaction   
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (5) :: amp2,amp3
  real*8, intent(in), dimension (4) :: p1,p2,p3
  real*8, intent(in), dimension (3) :: qrt
  complex*16, intent(out), dimension (4) :: nwam
!
  complex*16 :: a2a3,a2p,a3p,ax2,ax3
!
  integer j1
!
  a2a3=amp2(1)*amp3(1)-sum(amp2(2:4)*amp3(2:4))
  a2p=amp2(1)*(p3(1)-p1(1))-sum(amp2(2:4)*(p3(2:4)-p1(2:4)))
  a3p=amp3(1)*(p1(1)-p2(1))-sum(amp3(2:4)*(p1(2:4)-p2(2:4)))
  ax3=qrt(3)*amp3(5);   ax2=qrt(2)*amp2(5); 
  nwam(1:3)=a2a3*(p2(2:4)-p3(2:4))+(a2p+ax2)*amp3(2:4)+(a3p+ax3)*amp2(2:4)
  nwam(4)=qrt(1)*a2a3
!  nwam(1:4)=(ax2)*amp3(1:4)+(ax3)*amp2(1:4)
!  nwam(1)=-nwam(1)
!  nwam(5)=qrt(1)*a2a3
!
end subroutine triphw 
!
subroutine triw (qrt,p1,p2,amp2,p3,amp3,nwam)   !trilinear WWZ interaction   
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (5) :: amp2,amp3
  real*8, intent(in), dimension (4) :: p1,p2,p3
  real*8, intent(in), dimension (3) :: qrt
  complex*16, intent(out), dimension (5) :: nwam
!
  complex*16 :: a2a3,a2p,a3p,ax2,ax3
!
  integer j1
!
  a2a3=amp2(1)*amp3(1)-sum(amp2(2:4)*amp3(2:4))
  a2p=amp2(1)*(p3(1)-p1(1))-sum(amp2(2:4)*(p3(2:4)-p1(2:4)))
  a3p=amp3(1)*(p1(1)-p2(1))-sum(amp3(2:4)*(p1(2:4)-p2(2:4)))
  ax3=qrt(3)*amp3(5);   ax2=qrt(2)*amp2(5); 
  nwam(1:4)=a2a3*(p2(1:4)-p3(1:4))+(a2p+ax2)*amp3(1:4)+(a3p+ax3)*amp2(1:4)
  nwam(1)=-nwam(1)
  nwam(5)=qrt(1)*a2a3
!  nwam(1:4)=(ax2)*amp3(1:4)+(ax3)*amp2(1:4)
!  nwam(1)=-nwam(1)
!  nwam(5)=qrt(1)*a2a3
!
end subroutine triw 
!
subroutine triwph (qrt,p1,p2,amp2,p3,amp3,nwam)   !trilinear WWZ interaction   
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (5) :: amp2
  complex*16, intent(in), dimension (4) :: amp3
  real*8, intent(in), dimension (4) :: p1,p2,p3
  real*8, intent(in), dimension (3) :: qrt
  complex*16, intent(out), dimension (5) :: nwam
!
  complex*16 :: a2a3,a2p,a3p,ax2,ax3
!
  integer j1
!
  a2a3=-sum(amp2(2:4)*amp3(1:3))
  a2p=amp2(1)*(p3(1)-p1(1))-sum(amp2(2:4)*(p3(2:4)-p1(2:4)))
  a3p=-sum(amp3(1:3)*(p1(2:4)-p2(2:4)))
  ax3=qrt(3)*amp3(4);   ax2=qrt(2)*amp2(5); 
  nwam(1)=a2a3*(p2(1)-p3(1))+(a3p+ax3)*amp2(1)
  nwam(2:4)=a2a3*(p2(2:4)-p3(2:4))+(a2p+ax2)*amp3(1:3)+(a3p+ax3)*amp2(2:4)
  nwam(1)=-nwam(1)
  nwam(5)=qrt(1)*a2a3
!  nwam(1:4)=(ax2)*amp3(1:4)+(ax3)*amp2(1:4)
!  nwam(1)=-nwam(1)
!  nwam(5)=qrt(1)*a2a3
!
end subroutine triwph 
!
subroutine waxtow (amp2,amp3,nwam)   !fusion of W+ Y++ into W+
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1:4)=amp3(1)*amp2(1:4)
  nwam(2:4)=-nwam(2:4)
  nwam(5)=0.d0
!
end subroutine waxtow
!
subroutine pax2tow (qrt,amp2,amp3,nwam)   !fusion of W+ Y++ into W+
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in) :: qrt(2)
  complex*16, intent(out), dimension (nlor) :: nwam
!
  complex*16 :: tmp
!
  tmp=qrt(1)*amp3(1)+qrt(2)*amp3(2)
  nwam(2:4)=-tmp*amp2(1:3)
!  nwam(1)=-nwam(1)
  nwam(1)=0.d0
  nwam(5)=0.d0
!
end subroutine pax2tow
!
subroutine pax2top (qrt,amp2,amp3,nwam)   !fusion of W+ Y++ into W+
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in) :: qrt(2)
  complex*16, intent(out), dimension (nlor) :: nwam
!
  complex*16 :: tmp
!
  tmp=qrt(1)*amp3(1)+qrt(2)*amp3(2)
  nwam(1:3)=-tmp*amp2(1:3)
!  nwam(1)=-nwam(1)
  nwam(4)=0.d0
!
end subroutine pax2top
!
subroutine wax2top (qrt,amp2,amp3,nwam)   !fusion of W+ Y++ into W+
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in) :: qrt(2)
  complex*16, intent(out), dimension (nlor) :: nwam
!
  complex*16 :: tmp
!
  tmp=qrt(1)*amp3(1)+qrt(2)*amp3(2)
  nwam(1:3)=-tmp*amp2(2:4)
  nwam(4)=0.d0
!
end subroutine wax2top
!
subroutine wax2tow (qrt,amp2,amp3,nwam)   !fusion of W+ Y++ into W+
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in) :: qrt(2)
  complex*16, intent(out), dimension (nlor) :: nwam
!
  complex*16 :: tmp
!
  tmp=qrt(1)*amp3(1)+qrt(2)*amp3(2)
  nwam(1:4)=tmp*amp2(1:4)
!  nwam(1)=-nwam(1)
  nwam(2:4)=-nwam(2:4)
  nwam(5)=0.d0
!
end subroutine wax2tow
!
subroutine wwtoax (amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=amp3(1)*amp2(1)-sum(amp2(2:4)*amp3(2:4))
!
end subroutine wwtoax
!
subroutine wwtoh (cp,amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  use incnst
  implicit none
!
  real*8, intent(in) :: cp
  complex*16, intent(in), dimension (4) :: amp2,amp3
  complex*16, intent(out), dimension (3) :: nwam
!
  nwam(1)=amp3(1)*amp2(1)-sum(amp2(2:4)*amp3(2:4))
  nwam(2)=cp*nwam(1)
  nwam(3)=0.d0
!
end subroutine wwtoh
!
subroutine whtow (cp,amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  use incnst
  implicit none
!
  real*8, intent(in) :: cp
  complex*16, intent(in) :: amp2(4),amp3(3)
  complex*16, intent(out), dimension (5) :: nwam
!
  complex*16 :: ax
!
  ax=amp3(1)+cp*amp3(2)
  nwam(1)=ax*amp2(1)
  nwam(2:4)=-ax*amp2(2:4)
  nwam(5)=0.d0
!
end subroutine whtow
!
subroutine trih (amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  use couplings, only : qrtchh
  use incnst
  implicit none
!
  complex*16, intent(in) :: amp2(3),amp3(3)
  complex*16, intent(out), dimension (3) :: nwam
!
  complex*16 :: ax
!
  ax=amp3(1)*amp2(1); nwam=ax*qrtchh
! nwam(1)=ax*qrtchh(1); nwam(2)=ax*qrtchh(3);  nwam(3)=ax*qrtchh(2);
  nwam(1)=nwam(1)+amp2(1)*(qrtchh(2)*amp3(2)+qrtchh(3)*amp3(3)) &
                 +amp3(1)*(qrtchh(2)*amp2(2)+qrtchh(3)*amp2(3))
!
end subroutine trih
!
subroutine fbftoh (amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  implicit none
!
  complex*16, intent(in) :: amp2(4),amp3(4)
  complex*16, intent(out), dimension (3) :: nwam
!
  nwam(1)=sum(amp2*amp3); nwam(2:3)=(0.d0,0.d0)
!
end subroutine fbftoh
!
subroutine fhtof (amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  implicit none
!
  complex*16, intent(in) :: amp2(4),amp3(3)
  complex*16, intent(out), dimension (4) :: nwam
!
  nwam=amp3(1)*amp2
!
end subroutine fhtof
!
subroutine wwtoax2 (qrt,amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in) :: qrt(2)
  complex*16, intent(out), dimension (2) :: nwam
!
  complex*16 :: tmp
!
  tmp=amp3(1)*amp2(1)-sum(amp2(2:4)*amp3(2:4))
  nwam=tmp*qrt
!
end subroutine wwtoax2
!
subroutine zptoax2 (qrt,amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in) :: qrt(2)
  complex*16, intent(out), dimension (2) :: nwam
!
  complex*16 :: tmp
!
  tmp=-sum(amp2(2:4)*amp3(1:3))
  nwam=tmp*qrt
!
end subroutine zptoax2
!
subroutine pptoax2 (qrt,amp2,amp3,nwam)   !fusion of W+ W+ into X++
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  real*8, intent(in) :: qrt(2)
  complex*16, intent(out), dimension (2) :: nwam
!
  complex*16 :: tmp
!
  tmp=-sum(amp2(1:3)*amp3(1:3))
  nwam=tmp*qrt
!
end subroutine pptoax2
!
subroutine fbftog (amp2,amp3,nwam)                       !fusion of fermion and fermion-bar into gluon
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=amp2(1)*amp3(4)+amp2(2)*amp3(3)-amp2(3)*amp3(2)-amp2(4)*amp3(1)
  nwam(2)=-amp2(1)*amp3(4)+amp2(2)*amp3(3)+amp2(3)*amp3(2)-amp2(4)*amp3(1)
  nwam(2)=(0.d0,1.d0)*nwam(2)
  nwam(3)=amp2(1)*amp3(3)-amp2(2)*amp3(4)-amp2(3)*amp3(1)+amp2(4)*amp3(2)
  nwam(4:6)=0
!
end subroutine fbftog
!
subroutine fbftow (amp2,amp3,nwam)                       !fusion of fermion and fermion-bar into gluon
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
  complex*16 :: c11,c12,c21,c22
!
  c11=amp2(3)*amp3(1)
  c12=amp2(3)*amp3(2)
  c21=amp2(4)*amp3(1)
  c22=amp2(4)*amp3(2)
!
  nwam(1)=c11+c22
  nwam(2)=c12+c21
  nwam(3)=(0.d0,1.d0)*(c21-c12)
  nwam(4)=c11-c22
  nwam(5)=(0.d0,0.d0)
!
end subroutine fbftow
!
subroutine fbftogr (amp2,amp3,nwam)                       !fusion of fermion and fermion-bar into gluon
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=amp2(1)*amp3(2)+amp2(2)*amp3(1)
  nwam(2)=-amp2(1)*amp3(2)+amp2(2)*amp3(1)
  nwam(2)=(0.d0,1.d0)*nwam(2)
  nwam(3)=amp2(1)*amp3(1)-amp2(2)*amp3(2)
  nwam(4:6)=0
!
end subroutine fbftogr
!
subroutine fbftogl (amp2,amp3,nwam)                       !fusion of fermion and fermion-bar into gluon
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=-amp2(1)*amp3(2)-amp2(2)*amp3(1)
  nwam(2)=amp2(1)*amp3(2)-amp2(2)*amp3(1)
  nwam(2)=(0.d0,1.d0)*nwam(2)
  nwam(3)=-amp2(1)*amp3(1)+amp2(2)*amp3(2)
  nwam(4:6)=0
!
end subroutine fbftogl
!
subroutine fbftowl (amp2,amp3,nwam)                       !fusion of fermion and fermion-bar into gluon
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
  complex*16 :: c11,c12,c21,c22
!
  c11=amp2(1)*amp3(1)
  c12=amp2(1)*amp3(2)
  c21=amp2(2)*amp3(1)
  c22=amp2(2)*amp3(2)
!
  nwam(1)=c11+c22
  nwam(2)=(c12+c21)
  nwam(3)=(0.d0,1.d0)*(c21-c12)
  nwam(4)=c11-c22
  nwam(5)=(0.d0,0.d0)
!
end subroutine fbftowl
!
subroutine fbftowr (amp2,amp3,nwam)                       !fusion of fermion and fermion-bar into gluon
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
  complex*16 :: c11,c12,c21,c22
!
  c11=amp2(1)*amp3(1)
  c12=amp2(1)*amp3(2)
  c21=amp2(2)*amp3(1)
  c22=amp2(2)*amp3(2)
!
  nwam(1)=c11+c22
  nwam(2)=-(c12+c21)
  nwam(3)=(0.d0,-1.d0)*(c21-c12)
  nwam(4)=c22-c11
  nwam(5)=(0.d0,0.d0)
!
end subroutine fbftowr
!
subroutine glgr(fl,gl,gr)
!
  use incnst
  use naming
  use couplings, only :  gldn, gllp, grup, grdn, grlp, glup
  implicit none
!
  integer, intent(in) :: fl
  real*8, intent (out) :: gl,gr
!
  vma : select case (conv(fl))
    case (up)
      gl = glup; gr = grup
    case (dn)
      gl = gldn; gr = grdn
    case (lpt)
      gl = gllp; gr = grlp
    case default
      write(*,*) 'something wrong in vma neutral interaction'
      stop
  end select vma
!
end subroutine glgr
!
subroutine fbftoz (flv,amp2,amp3,nwam)                       !fusion of fermion and fermion-bar into gluon
!
  use incnst
  implicit none
!
  integer, intent(in) :: flv
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
  complex*16 :: c11,c12,c21,c22
  real*8 :: gl,gr
!
  call glgr(flv,gl,gr)
!
  c11=amp2(3)*amp3(1)
  c12=amp2(3)*amp3(2)
  c21=amp2(4)*amp3(1)
  c22=amp2(4)*amp3(2)
!
  nwam(1)=c11+c22
  nwam(2)=(c12+c21)
  nwam(3)=(0.d0,1.d0)*(c21-c12)
  nwam(4)=c11-c22
  nwam=gl*nwam
!
  c11=amp2(1)*amp3(3)
  c12=amp2(1)*amp3(4)
  c21=amp2(2)*amp3(3)
  c22=amp2(2)*amp3(4)
!
  nwam(1)=nwam(1)+gr*(c11+c22)
  nwam(2)=nwam(2)-gr*(c12+c21)
  nwam(3)=nwam(3)-gr*(0.d0,1.d0)*(c21-c12)
  nwam(4)=nwam(4)-gr*(c11-c22)
  nwam(5)=0.d0
!
end subroutine fbftoz
!
subroutine gfbtofb (amp2,amp3,nwam)                       !fusion of fermionbar and gluon into fermionbar
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=-amp3(3)*amp2(3)-amp3(4)*(amp2(1)+(0.d0,1.d0)*amp2(2))
  nwam(2)=-amp3(3)*(amp2(1)-(0.d0,1.d0)*amp2(2))+amp3(4)*amp2(3)
  nwam(3)=amp3(1)*amp2(3)+amp3(2)*(amp2(1)+(0.d0,1.d0)*amp2(2))
  nwam(4)=amp3(1)*(amp2(1)-(0.d0,1.d0)*amp2(2))-amp3(2)*amp2(3)
!
end subroutine gfbtofb                                   !fusion of fermion and gluon into fermion
!
subroutine wfbtofb (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=(amp2(1)+amp2(4))*amp3(3)+(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(4)
  nwam(2)=(amp2(1)-amp2(4))*amp3(4)+(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(3)
  nwam(3)=0.d0
  nwam(4)=0.d0
!
end subroutine wfbtofb
!
subroutine gfbtofbl (amp2,amp3,nwam)                       !fusion of fermionbar and gluon into fermionbar
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=-amp3(1)*amp2(3)-amp3(2)*(amp2(1)+(0.d0,1.d0)*amp2(2))
  nwam(2)=-amp3(1)*(amp2(1)-(0.d0,1.d0)*amp2(2))+amp3(2)*amp2(3)
!
end subroutine gfbtofbl                                   !fusion of fermion and gluon into fermion
!
subroutine gfbtofbr (amp2,amp3,nwam)                       !fusion of fermionbar and gluon into fermionbar
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=amp3(1)*amp2(3)+amp3(2)*(amp2(1)+(0.d0,1.d0)*amp2(2))
  nwam(2)=amp3(1)*(amp2(1)-(0.d0,1.d0)*amp2(2))-amp3(2)*amp2(3)
!
end subroutine gfbtofbr                                   !fusion of fermion and gluon into fermion
!
subroutine wfbtofbl (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=(amp2(1)+amp2(4))*amp3(1)+(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(2)
  nwam(2)=(amp2(1)-amp2(4))*amp3(2)+(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(1)
!
end subroutine wfbtofbl
!
subroutine wfbtofbr (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=(amp2(1)-amp2(4))*amp3(1)-(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(2)
  nwam(2)=(amp2(1)+amp2(4))*amp3(2)-(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(1)
!
end subroutine wfbtofbr
!
subroutine zfbtofb (flv,amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  integer, intent(in) :: flv
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
  real*8 :: gl,gr
!
  call glgr(flv,gl,gr)
!
  nwam(3)=(amp2(1)-amp2(4))*amp3(1)-(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(2)
  nwam(4)=(amp2(1)+amp2(4))*amp3(2)-(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(1)
  nwam(3:4)=gr*nwam(3:4)
!
  nwam(1)=(amp2(1)+amp2(4))*amp3(3)+(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(4)
  nwam(2)=(amp2(1)-amp2(4))*amp3(4)+(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(3)
  nwam(1:2)=gl*nwam(1:2)
!
end subroutine zfbtofb
!
subroutine gftof (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=amp3(3)*amp2(3)+amp3(4)*(amp2(1)-(0.d0,1.d0)*amp2(2))
  nwam(2)=amp3(3)*(amp2(1)+(0.d0,1.d0)*amp2(2))-amp3(4)*amp2(3)
  nwam(3)=-amp3(1)*amp2(3)-amp3(2)*(amp2(1)-(0.d0,1.d0)*amp2(2))
  nwam(4)=-amp3(1)*(amp2(1)+(0.d0,1.d0)*amp2(2))+amp3(2)*amp2(3)
!
end subroutine gftof
!
subroutine wftof (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=0.d0
  nwam(2)=0.d0
  nwam(3)=(amp2(1)+amp2(4))*amp3(1)+(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(2)
  nwam(4)=(amp2(1)-amp2(4))*amp3(2)+(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(1)
!
end subroutine wftof
!
subroutine gftofl (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=-amp3(1)*amp2(3)-amp3(2)*(amp2(1)-(0.d0,1.d0)*amp2(2))
  nwam(2)=-amp3(1)*(amp2(1)+(0.d0,1.d0)*amp2(2))+amp3(2)*amp2(3)
!
end subroutine gftofl
!
subroutine gftofr (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in), dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=amp3(1)*amp2(3)+amp3(2)*(amp2(1)-(0.d0,1.d0)*amp2(2))
  nwam(2)=amp3(1)*(amp2(1)+(0.d0,1.d0)*amp2(2))-amp3(2)*amp2(3)
!
end subroutine gftofr
!
subroutine wftofl (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in),  dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=(amp2(1)+amp2(4))*amp3(1)+(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(2)
  nwam(2)=(amp2(1)-amp2(4))*amp3(2)+(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(1)
!
end subroutine wftofl
!
subroutine wftofr (amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  complex*16, intent(in),  dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
!
  nwam(1)=(amp2(1)-amp2(4))*amp3(1)-(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(2)
  nwam(2)=(amp2(1)+amp2(4))*amp3(2)-(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(1)
!
end subroutine wftofr
!
subroutine zftof (flv,amp2,amp3,nwam)                            
!
  use incnst
  implicit none
!
  integer, intent(in) :: flv
  complex*16, intent(in),  dimension (nlor) :: amp2,amp3
  complex*16, intent(out), dimension (nlor) :: nwam
  real*8 :: gl,gr
!
  call glgr(flv,gl,gr)
!
  nwam(3)=(amp2(1)+amp2(4))*amp3(1)+(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(2)
  nwam(4)=(amp2(1)-amp2(4))*amp3(2)+(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(1)
  nwam(3:4)=gl*nwam(3:4)
!
  nwam(1)=(amp2(1)-amp2(4))*amp3(3)-(amp2(2)-(0.d0,1.d0)*amp2(3))*amp3(4)
  nwam(2)=(amp2(1)+amp2(4))*amp3(4)-(amp2(2)+(0.d0,1.d0)*amp2(3))*amp3(3)
  nwam(1:2)=gr*nwam(1:2)
!
end subroutine zftof
!
subroutine prpgt (fl,imp,fus,pol)                         
!
  use naming
  use incnst, only : nlor                             !number of lorentz degrees of freedom
  use couplings, only : masses,width                        !particle mass
  use strfld, only : lortyp                           !labelling lorentz representation
  implicit none
!
  integer, intent(in) :: fl                           !flavour
  real*8, intent(in), dimension (4) :: imp            !momenta
  complex*16, intent(in), dimension(nlor) :: fus      !amplitude prior to matching with the propagator 
  complex*16, intent(out), dimension(nlor) :: pol     !propagator times fusion
!
  propagator: select case (lortyp(fl))
    case(glu)                                       !massless gauge boson, F_munu
      call gbzms (imp,fus,pol)
    case(pht)                                       !massless gauge boson, F_munu
      call gbph (imp,fus,pol)
    case(msgbs)                                       !massless gauge boson
      call gbms (imp,fus,masses(fl),width(fl),pol)
    case(frm)                                       !fermion
      call psl (imp,fus,masses(fl),pol)
    case(frmbr)                                       !fermion-bar
      call psltr (imp,fus,masses(fl),pol)
    case(frmhl)                                       !fermion
      call psll (imp,fus,pol)
    case(frmbrhl)                                       !fermion-bar
      call psltrl (imp,fus,pol)
    case(frmhr)                                         !fermion
      call pslr (imp,fus,pol)
    case(frmbrhr)                                       !fermion-bar
      call psltrr (imp,fus,pol)
    case(axww)                                          !aux X++
      pol(1)=fus(1)
    case(axzzww)                                          !aux X++
      pol(1:2)=fus(1:2)
    case(scal)
      call scl(imp,fus,masses(fl),width(fl),pol)
    case(xgl)
      pol(1:6)=fus(1:6)
    case default
      write(*,*)'Wrong choice of propagator'
      stop
  end select propagator
!
end subroutine prpgt
!
subroutine gbzms (imp,fus,pol)                   !massless gauge boson propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in), dimension (4) :: imp
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 pq,ip0
  complex*16 pa
!
  logical flg90
  common/f90/flg90           ! to avoid spurious divergences in temporal gauge
!
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
  pq=1/pq
  ip0=imp(1)*imp(1)
  if(ip0.eq.0.d0) then
    write(*,*)'singular coulomb propagator'
    stop
  endif
  pa=sum(imp(2:)*fus(1:3))
  pa=pa/ip0
  pol(1:3)=(fus(1:3)-imp(2:4)*pa)*pq
  pol(4:6)=-fus(4:6)
!
  if(.not.flg90) then
    ip0=sum(abs(imp(2:4)))
    if(ip0.ne.0.d0.and.(abs(imp(1)/ip0).lt.1.d-4)) flg90=.true.
  endif
!
end subroutine gbzms
!
subroutine gbph (imp,fus,pol)                   !massless gauge boson propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in), dimension (4) :: imp
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 pq,ip0
  complex*16 pa
!
  logical flg90
  common/f90/flg90           ! to avoid spurious divergences in temporal gauge
!
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
  pq=1/pq
  ip0=imp(1)*imp(1)
  if(ip0.eq.0.d0) then
    write(*,*)'singular coulomb propagator'
    stop
  endif
  pa=sum(imp(2:)*fus(1:3))
  pa=pa/ip0
  pol(1:3)=(fus(1:3)-imp(2:4)*pa)*pq
  pol(4)=fus(4)
!
  if(.not.flg90) then
    ip0=sum(abs(imp(2:4)))
    if(ip0.ne.0.d0.and.(abs(imp(1)/ip0).lt.1.d-4)) flg90=.true.
  endif
!
end subroutine gbph
!
subroutine gbms (imp,fus,mass,width,pol)                   !massless gauge boson propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in) :: imp(4),mass,width
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 :: mq,theta
  complex*16 :: pa,pq
!
  character*2 wmode
  character resonance
  real*8 winsize
  common/gauinv/winsize,resonance,wmode
!
  if (resonance.eq.'y') then; pol=(0.d0,0.d0); return; endif
!
  mq=mass*mass
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))-mq             ! +(0.d0,1.d0)*mass*width
  if(wmode.eq.'yy') then
    pq=1.d0/(pq+(0.d0,1.d0)*mass*width); resonance='n'
  elseif(wmode.eq.'yn') then
    if(abs(pq).lt.winsize*mass*width) then; resonance='y'; else; resonance='n'; endif
    if(pq.eq.0.d0) pq=1.d0; pq=1/pq
  elseif(wmode.eq.'nn') then
    pq=1/pq
  endif
!!$  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
!!$  if(real(pq).gt.0) then; theta=1.d0; else; theta=0.d0; endif 
!!$  pq=pq-mq+theta*(0.d0,1.d0)*pq*width/mass
  pa=sum(imp(1:4)*fus(1:4))
  pa=pa/mq
  pol(1:4)=-fus(1)+imp(1)*pa
  pol(2:4)=fus(2:4)+imp(2:4)*pa
  pol(1:4)=pol(1:4)*pq
!  pol(5)=fus(5)
!  pol(1)=-pol(1)
!
end subroutine gbms
!
subroutine scl (imp,fus,mass,width,pol)                   !massless gauge boson propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in) :: imp(4),mass,width
  complex*16, intent(in), dimension(3) :: fus
  complex*16, intent(out), dimension(3) :: pol
!
  real*8 :: mq,theta
  complex*16 :: pa,pq
!
  character*2 wmode
  character resonance
  real*8 winsize
  common/gauinv/winsize,resonance,wmode
!
  if (resonance.eq.'y') then; pol=(0.d0,0.d0); return; endif
!
  mq=mass*mass
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))-mq
  if(wmode.eq.'yy') then
    pq=1.d0/(pq+(0.d0,1.d0)*mass*width); resonance='n'
  elseif(wmode.eq.'yn') then
    if(abs(pq).lt.winsize*mass*width) then; resonance='y'; else; resonance='n'; endif
    if(pq.eq.0.d0) pq=1.d0; pq=1/pq
  elseif(wmode.eq.'nn') then
    pq=1/pq
  endif

!  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))-mq+(0.d0,1.d0)*mass*width
!!$  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
!!$  if(real(pq).gt.0) then; theta=1.d0; else; theta=0.d0; endif 
!!$  pq=pq-mq+theta*(0.d0,1.d0)*pq*width/mass
!  pq=1.d0/pq
  pol(1)=pq*fus(1); pol(2)=fus(3); pol(3)=fus(2);
!
end subroutine scl
!
subroutine psl(imp,fus,mass,pol)                !fermion propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in), dimension (4) :: imp
  real*8, intent(in) :: mass
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 pq
!
  pol(1)=fus(1)*mass+fus(3)*(imp(1)-imp(4))-fus(4)*(imp(2)-(0.d0,1.d0)*imp(3))
  pol(2)=fus(2)*mass-fus(3)*(imp(2)+(0.d0,1.d0)*imp(3))+fus(4)*(imp(1)+imp(4))
  pol(3)=fus(1)*(imp(4)+imp(1))+fus(2)*(imp(2)-(0.d0,1.d0)*imp(3))+fus(3)*mass
  pol(4)=fus(1)*(imp(2)+(0.d0,1.d0)*imp(3))+fus(2)*(imp(1)-imp(4))+fus(4)*mass
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
  pq=pq-mass**2
  pq=1/pq
  pol(1:4)=pol(1:4)*pq
!
end subroutine psl
!
subroutine pslr(imp,fus,pol)                !fermion propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in), dimension (4) :: imp
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 pq
!
  pol(1)=fus(1)*(imp(4)+imp(1))+fus(2)*(imp(2)-(0.d0,1.d0)*imp(3))
  pol(2)=fus(1)*(imp(2)+(0.d0,1.d0)*imp(3))+fus(2)*(imp(1)-imp(4))
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
  pq=1/pq
  pol(1:2)=pol(1:2)*pq
!
end subroutine pslr
!
subroutine psll(imp,fus,pol)                !fermion propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in), dimension (4) :: imp
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 pq
!
  pol(1)=fus(1)*(imp(1)-imp(4))-fus(2)*(imp(2)-(0.d0,1.d0)*imp(3))
  pol(2)=-fus(1)*(imp(2)+(0.d0,1.d0)*imp(3))+fus(2)*(imp(1)+imp(4))
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
  pq=1/pq
  pol(1:2)=pol(1:2)*pq
!
end subroutine psll
!
subroutine psltr(imp,fus,mass,pol)                  !fermion bar propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in), dimension (4) :: imp
  real*8, intent(in) :: mass
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 pq
!
  pol(1)=fus(1)*mass-fus(3)*(imp(4)+imp(1))-fus(4)*(imp(2)+(0.d0,1.d0)*imp(3))
  pol(2)=fus(2)*mass-fus(3)*(imp(2)-(0.d0,1.d0)*imp(3))+fus(4)*(imp(4)-imp(1))
  pol(3)=fus(1)*(imp(4)-imp(1))+fus(2)*(imp(2)+(0.d0,1.d0)*imp(3))+fus(3)*mass
  pol(4)=fus(1)*(imp(2)-(0.d0,1.d0)*imp(3))-fus(2)*(imp(1)+imp(4))+fus(4)*mass
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
  pq=pq-mass**2
  pq=1/pq
  pol(1:4)=pol(1:4)*pq
!
end subroutine psltr
!
subroutine psltrr(imp,fus,pol)                  !fermion bar propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in), dimension (4) :: imp
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 pq
!
  pol(1)=-fus(1)*(imp(4)+imp(1))-fus(2)*(imp(2)+(0.d0,1.d0)*imp(3))
  pol(2)=-fus(1)*(imp(2)-(0.d0,1.d0)*imp(3))+fus(2)*(imp(4)-imp(1))
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
  pq=1/pq
  pol(1:2)=pol(1:2)*pq
!
end subroutine psltrr
!
subroutine psltrl(imp,fus,pol)                  !fermion bar propagator
!
  use incnst, only : nlor
  implicit none
!
  real*8, intent(in), dimension (4) :: imp
  complex*16, intent(in), dimension(nlor) :: fus
  complex*16, intent(out), dimension(nlor) :: pol
!
  real*8 pq
!
  pol(1)=fus(1)*(imp(4)-imp(1))+fus(2)*(imp(2)+(0.d0,1.d0)*imp(3))
  pol(2)=fus(1)*(imp(2)-(0.d0,1.d0)*imp(3))-fus(2)*(imp(1)+imp(4))
  pq=imp(1)*imp(1)-sum(imp(2:4)*imp(2:4))
  pq=1/pq
  pol(1:2)=pol(1:2)*pq
!
end subroutine psltrl
!
subroutine colprdsu3
!
  use incol                                                          
  implicit none                                                     
  integer j1,j2,j3,i2,k1,k2,i1
  integer*4, dimension(2,-11:11,-11:11) :: tabmom     !tabmom(i,j,k) = coulor(j) * coulor(k) (i=1,2 because of
!                                                      neutral gluons 
  integer*4, dimension(-11:11,-11:11) :: tabnew       ! number of coulors of the new object, == 0,1,2
  real*8,  dimension(2,-11:11,-11:11) :: tabcoeff     !coulor clebsch
!
  tabnew=0
  tabmom=0
  tabcoeff=0
!                                                                               
  tabnew(1,1)=0                                                
!                                                                               
  tabnew(1,3)=1; tabcoeff(1,1,3)=1./sqrt(2.d0); tabmom(1,1,3)=3
!                                                                               
  tabnew(1,4)=1; tabcoeff(1,1,4)=sqrt(2.d0);  tabmom(1,1,4)=4 
!                                                                               
  tabnew(1,6)=1; tabcoeff(1,1,6)=-1./sqrt(2.d0); tabmom(1,1,6)=6
!                                                                               
  tabnew(1,2)=0                                                
!                                                                               
  tabnew(1,5)=1; tabcoeff(1,1,5)=1./sqrt(2.d0); tabmom(1,1,5)=5
!                                                                               
  tabnew(1,7)=1; tabcoeff(1,1,7)=-sqrt(2.d0); tabmom(1,1,7)=7                                              
!                                                                               
  tabnew(1,8)=1; tabcoeff(1,1,8)=-1./sqrt(2.d0); tabmom(1,1,8)=8                                              
!                                                                               
  tabnew(3,1)=1; tabcoeff(1,3,1)=-1./sqrt(2.d0); tabmom(1,3,1)=3                                              
!                                                                               
  tabnew(3,3)=0                                                
!                                                                               
  tabnew(3,4)=0                                                
!                                                                               
  tabnew(3,6)=2; tabcoeff(1:2,3,6)=(/1./sqrt(2.d0),-sqrt(3.d0)/sqrt(2.d0)/); tabmom(1:2,3,6)=(/1,2/)
!                                                                               
  tabnew(3,2)=1; tabcoeff(1,3,2)=sqrt(3.d0)/sqrt(2.d0); tabmom(1,3,2)=3                                              
!                                                                               
  tabnew(3,5)=1; tabcoeff(1,3,5)=1.; tabmom(1,3,5)=4                                              
!                                                                               
  tabnew(3,7)=1; tabcoeff(1,3,7)=-1.; tabmom(1,3,7)=8                                              
!                                                                               
  tabnew(3,8)=0                                                
!                                                                               
  tabnew(4,1)=1; tabcoeff(1,4,1)=-sqrt(2.d0); tabmom(1,4,1)=4                                              
!                                                                               
  tabnew(4,3)=0                                                
!                                                                               
  tabnew(4,4)=0                                                
!                                                                               
  tabnew(4,6)=1; tabcoeff(1,4,6)=-1.; tabmom(1,4,6)=5                                              
!                                                                               
  tabnew(4,2)=0                                                
!                                                                               
  tabnew(4,5)=0                                                 
!                                                                               
  tabnew(4,7)=2; tabcoeff(1:2,4,7)=(/sqrt(2.d0),0.d0/); tabmom(1:2,4,7)=(/1,2/)
!                                                                               
  tabnew(4,8)=1; tabcoeff(1,4,8)=1.; tabmom(1,4,8)=3                                              
!                                                                               
  tabnew(6,1)=1; tabcoeff(1,6,1)=1./sqrt(2.d0); tabmom(1,6,1)=6                                              
!                                                                               
  tabnew(6,3)=2; tabcoeff(1:2,6,3)=(/-1./sqrt(2.d0),sqrt(3.d0)/sqrt(2.d0)/); tabmom(1:2,6,3)=(/1,2/)
!                                                                               
  tabnew(6,4)=1; tabcoeff(1,6,4)=1.; tabmom(1,6,4)=5                                              
!                                                                               
  tabnew(6,6)=0                                                
!                                                                               
  tabnew(6,2)=1; tabcoeff(1,6,2)=-sqrt(3.d0)/sqrt(2.d0); tabmom(1,6,2)=6                                              
!                                                                               
  tabnew(6,5)=0                                                
!                                                                               
  tabnew(6,7)=0                                                
!                                                                               
  tabnew(6,8)=1; tabcoeff(1,6,8)=-1.; tabmom(1,6,8)=7                                              
!                                                                               
  tabnew(2,1)=0                                                
!                                                                               
  tabnew(2,3)=1; tabcoeff(1,2,3)=-sqrt(3.d0)/sqrt(2.d0); tabmom(1,2,3)=3                                              
!                                                                               
  tabnew(2,4)=0                                                
!                                                                               
  tabnew(2,6)=1; tabcoeff(1,2,6)=sqrt(3.d0)/sqrt(2.d0); tabmom(1,2,6)=6                                              
!                                                                               
  tabnew(2,2)=0                                                
!                                                                               
  tabnew(2,5)=1; tabcoeff(1,2,5)=sqrt(3.d0)/sqrt(2.d0); tabmom(1,2,5)=5                                              
!                                                                               
  tabnew(2,7)=0                                                
!                                                                               
  tabnew(2,8)=1; tabcoeff(1,2,8)=-sqrt(3.d0)/sqrt(2.d0); tabmom(1,2,8)=8                                              
!                                                                               
  tabnew(5,1)=1; tabcoeff(1,5,1)=-1./sqrt(2.d0); tabmom(1,5,1)=5                                              
!                                                                               
  tabnew(5,3)=1; tabcoeff(1,5,3)=-1; tabmom(1,5,3)=4                                              
!                                                                               
  tabnew(5,4)=0                                                
!                                                                               
  tabnew(5,6)=0                                                
!                                                                               
  tabnew(5,2)=1; tabcoeff(1,5,2)=-sqrt(3.d0)/sqrt(2.d0); tabmom(1,5,2)=5                                              
!                                                                               
  tabnew(5,5)=0                                                
!                                                                               
  tabnew(5,7)=1; tabcoeff(1,5,7)=1.; tabmom(1,5,7)=6                                              
!                                                                               
  tabnew(5,8)=2; tabcoeff(1:2,5,8)=(/1./sqrt(2.d0), sqrt(3.d0)/sqrt(2.d0)/); tabmom(1:2,5,8)=(/1,2/)
!                                                                               
  tabnew(7,1)=1;   tabcoeff(1,7,1)=sqrt(2.d0);  tabmom(1,7,1)=7                                              
!                                                                               
  tabnew(7,3)=1;  tabcoeff(1,7,3)=1.;  tabmom(1,7,3)=8                                              
!                                                                               
  tabnew(7,4)=2;  tabcoeff(1:2,7,4)=(/-sqrt(2.d0), 0.d0 /);   tabmom(1:2,7,4)=(/1,2/)
!                                                                               
  tabnew(7,6)=0                                                
!                                                                               
  tabnew(7,2)=0                                                
!                                                                               
  tabnew(7,5)=1;  tabcoeff(1,7,5)=-1.;  tabmom(1,7,5)=6                                              
!                                                                               
  tabnew(7,7)=0                                                
!                                                                               
  tabnew(7,8)=0                                                
!                                                                               
  tabnew(8,1)=1;  tabcoeff(1,8,1)=1./sqrt(2.d0);  tabmom(1,8,1)=8                                              
!                                                                               
  tabnew(8,3)=0                                                
!                                                                               
  tabnew(8,4)=1;  tabcoeff(1,8,4)=-1.;  tabmom(1,8,4)=3                                              
!                                                                               
  tabnew(8,6)=1;  tabcoeff(1,8,6)=1.;  tabmom(1,8,6)=7                                              
!                                                                               
  tabnew(8,2)=1;  tabcoeff(1,8,2)=sqrt(3.d0/2.d0);  tabmom(1,8,2)=8                                              
!                                                                               
  tabnew(8,5)=2;  tabcoeff(1:2,8,5)=(/-1./sqrt(2.d0), -sqrt(3.d0)/sqrt(2.d0)/);  tabmom(1:2,8,5)=(/1,2/)
!                                                                               
  tabnew(8,7)=0                                                
!                                                                               
  tabnew(8,8)=0                                                
!                                                                               
  tabnew(9,-9)=2;  tabcoeff(1:2,9,-9)=(/1./sqrt(2.d0), -1./sqrt(6.d0)/);  tabmom(1:2,9,-9)=(/1,2/)
!                                                                               
  tabnew(9,-10)=1;  tabcoeff(1,9,-10)=1.;  tabmom(1,9,-10)=6                                              
!                                                                               
  tabnew(9,-11)=1;  tabcoeff(1,9,-11)=1.;  tabmom(1,9,-11)=7                                              
!                                                                               
  tabnew(10,-9)=1;  tabcoeff(1,10,-9)=1.;  tabmom(1,10,-9)=3                                              
!                                                                               
  tabnew(10,-10)=2;  tabcoeff(1:2,10,-10)=(/0.d0, sqrt(2.d0/3.d0)/);  tabmom(1:2,10,-10)=(/1,2/)
!                                                                               
  tabnew(10,-11)=1;  tabcoeff(1,10,-11)=1.;  tabmom(1,10,-11)=8                                              
!                                                                               
  tabnew(11,-9)=1;  tabcoeff(1,11,-9)=1.;  tabmom(1,11,-9)=4                                              
!                                                                               
  tabnew(11,-10)=1;  tabcoeff(1,11,-10)=1.;  tabmom(1,11,-10)=5                                              
!                                                                               
  tabnew(11,-11)=2;  tabcoeff(1:2,11,-11)=(/-1./sqrt(2.d0),-1./sqrt(6.d0)/);  tabmom(1:2,11,-11)=(/1,2/)
!                                                                               
  tabnew(1,9)=1;  tabcoeff(1,1,9)=1./sqrt(2.d0);  tabmom(1,1,9)=9                                              
!                                                                               
  tabnew(1,11)=1;  tabcoeff(1,1,11)=-1./sqrt(2.d0);  tabmom(1,1,11)=11                                             
!                                                                               
  tabnew(3,9)=1;  tabcoeff(1,3,9)=1.;  tabmom(1,3,9)=10                                              
!                                                                               
  tabnew(4,9)=1;  tabcoeff(1,4,9)=1.;  tabmom(1,4,9)=11                                              
!                                                                               
  tabnew(6,10)=1;  tabcoeff(1,6,10)=1.;  tabmom(1,6,10)=9                                              
!                                                                               
  tabnew(2,9)=1;  tabcoeff(1,2,9)=-1./sqrt(6.d0);  tabmom(1,2,9)=9                                              
!                                                                               
  tabnew(2,10)=1;  tabcoeff(1,2,10)=sqrt(2.d0/3.d0);  tabmom(1,2,10)=10                                              
!                                                                               
  tabnew(2,11)=1;  tabcoeff(1,2,11)=-1./sqrt(6.d0);  tabmom(1,2,11)=11                                              
!                                                                               
  tabnew(5,10)=1;  tabcoeff(1,5,10)=1.;  tabmom(1,5,10)=11                                              
!                                                                               
  tabnew(7,11)=1;  tabcoeff(1,7,11)=1.;  tabmom(1,7,11)=9                                              
!                                                                               
  tabnew(8,11)=1;  tabcoeff(1,8,11)=1.;  tabmom(1,8,11)=10                                              
!                                                                               
  tabnew(1,-9)=1;  tabcoeff(1,1,-9)=1./sqrt(2.d0);  tabmom(1,1,-9)=-9                                              
!                                                                               
  tabnew(1,-11)=1;  tabcoeff(1,1,-11)=-1./sqrt(2.d0);  tabmom(1,1,-11)=-11
!                                                                               
  tabnew(3,-10)=1;  tabcoeff(1,3,-10)=1.;  tabmom(1,3,-10)=-9                                              
!                                                                               
  tabnew(4,-11)=1;  tabcoeff(1,4,-11)=1.;  tabmom(1,4,-11)=-9                                              
!                                                                               
  tabnew(6,-9)=1;  tabcoeff(1,6,-9)=1.;  tabmom(1,6,-9)=-10                                              
!                                                                               
  tabnew(2,-9)=1;  tabcoeff(1,2,-9)=-1./sqrt(6.d0);  tabmom(1,2,-9)=-9                                              
!                                                                               
  tabnew(2,-10)=1;  tabcoeff(1,2,-10)=sqrt(2.d0/3.d0);  tabmom(1,2,-10)=-10
!                                                                               
  tabnew(2,-11)=1;  tabcoeff(1,2,-11)=-1./sqrt(6.d0);  tabmom(1,2,-11)=-11
!                                                                               
  tabnew(5,-11)=1;  tabcoeff(1,5,-11)=1.;  tabmom(1,5,-11)=-10                                             
!                                                                               
  tabnew(7,-9)=1;  tabcoeff(1,7,-9)=1.;  tabmom(1,7,-9)=-11                                              
!                                                                               
  tabnew(8,-10)=1;  tabcoeff(1,8,-10)=1.;  tabmom(1,8,-10)=-11                                              
!
  do j1=9,11
    tabnew(0,j1)=1
    tabcoeff(1,0,j1)=1.
    tabmom(1,0,j1)=j1
    tabnew(j1,0)=1
    tabcoeff(1,j1,0)=1.
    tabmom(1,j1,0)=j1
  enddo
!                           
  do j1=-11,-9
    tabnew(0,j1)=1
    tabcoeff(1,0,j1)=1.
    tabmom(1,0,j1)=j1
    tabnew(j1,0)=1
    tabcoeff(1,j1,0)=1.
    tabmom(1,j1,0)=j1
  enddo
!
  k1=0
  k2=0
  do j1=1,23
    do j2=1,23
      i1=j1-12
      i2=j2-12
      k1=k1+1
      tbnew(k1)=tabnew(i2,i1)
      do j3=1,2
       k2=k2+1
       tbmom(k2)=tabmom(j3,i2,i1)
       tbcoeff(k2)=tabcoeff(j3,i2,i1)
      enddo
    enddo
  enddo
!
  tbnew(266:273)=1; 
  tbcoeff(531:545:2)=1.d0
  tbmom(531:545:2)=(/1,2,3,4,5,6,7,8/)
!
  tbnew(288:449:23)=1; 
  tbcoeff(575:897:46)=1.d0
  tbmom(575:897:46)=(/1,2,3,4,5,6,7,8/)
! 
end subroutine colprdsu3           
!
subroutine permtn(nfr,p1,p2,perm)
!
  implicit none
!
  integer*4, intent (in) :: nfr,p1,p2
  real*8, intent (out) :: perm
!
  integer :: odd,j3,j1
!
  j3=0
  odd=0
  do j1=0,nfr-1
    if(btest(p2,j1)) then
      odd=odd+1
    else
      if(btest(odd,0)) then
        if(btest(p1,j1)) j3=j3+1
      endif
    endif   
  enddo
!
  if(btest(j3,0)) then
    perm=-1.d0
  else
    perm=1.d0
  endif
!
end subroutine permtn
!
subroutine opt(j)
!
  implicit none
  integer j
!
  j=j+1
!
end subroutine opt
!
       function popcnt0(i)
!
       implicit none
       integer*4 i,popcnt0
       integer*4 j1
!
       popcnt0=0
       do j1=0,31
         if(btest(i,j1))popcnt0=popcnt0+1
       enddo
!
       return
       end
!
subroutine ampd (prcss,dim,mtel)
!
  use alptyp
  use strfld
  use prcint
  use couplings, only : oper, coup 
  use incol
  use dimen
  use utilities
  use running
  implicit none
!
  type (proc), intent(in) :: prcss
  integer, intent(in) :: dim(7)
  complex*16, intent(out) :: mtel             !amplitude
!
  integer*4, parameter :: nca=4
  integer*4 :: nprt,nfld,flvold,padd,nconf,j5,padd0
  integer*4 :: j1,j2,iter,steptm,stepex,ex2,ex3,nf2,nf3,psfl2,psfl3
  integer*4 :: pscfl2,pscfl3,nl1,nl2,nl3,j3,j4,beg3,nc,psfl3i,pscfl3i
  integer*4 :: p1,p2,p3,posfld,poscol,nlold
  integer*4 :: nf1,psfl1,pscfl1,pref,beg
  integer*4 :: cl(4),clnw(2)
  integer*4 :: psnew(dim(d_psnew))
  integer*4, save :: fststp=0
  integer*4, dimension (9*nintmax+1) :: tabmrgpr
  integer*4, dimension(nexcmax) :: exc,excsm,excns
  integer*4 :: offshcol(dim(d_offshcol)),ampcol(dim(d_ampcol)),ordtmp(dim(d_ordtmp))
  complex*16  ::  offsh(dim(d_offsh)),ampfld(dim(d_ampfld))
  complex*16, dimension (nlor) :: fusion,polar
  complex*16 :: amplitude,tmp
  real*8 :: colcff
  real*8 :: impul(4,dim(d_impul))
  real*8, dimension (4) :: imp
  real*8 :: perm
!
  integer*4 :: nphf,even,nmom,ot,npmit,tmpi,ci,cf,ptr,nfr2
  integer*4 :: operot,jop,op(3*nintmax),smold
  logical :: choose,choose2,choose3,choose4,cnsold(0:npmax)
  real*8, dimension(4) :: imp3,imp2
  real*8 :: coupot,cpot1,cp(3*nintmax)
  integer,save :: psfit=0,pscit=0,psfa=0,psca=0
  logical, parameter :: dbg=.false.
!
! once for all inizialization
  if (fststp.eq.0) then
    fststp=1
! inizialization of arrays
    call initcoup
    call initmom
    call inittab
    call colprdsu3
!
! inizialization of merging tables
    call initprc
!
  endif
!
  fldptr=-9999
  nprt=0                                   ! number of external particle
  do j1=1,nmax
    if (prcss%flv(j1).eq.1001) exit
    nprt=nprt+1
  enddo
!
  nphf=nprt/2              ! for optmization
  even=mod(nprt,2)         ! for optmization
  nmom=2**nprt-1           ! maximum number of different combinations of nprt 
!                            momenta, each momenta enters at most once
! inizialization step
!
  call rncpl         !loading running couplings (alpha-strong)
!
  posfld=1           ! each off shell amplitude will be collected into a single 
  poscol=1           ! array OFFSH, POSFLD and POSCOL allows to find the position
  flvold=-99          
  do j1=1,nprt
    if(prcss%flv(j1).ne.flvold) then
      if(flvold.gt.0) fldptr(1,flvold,1)=nfld 
      nfld=1
      flvold=prcss%flv(j1)
      fldptr(2,flvold,1)=poscol
    else
      nfld=nfld+1
    endif
    p1=2**(j1-1)
    call initpol(prcss,j1,offsh(posfld)) !OFFSH(posfld,....) contains the 
!                                        polarization of the jth particle
    offshcol(poscol:poscol+nca-1)=(/p1,posfld,prcss%col(2*j1-1),prcss%col(2*j1)/)
    impul(1:4,offshcol(poscol))=prcss%mom(4*(j1-1)+1:4*j1)
    posfld=posfld+nlordof(prcss%flv(j1))
    poscol=poscol+nca
  enddo
  fldptr(1,flvold,1)=nfld
!
! iteration 
! we merge on/off shell field2 and field3 into field1 via equation of motion up to the step required by the algorithm
  tabmrgpr=tabintmrg(1:,nprc)
  cp(1:tabintmrg(1,nprc))=cpfs(1:tabintmrg(1,nprc),nprc)
  op(1:tabintmrg(1,nprc))=opfs(1:tabintmrg(1,nprc),nprc)
  do j1=1,nvasfs(1,nprc)
    cp(nvasfs(j1+1,nprc))=cp(nvasfs(j1+1,nprc))*rnas
  enddo
  psnew(1:nmom)=0        
  do iter=2,nphf
    cnsold=.false.
    steptm=2                  !for optimization
    flvold=-99   
    excsm=nexcsm(1:,iter)
    excns=nexc(1:,iter)
    do j1=1,tabmrgpr(1)
      fl1=tabmrgpr(steptm)    !labels rthe flavour of the daughter offshell particle
      fl2=tabmrgpr(steptm+1)  !fl2 and fl3 label the flavour of the parents offshell
      fl3=tabmrgpr(steptm+2)  !particles to be merged
!
      coupot=cp(j1)                                    !coup(fl1,fl2,fl3)
      operot=op(j1)                                    !oper(fl1,fl2,fl3)
!
      nl1=nlordof(fl1)        !nl1,nl2,nl3 are the corresponding numbers of lorentz
      nl2=nlordof(fl2)        !degrees of freedom
      nl3=nlordof(fl3)        
      if(fl1.ne.flvold) then
        if(flvold.gt.0) then
          fldptr(1,chcg(flvold),iter)=nconf
        else
          smold=0; if(fl2.eq.fl3)smold=1
        endif
        flvold=fl1
        if (cnsold(cnschg(fl1))) then
          psnew(1:nmom)=0                      !for optimization
          cnsold=.false.
        endif
        cnsold(cnschg(fl1))=.true.
        nconf=0
        fldptr(2,chcg(fl1),iter)=poscol
      endif
      steptm=steptm+3
      if(fl2.ne.fl3) then
        if(smold.eq.0) exc=excns
        smold=1
      else
        if(smold.eq.1) exc=excsm
        smold=0
      endif
      stepex=2
      do j2=1,exc(1)
        ex2=exc(stepex)
        ex3=exc(stepex+1)
        stepex=stepex+2
        nf2=fldptr(1,fl2,ex2)                               !number of excitation of field2
        nf3=fldptr(1,fl3,ex3)                               !number of excitation of field3
!
        if(nf2.le.0.or.nf3.le.0) cycle
!
        pscfl2=fldptr(2,fl2,ex2)                            !starting position of field2 in array OFFSHCOL
        pscfl3=fldptr(2,fl3,ex3)                            !starting position of field2 in array OFFSHCOL
        psfl2=offshcol(pscfl2+1)                             !starting position of field2 in array OFFSH
        psfl3=offshcol(pscfl3+1)                             !starting position of field3 in array OFFSH
        do j3=1,nf2
          if (ex2.eq.ex3.and.fl2.eq.fl3) then               !if the flavour of field 2 and 3 are the same don't 
            beg3=j3+1                                       !repeat twice the same computation 
            psfl3i=psfl3+nl3*(beg3-2)
            pscfl3i=pscfl3+nca*(beg3-2)
          else
           beg3=1
           psfl3i=psfl3-nl3
           pscfl3i=pscfl3-nca
          endif
!
          p2=offshcol(pscfl2)                               !momentum of field2 off-shell amplitude 
!          cl(1:2)=(/offshcol(pscfl2+2:pscfl2+3)/)
          cl(1)=offshcol(pscfl2+2); cl(2)=offshcol(pscfl2+3)
          imp2=impul(1:4,p2)
!          imp2(1)=impul(1,p2); imp2(2)=impul(2,p2); imp2(3)=impul(3,p2); imp2(4)=impul(4,p2); 
!
          choose=iter.eq.nphf.and.even.eq.0          !for optimization: true for even number of external particle and
!                                                    !the resulting off shell field is the merging of half this number
!
          if(prcss%nfrm.gt.2) then
            nfr2=ibits(p2,0,prcss%nfrm)
          else
            nfr2=0
          endif
!
          do j4=beg3,nf3
!
            pscfl3i=pscfl3i+nca
            psfl3i=offshcol(pscfl3i+1)                             !starting position of field3 in array OFFSH

!
            p3=offshcol(pscfl3i)                         !momentum of field3
            p1=p2+p3                                     !momentum of field1: resulting out of the merger of field2,3
!            if(p1.gt.nmom) cycle
            if(momd(p1).ne.iter) cycle               !if this is not true some momenta enter more than once in p1
            if(choose.and.btest(p1+1,0)) cycle       !for even number of particle we accept merger of npart/2 
!                                                        !momenta only if momenta 1 contribute 
!            cl(3:4)=(/offshcol(pscfl3i+2:pscfl3i+3)/)  !contains the coulors of field2 and 3
            cl(3)=offshcol(pscfl3i+2); cl(4)=offshcol(pscfl3i+3)  !contains the coulors of field2 and 3
            colprdit : select case (colrep(fl1))         
              case(3)                                  !quarks and gluons
                if(cl(2).eq.0) cycle
                if (cl(2).eq.cl(3)) then
                  if(cl(1).ne.0) then
                    if(cl(1).eq.cl(4)) cycle
                    colcff=1.d0
                  else
                    colcff=-0.333333333333333d0 
                  endif
                  clnw(1)=cl(1); clnw(2)=cl(4)
!                  clnw=(/cl(1),cl(4)/)
                elseif (cl(1).eq.cl(4)) then
                  if(cl(1).ne.0) then
                    colcff=-1.d0
                  else
                    colcff=1.d0 
                  endif
                  clnw(1)=cl(3); clnw(2)=cl(2)
!                  clnw=(/cl(3),cl(2)/)
                elseif (colrep(fl2).eq.1.and.colrep(fl3).eq.3) then
                  colcff=1.d0; clnw=cl(3:4);
                elseif (colrep(fl2).eq.3.and.colrep(fl3).eq.1) then
                  colcff=1.d0; clnw=cl(1:2);
                else
                  cycle
                endif
              case(2)                                    !coulorless object
                if (cl(2).eq.cl(3)) then
                  clnw(1)=cl(1); clnw(2)=cl(4)
!                  clnw=(/cl(1),cl(4)/)
                elseif (cl(1).eq.cl(4)) then
                  clnw(1)=cl(3); clnw(2)=cl(2)
!                  clnw=(/cl(3),cl(2)/)
                else
                  cycle
                endif
                colcff=1.d0
              case(1)                                    !coulorless object
                if (cl(2).eq.cl(3).and.cl(1).eq.cl(4)) then
!                  clnw=(/0,0/)
                  clnw(1)=0; clnw(2)=0
                else
                  cycle
                endif
                colcff=1.d0
              case default
                write(*,*)'something wrong in color assignment'
            end select colprdit
            imp3=impul(1:4,p3)
            imp=imp2+imp3             !momentum of the new fields
            impul (1:4,p1)=imp                          !array to store momenta combinations
            call fuse (-imp,imp2,offsh(psfl2),imp3,  &     !merging the fields 2 and 3 into fusion
                         offsh(psfl3i),operot,fusion)         !via the proper trilinear interaction (oper)
            call prpgt (chcg(fl1),imp,fusion,polar)                      !multiplying by the propagator -> polar
            
            if(nfr2.gt.0) then
              call permtn(prcss%nfrm,p2,p3,perm)                            !returning proper fermi statistics in perm
            else
              perm=1.d0
            endif
            padd=psnew(p1)                    !returning the position, in OFFSHCOL, of field1
!
!!$            write(*,*)'fl',fl1,fl2,fl3
!!$            write(*,*)'cl',cl,clnw,colcff
!!$            write(*,*)'pr',perm,p2,p3,p1,'padd',padd,'cp',coupot
!!$!
            if (padd.eq.0) then             !new configuration
              psnew(p1)=poscol
              cpot1=coupot*colcff*perm
              offsh(posfld:posfld+nl1-1)=cpot1*polar(1:nl1)            !storing field1
              offshcol(poscol:poscol+nca-1)=(/p1,posfld,clnw(1),clnw(2)/)     !storing pointers to field1
              nconf=nconf+1                      !number of new excitation of the same flavour and number of momenta
              poscol=poscol+nca
              posfld=posfld+nl1
            else                                                         !new contribution to an old configuration
              cpot1=coupot*colcff*perm
              padd0=offshcol(padd+1)                                                !position of field1 in OFFSH
              offsh(padd0:padd0+nl1-1)=offsh(padd0:padd0+nl1-1)+cpot1*polar(1:nl1)      !adding the new contribution
            endif
!
          enddo
          pscfl2=pscfl2+nca
          psfl2=offshcol(pscfl2+1)
        enddo
      enddo
!
    enddo
!
    fldptr(1,chcg(fl1),iter)=nconf
  enddo
!
  if (fststp.eq.1.and.dbg) then
    fststp=2
    if (poscol.gt.pscit) then
      pscit=poscol
      write(*,*)'iteration,poscol,posfld',poscol,posfld
    endif
    if (posfld.gt.psfit) then
      psfit=posfld
      write(*,*)'iteration,poscol,posfld',poscol,posfld
    endif
  endif
!
  if(poscol.gt.dim(d_offshcol).or.psfit.gt.dim(d_offsh)) then
    write(*,*)'ill dimensioned array OFFSH and/or OFFSHCOL in AMP'
    write(*,*)'iteration,poscol,posfld',poscol,posfld,dim(d_offshcol),dim(d_offsh)
    write(*,*)'prcss',prcss
    stop
  endif
!
! amplitude 
!
  amplitude=(0.,0.)
!
  tabmrgpr(1:3*nintmax+1)=tabint(1:3*nintmax+1,nprc)      !most of the steps are the same as in iteration  only the list of acceptable fusions is different
  cp(1:tabint(1,nprc))=cpam(1:tabint(1,nprc),nprc)
  op(1:tabint(1,nprc))=opam(1:tabint(1,nprc),nprc)
  do j1=1,nvasam(1,nprc)
    cp(nvasam(j1+1,nprc))=cp(nvasam(j1+1,nprc))*rnas
  enddo  
!                                            we first merge field2 and field3 and finally contract with field1 
!
  psnew(1:nmom)=0             !none of these mergin has been computed before
  ordtmp(1:nmom)=0            !reporting the order of field1 momenta to avoid repeating the sam    
  cnsold=.false.
  steptm=2 
  flvold=-99   
  do j1=1,tabmrgpr(1)+1
!
    if(j1.eq.tabmrgpr(1)+1) then
      flvold=9999                                !last step compute the final contribution to the amplitude and exit
    else
      fl1=tabmrgpr(steptm)
      fl2=tabmrgpr(steptm+1)
      fl3=tabmrgpr(steptm+2)
      steptm=steptm+3    
      nl1=nlordof(fl1)
      nl2=nlordof(fl2)
      nl3=nlordof(fl3)
!
      coupot=cp(j1)          !coup(fl1,fl2,fl3)
      operot=op(j1)          !oper(fl1,fl2,fl3)
!
    endif
!
    if(fl1.ne.flvold) then             !we have computed all 2,3 merging compatible with contraction with old 1
!
      if (fststp.eq.2.and.j1.ne.1) then
        if(poscol.gt.dim(d_ampcol).or.posfld.gt.dim(d_ampfld)) then
          write(*,*)'ill dimensioned array OFFSH and/or OFFSHCOL in AMP'
          write(*,*)'amplitude,poscol,posfld',poscol,posfld,dim(d_ampcol),dim(d_ampfld),j1
          stop
        endif
!
        if(dbg) then
          fststp=2
          if (poscol.gt.psca) then
            psca=poscol
            write(*,*)'amplitude,poscol,posfld',poscol,posfld
          endif
          if (posfld.gt.psfa) then
            psfa=posfld
            write(*,*)'amplitude,poscol,posfld',poscol,posfld
          endif
        endif
      endif
!
      if (flvold.gt.0) then            !compute the contribution to the amplitude (flvold < 0 iteration just finished)
        if(flvold.eq.9999)flvold=fl1   !flvold=9999 last contribution
        nlold=nlordof(flvold)          !number of lorentz degrees of freedom          
        do iter=1,nphf
          nf1=fldptr(1,flvold,iter)       !same as in iteration
          if(nf1.le.0) cycle
          pscfl1=fldptr(2,flvold,iter)
          psfl1=offshcol(pscfl1+1)
          pref=0                          !for optimization
          do j2=1,nf1
            psfl1=offshcol(pscfl1+1)
            p1=offshcol(pscfl1)
            p1=nmom-p1
            padd=psnew(p1)                !locating position of field1
            if (padd.eq.0) then           !no 2,3 merging compatible with field1
              pscfl1=pscfl1+nca
              cycle
            endif
            if (p1.ne.pref) then
              pref=p1
            else
              padd=padd+nca               !there are two neutral gluons
            endif
            cl=(/offshcol(pscfl1+2),offshcol(pscfl1+3),ampcol(padd+2),ampcol(padd+3)/)
            if(abs(cl(1)-cl(4))+abs(cl(2)-cl(3)).ne.0) cycle
            padd=ampcol(padd+1)
            tmp=sum(offsh(psfl1:psfl1+nlold-1)*ampfld(padd:padd+nlold-1))
            if(prcss%nfrm.ge.2) then
              call permtn(prcss%nfrm,nmom-p1,p1,perm)                            !returning proper fermi statistics in perm
            else
              perm=1.d0
            endif
            amplitude=amplitude+tmp*perm
            pscfl1=pscfl1+nca
          enddo
        enddo
        if(j1.eq.tabmrgpr(1)+1) cycle                             !end of computation
!
      endif
! 
      flvold=fl1                  !a new class of fusion 2,3 beginning, compatible to match with new flavour fl1
      posfld=1                    !pointer to offsh
      poscol=1                    !pointer to offshcol    
      if (cnsold(cnschg(fl1))) then
        psnew(1:nmom)=0             !none of these mergin has been computed before
        ordtmp(1:nmom)=0            !reporting the order of field1 momenta to avoid repeating the same computation for
      endif
      cnsold(cnschg(flvold))=.true.
!                                  interaction quadratic or cubic in the same field
!
    endif
!
    do iter=1,nphf
!
      npmit=nprt-iter
      choose2=iter.eq.nphf.and.even.eq.0
      nf1=fldptr(1,fl1,iter)
      if (nf1.le.0) cycle
      pscfl1=fldptr(2,fl1,iter)
!
      if(fl1.eq.fl2) then                     !again to deal with equal flavour
        beg=iter
      else
        beg=1
      endif
!
      do ex2=beg,nphf
!
        choose=fl1.eq.fl2.and.iter.eq.ex2     !again to deal with equal flavour  
!
        ex3=nprt-iter-ex2                     !the total number of momenta must match the number of external particle
        if(ex3.lt.1.or.ex3.gt.nphf) cycle
        if (fl2.eq.fl3.and.ex3.lt.ex2) cycle  !again to deal with equal flavour  
        nf2=fldptr(1,fl2,ex2)
        nf3=fldptr(1,fl3,ex3)
        if(nf2.le.0.or.nf3.le.0) cycle
        pscfl2=fldptr(2,fl2,ex2)
        pscfl3=fldptr(2,fl3,ex3)
!
        choose3=ex2.eq.ex3.and.fl2.eq.fl3     !again to deal with equal flavour for optimization
!
          j5=pscfl1
          do j4=1,nf1
            ordtmp(offshcol(j5))=j4           !now ordtmp register the ordering of momenta of field1
            j5=j5+nca
          enddo
!
        do j3=1,nf2
!
!          cl(1:2)=offshcol(pscfl2+2:pscfl2+3)
          cl(1)=offshcol(pscfl2+2);           cl(2)=offshcol(pscfl2+3)
          psfl2=offshcol(pscfl2+1)
!
          if (choose3) then
            beg3=j3+1
            pscfl3i=pscfl3+nca*(beg3-2)
          else
            beg3=1
            pscfl3i=pscfl3-nca
          endif
!
          p2=offshcol(pscfl2)
          imp2=impul(1:4,p2)
          if(prcss%nfrm.gt.2) then
            nfr2=ibits(p2,0,prcss%nfrm)
          else
            nfr2=0
          endif
!
          do j4=beg3,nf3
            pscfl3i=pscfl3i+nca
!
            p3=offshcol(pscfl3i) 
            p1=p2+p3
!            if(p1.gt.nmom) cycle
            if(momd(p1).ne.npmit) cycle          !momenta appearing more than once
            ot=ordtmp(nmom-p1)
            if(ot.eq.0) cycle              !field1 does not contain the frequency to match the proposed 2,3 fusion
            if(choose.and.ordtmp(p2).le.ot)cycle     !to avoid repeating the same computation twice
!            if(choose2.and.momb%evn(p1).eq.1) cycle
            if(choose2.and.btest(p1,0)) cycle
!            cl(3:4)=offshcol(pscfl3i+2:pscfl3i+3)
            cl(3)=offshcol(pscfl3i+2); cl(4)=offshcol(pscfl3i+3)
            colprdam : select case (colrep(fl1))
              case(3)                                  !quarks and gluons
                if(cl(2).eq.0) cycle
                if (cl(2).eq.cl(3)) then
                  if(cl(1).ne.0) then
                    if(cl(1).eq.cl(4))cycle
                    colcff=1.d0
                  else
                    colcff=1.d0 
                  endif
!                  clnw=(/cl(1),cl(4)/)
                  clnw(1)=cl(1); clnw(2)=cl(4)
                elseif (cl(1).eq.cl(4)) then
                  if(cl(1).ne.0) then
                    colcff=-1.d0
                  else
                    colcff=1.d0 
                  endif
!                  clnw=(/cl(3),cl(2)/)
                  clnw(1)=cl(3); clnw(2)=cl(2)
                elseif(colrep(fl2).eq.1.and.colrep(fl3).eq.3) then
                  colcff=1.d0; clnw=cl(3:4);
                elseif(colrep(fl2).eq.3.and.colrep(fl3).eq.1) then
                  colcff=1.d0; clnw=cl(1:2);
                else
                  cycle
                endif
              case(2)                                    !coulorless object
                if (cl(2).eq.cl(3)) then
                  clnw(1)=cl(1); clnw(2)=cl(4)
!                  clnw=(/cl(1),cl(4)/)
                elseif (cl(1).eq.cl(4)) then
                  clnw(1)=cl(3); clnw(2)=cl(2)
!                  clnw=(/cl(3),cl(2)/)
                else
                  cycle
                endif
                colcff=1.d0
              case(1)                                    !coulorless object
                if (cl(2).eq.cl(3).and.cl(1).eq.cl(4)) then
                  clnw(1)=0; clnw(2)=0
!                  clnw=(/0,0/)
                else
                  cycle
                endif
                colcff=1.d0
              case default
                write(*,*)'something wrong in color assignment'
            end select colprdam
!
            imp3=impul(1:4,p3)
            imp=imp2+imp3
            psfl3i=offshcol(pscfl3i+1)
            call fuse (-imp,imp2,offsh(psfl2),imp3,offsh(psfl3i),operot,fusion)
!              call permt(prcss%nfrm,p2,p3,perm)
            if(nfr2.gt.0) then
!                if(ibits(p3,0,prcss%nfrm).ne.0) then
              call permtn(prcss%nfrm,p2,p3,perm)
!                else
!                  perm=1.d0
!                endif
            else
              perm=1.d0
            endif
            padd=psnew(p1)
!
!!$            write(*,*)'fl',fl1,fl2,fl3
!!$            write(*,*)'cl',cl,clnw,colcff
!!$            write(*,*)'pr',perm,p2,p3,p1,'padd',padd,'cp',coupot
!
            if (padd.eq.0) then             !new configuration
              psnew(p1)=poscol
              cpot1=coupot*colcff*perm
              ampfld(posfld:posfld+nl1-1)=cpot1*fusion(1:nl1)
              ampcol(poscol:poscol+nca-1)=(/p1,posfld,clnw(1),clnw(2)/)
              nconf=nconf+1
              posfld=posfld+nl1
              poscol=poscol+nca
            else                            !new contribution to an old configuration
              padd0=ampcol(padd+1)
              cpot1=coupot*colcff*perm
              ampfld(padd0:padd0+nl1-1)=ampfld(padd0:padd0+nl1-1)+cpot1*fusion(1:nl1)
            endif
!
          enddo
          pscfl2=pscfl2+nca
        enddo
      enddo
!
    enddo
!
  enddo
!
  mtel=amplitude
!
end subroutine ampd
!***************************************************************************              
          subroutine sourcemassboson(spinsour,fieldaux,mom,mass)                          
!***************************************************************************              
!                                                                                         
! This subroutine given the four momentum of a massive boson MOM the required             
! source polarization SPINSOUR and the boson mass MASS returns in                         
! FIELDAUX the source term                                                                
!                                                                                         
          implicit none                                                                   
!                                                                                         
          integer n               !loop index                                             
!                                                                                         
          real*8 knorm2           !squared modulus of the three momentum                  
          real*8 mass             !massive boson mass                                     
          real*8 massa2           !squared mass                                           
          real*8 mom(4)           !four momentum                                          
          real*8 nz(4),nx(4)                                                              
          real*8 sour1(4),sour2(4),sour3(4)   !longitudinal (1)                           
!                                and transverse (2,3) polarizations                       
          integer spinsour   !required polarization                                      
          real*8 scalar3                                                                  
          real*8 xnorm      !to normalize the sources                                     
          complex*16 fieldaux(4)   !to return the relevant polarization                  
!                                                                                         
          data nz/0.,0.,0.,1./                                                            
          data nx/0.,1.,0.,0./                                                            
          data sour1/4*0./                                                                
          data sour2/4*0./                                                                
          data sour3/4*0./                                                                
!                                                                                         
          save nx,nz                                                                      
!                                                                                         
!                                                                                         
!                                                                                         
          knorm2=mom(2)**2+mom(3)**2+mom(4)**2                                            
!                                                                                         
          if (knorm2.eq.0) then                                                           
           sour1(2)=1.                                                                    
           sour2(3)=1.                                                                    
           sour3(4)=1.                                                                    
           return                                                                         
          endif                                                                           
          massa2=mass**2                                                                  
!                                                                                         
!    pol. longitudinal                                                                    
!                                                                                         
          do n=2,4                                                                        
           sour1(n)=mom(n)                                                                
          enddo                                                                           
	  sour1(1)=knorm2/mom(1)                                                                 
          xnorm=sqrt(scalar3(sour1,sour1))                                                
          if(xnorm.ne.0.) then                                                            
           do n=1,4                                                                       
            sour1(n)=sour1(n)/xnorm                                                       
           enddo                                                                          
          endif                                                                           
!                                                                                         
!    pol. transverse +                                                                    
!                                                                                         
          if (mom(4)**2.gt.mom(2)**2) then                                                
           call vector3prod(mom,nx,sour2)                                                 
          else                                                                            
           call vector3prod(mom,nz,sour2)                                                 
          endif                                                                           
          xnorm=Sqrt(scalar3(sour2,sour2))                                                
          do n=1,4                                                                        
           sour2(n)=sour2(n)/xnorm                                                        
          enddo                                                                           
!                                                                                         
!    pol. transverse -                                                                    
!                                                                                         
          call vector3prod(mom,sour2,sour3)                                               
          xnorm=sqrt(scalar3(sour3,sour3))                                                
          do n=1,4                                                                        
           sour3(n)=sour3(n)/xnorm                                                        
          enddo                                                                           
!                                                                                         
          do n=1,4                                                                        
           fieldaux(n)=0.                                                                 
           if(spinsour.eq.0)fieldaux(n)=sour1(n)                                    
           if(spinsour.eq.1)fieldaux(n)=sour2(n)                                    
           if(spinsour.eq.-1)fieldaux(n)=sour3(n)                                   
          enddo                                                                           
!                                                                                         
          return                                                                          
          end                                                                             
!*********************************************************************                    
	  subroutine vector3prod(v1,v2,vfin)                                                     
!*********************************************************************                    
!                                                                                         
          implicit none                                                                   
!                                                                                         
           real*8 v1(4)                                                                   
           real*8 v2(4)                                                                   
           real*8 vfin(4)                                                                 
!                                                                                         
           vfin(2)=v1(3)*v2(4)-v1(4)*v2(3)                                                
           vfin(3)=v1(4)*v2(2)-v1(2)*v2(4)                                                
           vfin(4)=v1(2)*v2(3)-v1(3)*v2(2)                                                
           vfin(1)=0                                                                      
!                                                                                         
           return                                                                         
           end                                                                            
!***********************************************************************                  
           function scalar3(v1,v2)                                                        
!***********************************************************************                  
!                                                                                         
           implicit none                                                                  
!                                                                                         
           real*8 scalar3                                                                 
	   real*8 v1(4)                                                                          
           real*8 v2(4)                                                                   
!                                                                                         
           scalar3=abs(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))                   
!                                                                                         
            return                                                                        
            end                                                                           
!***********************************************************************                  
          subroutine fermionsourceshl(spinsour,fieldaux ,mom)                         
!***********************************************************************                  
!                                                                                         
! This subroutine return the initialized fermion field.                                   
! Elicity heigenstates. Massless fermions
!                                                                                         
          implicit none                                                                   
!                                                                                         
          real*8  mom(4)      !momenta of the external particle              
          real*8 p,p3p,p3m,coeffp,coeffm                                               
!                                                                                         
          complex*16 fieldaux(2)  !array returning the fermion field configurati          
          complex*16 im  !immaginary unity in complex representation                      
          complex*16 p1p                                                                  
          data im/(0.,1.)/                                                                
          integer spinsour
          real*8 :: usd 
!                                                                                         
! Static variables                                                                        
!                                                                                         
          save im                                                                         
!                                                                                         
!                                                                                         
!                                                                                         
          p=sqrt(mom(2)**2+mom(3)**2+mom(4)**2)                                           
          p3p=p+mom(4)                                                                    
          p3m=p-mom(4)                                                                    
          p1p=mom(2)+mom(3)*im                                                            
!                                                                                         
          if(abs(spinsour).eq.1) spinsour=-spinsour                                  
          if(abs(p3m).lt.1.d-10.or.abs(p3p).lt.1.d-10)then                                
           call fshl(spinsour,fieldaux ,mom)                                           
           return                                                                         
          endif                                                                           
!                                                                                         
          coeffp=1./Sqrt(p3p)                                                     
          coeffm=1./Sqrt(p3m)                                                     
!                                                                                         
! "spin up" ingoing fermion                                                               
!                                                                                         
           if (spinsour.eq.1) then                                                       
            fieldaux(1)=coeffp*p3p
            fieldaux(2)=coeffp*p1p
           endif                                                                          
!                                                                                         
! "spin down" ingoing fermion                                                             
!                                                                                         
           if (spinsour.eq.-1) then                                                      
            fieldaux(1)=coeffm*p3m                                                     
            fieldaux(2)=-coeffm*p1p                                                    
           endif                                                                          
!                                                                                         
! "spin up" outgoing antifermion                                                          
!                                                                                         
           if (spinsour.eq.2) then                                                       
            fieldaux(1)=-coeffm*p3m                                                     
            fieldaux(2)=coeffm*p1p                                                      
           endif                                                                          
!                                                                                         
! "spin down" outgoing antifermion                                                        
!                                                                                         
           if (spinsour.eq.-2) then                                                      
            fieldaux(1)=coeffp*p3p                                                     
            fieldaux(2)=coeffp*p1p                                                     
           endif                                                                          
!                                                                                         
           return                                                                         
           end                                                                            
                                                                                          
!***********************************************************************                  
          subroutine fermionbarsourceshl(spinsour,fieldaux ,mom)                      
!***********************************************************************                  
!                                                                                         
! This subroutine return the initialized fermionbar field.                                
! Elicity heigenstates. Massless fermion                                                  
!                                                                                         
          implicit none                                                                   
!                                                                                         
          real*8  mom(4)          !array containing the momenta of the external           
          real*8 p,p3p,p3m,coeffp,coeffm                                               
!                                                                                         
          complex*16 fieldaux(2)     !array returning the fermion field configur          
          complex*16 im  !immaginary unity in complex representation                      
          complex*16 p1p                                                                  
          data im/(0.,1.)/                                                                
          integer spinsour                                                                
!                                                                                         
! Static variables                                                                        
!                                                                                         
          save im                                                                         
!                                                                                         
!                                                                                         
          p=sqrt(mom(2)**2+mom(3)**2+mom(4)**2)                                           
          p3p=p+mom(4)                                                                    
          p3m=p-mom(4)                                                                    
          p1p=mom(2)-mom(3)*im                                                            
!                                                                                         
          if(abs(spinsour).eq.1) spinsour=-spinsour                                  
          if(abs(p3m).lt.1.d-10.or.abs(p3p).lt.1.d-10)then                                
           call fbshl(spinsour,fieldaux ,mom)                                          
           return                                                                         
          endif                                                                           
!                                                                                         
          coeffp=1./Sqrt(p3p)                                                     
          coeffm=1./Sqrt(p3m)
!                                                                                         
! "spin up" ingoing fermion                                                               
!                                                                                         
           if (spinsour.eq.1) then                                                       
            fieldaux(1)=coeffp*p3p                                                     
            fieldaux(2)=coeffp*p1p
           endif                                                                          
!                                                                                         
! "spin down" ingoing fermion                                                             
!                                                                                         
           if (spinsour.eq.-1) then                                                      
            fieldaux(1)=coeffm*p3m
            fieldaux(2)=-coeffm*p1p
           endif                                                                          
!                                                                                         
! "spin up" outgoing antifermion                                                          
!                                                                                         
           if (spinsour.eq.2) then                                                       
            fieldaux(1)=-coeffm*p3m
            fieldaux(2)=coeffm*p1p
           endif                                                                          
!                                                                                         
! "spin down" outgoing antifermion                                                        
!                                                                                         
           if (spinsour.eq.-2) then                                                      
            fieldaux(1)=coeffp*p3p
            fieldaux(2)=coeffp*p1p
           endif                                                                          
!                                                                                         
           return                                                                         
           end                                                                            
!***********************************************************************                  
          subroutine fshl(spinsour,fieldaux ,mom)                                      
!***********************************************************************                  
!                                                                                         
! This subroutine return the initialized fermion field. Massless fermion
!                                                                                         
          implicit none                                                                   
!                                                                                         
          real*8  mom(4)       !array containing the momenta of the external par          
          integer spinsour    !array containing the spin of the source                    
!                                                                                         
          complex*16 fieldaux(2)     !array returning the fermion field configur          
          complex*16 im  !immaginary unity in complex representation                      
          data im/(0.,1.)/                                                                
!                                                                                         
! Static variables                                                                        
!                                                                                         
          save im                                                                         
!                                                                                         
!                                                                                         
! "spin up" ingoing fermion                                                               
!                                                                                         
           if(mom(4).lt.0) spinsour=-spinsour                                             
           if (spinsour.eq.1) then                                                       
            fieldaux(1)=sqrt(2.d0*mom(1))
            fieldaux(2)=0.d0
!                                                                                         
! "spin down" ingoing fermion                                                             
!                                                                                         
           elseif (spinsour.eq.-1) then                                                  
            fieldaux(1)=0.d0
            fieldaux(2)=sqrt(2.d0*mom(1))
!                                                                                         
! "spin up" outgoing antifermion                                                          
!                                                                                         
           elseif (spinsour.eq.2) then                                                   
            fieldaux(1)=sqrt(2.d0*mom(1))*mom(4)/abs(mom(4))
            fieldaux(2)=0.d0
!                                                                                         
! "spin down" outgoing antifermion                                                        
!                                                                                         
           elseif (spinsour.eq.-2) then                                                  
            fieldaux(1)=0.d0
            fieldaux(2)=-sqrt(2.d0*mom(1))*mom(4)/abs(mom(4))
           endif                                                                          
!
           return                                                                         
           end                                                                            
!***********************************************************************                  
          subroutine fbshl(spinsour,fieldaux,mom)                                      
!***********************************************************************                  
!                                                                                         
! This subroutine return the initialized fermion field.                                   
!                                                                                         
          implicit none                                                                   
!                                                                                         
          real*8  mass        !containing the mass of the external particle               
          real*8  mom(4)      !array containing the momenta of the external part          
          integer spinsour    !array containing the spin of the source                    
!                                                                                         
          complex*16 fieldaux(2)   !array returning the fermionbar field configu          
          complex*16 im          !immaginary unity in complex representation             
          data im/(0.,1.)/                                                                
!                                                                                         
! Static variables                                                                        
!                                                                                         
          save im                                                                         
!                                                                                         
!                                                                                         
! "spin up" outgoing fermion                                                              
!                                                                                         
           if(mom(4).gt.0) spinsour=-spinsour                                             
           if (spinsour.eq.1) then                                                       
            fieldaux(1)=sqrt(2.d0*mom(1)) 
            fieldaux(2)=0.d0                                  
!                                                                                         
! "spin down" outgoing fermion                                                            
!                                                                                         
           elseif (spinsour.eq.-1) then                                                  
            fieldaux(2)=sqrt(2.d0*mom(1)) 
            fieldaux(1)=0.d0                                  
!                                                                                         
! "spin up" ingoing antifermion                                                           
!                                                                                         
           elseif (spinsour.eq.2) then                                                   
            fieldaux(1)=-sqrt(2.d0*mom(1))*mom(4)/abs(mom(4))
            fieldaux(2)=0.d0                                  
!                                                                                         
! "spin down" ingoing antifermion                                                         
!                                                                                         
           elseif (spinsour.eq.-2) then                                                  
            fieldaux(2)=sqrt(2.d0*mom(1))*mom(4)/abs(mom(4)) 
            fieldaux(1)=0.d0                                  
           endif                                                                          
!                                                                                         
           return                                                                         
           end                                                                            
!***********************************************************************                  
          subroutine fermionsources(spinsour,fieldaux ,mom,mass)                         
!***********************************************************************                  
!                                                                                         
! This subroutine return the initialized fermion field.                                   
! Elicity heigenstates                                                                    
!                                                                                         
          implicit none                                                                   
!                                                                                         
          real*8  mass         !containing the mass of the external particle              
          real*8  mom(4)          !array containing the momenta of the external           
          real*8 p,p3p,p3m,mp,coeffp,coeffm                                               
!                                                                                         
          complex*16 fieldaux(4)  !array returning the fermion field configurati          
          complex*16 im  !immaginary unity in complex representation                      
          complex*16 p1p                                                                  
          data im/(0.,1.)/                                                                
          integer spinsour
          real*8 :: usd 
!                                                                                         
! Static variables                                                                        
!                                                                                         
          save im                                                                         
!                                                                                         
!                                                                                         
!                                                                                         
          usd=1.d0/Sqrt(2.d0)                                                                
          p=sqrt(mom(2)**2+mom(3)**2+mom(4)**2)                                           
          p3p=p+mom(4)                                                                    
          p3m=p-mom(4)                                                                    
          mp=mass+mom(1)                                                                  
          p1p=mom(2)+mom(3)*im                                                            
!                                                                                         
          if(abs(spinsour).eq.1) spinsour=-spinsour                                  
          if(abs(p3m).lt.1.d-10.or.abs(p3p).lt.1.d-10)then                                
           call fs(spinsour,fieldaux ,mom,mass)                                           
           return                                                                         
          endif                                                                           
!                                                                                         
          coeffp=1./Sqrt(2.*p*mp*p3p)                                                     
          coeffm=1./Sqrt(2.*p*mp*p3m)                                                     
!                                                                                         
! "spin up" ingoing fermion                                                               
!                                                                                         
           if (spinsour.eq.1) then                                                       
            fieldaux(1)=usd*coeffp*p3p*(mp-p)
            fieldaux(2)=usd*coeffp*p1p*(mp-p)
            fieldaux(3)=usd*coeffp*p3p*(p+mp)
            fieldaux(4)=usd*coeffp*p1p*(p+mp)
           endif                                                                          
!                                                                                         
! "spin down" ingoing fermion                                                             
!                                                                                         
           if (spinsour.eq.-1) then                                                      
            fieldaux(1)=usd*coeffm*p3m*(mp+p)
            fieldaux(2)=-usd*coeffm*p1p*(mp+p)
            fieldaux(3)=-usd*coeffm*p3m*(p-mp)
            fieldaux(4)=usd*coeffm*p1p*(p-mp)
           endif                                                                          
!                                                                                         
! "spin up" outgoing antifermion                                                          
!                                                                                         
           if (spinsour.eq.2) then                                                       
            fieldaux(1)=-usd*coeffm*p3m*(p+mp)
            fieldaux(2)=usd*coeffm*p1p*(p+mp)
            fieldaux(3)=usd*coeffm*p3m*(mp-p)
            fieldaux(4)=-usd*coeffm*p1p*(mp-p)
           endif                                                                          
!                                                                                         
! "spin down" outgoing antifermion                                                        
!                                                                                         
           if (spinsour.eq.-2) then                                                      
            fieldaux(1)=usd*coeffp*p3p*(p-mp)
            fieldaux(2)=usd*coeffp*p1p*(p-mp)
            fieldaux(3)=usd*coeffp*p3p*(mp+p)
            fieldaux(4)=usd*coeffp*p1p*(mp+p)
           endif                                                                          
!                                                                                         
           return                                                                         
           end                                                                            
                                                                                          
!***********************************************************************                  
          subroutine fermionbarsources(spinsour,fieldaux ,mom,mass)                      
!***********************************************************************                  
!                                                                                         
! This subroutine return the initialized fermionbar field.                                
! Elicity heigenstates                                                                    
!                                                                                         
          implicit none                                                                   
!                                                                                         
          real*8  mass         !containing the mass of the external particle              
          real*8  mom(4)          !array containing the momenta of the external           
          real*8 p,p3p,p3m,mp,coeffp,coeffm                                               
!                                                                                         
          complex*16 fieldaux(4)     !array returning the fermion field configur          
          complex*16 im  !immaginary unity in complex representation                      
          complex*16 p1p                                                                  
          data im/(0.,1.)/                                                                
          integer spinsour                                                                
!                                                                                         
! Static variables                                                                        
!                                                                                         
          save im                                                                         
!                                                                                         
!                                                                                         
          p=sqrt(mom(2)**2+mom(3)**2+mom(4)**2)                                           
          p3p=p+mom(4)                                                                    
          p3m=p-mom(4)                                                                    
          mp=mass+mom(1)                                                                  
          p1p=mom(2)-mom(3)*im                                                            
!                                                                                         
          if(abs(spinsour).eq.1) spinsour=-spinsour                                  
          if(abs(p3m).lt.1.d-10.or.abs(p3p).lt.1.d-10)then                                
           call fbs(spinsour,fieldaux ,mom,mass)                                          
           return                                                                         
          endif                                                                           
!                                                                                         
          coeffp=1./Sqrt(p*mp*p3p)/2.d0                                                     
          coeffm=1./Sqrt(p*mp*p3m)/2.d0                                                     
!                                                                                         
! "spin up" ingoing fermion                                                               
!                                                                                         
           if (spinsour.eq.1) then                                                       
            fieldaux(1)=coeffp*p3p*(mp+p)
            fieldaux(2)=coeffp*p1p*(mp+p)  
            fieldaux(3)=-coeffp*p3p*(p-mp)
            fieldaux(4)=-coeffp*p1p*(p-mp)
           endif                                                                          
!                                                                                         
! "spin down" ingoing fermion                                                             
!                                                                                         
           if (spinsour.eq.-1) then                                                      
            fieldaux(1)=coeffm*p3m*(mp-p)
            fieldaux(2)=-coeffm*p1p*(mp-p)
            fieldaux(3)=coeffm*p3m*(p+mp)
            fieldaux(4)=-coeffm*p1p*(p+mp)
           endif                                                                          
!                                                                                         
! "spin up" outgoing antifermion                                                          
!                                                                                         
           if (spinsour.eq.2) then                                                       
            fieldaux(1)=-coeffm*p3m*(p-mp)
            fieldaux(2)=coeffm*p1p*(p-mp)
            fieldaux(3)=-coeffm*p3m*(mp+p)
            fieldaux(4)=coeffm*p1p*(mp+p)
           endif                                                                          
!                                                                                         
! "spin down" outgoing antifermion                                                        
!                                                                                         
           if (spinsour.eq.-2) then                                                      
            fieldaux(1)=coeffp*p3p*(p+mp)
            fieldaux(2)=coeffp*p1p*(p+mp)
            fieldaux(3)=-coeffp*p3p*(mp-p)
            fieldaux(4)=-coeffp*p1p*(mp-p)
           endif                                                                          
!                                                                                         
           return                                                                         
           end                                                                            
!***********************************************************************                  
          subroutine fs(spinsour,fieldaux ,mom,mass)                                      
!***********************************************************************                  
!                                                                                         
! This subroutine return the initialized fermion field.                                   
!                                                                                         
          implicit none                                                                   
!                                                                                         
          real*8  mass         !containing the mass of the external particle              
          real*8  mom(4)       !array containing the momenta of the external par          
          integer spinsour    !array containing the spin of the source                    
!                                                                                         
          complex*16 fieldaux(4)     !array returning the fermion field configur          
          complex*16 im,f1,f2,f3,f4  !immaginary unity in complex representation                      
          data im/(0.,1.)/                                                                
!                                                                                         
! Static variables                                                                        
!                                                                                         
          save im                                                                         
!                                                                                         
!                                                                                         
! "spin up" ingoing fermion                                                               
!                                                                                         
           if(mom(4).lt.0) spinsour=-spinsour                                             
           if (spinsour.eq.1) then                                                       
            f1=sqrt((mom(1)+mass))                                               
            f2=0.                                                                
            f3=mom(4)/sqrt((mass+mom(1)))                                        
            f4=mom(2)/sqrt((mass+mom(1))) + mom(3)/sqrt((mass+mom(1)))*im        
!                                                                                         
! "spin down" ingoing fermion                                                             
!                                                                                         
           elseif (spinsour.eq.-1) then                                                  
            f1=0.                                                                
            f2=sqrt((mom(1)+mass))                                               
            f3= mom(2)/sqrt((mass+mom(1)))- mom(3)/sqrt((mass+mom(1)))*im        
            f4=-mom(4)/sqrt((mass+mom(1)))                                       
!                                                                                         
! "spin up" outgoing antifermion                                                          
!                                                                                         
           elseif (spinsour.eq.2) then                                                   
            f1=mom(4)/sqrt((mass+mom(1)))                                        
            f2=mom(2)/sqrt((mass+mom(1))) + mom(3)/sqrt((mass+mom(1)))*im        
            f3=sqrt((mom(1)+mass))                                               
            f4=0.                                                                
!                                                                                         
! "spin down" outgoing antifermion                                                        
!                                                                                         
           elseif (spinsour.eq.-2) then                                                  
            f1= mom(2)/sqrt((mass+mom(1))) - mom(3)/sqrt((mass+mom(1)))*im       
            f2=-mom(4)/sqrt((mass+mom(1)))                                       
            f3=0.                                                                
            f4=sqrt((mom(1)+mass))                                               
           endif                                                                          
!
           fieldaux(1)=sqrt(1.d0/2.d0)*(f1-f3)                                              
           fieldaux(2)=sqrt(1.d0/2.d0)*(f2-f4)                                              
           fieldaux(3)=sqrt(1.d0/2.d0)*(f1+f3)                                              
           fieldaux(4)=sqrt(1.d0/2.d0)*(f2+f4)                                              
!                                                                                         
           return                                                                         
           end                                                                            
!***********************************************************************                  
          subroutine fbs(spinsour,fieldaux,mom,mass)                                      
!***********************************************************************                  
!                                                                                         
! This subroutine return the initialized fermion field.                                   
!                                                                                         
          implicit none                                                                   
!                                                                                         
          real*8  mass        !containing the mass of the external particle               
          real*8  mom(4)      !array containing the momenta of the external part          
          integer spinsour    !array containing the spin of the source                    
!                                                                                         
          complex*16 fieldaux(4)   !array returning the fermionbar field configu          
          complex*16 im,f1,f2,f3,f4          !immaginary unity in complex representation             
          data im/(0.,1.)/                                                                
!                                                                                         
! Static variables                                                                        
!                                                                                         
          save im                                                                         
!                                                                                         
!                                                                                         
! "spin up" outgoing fermion                                                              
!                                                                                         
           if(mom(4).gt.0) spinsour=-spinsour                                             
           if (spinsour.eq.1) then                                                       
            f1=sqrt((mom(1)+mass))                                               
            f2=0.                                                                
            f3=-mom(4)/sqrt((mass+mom(1)))                                       
            f4=-mom(2)/sqrt((mass+mom(1)))+mom(3)/sqrt((mass+mom(1)))*im         
!                                                                                         
! "spin down" outgoing fermion                                                            
!                                                                                         
           elseif (spinsour.eq.-1) then                                                  
            f1=0.                                                                
            f2=sqrt((mom(1)+mass))                                               
            f3= -mom(2)/sqrt((mass+mom(1)))-mom(3)/sqrt((mass+mom(1)))*im        
            f4=mom(4)/sqrt((mass+mom(1)))                                        
!                                                                                         
! "spin up" ingoing antifermion                                                           
!                                                                                         
           elseif (spinsour.eq.2) then                                                   
            f1=-mom(4)/sqrt((mass+mom(1)))                                       
            f2=-mom(2)/sqrt((mass+mom(1)))+mom(3)/sqrt((mass+mom(1)))*im         
            f3=sqrt((mom(1)+mass))                                               
            f4=0.                                                                
!                                                                                         
! "spin down" ingoing antifermion                                                         
!                                                                                         
           elseif (spinsour.eq.-2) then                                                  
            f1= -mom(2)/sqrt((mass+mom(1))) - mom(3)/sqrt((mass+mom(1)))*im      
            f2=mom(4)/sqrt((mass+mom(1)))                                        
            f3=0.                                                                
            f4=sqrt((mom(1)+mass))                                               
           endif                                                                          
!                                                                                         
           fieldaux(1)=sqrt(1.d0/2.d0)*(f1-f3)                                              
           fieldaux(2)=sqrt(1.d0/2.d0)*(f2-f4)                                              
           fieldaux(3)=sqrt(1.d0/2.d0)*(f1+f3)                                              
           fieldaux(4)=sqrt(1.d0/2.d0)*(f2+f4)                                              
!                                                                                         
           return                                                                         
           end                                                                            



