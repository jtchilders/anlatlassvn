!
subroutine inital
!
  use couplings
  use prcint
!
  implicit none
  real*8 :: apar(100)
  common/alphapar/apar
!
!
! initializing on/off interactions
!!$  tabint(1:61,1)=(/20,qupl,qupbrl,gluon, qupl,qupbrl,photon,qupl,qupbrl,z, &  
!!$                      qupr,qupbrr,gluon, qupr,qupbrr,photon,qupr,qupbrr,z, &  
!!$                      qdr,qdbrr,gluon, qdr,qdbrr,photon,qdr,qdbrr,z, &  
!!$                      qdl,qdbrl,gluon, qdl,qdbrl,photon,qdl,qdbrl,z, &  
!!$                      qbt,qbtbr,gluon, qbt,qbtbr,photon,qbt,qbtbr,z, &  
!!$                      ebrr,z,er,ebrr,photon,er,ebrl,z,el,ebrl,photon,el, &
!!$                                                       gluon,gluon,gluon /)       
  tabint(1:43,1)=(/14,qupl,qupbrl,gluon, qupr,qupbrr,gluon, qdr,qdbrr,gluon, &
                      qdl,qdbrl,gluon,qcr,qcbrr,gluon,qcl,qcbrl,gluon,       &
                      qsr,qsbrr,gluon,qsl,qsbrl,gluon,qbt,qbtbr,gluon,       &
                      qtop,qtopbr,gluon, gluon,gluon,gluon,                  & 
                      higgs,qtop,qtopbr,higgs,qbt,qbtbr,higgs,higgs,higgs /)       ! the first entry is

!    the number of relevant interaction, the following entries in group of three label the particles
!    entering the lagrangian interaction term (the ordering is immaterial)
!
!
! initializing coupling constants
  gstrong=0.d0; g2weak=0.d0; masses=0.d0; width=0.d0
!
  gstrong=1.d0
  g2weak= apar(51)*sqrt(2.d0); gbar=apar(55)*sqrt(2.d0); ctw=apar(53) ; stw= apar(54);
  masses(wp)=apar(3); masses(wm)=masses(wp); width(wp)=apar(4);  width(wm)=width(wp)
  masses(z)=apar(1) ; width(z)=apar(2) ; masses(higgs)=apar(5); width(higgs)=apar(6);
  masses(qbt)= apar(15);   masses(qbtbr)= masses(qbt);  masses(qtop)=apar(16); 
  masses(qtopbr)=masses(qtop); 
!
end subroutine inital
!
!
subroutine matrix0(flvmlm,posi,imp,hel,result)
! 
! translating old alpha input into new alpha one
!
  use incnst
  use alptyp
  use naming
  use dimen
  use extsrc
  use prcint, only : nprc
!
  implicit none
!
  integer, intent (in) :: flvmlm(nmax),posi(2), hel(nmax)
  real*8, intent (in) :: imp(4*nmax)
  complex*16, intent (out) :: result
!
  integer :: j1,nfr,flv(nmax),npr,nlp,nlepbr,spn,hl(nmax)
  real*8 :: impul(4*nmax),plep(4),plepbr(4),mw,ww
  complex*16 :: lp(4),lpbr(4),w(4),wpol(4)
  integer :: clr(2*nmax)
  integer color(2*nmax)     !coulor string
  common/colore/color
  type (proc) :: prcss
  integer,save :: convert (-100:100), cnvmlsl(-16:16), cnvmlsr(-16:16),dim(7)
  logical,save :: massless (-16:16)
  real*8 :: apar(100)
  common/alphapar/apar
  integer flgdual          !dual (0) or su3 (1) amplitudes
  common/dual/flgdual
  integer rep(nmax),nprt,nglu,nqrk,nlep,ngb,nphot
!  common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot     
! 
  data convert/75*-101,higgs,wm,z,photon,-101, 8*-101,nuebr, 5*-101,qtopbr,qbtbr,qcbr,-101,qupbr,qdbr, gluon, &
           qd,qup,-101,qc,qbt,qtop,5*-101,nue,8*-101, -101,photon,z,wp,higgs,75*-101/ 
  data cnvmlsl /4*-101,nuebr,ebrl,6*-101,qcbrl,qsbrl,qupbrl,qdbrl, gluon, &
               qdl,qupl,qsl,qcl,6*-101, el,nue, 4*-101/
  data cnvmlsr /5*-101,ebrr,6*-101,qcbrr,qsbrr,qupbrr,qdbrr, gluon, &
                qdr,qupr,qsr,qcr,6*-101,er,  5*-101/
  data massless /4*.false.,2*.true., 6*.false.,4*.true., .false., &
                  4*.true., 6*.false., 2*.true.,4*.false./
!
  nprc=1
!
  impul=-imp
  impul(4*posi(1)-3:4*posi(1))=-impul(4*posi(1)-3:4*posi(1))
  impul(4*posi(2)-3:4*posi(2))=-impul(4*posi(2)-3:4*posi(2))
  nfr=0
  clr=color
  hl=hel
  flv=-flvmlm
  flv(posi(1))=-flv(posi(1))
  flv(posi(2))=-flv(posi(2))
  nfr=0
  npr=0
  do j1=1,nmax
    if(abs(flvmlm(j1)).lt.17) then
      if(flvmlm(j1).ne.0) then
        nfr=nfr+1
        if(flvmlm(j1).lt.0) hl(j1)=2*hl(j1)
      endif
    endif
    if(flvmlm(j1).ne.1001) npr=npr+1
  enddo
!
  flv(npr+1:nmax)=-1001
!!$  write(*,*) 'fl',flv
!!$  write(*,*) 'hl',hl
!
  exflg=.false.
  do j1=1,nmax
    if(flv(j1).eq.-1001) then
      flv(j1)=1001
    else
      if(abs(flv(j1)).lt.17) then
        if(massless(flv(j1))) then
          if(hl(j1).gt.0) then
            flv(j1)=cnvmlsl(flv(j1))
          else
            flv(j1)=cnvmlsr(flv(j1))
          endif
        else
          flv(j1)=convert(flv(j1))
        endif
      else
        flv(j1)=convert(flv(j1))
      endif
    endif
  enddo
!
  prcss%flv=flv; prcss%hel=hl; prcss%col=clr; prcss%nfrm=nfr; prcss%mom=impul
!
!  write(*,*) prcss
!
  chdim : select case (npr)
    case (5)
      dim(d_offsh)=80; dim(d_offshcol)=50; dim(d_ampfld)=30; dim(d_ampcol)=20
    case (6)
      dim(d_offsh)=180; dim(d_offshcol)=120; dim(d_ampfld)=50; dim(d_ampcol)=30
    case (7)
      dim(d_offsh)=400; dim(d_offshcol)=260; dim(d_ampfld)=220; dim(d_ampcol)=110
    case (8)
      dim(d_offsh)=800; dim(d_offshcol)=520; dim(d_ampfld)=430; dim(d_ampcol)=220
    case (9)
      dim(d_offsh)=1540; dim(d_offshcol)=1000; dim(d_ampfld)=820; dim(d_ampcol)=410
    case default
      write(*,*)'too large/small number of particles', npr
  end select chdim
  dim(d_psnew)=2**npr
  dim(d_ordtmp)=2**npr
  dim(d_impul)=2**npr
!
  if(flgdual.eq.1) then
    call amp(prcss,dim,result)
  elseif(flgdual.eq.0.or.flgdual.eq.2) then
    call ampd(prcss,dim,result)
  else
    write(*,*)'failure in INITALP'
    stop
  endif
!     
end subroutine matrix0






