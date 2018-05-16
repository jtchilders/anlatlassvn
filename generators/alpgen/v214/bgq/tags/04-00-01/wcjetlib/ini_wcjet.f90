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
!  tabint(1:19,1)=(/6,qupl,qupbrl,gluon,qdl,qdbrl,gluon,qbt,qbtbr,gluon,qdbrl,wm,qupl,nuebr,wp,el,gluon,gluon,gluon/)       ! the first entry is
  tabint(1:34,1)=(/11,qupl,qupbrl,gluon,qdl,qdbrl,gluon,qbt,qbtbr,gluon,qupbrl,wp,qdl, &
                   ebrl,wm,nue,gluon,gluon,gluon,qcbrr,gluon,qcr,qcbrl,gluon,qcl,     &
                   qupbrr,gluon,qupr,qdbrr,gluon,qdr,gluon,qtop,qtopbr/)       ! the first entry is
!    the number of relevant interaction, the following entries in group of three label the particles
!    entering the lagrangian interaction term (the ordering is immaterial)
!
!
! initializing coupling constants
  gstrong=0.d0; g2weak=0.d0; masses=0.d0; width=0.d0
!
  gstrong=1.d0
  g2weak= apar(51)*sqrt(2.d0); gbar=apar(55)*sqrt(2.d0); ctw=apar(53) ; stw= apar(54);
  masses(wp)=apar(3); masses(wm)=masses(wp); width(wp)=apar(4);  width(wm)=width(wp);
  masses(z)=apar(1) ; width(z)=apar(2); 
  masses(qbt)= apar(15);  masses(higgs)=apar(5); width(higgs)=apar(6);
  masses(qbtbr)= masses(qbt); masses(qtop)=apar(16);
  masses(qtopbr)=masses(qtop);
!  ctw=1.d0           !not used, but coupling initialization has 1/ctw in it
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
  integer :: j1,nfr,flv(nmax),dim(7),npr,nlp,nlepbr,spn,hl(nmax)
  real*8 :: impul(4*nmax),plep(4),plepbr(4),mw,ww
  complex*16 :: lp(4),lpbr(4),w(4),wpol(4)
  integer :: clr(2*nmax)
  integer color(2*nmax)     !coulor string
  common/colore/color
  type (proc) :: prcss
  integer,save :: convert (-100:100), cnvmlsl(-16:16), cnvmlsr(-16:16)
  logical,save :: massless (-16:16)
  real*8 :: apar(100)
  common/alphapar/apar
  integer flgdual          !dual (0) or su3 (1) amplitudes
  common/dual/flgdual
  integer rep(nmax),nprt,nglu,nqrk,nlep,ngb,nphot
!  common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot     
! 
  data convert/76*-101,wm,3*-101, 8*-101,nuebr, 5*-101,qtopbr,qbtbr,qcbr,-101,qupbr,qdbr, gluon, &
           qd,qup,-101,qc,qbt,qtop,5*-101,nue,8*-101, 3*-101,wp,76*-101/ 
  data cnvmlsl /4*-101,nuebr,ebrl,6*-101,qcbrl,qsbrl,qupbrl,qdbrl, gluon, &
               qdl,qupl,qsl,qcl,6*-101, el,nue, 4*-101/
  data cnvmlsr /6*-101,6*-101,qcbrr,qsbrr,qupbrr,qdbrr, gluon, &
                qdr,qupr,qsr,qcr,6*-101,  6*-101/
  data massless /4*.false.,2*.true., 6*.false.,4*.true., .false., &
                  4*.true.,6*.false., 2*.true.,4*.false./
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
    if(flvmlm(j1).eq.-11.or.flvmlm(j1).eq.-12)nlp=j1
    if(flvmlm(j1).eq.11.or.flvmlm(j1).eq.12)nlepbr=j1
    if(flvmlm(j1).ne.1001) npr=npr+1
  enddo
!
!!$  write(*,*)flv
!!$  write(*,*)hel
!!$  write(*,*)color
!!$  write(*,*)impul(1:4*npr)
!!$  write(*,*)
!!$!
  plep=impul(4*nlp-3:4*nlp)
  plepbr=impul(4*nlepbr-3:4*nlepbr)
  spn=2
  call fermionsourceshl(spn,lp,-plep)
  spn=1
  call fermionbarsourceshl(spn,lpbr,-plepbr)
  call fbftowl (lpbr,lp,w)
!!$  write(*,*)lpbr
!!$  write(*,*)lp
!!$  write(*,*)w
  mw=apar(3)
  ww=apar(4)
  call gbms (plep+plepbr,w,mw,ww,wpol)
!!$  write(*,*)wpol,' ',mw,' ',ww
  wpol=wpol*apar(51)/sqrt(2.d0)
!!$  write(*,*)wpol,' ',mw,' ',ww
  do j1=nlp+1,nlepbr-1
    flv(j1-1)=flv(j1)
    hl(j1-1)=hl(j1)
    clr(2*j1-3:2*j1-2)=clr(2*j1-1:2*j1)
    impul(4*j1-7:4*j1-4)=impul(4*j1-3:4*j1)
!    rep(j1-1)=rep(j1)
  enddo
  do j1=nlepbr+1,npr
    flv(j1-2)=flv(j1)
    hl(j1-2)=hl(j1)
    clr(2*j1-5:2*j1-4)=clr(2*j1-1:2*j1)
    impul(4*j1-11:4*j1-8)=impul(4*j1-3:4*j1)
!    rep(j1-2)=rep(j1)
  enddo
  if (flvmlm(nlp).eq.-11) then
    flv(npr-1)=-24
  else
    flv(npr-1)=24
  endif
!  rep(npr-1)=0
!  rep(npr:nmax)=99
  flv(npr:nmax)=-1001
  npr=npr-1
  nfr=nfr-2
!  ngb=ngb+1
  nlep=nlep-2
  clr(2*npr-1:2*npr)=0
  exflg=.true.
  src(1:4)=wpol
  impul(4*npr-3:4*npr)=plep+plepbr
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
!!$  write(*,*)flv
!!$  write(*,*)hel
!!$  write(*,*)color
!!$  write(*,*)impul(1:4*npr)
!!$  write(*,*) wp,wm,'    ',plep,plepbr
!!$  write(*,*)
!
     if(npr.eq.9) then 
      dim(d_offshcol)=650; dim(d_offsh)=780; dim(d_ampcol)=220; dim(d_ampfld)=  433 
     elseif(npr.eq.8) then 
      dim(d_offshcol)=330; dim(d_offsh)=410; dim(d_ampcol)=120; dim(d_ampfld)=  230 
     elseif(npr.eq.7) then 
      dim(d_offshcol)=170; dim(d_offsh)=220; dim(d_ampcol)=70; dim(d_ampfld)=  130 
     elseif(npr.eq.6) then 
      dim(d_offshcol)=90; dim(d_offsh)=110; dim(d_ampcol)=40; dim(d_ampfld)=  70 
     elseif(npr.eq.5) then 
      dim(d_offshcol)=90; dim(d_offsh)=110; dim(d_ampcol)=40; dim(d_ampfld)=  40 
     elseif(npr.eq.4) then 
      dim(d_offshcol)=30; dim(d_offsh)=40; dim(d_ampcol)=10; dim(d_ampfld)=  20
     elseif(npr.eq.3) then 
      dim(d_offshcol)=20; dim(d_offsh)=10; dim(d_ampcol)=10; dim(d_ampfld)=  10
     else 
      write(*,*)'wrong particle number in MATRIX0' 
      stop 
     endif 
! 
      dim(d_impul)= 2**npr-1 
      dim(d_psnew)= 2**npr-1 
      dim(d_ordtmp)=2**npr-1 
! 
  prcss%flv=flv; prcss%hel=hl; prcss%col=clr; prcss%nfrm=nfr; prcss%mom=impul
!!$  write(*,*) prcss
!!$  write(7,*)posi,'    ',flv
!!$  write(7,*)color
!!$!
!!$  write(*,*) 'flg',flgdual

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






