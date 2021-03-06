* ffS3n.F
* calculate S3n = \int_0^1 dy (ra y^2 + cb y + cc + I signc)/(y - cy)
* where ra can be zero.
* input: cy=y0, ra=a (real), cb=b, cc=c
* signc=sign(im(c)), signy=sign(im(cy)) in case they are real.
* cza and czb are the 2 roots of: a y^2 + b y + c == 0
* remarks: ieps is needed for cza, czb and y0.
* this file is part of LoopTools
* last modified 8 Dec 10 th

* Written by Le Duc Ninh, MPI, Munich (2008).
* Spence, log and eta functions are taken from FF.
* Oct 27 2008

#include "externals.h"


	double complex function ffS3n(cy, signy, ra, cb, cc, signc,
     &    ier)
	implicit none
	double precision ra, signy, signc
	double complex cy, cb, cc
	integer ier

#include "ff.h"

	double complex cl, crdisc, cza, czb
	double precision sza, szb, sy1, sy2, sc

	double complex ffRn, zfflog
	integer nffet1
	external ffRn, zfflog, nffet1
	
* check for end-point sing.
	if( abs(cy) .lt. precx .or. abs(cy - 1) .lt. precx ) then
	  call fferr(90, ier)
	  ffS3n = 0
	  return
	endif

	cl = zfflog((cy - 1)/cy, 1, DCMPLX(signy), ier)

	sc = DIMAG(cc)
	if( sc .eq. 0 ) sc = signc

	if( abs(ra) .lt. precx ) then
	  if( abs(cb) .lt. precx ) then
* 0 roots:
	    if( abs(cc) .lt. precx ) then
	       call fferr(91, ier)
	       ffS3n = 0
	       return
	    endif
	    ffS3n = cl*zfflog(cc, 1, DCMPLX(signc), ier)
	    return
	  endif

* 1 root:
* eq.: b y + c == 0
	  cza = -cc/cb
	  sza = -signc*DBLE(cb)
	  if( sza .eq. 0 ) sza = -signc

	  ffS3n = cl*zfflog(cb, 1, DCMPLX(signc), ier) +
     &      ffRn(cy, signy, cza, sza, ier)

	  if( abs(DIMAG(cb)) .gt. precx ) then
	    szb = DIMAG(cza)
	    if( szb .eq. 0 ) szb = sza
	    ffS3n = ffS3n + cl*c2ipi*
     &        nffet1(cb, DCMPLX(0D0, -szb), DCMPLX(0D0, sc), ier)
	  endif
	  return
	endif

* 2 roots: cza = y1, czb = y2
* eq.: y**2 + (b/a) y + (c/a) = 0
* the ieps is irrelevant here since we take into account
* the contributions of both roots

*** Ninh changed: 14 Aug 2009
	crdisc = sqrt(cb**2/ra**2 - 4*cc/ra)
	cza = -.5D0*(cb/ra + crdisc)
	czb = -.5D0*(cb/ra - crdisc)
	if( abs(cza) .gt. abs(czb) ) then
	  czb = cc/(ra*cza)
	else if( abs(czb) .gt. 1D-13 ) then
	  cza = cc/(ra*czb)
	endif

* calculate the sign of im(cza) and im(czb) which are related to ieps
	sza = sc/ra
	if( abs(DBLE(crdisc)) .gt. precx ) sza = sza/DBLE(crdisc)
	szb = -sza

	sy1 = DIMAG(cza)
	if( sy1 .eq. 0 ) sy1 = sza
	sy2 = DIMAG(czb)
	if( sy2 .eq. 0 ) sy2 = szb

* calculate the log, etas, and the 2 R-functions
* ieps=1 to choose the cut along the real axis

	ffS3n =
     &    cl*( zfflog(DCMPLX(ra), 1, DCMPLX(sc), ier) +
     &      c2ipi*nffet1(DCMPLX(0D0, -sy1), DCMPLX(0D0, -sy2),
     &        DCMPLX(0D0, sc/ra), ier) ) +
     &    ffRn(cy, signy, cza, sza, ier) +
     &    ffRn(cy, signy, czb, szb, ier)
	end

