!*****************************************************************************

!   MPFUN: A MULTIPLE PRECISION FLOATING POINT COMPUTATION PACKAGE

!   IEEE Fortran-90 version
!   Version Date:  2001-08-02

!   Author:

!      David H. Bailey
!      NERSC, Lawrence Berkeley Lab
!      Mail Stop 50B-2239
!      Berkeley, CA 94720
!      Email: dhbailey@lbl.gov

!   Restrictions:

!   This software was written while the author was an employee of NASA.
!   This software has been approved by NASA for unrestricted distribution.
!   However, usage of this software is subject to the following:

!   1. This software is offered without warranty of any kind, either expressed
!      or implied.  The author would appreciate, however, any reports of bugs
!      or other difficulties that may be encountered.
!   2. If modifications or enhancements to this software are made to this
!      software by others, NASA Ames reserves the right to obtain this enhanced
!      software at no cost and with no restrictions on its usage.
!   3. The author and NASA Ames are to be acknowledged in any published paper
!      based on computations using this software.  Accounts of practical
!      applications or other benefits resulting from this software are of
!      particular interest.  Please send a copy of such papers to the author.

!   Description:

!   The following information is a brief description of this program.  For
!   full details and instructions for usage, see the paper "A Portable High
!   Performance Multiprecision Package", available from the author.

!   This package of Fortran subroutines performs multiprecision floating point
!   arithmetic.  If sufficient main memory is available, the maximum precision
!   level is at least 16 million digits.  The maximum dynamic range is at
!   least 10^(+-14,000,000).  It employs advanced algorithms, including an
!   FFT-based multiplication routine and some recently discovered
!   quadratically convergent algorithms for pi, exp and log.  The package also
!   features extensive debug and self-checking facilities, so that it can be
!   used as a rigorous system integrity test.  All of the routines in this
!   package have been written to facilitate vector, parallel processing and/or
!   RISC processing.

!   Most users will not wish to manually write code that directly calls these
!   routines, since this is a tedious and error-prone process.  To assist
!   such users, the author has prepared Fortran-90 modules that permit one
!   to write ordinary Fortran program, yet have the required MPFUN routines
!   automatically called.  Contact the author for details.

!   Machine-specific tuning notes may be located by searching for the text
!   string !> in this program file.  It is highly recommended that these notes
!   be read before running this package on a specific system.  Certain
!   vectorizable DO loops that are often not recognized as such by vectorizing
!   compilers are prefaced with Cray !dir$ ivdep directives.  On other vector
!   systems these directives should be replaced by the appropriate equivalents.
!   Double precision should be disabled when compiling this code on Cray
!   systems (use the -dp flag).

!   Instructions for compiling and testing this program on various specific
!   systems are included in the readme file that accompanies this file.

!*****************************************************************************

module mpfuna

!   This section initializes global parameters and arrays default values.
!   Note that this Fortran-90 version completely dispenses with common blocks.
!   The global parameters previously in common block MPCOM1 are now global
!   parameters in this module and can be accessed in any user program with a
!   USE MPFUNMOD statement.  Note however the names now start with 'MP'.  The
!   user also no longer needs to allocate scratch space in common MPCOM3,
!   MPCOM4 or MPCOM5, since the required space is now allocated dynamically
!   by each routine as required.

double precision:: mpbbx, mpbdx, mpbx2, mprbx, mprdx, mprx2, mprxx
parameter (mpkdp = kind (0.d0))
complex (mpkdp), allocatable:: mpuu1, mpuu2
!>
!   The parameters MPBBX, MPNBT, MPNPR and MPMCRX depend on system word size.
!   On IEEE systems and most other 32 bit systems, set MPBBX = 4096.D0,
!   MPNBT = 24, MPNPR = 32, and MPMCRX = 7.  On Cray (non-IEEE) systems, set 
!   MPBBX = 2048.D0, MPNBT = 22, MPNPR = 16, and MPMCRX = 8.  The parameters 
!   MPNROW, MPNSP1 and MPNSP2 are spacing parameters to avoid bank and cache 
!   conflict performance problems in the FFT routines.  MPNROW = 16, 
!   MPNSP1 = 2 and MPNSP2 = 9 appear to work well on most systems.  Set 
!   MPNROW = 64 on Cray vector systems.

parameter (mpbbx = 4096.d0, mpnbt = 24, mpnpr = 32, mpmcrx = 7, &
  mpbdx = mpbbx ** 2, mpbx2 = mpbdx ** 2, mprbx = 1.d0 / mpbbx, &
  mprdx = mprbx ** 2, mprx2 = mprdx ** 2, mprxx = 16.d0 * mprx2, &
  mpnrow = 16, mpnsp1 = 2, mpnsp2 = 9)
dimension mpker(73), mpuu1(:), mpuu2(:)

data mpnw, mpidb, mpldb, mpndb, mpier, mpmcr, mpird / 16, 0, 6, 22, 99, &
  mpmcrx, 1/
data mpker /73 * 2/

contains

subroutine mpabrt
!>
!   This routine terminates execution.  Users may wish to replace the
!   default STOP with a call to a system routine that provides a traceback.
!   Examples of code that produce traceback are included here (commented out)
!   for some systems.

if (mpier .eq. 99) then
  write (mpldb, 1)
1 format ('*** The MPFUN library has not been initialized.  If you are using'/&
  'the Fortran-90 translation modules, you must insert the following line'/ &
  'at the start of execution in your main program:'/ ' '/ 'CALL MPINIT')
else
  write (mpldb, 2) mpier
2 format ('*** MPABRT: Execution terminated, error code =',i4)
endif

!   Use this line on Cray systems.

! call abort

!   On other systems, merely terminate execution.

stop
end subroutine

end module

module mpfunb

!   This module defines the 'dp' routines, i.e. routines that operate on 
!   numbers of the DPE (double precision plus exponent) datatype.

use mpfuna
contains

subroutine dpadd (a, na, b, nb, c, nc)

!   This adds the DPE numbers (A, NA) and (B, NB) to yield the sum (C, NC).

implicit double precision (a-h, o-z)
dimension pt(64)
save pt
data pt/ 64 * 0.d0/

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c = 0.d0
  nc = 0
  return
endif

!   If this is the first call to DPADD, initialize the PT table.

if (pt(1) .eq. 0.d0) then
  pt(1) = 0.5d0

  do i = 2, 64
    pt(i) = 0.5d0 * pt(i-1)
  enddo
endif

!   This operation reduces to five cases.

if (b .eq. 0.d0) then
  c = a
  nc = na
else if (a .eq. 0.d0) then
  c = b
  nc = nb
else if (na .eq. nb) then
  c = a + b
  nc = na
else if (na .gt. nb) then
  k = na - nb
  nc = na
  if (k .gt. 64) then
    c = a
  else
    c = a + b * pt(k)
  endif
else
  k = nb - na
  nc = nb
  if (k .gt. 64) then
    c = b
  else
    c = b + a * pt(k)
  endif
endif
if (c .eq. 0.d0) then
  nc = 0
  goto 130
endif

!   Normalize the result to a decent range if it is not.

110  if (abs (c) .ge. mpbdx) then
  c = mprdx * c
  nc = nc + mpnbt
  goto 110
endif

120  if (abs (c) .lt. 1.d0) then
  c = mpbdx * c
  nc = nc - mpnbt
  goto 120
endif

130  return
end subroutine

subroutine dpdec (a, na, b, nb)

!   This converts the DPE number (A, NA) to decimal form, i.e. B * 10^NB,
!   where |B| is between 1 and 10.

implicit double precision (a-h, o-z)
parameter (xlt = 0.3010299956639812d0)

if (a .ne. 0.d0) then
  t1 = xlt * na + log10 (abs (a))
  nb = t1
  if (t1 .lt. 0.d0) nb = nb - 1
  b = sign (10.d0 ** (t1 - nb), a)
else
  b = 0.d0
  nb = 0
endif

return
end subroutine

subroutine dpdiv (a, na, b, nb, c, nc)

!   This divides the DPE number (A, NA) by (B, NB) to yield the quotient
!   (C, NC).

implicit double precision (a-h, o-z)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c = 0.d0
  nc = 0
  return
endif
if (b .eq. 0.d0) then
  if (mpker(1) .ne. 0) then
    write (mpldb, 1)
1   format ('*** DPDIV: Divisor is zero.')
    mpier = 1
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Divide A by B and subtract exponents, unless A is zero.

if (a .eq. 0.d0) then
  c = 0.d0
  nc = 0
  goto 120
else
  c = a / b
  nc = na - nb
endif

!   Normalize the result to a decent range if it is not.

100  if (abs (c) .ge. mpbdx) then
  c = mprdx * c
  nc = nc + mpnbt
  goto 100
endif

110  if (abs (c) .lt. 1.d0) then
  c = mpbdx * c
  nc = nc - mpnbt
  goto 110
endif

120  return
end subroutine

subroutine dpmul (a, na, b, nb, c, nc)

!   This multiplies the DPE number (A, NA) by (B, NB) to yield the product
!   (C, NC).

implicit double precision (a-h, o-z)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c = 0.d0
  nc = 0
  return
endif

!   Multiply A by B and add exponents, unless either is zero.

if (a .eq. 0.d0 .or. b .eq. 0.d0) then
  c = 0.d0
  nc = 0
  goto 120
else
  c = a * b
  nc = na + nb
endif

!   Normalize the result to a decent range if it is not.

100  if (abs (c) .ge. mpbdx) then
  c = mprdx * c
  nc = nc + mpnbt
  goto 100
endif

110  if (abs (c) .lt. 1.d0) then
  c = mpbdx * c
  nc = nc - mpnbt
  goto 110
endif

120  return
end subroutine

subroutine dppwr (a, na, b, nb, c, nc)

!   This raises the DPE number (A, NA) to the (B, NB) power and places the
!   result in (C, NC).

implicit double precision (a-h, o-z)
parameter (cl2 = 1.4426950408889633d0)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c = 0.d0
  nc = 0
  return
endif
if (a .le. 0.d0) then
  if (mpker(2) .ne. 0) then
    write (mpldb, 1)
1   format ('*** DPPWR: Argument is less than or equal to zero.')
    mpier = 2
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

if (b .eq. 0.d0) then
  c = 1.d0
  nc = 0
  goto 120
endif

if (b .eq. 1.d0 .and. nb .eq. 0) then
  c = a
  nc = na
  goto 120
endif

!   Compute the base 2 logarithm of A and multiply by B.

al = cl2 * log (a) + na
call dpmul (al, 0, b, nb, t1, n1)

!   Check for possible overflow or underflow.

if (n1 .gt. 6) then
  if (t1 .gt. 0.d0) then
    if (mpker(3) .ne. 0) then
      write (mpldb, 2)
2     format ('*** DPPWR: Overflow')
      mpier = 3
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  else
    c = 0.d0
    nc = 0
    goto 120
  endif
endif

!   Compute 2 raised to the power B * Log_2 (A).

t1 = t1 * 2.d0 ** n1
nc = int (t1)
c = 2.d0 ** (t1 - nc)

!   Normalize the result to a decent range if it is not.

100  if (abs (c) .ge. mpbdx) then
  c = mprdx * c
  nc = nc + mpnbt
  goto 100
endif

110  if (abs (c) .lt. 1.d0) then
  c = mpbdx * c
  nc = nc - mpnbt
  goto 110
endif

120  return
end subroutine

subroutine dpsqrt (a, na, b, nb)

!   This computes the square root of the DPE number (A, NA) and places the
!   result in (B, NB).

implicit double precision (a-h, o-z)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b = 0.d0
  nb = 0
  return
endif
if (a .lt. 0.d0) then
  if (mpker(4) .ne. 0) then
    write (mpldb, 1)
1   format ('*** DPSQRT: Argument is negative.')
    mpier = 4
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

if (a .eq. 0.d0) then
  b = 0.d0
  nb = 0
  goto 120
endif

!   Divide the exponent of A by two and then take the square root of A.  If
!   NA is not an even number, then we have to multiply A by 10 before taking
!   the square root.

nb = na / 2
if (na .eq. 2 * nb) then
  b = sqrt (a)
else
  b = sqrt (2.d0 * a)
  if (na .lt. 0) nb = nb - 1
endif

!   Normalize the result to a decent range if it is not.

100  if (abs (b) .ge. mpbdx) then
  b = mprdx * b
  nb = nb + mpnbt
  goto 100
endif

110  if (abs (b) .lt. 1.d0) then
  b = mpbdx * b
  nb = nb - mpnbt
  goto 110
endif

120  return
end subroutine

subroutine dpsub (a, na, b, nb, c, nc)

!   This subtracts the DPE number (B, NB) from (A, NA) to yield the difference
!   (C, NC).

implicit double precision (a-h, o-z)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c = 0.d0
  nc = 0
  return
endif

bb = -b
call dpadd (a, na, bb, nb, c, nc)

return
end subroutine

end module

module mpfunc

!   This module defines basic arithmetic routines.

use mpfuna
contains

subroutine mpadd (a, b, c)

!   This routine adds MP numbers A and B to yield the MP sum C.  It attempts
!   to include all significance of A and B in the result, up to the maximum
!   mantissa length MPNW.  Debug output starts with MPIDB = 9.  This is a new
!   simplified version.

!   Max SP space for C: MPNW + 4 cells.

double precision d(mpnw+5), db
dimension a(mpnw+2), b(mpnw+2), c(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 9) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPADD I'/(6f12.0))
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 1) (b(i), i = 1, no)
endif

ia = sign (1., a(1))
ib = sign (1., b(1))
na = min (int (abs (a(1))), mpnw)
nb = min (int (abs (b(1))), mpnw)

!   Check for zero inputs.

if (na .eq. 0) then

!   A is zero -- the result is B.

  c(1) = sign (nb, ib)

  do i = 2, nb + 2
    c(i) = b(i)
  enddo

  goto 100
elseif (nb .eq. 0) then

!   B is zero -- the result is A.

  c(1) = sign (na, ia)

  do i = 2, na + 2
    c(i) = a(i)
  enddo

  goto 100
endif

if (ia .eq. ib) then
  db = 1.d0
else
  db = -1.d0
endif
ixa = a(2)
ixb = b(2)
ish = ixa - ixb

if (mpidb .ge. 9) write (6, *) 'ish =', ish

if (ish .ge. 0) then

!   A has greater exponent than B, so B must be shifted to the right.

  m1 = min (na, ish)
  m2 = min (na, nb + ish)
  m3 = na
  m4 = min (max (na, ish), mpnw + 1)
  m5 = min (max (na, nb + ish), mpnw + 1)
  d(1) = 0.d0
  d(2) = 0.d0

  do i = 1, m1
    d(i+2) = a(i+2)
  enddo

  do i = m1 + 1, m2
    d(i+2) = a(i+2) + db * b(i+2-ish)
  enddo

  do i = m2 + 1, m3
    d(i+2) = a(i+2)
  enddo

  do i = m3 + 1, m4
    d(i+2) = 0.d0
  enddo

  do i = m4 + 1, m5
    d(i+2) = db * b(i+2-ish)
  enddo

  nd = m5
  ixd = ixa
  d(nd+3) = 0.d0
  d(nd+4) = 0.d0
else

!   B has greater exponent than A, so A must be shifted to the right.

  nsh = - ish
  m1 = min (nb, nsh)
  m2 = min (nb, na + nsh)
  m3 = nb
  m4 = min (max (nb, nsh), mpnw + 1)
  m5 = min (max (nb, na + nsh), mpnw + 1)
  d(1) = 0.d0
  d(2) = 0.d0

  do i = 1, m1
    d(i+2) = db * b(i+2)
  enddo

  do i = m1 + 1, m2
    d(i+2) = a(i+2-nsh) + db * b(i+2)
  enddo

  do i = m2 + 1, m3
    d(i+2) = db * b(i+2)
  enddo

  do i = m3 + 1, m4
    d(i+2) = 0.d0
  enddo

  do i = m4 + 1, m5
    d(i+2) = a(i+2-nsh)
  enddo

  nd = m5
  ixd = ixb
  d(nd+3) = 0.d0
  d(nd+4) = 0.d0
endif

!   Call mpnorm to fix up result and store in c.

d(1) = sign (nd, ia)
d(2) = ixd
call mpnorm (d, c)
  
100 continue
if (mpidb .ge. 9) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 2) (c(i), i = 1, no)
2 format ('MPADD O'/(6f12.0))
endif
return
end subroutine

subroutine mpcbrt (a, b)

!   This computes the cube root of the MP number A and returns the MP result
!   in B.  For extra high levels of precision, use MPCBRX.  Debug output
!   starts with MPIDB = 7.

!   Max SP space for B: MPNW + 4 cells.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to A ^ (-2/3):

!    X_{k+1} = X_k + (1 - X_k^3 * A^2) * X_k / 3

!   where the muliplication () * X_k is performed with only half of the
!   normal level of precision.  These iterations are performed with a
!   maximum precision level MPNW that is dynamically changed, doubling with
!   each iteration.  The final iteration is performed as follows (this is
!   due to A. Karp):

!    Cbrt(A) = (A * X_n) + [A - (A * X_n)^3] * X_n / 3 (approx.)

!   where the multiplications A * X_n and [] * X_n are performed with only
!   half of the final level of precision.  See the comment about the parameter
!   NIT in MPDIVX.

double precision cl2, t1, t2
parameter (cl2 = 1.4426950408889633d0, nit = 3)
dimension a(mpnw+2), b(mpnw+4), f(8), s(3*mpnw+15)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 7) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPCBRT I'/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)

if (na .eq. 0) then
  b(1) = 0.
  b(2) = 0.
  goto 120
endif
if (ia .lt. 0.d0) then
  if (mpker(13) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPCBRT: Argument is negative.')
    mpier = 13
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
nws = mpnw

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Compute A^2 outside of the iteration loop.

mpnw = nws + 1
call mpmul (a, a, s(k0))

!   Compute the initial approximation of A ^ (-2/3).

call mpmdc (a, t1, n)
n3 = - 2 * n / 3
t2 = (t1 * 2.d0 ** (n + 3.d0 * n3 / 2.d0)) ** (-2.d0 / 3.d0)
call mpdmc (t2, n3, b)
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.
mpnw = 3
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 2, mq - 1
  nw1 = mpnw
  mpnw = min (2 * mpnw - 2, nws) + 1
  nw2 = mpnw
100  continue
  call mpmul (b, b, s(k1))
  call mpmul (b, s(k1), s(k2))
  call mpmul (s(k0), s(k2), s(k1))
  call mpsub (f, s(k1), s(k2))
  mpnw = nw1
  call mpmul (b, s(k2), s(k1))
  call mpdivd (s(k1), 3.d0, 0, s(k2))
  mpnw = nw2
  call mpadd (b, s(k2), s(k1))
  call mpeq (s(k1), b)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
 enddo

!   Perform last iteration using Karp's trick.

call mpmul (a, b, s(k0))
nw1 = mpnw
mpnw = min (2 * mpnw - 2, nws) + 1
nw2 = mpnw
call mpmul (s(k0), s(k0), s(k1))
call mpmul (s(k0), s(k1), s(k2))
call mpsub (a, s(k2), s(k1))
mpnw = nw1
call mpmul (s(k1), b, s(k2))
call mpdivd (s(k2), 3.d0, 0, s(k1))
mpnw = nw2
call mpadd (s(k0), s(k1), s(k2))
call mpeq (s(k2), b)

!   Restore original precision level.

mpnw = nws
call mproun (b)

120  if (mpidb .ge. 7) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPCBRT O'/(6f12.0))
endif
return
end subroutine

subroutine mpcpr (a, b, ic)

!   This routine compares the MP numbers A and B and returns in IC the value
!   -1, 0, or 1 depending on whether A < B, A = B, or A > B.  It is faster
!   than merely subtracting A and B and looking at the sign of the result.
!   Debug output begins with MPIDB = 9.

dimension a(mpnw+4), b(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  ic = 0
  return
endif
if (mpidb .ge. 9) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPCPR I'/(6f12.0))
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 1) (b(i), i = 1, no)
endif
ia = sign (1., a(1))
if (a(1) .eq. 0.) ia = 0
ib = sign (1., b(1))
if (b(1) .eq. 0.) ib = 0

!   Compare signs.

if (ia .ne. ib) then
  ic = sign (1, ia - ib)
  goto 110
endif

!   The signs are the same.  Compare exponents.

ma = a(2)
mb = b(2)
if (ma .ne. mb) then
  ic = ia * sign (1, ma - mb)
  goto 110
endif

!   The signs and the exponents are the same.  Compare mantissas.

na = min (int (abs (a(1))), mpnw)
nb = min (int (abs (b(1))), mpnw)

do i = 3, min (na, nb) + 2
  if (a(i) .ne. b(i)) then
    ic = ia * sign (1., a(i) - b(i))
    goto 110
  endif
enddo

!   The mantissas are the same to the common length.  Compare lengths.

if (na .ne. nb) then
  ic = ia * sign (1, na - nb)
  goto 110
endif

!   The signs, exponents, mantissas and lengths are the same.  Thus A = B.

ic = 0

110  if (mpidb .ge. 9) write (mpldb, 2) ic
2 format ('MPCPR O',i4)
return
end subroutine

subroutine mpdotd (n, isa, a, isb, db, c)

!   This routine computes the dot product of the MP vector A with the DP
!   vector DB, returning the MP result in C.  This routine is used in the
!   author's customized PSLQ routine, resulting in substantial speedup.
!   The length of both the A and DB vectors is N, and ISA and ISB are the 
!   skip distances between successive elements of A and DB, measured in 
!   single precision words, and in DP words, respectively.  ISA must be
!   at least MPNW + 4.  The DP values in DB must be whole numbers, so for
!   example they cannot be larger than 2^53 on IEEE systems.  Debug output
!   begins with MPIDB = 8.

double precision d1(mpnw+5), d2(mpnw+5), db(isb*n), dmax, &
  dt0, dt1, dt2, dt3, dt4
!>
!   Set DMAX to the largest DP whole number that can be represented exactly:
!   2.d0 ** 53 on IEEE systems, or 2.d0 ** 48 on Cray non-IEEE systems.

parameter (dmax = 2.d0 ** 53)
real a(isa*n), c(mpnw+4), s1(mpnw+5)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 8) then
  write (6, 1) n, isa, isb
1 format ('mpdotd input:',3i8)

  do k = 1, n
    ka = (k - 1) * isa
8    kb = (k - 1) * isb + 1
    call mpmdc (a(1+ka), dt1, n1)
    dt1 = dt1 * 2.d0 ** n1
    write (6, '(i4,1p,2d25.16)') k, dt1, db(kb)
  enddo
endif

do i = 1, 4
  d1(i) = 0.d0
  d2(i) = 0.d0
enddo

!   ND is the length of D1, IXD is the exponent (as in ordinary mpfun format).
!   In the code below ND + 1 mantissa words are maintained whenever possible.

nd = 0
ixd = 0
nrel = 0

!   Loop over the n input data pairs.

do k = 1, n
  ka = (k - 1) * isa
  kb = (k - 1) * isb + 1

  do kk = 1, mpnw + 4
    s1(kk) = a(ka+kk)
  enddo

  na = min (int (abs (s1(1))), mpnw)
  dt0 = db(kb)

  if (na .eq. 0 .or. dt0 .eq. 0.d0) goto 100

!   Check to make sure the input DP value satisfies the requirements.

  if (abs (dt0) .ge. dmax .or. mod (dt0, 1.d0) .ne. 0.d0) then
    if (mpker(73) .ne. 0) then
      write (6, 2) k, dt0
2     format ('mpdotd: improper dp value:',i4,1p,d25.15)
      mpier = 73
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  endif

!   Save the two initial and the two final words of A, then zero these words.

  ia1 = s1(1)
  ia2 = s1(2)
  ia3 = s1(na+3)
  ia4 = s1(na+4)
  s1(1) = 0.
  s1(2) = 0.
  s1(na+3) = 0.
  s1(na+4) = 0.
  if (ia1 .lt. 0) dt0 = - dt0

!   Split the input DP value into high-order, mid-order and a low-order values.

  dt1 = int (mprx2 * dt0)
  dt0 = dt0 - mpbx2 * dt1
  dt2 = int (mprdx * dt0)
  dt3 = dt0 - mpbdx * dt2

  if (dt1 .eq. 0.d0 .and. dt2 .eq. 0.d0) then

!   Only the low-order part of the input DP value is nonzero.

    ish = ia2 - ixd
    if (nd .eq. 0) ish = 0

    if (ish .ge. 0) then

!   The product a(k) * db(k) has greater exponent than the cumulative sum.
!   Thus the cumulative sum must be shifted to the right by ish words.

      m1 = min (na, ish)
      m2 = min (na, nd + ish)
      m3 = na
      m4 = min (max (na, ish), mpnw + 1)
      m5 = min (max (na, nd + ish), mpnw + 1)
      d2(1) = 0.d0
      d2(2) = 0.d0

      do i = 1, m1
        d2(i+2) = dt3 * s1(i+2)
      enddo

      do i = m1 + 1, m2
        d2(i+2) = d1(i+2-ish) + dt3 * s1(i+2)
      enddo

      do i = m2 + 1, m3
        d2(i+2) = dt3 * s1(i+2)
      enddo

      do i = m3 + 1, m4
        d2(i+2) = 0.d0
      enddo

      do i = m4 + 1, m5
        d2(i+2) = d1(i+2-ish)
      enddo

!   Copy d2 back to d1.

      do i = 1, m5 + 2
        d1(i) = d2(i)
      enddo

      nd = m5
      ixd = ia2
      d1(nd+3) = 0.d0
      d1(nd+4) = 0.d0
    else

!   The product a(k) * db(k) has smaller exponent than the cumulative sum.
!   Thus the product must be shifted to the right by -ish words.

      nsh = - ish
      m1 = min (nd, nsh)
      m2 = min (nd, na + nsh)
      m3 = nd
      m4 = min (max (nd, nsh), mpnw + 1)
      m5 = min (max (nd, na + nsh), mpnw + 1)

      do i = m1 + 1, m2
        d1(i+2) = d1(i+2) + dt3 * s1(i+2-nsh)
      enddo

      do i = m3 + 1, m4
        d1(i+2) = 0.d0
      enddo
        
      do i = m4 + 1, m5
        d1(i+2) = dt3 * s1(i+2-nsh)
      enddo

      nd = m5
      d1(nd+3) = 0.d0
      d1(nd+4) = 0.d0
    endif
    nrel = nrel + 1
  elseif (dt1 .eq. 0.d0) then

!   Only the lower two parts of the input DP value are nonzero.

    ish = ia2 + 1 - ixd
    if (nd .eq. 0) ish = 0

    if (ish .ge. 0) then

!   The product a(k) * db(k) has greater exponent than the cumulative sum.
!   Thus the cumulative sum must be shifted to the right by ish words.

      m1 = min (na + 1, ish)
      m2 = min (na + 1, nd + ish)
      m3 = na + 1
      m4 = min (max (na + 1, ish), mpnw + 1)
      m5 = min (max (na + 1, nd + ish), mpnw + 1)
      d2(1) = 0.d0
      d2(2) = 0.d0

      do i = 1, m1
        d2(i+2) = dt2 * s1(i+2) + dt3 * s1(i+1)
      enddo

      do i = m1 + 1, m2
        d2(i+2) = d1(i+2-ish) + dt2 * s1(i+2) + dt3 * s1(i+1)
      enddo

      do i = m2 + 1, m3
        d2(i+2) = dt2 * s1(i+2) + dt3 * s1(i+1)
      enddo

      do i = m3 + 1, m4
        d2(i+2) = 0.d0
      enddo

      do i = m4 + 1, m5
        d2(i+2) = d1(i+2-ish)
      enddo

!   Copy d2 back to d1.

      do i = 1, m5 + 2
        d1(i) = d2(i)
      enddo

      nd = m5
      ixd = ia2 + 1
      d1(nd+3) = 0.d0
      d1(nd+4) = 0.d0
    else

!   The product a(k) * db(k) has smaller exponent than the cumulative sum.
!   Thus the product must be shifted to the right by -ish words.

      nsh = - ish
      m1 = min (nd, nsh)
      m2 = min (nd, na + 1 + nsh)
      m3 = nd
      m4 = min (max (nd, nsh), mpnw + 1)
      m5 = min (max (nd, na + 1 + nsh), mpnw + 1)

      do i = m1 + 1, m2
        d1(i+2) = d1(i+2) + dt2 * s1(i+2-nsh) + dt3 * s1(i+1-nsh)
      enddo

      do i = m3 + 1, m4
        d1(i+2) = 0.d0
      enddo
        
      do i = m4 + 1, m5
        d1(i+2) = dt2 * s1(i+2-nsh) + dt3 * s1(i+1-nsh)
      enddo

      nd = m5
      d1(nd+3) = 0.d0
      d1(nd+4) = 0.d0
    endif
    nrel = nrel + 2
  else

!   All three parts of the input DP value are nonzero.

    ish = ia2 + 2 - ixd
    if (nd .eq. 0) ish = 0

    if (ish .ge. 0) then

!   The product a(k) * db(k) has greater exponent than the cumulative sum.
!   Thus the cumulative sum must be shifted to the right by ish words.

      m1 = min (na + 2, ish)
      m2 = min (na + 2, nd + ish)
      m3 = na + 2
      m4 = min (max (na + 2, ish), mpnw + 1)
      m5 = min (max (na + 2, nd + ish), mpnw + 1)
      d2(1) = 0.d0
      d2(2) = 0.d0

      do i = 1, m1
        d2(i+2) = dt1 * s1(i+2) + dt2 * s1(i+1) + dt3 * s1(i)
      enddo

      do i = m1 + 1, m2
        d2(i+2) = d1(i+2-ish) + dt1 * s1(i+2) + dt2 * s1(i+1) &
          + dt3 * s1(i)
      enddo

      do i = m2 + 1, m3
        d2(i+2) = dt1 * s1(i+2) + dt2 * s1(i+1) + dt3 * s1(i)
      enddo

      do i = m3 + 1, m4
        d2(i+2) = 0.d0
      enddo

      do i = m4 + 1, m5
        d2(i+2) = d1(i+2-ish)
      enddo

!   Copy d2 back to d1.

      do i = 1, m5 + 2
        d1(i) = d2(i)
      enddo

      nd = m5
      ixd = ia2 + 2
      d1(nd+3) = 0.d0
      d1(nd+4) = 0.d0
    else

!   The product a(k) * db(k) has smaller exponent than the cumulative sum.
!   Thus the product must be shifted to the right by -ish words.

      nsh = - ish
      m1 = min (nd, nsh)
      m2 = min (nd, na + 2 + nsh)
      m3 = nd
      m4 = min (max (nd, nsh), mpnw + 1)
      m5 = min (max (nd, na + 2 + nsh), mpnw + 1)

      do i = m1 + 1, m2
        d1(i+2) = d1(i+2) + dt1 * s1(i+2-nsh) + dt2 * s1(i+1-nsh) &
          + dt3 * s1(i-nsh)
      enddo

      do i = m3 + 1, m4
        d1(i+2) = 0.d0
      enddo
        
      do i = m4 + 1, m5
        d1(i+2) = dt1 * s1(i+2-nsh) + dt2 * s1(i+1-nsh) &
          + dt3 * s1(i-nsh)
      enddo

      nd = m5
      d1(nd+3) = 0.d0
      d1(nd+4) = 0.d0
    endif
    nrel = nrel + 3
  endif

100 continue

  if (nd .eq. 0) goto 120

  if (nrel .ge. mpnpr - 1 .or. k .eq. n) then

!   Release carries using a vectorizable scheme.  Results may be negative, but
!   that is not a problem -- these will be fixed in the final call to mpnorm.

    nrel = 0

!dir$ ivdep
    do i = 3, nd + 2
      dt1 = d1(i)
      dt2 = int (mprdx * dt1)
      d1(i) = dt1 - mpbdx * dt2
      d1(i-1) = d1(i-1) + dt2
    enddo

!   If d1(2) is nonzero due to carry release, shift result to right.

    if (d1(2) .ne. 0.d0) then
      ish = 1
      ixd = ixd + 1
      nd = min (nd + 1, mpnw + 1)
    else
      ish = 0
    endif

    if (ish .ne. 0) then
      do i = nd + 2, 3, -1
        d1(i) = d1(i-ish)
      enddo

      d1(1) = 0.d0
      d1(2) = 0.d0
    endif
    d1(nd+3) = 0.d0
    d1(nd+4) = 0.d0
  endif

!   Check to see if there are leading zeros.

  do i = 1, nd + 1
    if (d1(i+2) .ne. 0.d0) goto 110
  enddo

!   The cumulative sum is now zero.

  nd = 0
  ixd = 0
  nrel = 0
  d1(1) = 0.d0
  d1(2) = 0.d0
  goto 120

110 kz = i - 1

  if (kz .gt. 0) then

!   Leading zeroes -- shift cumulative sum to left.

    do i = 1, nd - kz + 1
      d1(i+2) = d1(i+2+kz)
    enddo

    nd = nd - kz
    ixd = ixd - kz
    d1(nd+3) = 0.d0
    d1(nd+4) = 0.d0
  endif

120 continue

enddo

!   Call mpnorm to fix up result and store in c.

d1(1) = nd
d1(2) = ixd
call mpnorm (d1, c)
  
if (mpidb .ge. 8) then
  call mpmdc (c, dt1, n1)
  dt1 = dt1 * 2.d0 ** n1
  write (6, 4) dt1
4 format ('mpdotd output:',1p,d25.15)
  write (6, '(6f13.0)') (c(i), i = 1, min (int (abs (c(1))) + 2, mpndb))
endif

return
end subroutine

subroutine mpdiv (a, b, c)

!   This divides the MP number A by the MP number B to yield the MP quotient C.
!   For extra high levels of precision, use MPDIVX.  Debug output starts with
!   MPIDB = 8.

!   Max SP space for C: MPNW + 4 cells.

double precision d, rb, ss, t0, t1, t2
dimension a(mpnw+2), b(mpnw+2), c(mpnw+4), d(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 8) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPDIV I'/(6f12.0))
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 1) (b(i), i = 1, no)
endif

ia = sign (1., a(1))
ib = sign (1., b(1))
na = min (int (abs (a(1))), mpnw)
nb = min (int (abs (b(1))), mpnw)

!   Check if dividend is zero.

if (na .eq. 0) then
  c(1) = 0.
  c(2) = 0.
  goto 190
endif
if (nb .eq. 1 .and. b(3) .eq. 1.) then

!   Divisor is 1 or -1 -- result is A or -A.

  c(1) = sign (na, ia * ib)
  c(2) = a(2) - b(2)

  do i = 3, na + 2
    c(i) = a(i)
  enddo

  goto 190
endif

!   Check if divisor is zero.

if (nb .eq. 0) then
  if (mpker(31) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPDIV: Divisor is zero.')
    mpier = 31
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Initialize trial divisor and trial dividend.

t0 = mpbdx * b(3)
if (nb .ge. 2) t0 = t0 + b(4)
if (nb .ge. 3) t0 = t0 + mprdx * b(5)
if (nb .ge. 4) t0 = t0 + mprx2 * b(6)
rb = 1.d0 / t0
md = min (na + nb, mpnw)
d(1) = 0.d0

do i = 2, na + 1
  d(i) = a(i+1)
enddo

do i = na + 2, md + 4
  d(i) = 0.d0
enddo

!   Perform ordinary long division algorithm.  First compute only the first
!   NA words of the quotient.

do j = 2, na + 1
  t1 = mpbx2 * d(j-1) + mpbdx * d(j) + d(j+1) + mprdx * d(j+2)
  t0 = int (rb * t1)
  j3 = j - 3
  i2 = min (nb, mpnw + 2 - j3) + 2
  ij = i2 + j3

  do i = 3, i2
    i3 = i + j3
    d(i3) = d(i3) - t0 * b(i)
  enddo

!   Release carries periodically to avoid overflowing the exact integer
!   capacity of double precision floating point words in D.

  if (mod (j - 1, mpnpr) .eq. 0) then
!dir$ ivdep
    do i = j + 1, ij
      t1 = d(i)
      t2 = int (mprdx * t1)
      d(i) = t1 - mpbdx * t2
      d(i-1) = d(i-1) + t2
    enddo
  endif
  d(j) = d(j) + mpbdx * d(j-1)
  d(j-1) = t0
enddo

!   Compute additional words of the quotient, as long as the remainder
!   is nonzero.

do j = na + 2, mpnw + 3
  t1 = mpbx2 * d(j-1) + mpbdx * d(j) + d(j+1)
  if (j .le. mpnw + 2) t1 = t1 + mprdx * d(j+2)
  t0 = int (rb * t1)
  j3 = j - 3
  i2 = min (nb, mpnw + 2 - j3) + 2
  ij = i2 + j3
  ss = 0.d0

  do i = 3, i2
    i3 = i + j3
    d(i3) = d(i3) - t0 * b(i)
    ss = ss + abs (d(i3))
  enddo

  if (mod (j - 1, mpnpr) .eq. 0) then
!dir$ ivdep
    do i = j + 1, ij
      t1 = d(i)
      t2 = int (mprdx * t1)
      d(i) = t1 - mpbdx * t2
      d(i-1) = d(i-1) + t2
    enddo
  endif
  d(j) = d(j) + mpbdx * d(j-1)
  d(j-1) = t0
  if (ss .eq. 0.d0) goto 170
  if (ij .le. mpnw + 1) d(ij+3) = 0.d0
enddo

!   Set sign and exponent, and fix up result.

j = mpnw + 3

170  d(j) = 0.d0
if (d(1) .eq. 0.d0) then
  is = 1
else
  is = 2
endif
nc = min (j - 1, mpnw)
d(nc+3) = 0.d0
d(nc+4) = 0.d0

do i = j + 1, 3, -1
  d(i) = d(i-is)
enddo

d(1) = sign (nc, ia * ib)
d(2) = a(2) - b(2) + is - 2

call mpnorm (d, c)

190  if (mpidb .ge. 8) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 3) (c(i), i = 1, no)
3 format ('MPDIV O'/(6f12.0))
endif
return
end subroutine

subroutine mpdivd (a, b, n, c)

!   This routine divides the MP number A by the DPE number (B, N) to yield
!   the MP quotient C.   Debug output starts with MPIDB = 9.

!   Max SP space for C: MPNW + 4 cells.  Max DP space: MPNW + 4 cells.

double precision b, bb, br, d, dd, t1
dimension a(mpnw+2), c(mpnw+4), f(8), d(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 9) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPDIVD I'/(6f12.0))
  write (mpldb, 2) b, n
2 format ('MPDIVD I',1pd25.15,i10)
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ib = sign (1.d0, b)

!   Check if dividend is zero.

if (na .eq. 0) then
  c(1) = 0.
  c(2) = 0.
  goto 150
endif

!   Check if divisor is zero.

if (b .eq. 0.d0) then
  if (mpker(32) .ne. 0) then
    write (mpldb, 3)
3   format ('*** MPDIVD: Divisor is zero.')
    mpier = 32
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif
n1 = n / mpnbt
n2 = n - mpnbt * n1
bb = abs (b) * 2.d0 ** n2

!   Reduce BB to within 1 and MPBDX.

if (bb .ge. mpbdx) then

  do k = 1, 100
    bb = mprdx * bb
    if (bb .lt. mpbdx) then
      n1 = n1 + k
      goto 120
    endif
 enddo

elseif (bb .lt. 1.d0) then

  do k = 1, 100
    bb = mpbdx * bb
    if (bb .ge. 1.d0) then
      n1 = n1 - k
      goto 120
    endif
 enddo

endif

!   If B cannot be represented exactly in a single mantissa word, use MPDIV.

120  if (bb .ne. aint (bb)) then
  bb = sign (bb, b)
  call mpdmc (bb, n1 * mpnbt, f)
  call mpdiv (a, f, c)
  goto 150
endif

br = 1.d0 / bb
dd = a(3)

!   Perform short division (not vectorizable at present).  Continue as long as
!   the remainder remains nonzero.

do j = 2, mpnw + 3
  t1 = int (br * dd)
  d(j+1) = t1
  dd = mpbdx * (dd - t1 * bb)
  if (j .le. na) then
    dd = dd + a(j+2)
  else
    if (dd .eq. 0.d0) goto 140
  endif
enddo

!   Set sign and exponent of result.

j = mpnw + 3

140  nc = min (j - 1, mpnw)
d(1) = sign (nc, ia * ib)
d(2) = a(2) - n1
if (j .le. mpnw + 2) d(j+2) = 0.d0
if (j .le. mpnw + 1) d(j+3) = 0.d0

call mpnorm (d, c)

150  if (mpidb .ge. 9) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 4) (c(i), i = 1, no)
4 format ('MPDIVD O'/(6f12.0))
endif
return
end subroutine

subroutine mpdmc (a, n, b)

!   This routine converts the DPE number (A, N) to MP form in B.  All bits of
!   A are recovered in B.  However, note for example that if A = 0.1D0 and N
!   is 0, then B will NOT be the multiprecision equivalent of 1/10.  Debug
!   output starts with MPIDB = 9.

!   Max SP space for B:  8 cells.

double precision a, aa
dimension b(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 9) write (mpldb, 1) a, n
1 format ('MPDMC I',1pd25.15,i10)

!   Check for zero.

if (a .eq. 0.d0) then
  b(1) = 0.
  b(2) = 0.
  goto 150
endif
n1 = n / mpnbt
n2 = n - mpnbt * n1
aa = abs (a) * 2.d0 ** n2

!   Reduce AA to within 1 and MPBDX.

if (aa .ge. mpbdx) then

  do k = 1, 100
    aa = mprdx * aa
    if (aa .lt. mpbdx) then
      n1 = n1 + k
      goto 120
    endif
 enddo

elseif (aa .lt. 1.d0) then

  do k = 1, 100
    aa = mpbdx * aa
    if (aa .ge. 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo

endif

!   Store successive sections of AA into B.

120  b(2) = n1
b(3) = aint (aa)
aa = mpbdx * (aa - b(3))
b(4) = aint (aa)
aa = mpbdx * (aa - b(4))
b(5) = aint (aa)
aa = mpbdx * (aa - b(5))
b(6) = aint (aa)
b(7) = 0.
b(8) = 0.

do i = 6, 3, -1
  if (b(i) .ne. 0.) goto 140
enddo

140  aa = i - 2
b(1) = sign (aa, a)

150  if (mpidb .ge. 9) then
  no = abs (b(1)) + 2.
  write (mpldb, 2) (b(i), i = 1, no)
2 format ('MPDMC O'/(6f12.0))
endif
return
end subroutine

subroutine mpeq (a, b)

!   This routine sets the MP number B equal to the MP number A.  Debug output
!   starts with MPIDB = 10.

!   Max SP space for B: MPNW + 3 cells.

!   The fact that only MPNW + 3 cells, and not MPNW + 4 cells, are copied is
!   important in some routines that increase the precision level by one.

dimension a(mpnw+3), b(mpnw+3)

! GZ 
!write(*,*) 'mpnw',mpnw
!write(*,*) 'a',a
!write(*,*) 'b',b

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 10) write (mpldb, 1)
1 format ('MPEQ')

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
if (na .eq. 0)  then
  b(1) = 0.
  b(2) = 0.
  goto 110
endif
b(1) = sign (na, ia)

do i = 2, na + 3
  b(i) = a(i)
enddo

!write(*,*) 'a2',a
!write(*,*) 'b2',b

110 return
end subroutine

subroutine mpinfr (a, b, c)

!   Sets B to the integer part of the MP number A and sets C equal to the
!   fractional part of A.  Note that if A = -3.3, then B = -3 and C = -0.3.
!   Debug output starts with MPIDB = 9.

!   Max SP space for B and C: MPNW + 4 cells.

dimension a(mpnw+2), b(mpnw+2), c(mpnw+2)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 9) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPINFR I'/(6f12.0))
endif

!   Check if  A  is zero.

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ma = a(2)
if (na .eq. 0)  then
  b(1) = 0.
  b(2) = 0.
  c(1) = 0.
  c(2) = 0.
  goto 120
endif

if (ma .ge. mpnw - 1) then
  if (mpker(40) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPINFR: Argument is too large.')
    mpier = 40
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Place integer part in  B.

nb = min (max (ma + 1, 0), na)
if (nb .eq. 0) then
  b(1) = 0.
  b(2) = 0.
else
  b(1) = sign (nb, ia)
  b(2) = ma
  b(nb+3) = 0.
  b(nb+4) = 0.

  do i = 3, nb + 2
    b(i) = a(i)
  enddo
endif

!   Place fractional part in C.

nc = na - nb
if (nc .le. 0) then
  c(1) = 0.
  c(2) = 0.
else
  c(1) = sign (nc, ia)
  c(2) = ma - nb
  c(nc+3) = 0.
  c(nc+4) = 0.

  do i = 3, nc + 2
    c(i) = a(i+nb)
  enddo
endif

!   Fix up results.  B may have trailing zeros and C may have leading zeros.

call mproun (b)
call mproun (c)

120  if (mpidb .ge. 9)  then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPINFR O'/(6f12.0))
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 3) (c(i), i = 1, no)
endif
return
end subroutine

subroutine mpmdc (a, b, n)

!   This converts the MP number A to the DPE form (B, N), accurate to between
!   14 and 17 digits, depending on system.  B will be between 1 and MPBDX.
!   Debug output starts with MPIDB = 9.

double precision aa, b
dimension a(mpnw+2)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b = 0.d0
  n = 0
  return
endif
if (mpidb .ge. 9) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPMDC I'/(6f12.0))
endif

if (a(1) .eq. 0.)  then
  b = 0.d0
  n = 0
  goto 100
endif

na = abs (a(1))
aa = a(3)
if (na .ge. 2) aa = aa + mprdx * a(4)
if (na .ge. 3) aa = aa + mprx2 * a(5)
if (na .ge. 4) aa = aa + mprdx * mprx2 * a(6)

n = mpnbt * a(2)
b = sign (aa, dble (a(1)))

100  if (mpidb .ge. 9) write (mpldb, 2) b, n
2 format ('MPMDC O',f10.0,i10)
return
end subroutine

subroutine mpmul (a, b, c)

!   This routine multiplies MP numbers A and B to yield the MP product C.
!   When one of the arguments has a much higher level of precision than the
!   other, this routine is slightly more efficient if A has the lower level of
!   precision.  For extra high levels of precision, use MPMULX.  Debug output
!   starts with MPIDB = 8.

!   Max SP space for C: MPNW + 4 cells.

!   This routine returns up to MPNW mantissa words of the product.  If the
!   complete double-long product of A and B is desired (for example in large
!   integer applications), then MPNW must be at least as large as the sum of
!   the mantissa lengths of A and B.  In other words, if the precision levels
!   of A and B are both 64 words, then MPNW must be at least 128 words to
!   obtain the complete double-long product in C.

double precision d, t1, t2, t3
dimension a(mpnw+2), b(mpnw+2), c(mpnw+4), d(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 8) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPMUL I'/(6f12.0))
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 1) (b(i), i = 1, no)
endif

ia = sign (1., a(1))
ib = sign (1., b(1))
na = min (int (abs (a(1))), mpnw)
nb = min (int (abs (b(1))), mpnw)
if (na .eq. 0 .or. nb .eq. 0) then

!   One of the inputs is zero -- result is zero.

  c(1) = 0.
  c(2) = 0.
  goto 170
endif
if (na .eq. 1 .and. a(3) .eq. 1.) then

!   A is 1 or -1 -- result is B or -B.

  c(1) = sign (nb, ia * ib)
  c(2) = a(2) + b(2)

  do i = 3, nb + 2
    c(i) = b(i)
  enddo

  goto 170
elseif (nb .eq. 1 .and. b(3) .eq. 1.) then

!   B is 1 or -1 -- result is A or -A.

  c(1) = sign (na, ia * ib)
  c(2) = a(2) + b(2)

  do i = 3, na + 2
    c(i) = a(i)
  enddo

  goto 170
endif

nc = min (na + nb, mpnw)
d2 = a(2) + b(2)

do i = 1, nc + 4
  d(i) = 0.d0
enddo

!   Perform ordinary long multiplication algorithm.  Accumulate at most MPNW+4
!   mantissa words of the product.

do j = 3, na + 2
  t1 = a(j)
  j3 = j - 3
  n2 = min (nb + 2, mpnw + 4 - j3)

  do i = 3, n2
    d(i+j3) = d(i+j3) + t1 * b(i)
  enddo

!   Release carries periodically to avoid overflowing the exact integer
!   capacity of double precision floating point words in D.

  if (mod (j - 2, mpnpr) .eq. 0) then
    i1 = max (3, j - mpnpr)
    i2 = n2 + j3

!dir$ ivdep
    do i = i1, i2
      t1 = d(i)
      t2 = int (mprdx * t1)
      d(i) = t1 - mpbdx * t2
      d(i-1) = d(i-1) + t2
    enddo
  endif
enddo

!   If D(2) is nonzero, shift the result one cell right.

if (d(2) .ne. 0.d0) then
  d2 = d2 + 1.

!dir$ ivdep
  do i = nc + 4, 3, -1
    d(i) = d(i-1)
  enddo
endif
d(1) = sign (nc, ia * ib)
d(2) = d2

!   Fix up result, since some words may be negative or exceed MPBDX.

call mpnorm (d, c)

170  if (mpidb .ge. 8) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 2) (c(i), i = 1, no)
2 format ('MPMUL O'/(6f12.0))
endif
return
end subroutine

subroutine mpmuld (a, b, n, c)

!   This routine multiplies the MP number A by the DPE number (B, N) to yield
!   the MP product C.  Debug output starts with MPIDB = 9.

!   Max SP space for C: MPNW + 4 cells.  Max DP space: MPNW + 4 cells.

double precision b, bb, d
dimension a(mpnw+2), c(mpnw+4), f(8), d(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 9) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPMULD I'/(6f12.0))
  write (mpldb, 2) b, n
2 format ('MPMULD I',1pd25.15,i10)
endif

!   Check for zero inputs.

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ib = sign (1.d0, b)
if (na .eq. 0 .or. b .eq. 0.d0) then
  c(1) = 0.
  c(2) = 0.
  goto 140
endif
n1 = n / mpnbt
n2 = n - mpnbt * n1
bb = abs (b) * 2.d0 ** n2

!   Reduce BB to within 1 and MPBDX.

if (bb .ge. mpbdx) then

  do k = 1, 100
    bb = mprdx * bb
    if (bb .lt. mpbdx) then
      n1 = n1 + k
      goto 120
    endif
  enddo
elseif (bb .lt. 1.d0) then
  do k = 1, 100
    bb = mpbdx * bb
    if (bb .ge. 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo
endif

!   If B cannot be represented exactly in a single mantissa word, use MPMUL.

120  if (bb .ne. aint (bb)) then
  bb = sign (bb, b)
  call mpdmc (bb, n1 * mpnbt, f)
  call mpmul (f, a, c)
  goto 140
endif

!   Perform short multiply operation.

!dir$ ivdep
do i = 3, na + 2
  d(i) = bb * a(i)
enddo

!   Set the exponent and fix up the result.

d(1) = sign (na, ia * ib)
d(2) = a(2) + n1
d(na+3) = 0.d0
d(na+4) = 0.d0

call mpnorm (d, c)

140  if (mpidb .ge. 9) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 3) (c(i), i = 1, no)
3 format ('MPMULD O'/(6f12.0))
endif
return
end subroutine

subroutine mpnint (a, b)

!   This sets B equal to the integer nearest to the MP number A.  Debug output
!   starts with MPIDB = 8.

!   Max SP space for B: MPNW + 4 cells.

dimension a(mpnw+2), b(mpnw+2), f(8), s(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 8) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPNINT I'/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ma = a(2)
if (na .eq. 0)  then

!   A is zero -- result is zero.

  b(1) = 0.
  b(2) = 0.
  goto 110
endif
if (ma .ge. mpnw) then

!   A cannot be represented exactly as an integer.

  if (mpker(56) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPNINT: Argument is too large.')
    mpier = 56
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

k0 = 1
f(1) = 1.
f(2) = -1.
f(3) = 0.5d0 * mpbdx
f(4) = 0.

!   Add or subtract 1/2 from the input, depending on its sign.

if (ia .eq. 1) then
  call mpadd (a, f, s(k0))
else
  call mpsub (a, f, s(k0))
endif
ic = sign (1., s(k0))
nc = abs (s(k0))
mc = s(k0+1)

!   Place integer part of S in B.

nb = min (max (mc + 1, 0), nc)
if (nb .eq. 0) then
  b(1) = 0.
  b(2) = 0.
else
  b(1) = sign (nb, ic)
  b(2) = mc
  b(nb+3) = 0.
  b(nb+4) = 0.

  do i = 3, nb + 2
    b(i) = s(i+k0-1)
  enddo
endif

110  if (mpidb .ge. 8) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPNINT O'/(6f12.0))
endif
return
end subroutine

subroutine mpnorm (d, a)

!   This converts the MP number in array D of MPCOM4 to the standard
!   normalized form in A.  The MP routines often leave negative numbers or
!   values exceeding the radix MPBDX in result arrays, and this fixes them.
!   MPNORM assumes that two extra mantissa words are input at the end of D.
!   This reduces precision loss when it is necessary to shift the result to
!   the left.  This routine is not intended to be called directly by the user.
!   Debug output starts with MPIDB = 10.

!   Max SP space for A: MPNW + 4 cells.

double precision d, s1, t1, t2, t3
dimension d(mpnw+4), a(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  a(1) = 0.
  a(2) = 0.
  return
endif
if (mpidb .ge. 9) then
  no = min (int (abs (d(1))), mpndb) + 4
  write (mpldb, 1) (d(i), i = 1, no)
1 format ('MPNORM I'/(4f18.0))
endif

ia = sign (1.d0, d(1))
na = min (int (abs (d(1))), mpnw)
if (na .eq. 0)  then
  a(1) = 0.
  a(2) = 0.
  goto 170
endif
n4 = na + 4
a2 = d(2)
d(2) = 0.d0

110 continue
!>
!   Try a vectorized fixup loop three times, unless A is very short.  This
!   should handle 99% of the inputs.  On RISC computers, it is more
!   efficient to completely bypass this loop, by uncommenting the next line.

goto 120
if (na .le. 8) goto 120

do k = 1, 3
  s1 = 0.d0

!dir$ ivdep
  do i = 3, n4
    t2 = mprdx * d(i)
    t1 = int (t2)
    if (t2 .lt. 0.d0 .and. t1 .ne. t2) t1 = t1 - 1.d0
    d(i) = d(i) - t1 * mpbdx
    d(i-1) = d(i-1) + t1
    s1 = s1 + abs (t1)
  enddo

  if (s1 .eq. 0.d0) goto 140
enddo

!   Still not fixed - use recursive loop.  This loop is not vectorizable,
!   but it is guaranteed to complete the job in one pass.

120  t1 = 0.d0

do i = n4, 3, -1
  t3 = t1 + d(i)
  t2 = mprdx * (t3)
  t1 = int (t2)
  if (t2 .lt. 0.d0 .and. t1 .ne. t2) t1 = t1 - 1.d0
  d(i) = t3 - t1 * mpbdx
enddo

d(2) = d(2) + t1

140  continue

if (d(2) .lt. 0.d0) then

!   D(2) is negative -- negate all words and re-normalize.

  ia = - ia
  d(3) = d(3) + mpbdx * d(2)
  d(2) = 0.d0

  do i = 2, n4
    d(i) = - d(i)
  enddo

  goto 110
elseif (d(2) .gt. 0.d0) then

!   The fixup loops above "spilled" a nonzero number into D(2).  Shift the
!   entire number right one cell.  The exponent and length of the result
!   are increased by one.

  do i = n4, 3, -1
    a(i) = d(i-1)
  enddo

  na = min (na + 1, mpnw)
  a2 = a2 + 1.
else
  do i = 3, n4
    a(i) = d(i)
  enddo
endif

!   Perform rounding and truncation.

a(1) = sign (na, ia)
a(2) = a2

call mproun (a)

170  if (mpidb .ge. 9) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPNORM O'/(6f12.0))
endif
return
end subroutine

subroutine mpnpwr (a, n, b)

!   This computes the N-th power of the MP number A and returns the MP result
!   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
!   of A ^ |N| is returned.  For extra high levels of precision, use MPNPWX.
!   Debug output starts with MPIDB = 7.

!   Max SP space for B: MPNW + 4 cells.

!   This routine employs the binary method for exponentiation.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0)
dimension a(mpnw+2), b(mpnw+4), f1(8), s(2*mpnw+10)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 7) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) n, (a(i), i = 1, no)
1 format ('MPNPWR I',i5/(6f12.0))
endif

na = min (int (abs (a(1))), mpnw)
if (na .eq. 0) then
  if (n .ge. 0) then
    b(1) = 0.
    b(2) = 0.
    goto 120
  else
    if (mpker(57) .ne. 0) then
      write (mpldb, 2)
2     format ('*** MPNPWR: Argument is zero and N is negative or zero.')
      mpier = 57
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  endif
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
nws = mpnw
mpnw = mpnw + 1
nn = abs (n)
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
if (nn .eq. 0) then
  call mpeq (f1, b)
  mpnw = nws
  goto 120
elseif (nn .eq. 1) then
  call mpeq (a, b)
  goto 110
elseif (nn .eq. 2) then
  call mpmul (a, a, s(k0))
  call mpeq (s(k0), b)
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN .GT. NN.

t1 = nn
mn = cl2 * log (t1) + 1.d0 + mprxx
call mpeq (f1, b)
call mpeq (a, s(k0))
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn .ne. 2 * kk) then
    call mpmul (b, s(k0), s(k1))
    call mpeq (s(k1), b)
  endif
  kn = kk
  if (j .lt. mn) then
    call mpmul (s(k0), s(k0), s(k1))
    call mpeq (s(k1), s(k0))
  endif
enddo

!   Compute reciprocal if N is negative.

110  if (n .lt. 0) then
  call mpdiv (f1, b, s(k0))
  call mpeq (s(k0), b)
endif

!   Restore original precision level.

mpnw = nws
call mproun (b)

120  if (mpidb .ge. 7) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPNPWR O'/(6f12.0))
endif
return
end subroutine

subroutine mpnrt (a, n, b)

!   This computes the N-th root of the MP number A and returns the MP result
!   in B.  N must be at least one and must not exceed 2 ^ 30.  For extra high
!   levels of precision, use MPNRTX.  Debug output starts with MPIDB = 7.

!   Max SP space for B: MPNW + 4 cells.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to A ^ (-1/N):

!    X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)

!   The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.
!   These iterations are performed with a maximum precision level MPNW that
!   is dynamically changed, approximately doubling with each iteration.
!   See the comment about the parameter NIT in MPDIVX.

!   When N is large and A is very near one, the following binomial series is
!   employed instead of the Newton scheme:

!   (1 + x)^(1/N)  =  1  +  x / N  +  x^2 * (1 - N) / (2! N^2)  +  ...

!   See the comment about the parameter NIT in MPDIVX.

double precision alt, cl2, t1, t2, tn
parameter (alt = 0.693147180559945309d0, cl2 = 1.4426950408889633d0, &
  nit = 3, n30 = 2 ** 30)
dimension a(mpnw+2), b(mpnw+4), f1(8), f2(8), s(4*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 7) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) n, (a(i), i = 1, no)
1 format ('MPNRT I',i5/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)

if (na .eq. 0) then
  b(1) = 0.
  b(2) = 0.
  goto 140
endif
if (ia .lt. 0) then
  if (mpker(59) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPNRT: Argument is negative.')
    mpier = 59
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif
if (n .le. 0 .or. n .gt. n30) then
  if (mpker(60) .ne. 0) then
    write (mpldb, 3) n
3   format ('*** MPNRT: Improper value of N',i10)
    mpier = 60
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   If N = 1, 2 or 3, call MPEQ, MPSQRT or MPCBRT.  These are faster.

if (n .eq. 1) then
  call mpeq (a, b)
  goto 140
elseif (n .eq. 2) then
  call mpsqrt (a, b)
  goto 140
elseif (n .eq. 3) then
  call mpcbrt (a, b)
  goto 140
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
k3 = k2 + n5
nws = mpnw
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Check how close A is to 1.

call mpsub (a, f1, s(k0))
if (s(k0) .eq. 0.) then
  call mpeq (f1, b)
  goto 140
endif
call mpmdc (s(k0), t1, n1)
n2 = cl2 * log (abs (t1))
t1 = t1 * 0.5d0 ** n2
n1 = n1 + n2
if (n1 .le. -30) then
  t2 = n
  n2 = cl2 * log (t2) + 1.d0 + mprxx
  n3 = - mpnbt * mpnw / n1
  if (n3 .lt. 1.25 * n2) then

!   A is so close to 1 that it is cheaper to use the binomial series.

    mpnw = mpnw + 1
    call mpdivd (s(k0), t2, 0, s(k1))
    call mpadd (f1, s(k1), s(k2))
    k = 0

100 k = k + 1
    t1 = 1 - k * n
    t2 = (k + 1) * n
    call mpmuld (s(k1), t1, 0, s(k3))
    call mpdivd (s(k3), t2, 0, s(k1))
    call mpmul (s(k0), s(k1), s(k3))
    call mpeq (s(k3), s(k1))
    call mpadd (s(k1), s(k2), s(k3))
    call mpeq (s(k3), s(k2))
    if (s(k1) .ne. 0. .and. s(k1+1) .ge. - mpnw) goto 100

    call mpeq (s(k2), b)
    call mpdiv (f1, s(k2), s(k0))
    goto 130
  endif
endif

!   Compute the initial approximation of A ^ (-1/N).

tn = n
call mpmdc (a, t1, n1)
n2 = - n1 / tn
t2 = exp (-1.d0 / tn * (log (t1) + (n1 + tn * n2) * alt))
call mpdmc (t2, n2, b)
call mpdmc (tn, 0, f2)
mpnw = 3
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 2, mq
  mpnw = min (2 * mpnw - 2, nws) + 1
110  continue
  call mpnpwr (b, n, s(k0))
  call mpmul (a, s(k0), s(k1))
  call mpsub (f1, s(k1), s(k0))
  call mpmul (b, s(k0), s(k1))
  call mpdivd (s(k1), tn, 0, s(k0))
  call mpadd (b, s(k0), s(k1))
  call mpeq (s(k1), b)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 110
  endif
enddo

!   Take the reciprocal to give final result.

call mpdiv (f1, b, s(k1))
call mpeq (s(k1), b)

!   Restore original precision level.

130  mpnw = nws
call mproun (b)

140  if (mpidb .ge. 7) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 4) (b(i), i = 1, no)
4 format ('MPNRT O'/(6f12.0))
endif
return
end subroutine

subroutine mprand (a)

!   This returns a pseudo-random MP number A between 0 and 1.  Debug output 
!   starts with MPIDB = 9.

!   Max SP space for A: MPNW + 4 cells.

double precision f7, r30, s0, sd, t1, t2, t30
parameter (f7 = 78125.d0, s0 = 314159265.d0)
dimension a(mpnw+4)
save r30, t30, sd
data sd/s0/, r30/0.d0/

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  a(1) = 0.
  a(2) = 0.
  return
endif
if (r30 .eq. 0.d0) then
  r30 = 1.d0
  t30 = 1.d0

  do i = 1, 30
    r30 = 0.5d0 * r30
    t30 = 2.d0 * t30
  enddo
endif

a(1) = mpnw
a(2) = -1.

do i = 3, mpnw + 4
  t1 = f7 * sd
  t2 = aint (r30 * t1)
  sd = t1 - t30 * t2
  a(i) = aint (mpbdx * r30 * sd)
enddo

call mproun (a)

if (mpidb .ge. 9) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPRAND O'/(6f12.0))
endif
return
end subroutine

subroutine mproun (a)

!   This performs rounding and truncation of the MP number A.  It is called
!   by MPNORM, and also by other subroutines when the precision level is
!   reduced by one.  It is not intended to be directly called by the user.

!   Maximum SP space for A:  MPNW + 4 cells.

!   The parameter AMX is the absolute value of the largest exponent word
!   allowed for MP numbers.

parameter (amx = 2.e6)
dimension a(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  a(1) = 0.
  a(2) = 0.
  return
endif

!   Check for initial zeroes.

a2 = a(2)
a(2) = 0.
ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
n4 = na + 4
if (a(3) .eq. 0.) then

!   Find the first nonzero word and shift the entire number left.  The length
!   of the result is reduced by the length of the shift.

  do i = 4, n4
    if (a(i) .ne. 0.) goto 110
  enddo

  a(1) = 0.
  a(2) = 0.
  goto 170

110  k = i - 3

!dir$ ivdep
  do i = 3, n4 - k
    a(i) = a(i+k)
  enddo

  a2 = a2 - k
  na = na - max (k - 2, 0)
  if (k .eq. 2) a(na+3) = 0.
endif

!   Perform rounding depending on MPIRD.

if (na .eq. mpnw .and. mpird .ge. 1) then
  if (mpird .eq. 1 .and. a(na+3) .ge. 0.5d0 * mpbdx .or. mpird .eq. 2 &
    .and. a(na+3) .ge. 1.) a(na+2) = a(na+2) + 1.

!   Release carries as far as necessary due to rounding.

  do i = na + 2, 3, -1
    if (a(i) .lt. mpbdx) goto 140
    a(i) = a(i) - mpbdx
    a(i-1) = a(i-1) + 1.
  enddo

!   Release of carries due to rounding continued all the way to the start --
!   i.e. number was entirely 9's.

  a(3) = a(2)
  na = 1
  a2 = a2 + 1.
endif

140  if (a(na+2) .eq. 0.) then

!   At least the last mantissa word is zero.  Find the last nonzero word
!   and adjust the length of the result accordingly.

  do i = na + 2, 3, -1
    if (a(i) .ne. 0.) goto 160
  enddo

  a(1) = 0.
  a(2) = 0.
  goto 170

160  na = i - 2
endif

!   Check for overflow and underflow.

if (a2 .lt. - amx) then
  if (mpker(68) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPROUN: Exponent underflow.')
    mpier = 68
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
elseif (a2 .gt. amx) then
  if (mpker(69) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPROUN: Exponent overflow.')
    mpier = 69
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
endif

!   Check for zero.

if (a(3) .eq. 0.) then
  a(1) = 0.
  a(2) = 0.
else
  a(1) = sign (na, ia)
  a(2) = a2
  a(na+3) = 0.
  a(na+4) = 0.
endif

170  return
end subroutine

subroutine mpsqrt (a, b)

!   This computes the square root of the MP number A and returns the MP result
!   in B.  For extra high levels of precision, use MPSQRX.  Debug output
!   starts with MPIDB = 7.

!   Max SP space for B: MPNW + 4 cells.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to 1 / Sqrt(A):

!    X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k

!   where the muliplication () * X_k is performed with only half of the
!   normal level of precision.  These iterations are performed with a
!   maximum precision level MPNW that is dynamically changed, doubling with
!   each iteration.  The final iteration is performed as follows (this is
!   due to A. Karp):

!    Sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n  (approx.)

!   where the multiplications A * X_n and [] * X_n are performed with only
!   half of the final level of precision.  See the comment about the parameter
!   NIT is MPDIVX.

double precision cl2, t1, t2
parameter (cl2 = 1.4426950408889633d0, nit = 3)
dimension a(mpnw+2), b(mpnw+4), f(8), s(3*mpnw+15)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 7) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPSQRT I'/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)

if (na .eq. 0) then
  b(1) = 0.
  b(2) = 0.
  goto 120
endif
if (ia .lt. 0.d0) then
  if (mpker(70) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPSQRT: Argument is negative.')
    mpier = 70
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
nws = mpnw

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Compute the initial approximation of 1 / Sqrt(A).

call mpmdc (a, t1, n)
n2 = - n / 2
t2 = sqrt (t1 * 2.d0 ** (n + 2 * n2))
t1 = 1.d0 / t2
call mpdmc (t1, n2, b)
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.
mpnw = 3
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 2, mq - 1
  nw1 = mpnw
  mpnw = min (2 * mpnw - 2, nws) + 1
  nw2 = mpnw
100  continue
  call mpmul (b, b, s(k0))
  call mpmul (a, s(k0), s(k1))
  call mpsub (f, s(k1), s(k0))
  mpnw = nw1
  call mpmul (b, s(k0), s(k1))
  call mpmuld (s(k1), 0.5d0, 0, s(k0))
  mpnw = nw2
  call mpadd (b, s(k0), s(k1))
  call mpeq (s(k1), b)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
enddo

!   Perform last iteration using Karp's trick.

call mpmul (a, b, s(k0))
nw1 = mpnw
mpnw = min (2 * mpnw - 2, nws) + 1
nw2 = mpnw
call mpmul (s(k0), s(k0), s(k1))
call mpsub (a, s(k1), s(k2))
mpnw = nw1
call mpmul (s(k2), b, s(k1))
call mpmuld (s(k1), 0.5d0, 0, s(k2))
mpnw = nw2
call mpadd (s(k0), s(k2), s(k1))
call mpeq (s(k1), b)

!   Restore original precision level.

mpnw = nws
call mproun (b)

120  if (mpidb .ge. 7) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPSQRT O'/(6f12.0))
endif

return
end subroutine

subroutine mpsub (a, b, c)

!   This routine subtracts MP numbers A and B to yield the MP difference C,
!   by negating B and adding.  Debug output starts with MPIDB = 9.

!   Max SP space for C: MPNW + 4 cells.

dimension a(mpnw+2), b(mpnw+2), c(mpnw+4)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 9) write (mpldb, 1)
1  format ('MPSUB')

!   Check if A = B.  This is necessary because A and B might be same array,
!   in which case negating B below won't work.

if (a(1) .ne. b(1)) goto 110

do i = 2, int (abs (a(1))) + 2
  if (a(i) .ne. b(i)) goto 110
enddo

!   A = B.  Result is zero.

c(1) = 0.
c(2) = 0.
if (mpidb .ge. 9) write (mpldb, 2) (c(i), i = 1, 2)
2  format ('MPSUB O'/2f9.0)
goto 120

!   Save the sign of B, and then negate B.

110  b1 = b(1)
b(1) = - b1

!   Perform addition and restore the sign of B.

call mpadd (a, b, c)
b(1) = b1

120  return
end subroutine

end module

module mpfund

!   This module defines some I/O and debug routines.

use mpfuna
use mpfunb
use mpfunc
contains

subroutine mpdeb (cs, a)

!   This outputs the character string CS, the exponent of the MP number A, and
!   the first 50 digits of A, all on one line.  CS must either be a literal
!   string not exceeding 12 characters in length or a variable of type
!   CHARACTER*n, where n does not exceed 12.

character*(*) cs
character*1 b(160)
dimension a(mpnw+2)

if (mpier .ne. 0) return
ids = mpidb
mpidb = 0
nws = mpnw
mpnw = min (mpnw, 10)
call mpoutc (a, b, n)
n = min (n, 70)
write (mpldb, 1) cs, ' ', (b(k), k = 1, 4), (b(k), k = 9, n)
1 format (a12,67a1:/(79a1))
mpidb = ids
mpnw = nws
return
end subroutine

subroutine mpinp (iu, a, cs)

!   This routine reads the MP number A from logical unit IU.  CS is a scratch
!   array of type CHARACTER*1.  CS must be dimensioned at least 7.225*MPNW
!   + 100.   The digits of A may span more than one line.  A comma at the end
!   of the last line denotes the end of the MP number.  The input lines may not
!   exceed 120 characters in length.  Embedded blanks are allowed anywhere.
!   However, if the input number contains more than 80 embedded blanks, then
!   the dimension of CS must be increased by a corresponding amount.  The
!   exponent is optional in the input number, but if present it must appear
!   first.  Two examples:

!   1073741824.,
!   10 ^  -4 x  3.14159 26535 89793 23846 26433 83279
!     50288 41971 69399 37510,

!   Max SP space for A: MPNW + 4 cells.

character*120 lin
character*1 cs(*)
dimension a(mpnw+2)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  a(1) = 0.
  a(2) = 0.
  return
endif
l = 0
nd = 7.225d0 * mpnw + 100.d0

100  read (iu, '(A)', end = 150) lin

do i = 120, 1, -1
  if (lin(i:i) .ne. ' ') goto 120
enddo

goto 100

120  k = i
if (l .gt. nd) goto 140

do i = 1, k
  l = l + 1
  if (l .gt. nd) goto 140
  cs(l)= lin(i:i)
enddo

140  if (lin(k:k) .ne. ',') goto 100
l = l - 1

call mpinpc (cs, l, a)
goto 160

150  if (mpker(72) .ne. 0) then
  write (mpldb, 1)
1 format ('*** MPINP: End-of-file encountered.')
  mpier = 72
  if (mpker(mpier) .eq. 2) call mpabrt
endif

160  return
end subroutine

subroutine mpinpc (a, n, b)

!   Converts the CHARACTER*1 array A of length N into the MP number B.  The
!   string A must be in the format '10^s a x tb.c' where a, b and c are digit
!   strings; s and t are '-', '+' or blank; x is either 'x' or '*'.  Blanks may
!   be embedded anywhere.  The digit string a is limited to nine digits and
!   80 total characters, including blanks.  The exponent portion (i.e. the
!   portion up to and including x) and the period may optionally be omitted.
!   Debug output starts with MPIDB = 7.

!   Max SP space for B: MPNW + 4 cells.

!   The following example shows how this routine may be used to input a MP
!   number:

!   CHARACTER*1 CX(800)
!   READ (1, '(80A1)') (CX(I), I = 1, ND)
!   CALL MPINPC (CX, ND, B)

double precision bi
character*1 a, ai
character*10 dig
character*80 ca
parameter (dig = '0123456789')
dimension a(n), b(mpnw+4), f(8), s(3*mpnw+15)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 7) then
  no = min (n, int (7.225 * mpndb) + 20)
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPINPC I'/(78a1))
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
nws = mpnw
mpnw = mpnw + 1
i1 = 1
nn = 0

!   Find the carat, period, plus or minus sign, whichever comes first.

do i = 1, n
  ai = a(i)
  if (ai .eq. '^') goto 110
  if (ai .eq. '.' .or. ai .eq. '+' .or. ai .eq. '-') goto 160
enddo

goto 160

!   Make sure number preceding the carat is 10.

110  i2 = i - 1
if (i2 .gt. 80) goto 210
ca = ' '

do i = 1, i2
  ai = a(i)
  if (ai .eq. ' ') then
    goto 120
  elseif (index (dig, ai) .eq. 0) then
    goto 210
  endif
  ca(i:i) = ai
120  continue
enddo

read (ca, '(BN,I80)') nn
if (nn .ne. 10) goto 210
i1 = i2 + 2

!   Find the x or *.

do i = i1, n
  ai = a(i)
  if (ai .eq. 'x' .or. ai .eq. '*') goto 140
enddo

goto 210

!   Convert the exponent.

140  i2 = i - 1
l1 = i2 - i1 + 1
if (l1 .gt. 80) goto 210
ca = ' '
id = 0
is = 1

do i = 1, l1
  ai = a(i+i1-1)
  if (ai .eq. ' ' .or. ai .eq. '+') then
    goto 150
  elseif (ai .eq. '-' .and. id .eq. 0) then
    id = 1
    is = -1
    ca(i:i) = ' '
  else
    if (index (dig, ai) .eq. 0) goto 210
    id = 1
    ca(i:i) = ai
  endif
150  continue
enddo

read (ca, '(BN,I80)') nn
nn = is * nn
i1 = i2 + 2

!   Find the next nonblank character.

160  do i = i1, n
  if (a(i) .ne. ' ') goto 180
enddo

goto 210

!   Check if the nonblank character is a plus or minus sign.

180  i1 = i
if (a(i1) .eq. '+') then
  i1 = i1 + 1
  is = 1
elseif (a(i1) .eq. '-') then
  i1 = i1 + 1
  is = -1
else
  is = 1
endif
nb = 0
ib = 0
id = 0
ip = 0
s(k2) = 0.
s(k2+1) = 0.
f(1) = 1.
f(2) = 0.
it = 0

190   ip = 0
ca(1:6) = '000000'

!   Scan for digits, looking for the period also.  On the first pass we just
!   count, so that on the second pass it will come out right.

do i = i1, n
  ai = a(i)
  if (ai .eq. ' ') then
  elseif (ai .eq. '.') then
    if (ip .ne. 0) goto 210
    ip = id
  elseif (index (dig, ai) .eq. 0) then
    goto 210
  else
    ib = ib + 1
    id = id + 1
    ca(ib:ib) = ai
  endif
  if (ib .eq. 6 .or. i .eq. n .and. ib .ne. 0) then
    if (it .ne. 0) then
      nb = nb + 1
      read (ca(1:6), '(F6.0)') bi
      call mpmuld (s(k2), 1.d6, 0, s(k0))
      if (bi .ne. 0) then
        f(1) = 1.
        f(3) = bi
      else
        f(1) = 0.
      endif
      call mpadd (s(k0), f, s(k2))
      ca(1:6) = '000000'
    endif
    if (i .ne. n) ib = 0
  endif
enddo

if (it .eq. 0) then
  ib = 6 - ib
  if (ib .eq. 6) ib = 0
  it = 1
  goto 190
endif
if (is .eq. -1) s(k2) = - s(k2)
if (ip .eq. 0) ip = id
nn = nn + ip - id
f(1) = 1.
f(3) = 10.
call mpnpwr (f, nn, s(k0))
call mpmul (s(k2), s(k0), s(k1))
call mpeq (s(k1), b)
mpnw = nws
call mproun (b)

if (mpidb .ge. 7) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 2) (b(i), i = 1, no)
2 format ('MPINPC O'/(6f12.0))
endif
goto 220

210  if (mpker(41) .ne. 0) then
  write (mpldb, 3)
3 format ('*** MPINPC: Syntax error in literal string.')
  mpier = 41
  if (mpker(mpier) .eq. 2) call mpabrt
  mpnw = nws
endif

220  return
end subroutine

subroutine mpmout (n1, n2, la, a)

!  This outputs the MP matrix a as a DP matrix.

integer  na(n2)
double precision da(n2), d1
real a(la,n1,n2)

do i = 1, n1
  write (mpldb, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpmdc (a(1,i,j), d1, n1)
    call dpdec (d1, n1, da(j), na(j))
  enddo

  write (mpldb, 2) (da(j), na(j), j = 1, n2)
2 format (4(f13.8,'D',i5))
enddo

return
end subroutine

subroutine mpout (iu, a, la, cs)

!   This routine writes the exponent plus LA mantissa digits of the MP number
!   A to logical unit IU.  CS is a scratch array of type CHARACTER*1.  CS must
!   be dimensioned at least LA + 25.  The digits of A may span more than one
!   line.  A comma is placed at the end of the last line to denote the end of
!   the MP number.  Here is an example of the output:

!   10 ^  -4 x  3.14159265358979323846264338327950288419716939937510,

character*1 cs(la+25)
dimension a(mpnw+2)

if (mpier .ne. 0) return

nws = mpnw
ll = la / log10 (mpbdx) + 2.d0
mpnw = min (mpnw, ll)
call mpoutc (a, cs, l)
mpnw = nws
l = min (l, la + 20) + 1
cs(l) = ','
write (iu, '(78A1)') (cs(i), i = 1, l)

return
end subroutine

subroutine mpoutc (a, b, n)

!   Converts the MP number A into character form in the CHARACTER*1 array B.
!   N (an output parameter) is the length of the output.  In other words, B is
!   contained in B(1), ..., B(N).  The format is analogous to the Fortran
!   exponential format (E format), except that the exponent is placed first.
!   Debug output starts with MPIDB = 7.

!   Max CHARACTER*1 space for B: 7.225 * MPNW + 30 cells.

!   This routine is called by MPOUT, but it may be directly called by the user
!   if desired for custom output.  Example:

!   CHARACTER*1 CX(800)
!   CALL MPOUTC (A, CX, ND)
!   WRITE (1, '(20A1/(72A1))') (CX(I), I = 1, ND)

double precision aa, al2, t1
character*1 b
character*16 ca
parameter (al2 = 0.301029995663981195d0, con = 0.8304820235d0)
dimension a(mpnw+2), b(*), f(8), s(2*mpnw+10)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = ' '
  n = 0
  return
endif
if (mpidb .ge. 7) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPOUTC I'/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
nws = mpnw
mpnw = mpnw + 1
f(1) = 1.
f(2) = 0.
f(3) = 10.

!   Determine exact power of ten for exponent.

if (na .ne. 0) then
  aa = a(3)
  if (na .ge. 2) aa = aa + mprdx * a(4)
  if (na .ge. 3) aa = aa + mprx2 * a(5)
  if (na .ge. 4) aa = aa + mprdx * mprx2 * a(6)
  t1 = al2 * mpnbt * a(2) + log10 (aa)
  if (t1 .ge. 0.d0) then
    nx = t1
  else
    nx = t1 - 1.d0
  endif
  call mpnpwr (f, nx, s(k0))
  call mpdiv (a, s(k0), s(k1))

!   If we didn't quite get it exactly right, multiply or divide by 10 to fix.

100 if (s(k1+1) .lt. 0.) then
    nx = nx - 1
    call mpmuld (s(k1), 10.d0, 0, s(k0))
    call mpeq (s(k0), s(k1))
    goto 100
  elseif (s(k1+2) .ge. 10.) then
    nx = nx + 1
    call mpdivd (s(k1), 10.d0, 0, s(k0))
    call mpeq (s(k0), s(k1))
    goto 100
  endif
  s(k1) = abs (s(k1))
else
  nx = 0
endif

!   Place exponent first instead of at the very end as in Fortran.

b(1) = '1'
b(2) = '0'
b(3) = ' '
b(4) = '^'
write (ca, '(I10)') nx

do i = 1, 10
  b(i+4) = ca(i:i)
enddo

b(15) = ' '
b(16) = 'x'
b(17) = ' '

!   Insert sign and first digit.

if (ia .eq. -1) then
  b(18) = '-'
else
  b(18) = ' '
endif
if (na .ne. 0) then
  nn = s(k1+2)
else
  nn = 0
endif
write (ca, '(I1)') nn
b(19) = ca(1:1)
b(20) = '.'
ix = 20
if (na .eq. 0) goto 190
f(3) = nn
call mpsub (s(k1), f, s(k0))
if (s(k0) .eq. 0) goto 190
call mpmuld (s(k0), 1.d6, 0, s(k1))
nl = max (mpnw * log10 (mpbdx) / 6.d0 - 1.d0, 1.d0)

!   Insert the digits of the remaining words.

do j = 1, nl
  if (s(k1+1) .eq. 0.) then
    nn = s(k1+2)
    f(1) = 1.
    f(3) = nn
  else
    f(1) = 0.
    nn = 0
  endif
  write (ca, '(I6.6)') nn

  do i = 1, 6
    b(i+ix) = ca(i:i)
  enddo

  ix = ix + 6
  call mpsub (s(k1), f, s(k0))
  call mpmuld (s(k0), 1.d6, 0, s(k1))
  if (s(k1) .eq. 0.) goto 140
enddo

!   Check if trailing zeroes should be trimmed.

j = nl + 1

140  l = ix
if (b(l) .eq. '0' .or. (j .gt. nl .and. b(l-1) .eq. '0' .and. &
  b(l-2) .eq. '0' .and. b(l-3) .eq. '0')) then
  b(l) = ' '

  do i = l - 1, 21, -1
    if (b(i) .ne. '0') then
      ix = i
      goto 190
    endif
    b(i) = ' '
  enddo

  ix = 20

!   Check if trailing nines should be rounded up.

elseif (j .gt. nl .and. b(l-1) .eq. '9' .and. b(l-2) .eq. '9' &
  .and. b(l-3) .eq. '9') then
  b(l) = ' '

  do i = l - 1, 21, -1
    if (b(i) .ne. '9') goto 180
    b(i) = ' '
  enddo

!   We have rounded away all digits to the right of the decimal point, and the
!   digit to the left of the digit is a 9.  Set the digit to 1 and increase
!   the exponent by one.

  ix = 20
  if (b(19) .eq. '9') then
    b(19) = '1'
    write (ca, '(I10)') nx + 1

    do i = 1, 10
      b(i+4) = ca(i:i)
    enddo
  else
    ca = b(19)
    read (ca, '(I1)') nn
    write (ca, '(I1)') nn + 1
    b(19) = ca(1:1)
  endif
  goto 190

180  ca = b(i)
  read (ca, '(I1)') nn
  write (ca, '(I1)') nn + 1
  b(i) = ca(1:1)
  ix = i
endif

190  n = ix
mpnw = nws
if (mpidb .ge. 7) then
  no = min (n, 6 * mpndb + 20)
  write (mpldb, 2) (b(i), i = 1, no)
2 format ('MPOUTC O'/(78a1))
endif
return
end subroutine

subroutine mpfform (a, n1, n2, b)

!   This routine converts the MP number A to F format, i.e. F N1.N2.
!   B is the output array (type CHARACTER*1) of size N1.

real a(mpnw+2)
character*1 b(n1), c(8*mpnw+50)
character*40 chr40

call mpoutc (a, c, n)
write (chr40, '(10a1)') (c(i), i = 5, 14)
read (chr40, '(i10)') ix
if (a(1) .ge. 0.) then
  ls = 0
else
  ls = 1
endif
if (ix .ge. 0 .and. a(1) .ne. 0.) then
  lz = 0
else
  lz = 1
endif
mx = max (ix, 0)

!   Check for overflow of field length.

if (ls + lz + mx + n2 + 2 .gt. n1) then
  do i = 1, n1
    b(i) = '*'
  enddo

  goto 200
endif

!   Check if a zero should be output.

if (a(1) .eq. 0 .or. -ix .gt. n2) then
  do i = 1, n1 - n2 - 2
    b(i) = ' '
  enddo

  b(n1-n2-1) = '0'
  b(n1-n2) = '.'

  do i = 1, n2
    b(i+n1-n2) = '0'
  enddo

  goto 200
endif

!   Process other cases.

do i = 1, n1 - n2 - mx - 2
  b(i) = ' '
enddo

if (a(1) .lt. 0.) b(n1-n2-mx-2) = '-'
if (ix .ge. 0) then
  b(n1-n2-ix-1) = c(19)
  kx = min (n - 20, ix)

  do i = 1, kx
    b(i+n1-n2-ix-1) = c(i+20)
  enddo

  do i = kx + 1, ix
    b(i+n1-n2-ix-1) = '0'
  enddo

  b(n1-n2) = '.'
  kx = max (min (n - ix - 20, n2), 0)

  do i = 1, kx
    b(i+n1-n2) = c(i+ix+20)
  enddo

  do i = kx + 1, n2
    b(i+n1-n2) = '0'
  enddo
else
  nx = - ix
  b(n1-n2-1) = '0'
  b(n1-n2) = '.'

  do i = 1, nx - 1
    b(i+n1-n2) = '0'
  enddo

  b(n1-n2+nx) = c(19)
  kx = min (n - 20, n2 - nx)

  do i = 1, kx
    b(i+n1-n2+nx) = c(i+20)
  enddo

  do i = kx + 1, n2 - nx
    b(i+n1-n2+nx) = '0'
  enddo
endif

200 continue

return
end subroutine

subroutine mpsort (n, la, a, ip)

!   This routine sorts the entries of the N-long MP vector A into ascending
!   order using the quicksort algorithm.  The entries of A are assumed to
!   start at A(1), A(LA+1), A(2*LA+1), etc. The permutation vector that would
!   sort the vector is returned in IP.  Debug output starts with MPIDB = 7.

!   Max integer space for IP: N cells.

character*8 cx
dimension a(la,n), ip(n), ik(50), jk(50), s(2*mpnw+8)

if (mpier .ne. 0) then

  do i = 1, n
    ip(i) = i
  enddo

  return
endif
if (mpidb .ge. 7) then
  write (mpldb, 1) n, la
1 format ('MPSORT I',2i6)

  do k = 1, n
    write (cx, '(I4)') k
    call mpdeb (cx, a(1,k))
  enddo
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4

do i = 1, n
  ip(i) = i
enddo

k = 1
ik(1) = 1
jk(1) = n

130  i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
call mpeq (a(1,ip(j)), s(k0))
j = j - 1

140  do l = i, j
  call mpcpr (s(k0), a(1,ip(l)), ic)
  if (ic .lt. 0) goto 160
enddo

i = j
goto 190

160  i = l

do l = j, i, -1
  call mpcpr (s(k0), a(1,ip(l)), ic)
  if (ic .gt. 0) goto 180
enddo

j = i
goto 190

180  j = l
if (i .ge. j)  goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190  call mpcpr (s(k0), a(1,ip(i)), ic)
if (ic .ge. 0) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200  k = k - 1
jz = 0
if (j .eq. iq)  goto 210
k = k + 1
jk(k) = j
jz = 1

210  i = i + 1
if (i .eq. jq)  goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz .eq. 0)  goto 220
if (j - iq .ge. jq - i)  goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220  if (k .gt. 0)  goto 130

if (mpidb .ge. 7) write (mpldb, 2) ip
2 format ('MPSORT O'/(8i9))
return
end subroutine

end module

module mpfune

!   This module defines algebraic and transcendental routines.

use mpfuna
use mpfunc
use mpfund
contains

subroutine mpang (x, y, pi, a)

!   This computes the MP angle A subtended by the MP pair (X, Y) considered as
!   a point in the x-y plane.  This is more useful than an arctan or arcsin
!   routine, since it places the result correctly in the full circle, i.e.
!   -Pi < A <= Pi.  PI is the MP value of Pi computed by a previous call to
!   MPPI.  For extra high levels of precision, use MPANGX.  The last word of
!   the result is not reliable.  Debug output starts with MPIDB = 5.

!   Max SP space for A: MPNW + 4 cells.

!   The Taylor series for Sin converges much more slowly than that of Arcsin.
!   Thus this routine does not employ Taylor series, but instead computes
!   Arccos or Arcsin by solving Cos (a) = x or Sin (a) = y using one of the
!   following Newton iterations, both of which converge to a:

!     z_{k+1} = z_k - [x - Cos (z_k)] / Sin (z_k)
!     z_{k+1} = z_k + [y - Sin (z_k)] / Cos (z_k)

!   The first is selected if Abs (x) <= Abs (y); otherwise the second is used.
!   These iterations are performed with a maximum precision level MPNW that
!   is dynamically changed, approximately doubling with each iteration.
!   See the comment about the parameter NIT in MPDIVX.

double precision cl2, cpi, t1, t2, t3
parameter (cl2 = 1.4426950408889633d0, cpi = 3.141592653589793d0, nit = 3)
dimension a(mpnw+4), pi(mpnw+2), x(mpnw+2), y(mpnw+2), s(5*mpnw+25)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  a(1) = 0.
  a(2) = 0.
  return
endif
if (mpidb .ge. 5) then
  call mpdeb ('MPANG I', x)
  call mpdeb ('MPANG I', y)
endif

ix = sign (1., x(1))
nx = min (int (abs (x(1))), mpnw)
iy = sign (1., y(1))
ny = min (int (abs (y(1))), mpnw)

!   Check if both X and Y are zero.

if (nx .eq. 0 .and. ny .eq. 0) then
  if (mpker(7) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPANG: Both arguments are zero.')
    mpier = 7
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if Pi has been precomputed.

call mpmdc (pi, t1, n1)
if (n1 .ne. 0 .or. abs (t1 - cpi) .gt. mprx2) then
  if (mpker(8) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPANG: PI must be precomputed.')
    mpier = 8
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if one of X or Y is zero.

if (nx .eq. 0) then
  if (iy .gt. 0) then
    call mpmuld (pi, 0.5d0, 0, a)
  else
    call mpmuld (pi, -0.5d0, 0, a)
  endif
  goto 120
elseif (ny .eq. 0) then
  if (ix .gt. 0) then
    a(1) = 0.
    a(2) = 0.
  else
    call mpeq (pi, a)
  endif
  goto 120
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
k3 = k2 + n5
k4 = k3 + n5
nws = mpnw
mpnw = mpnw + 1

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = nws
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Normalize x and y so that x^2 + y^2 = 1.

call mpmul (x, x, s(k0))
call mpmul (y, y, s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpsqrt (s(k2), s(k3))
call mpdiv (x, s(k3), s(k1))
call mpdiv (y, s(k3), s(k2))

!   Compute initial approximation of the angle.

call mpmdc (s(k1), t1, n1)
call mpmdc (s(k2), t2, n2)
n1 = max (n1, -66)
n2 = max (n2, -66)
t1 = t1 * 2.d0 ** n1
t2 = t2 * 2.d0 ** n2
t3 = atan2 (t2, t1)
call mpdmc (t3, 0, a)

!   The smaller of x or y will be used from now on to measure convergence.
!   This selects the Newton iteration (of the two listed above) that has the
!   largest denominator.

if (abs (t1) .le. abs (t2)) then
  kk = 1
  call mpeq (s(k1), s(k0))
else
  kk = 2
  call mpeq (s(k2), s(k0))
endif

mpnw = 3
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 2, mq
  mpnw = min (2 * mpnw - 2, nws) + 1

100  continue
  call mpcssn (a, pi, s(k1), s(k2))
  if (kk .eq. 1) then
    call mpsub (s(k0), s(k1), s(k3))
    call mpdiv (s(k3), s(k2), s(k4))
    call mpsub (a, s(k4), s(k1))
  else
    call mpsub (s(k0), s(k2), s(k3))
    call mpdiv (s(k3), s(k1), s(k4))
    call mpadd (a, s(k4), s(k1))
  endif
  call mpeq (s(k1), a)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
enddo

!   Restore original precision level.

mpnw = nws
call mproun (a)

120  if (mpidb .ge. 5) call mpdeb ('MPANG O', a)

return
end subroutine

subroutine mpcssh (a, al2, x, y)

!   This computes the hyperbolic cosine and sine of the MP number A and
!   returns the two MP results in X and Y, respectively.  AL2 is the MP value
!   of Log (10) computed by a previous call to MPLOG.  For extra high levels of
!   precision, use MPCSHX.  The last word of the result is not reliable.
!   Debug output starts with MPIDB = 5.

!   Max SP space for X and Y: MPNW + 4 cells.

dimension a(mpnw+2), f(8), al2(mpnw+2), x(mpnw+4), y(mpnw+4), s(4*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  x(1) = 0.
  x(2) = 0.
  y(1) = 0.
  y(2) = 0.
  return
endif
if (mpidb .ge. 5) call mpdeb ('MPCSSH I', a)

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
k3 = k2 + n5
nws = mpnw
mpnw = mpnw + 1
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.

call mpexp (a, al2, s(k0))
call mpdiv (f, s(k0), s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpmuld (s(k2), 0.5d0, 0, s(k3))
call mpeq (s(k3), x)
call mpsub (s(k0), s(k1), s(k2))
call mpmuld (s(k2), 0.5d0, 0, s(k3))
call mpeq (s(k3), y)

!   Restore original precision level.

mpnw = nws
call mproun (x)
call mproun (y)

if (mpidb .ge. 5) then
  call mpdeb ('MPCSSH O', x)
  call mpdeb ('MPCSSH O', y)
endif
return
end subroutine

subroutine mpcssn (a, pi, x, y)

!   This computes the cosine and sine of the MP number A and returns the two MP
!   results in X and Y, respectively.  PI is the MP value of Pi computed by a
!   previous call to MPPI.  For extra high levels of precision, use MPCSSX.
!   The last word of the result is not reliable.  Debug output starts with
!   MPIDB = 6.

!   Max SP space for X and Y: MPNW + 4 cells.

!   This routine uses the conventional Taylor's series for Sin (s):

!   Sin (s) =  s - s^3 / 3! + s^5 / 5! - s^7 / 7! ...

!   where s = t - a * pi / 2 - b * pi / 16 and the integers a and b are chosen
!   to minimize the absolute value of s.  We can then compute

!   Sin (t) = Sin (s + a * pi / 2 + b * pi / 16)
!   Cos (t) = Cos (s + a * pi / 2 + b * pi / 16)

!   by applying elementary trig identities for sums.  The sine and cosine of
!   b * pi / 16 are of the form 1/2 * Sqrt {2 +- Sqrt [2 +- Sqrt(2)]}.
!   Reducing t in this manner insures that -Pi / 32 < s <= Pi / 32, which
!   accelerates convergence in the above series.

double precision cpi, t1, t2
parameter (cpi = 3.141592653589793d0)
dimension a(mpnw+2), f(8), pi(mpnw+2), x(mpnw+4), y(mpnw+4), s(7*mpnw+35)

! GZ 
!y=0d0

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  x(1) = 0.
  x(2) = 0.
  y(1) = 0.
  y(2) = 0.
  return
endif
if (mpidb .ge. 6) call mpdeb ('MPCSSN I', a)

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
if (na .eq. 0) then
  x(1) = 1.
  x(2) = 0.
  x(3) = 1.
  y(1) = 0.
  y(2) = 0.
  l1 = 0
  goto 120
endif

!   Check if Pi has been precomputed.

call mpmdc (pi, t1, n1)
if (n1 .ne. 0 .or. abs (t1 - cpi) .gt. mprx2) then
  if (mpker(28) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPCSSN: PI must be precomputed.')
    mpier = 28
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
k3 = k2 + n5
k4 = k3 + n5
k5 = k4 + n5
k6 = k5 + n5
nws = mpnw
mpnw = mpnw + 1
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.

!   Reduce to between - Pi and Pi.

call mpmuld (pi, 2.d0, 0, s(k0))
call mpdiv (a, s(k0), s(k1))
call mpnint (s(k1), s(k2))
call mpsub (s(k1), s(k2), s(k3))

!   Determine nearest multiple of Pi / 2, and within a quadrant, the nearest
!   multiple of Pi / 16.  Through most of the rest of this subroutine, KA and
!   KB are the integers a and b of the algorithm above.

call mpmdc (s(k3), t1, n1)
if (n1 .ge. - mpnbt) then
  t1 = t1 * 2.d0 ** n1
  t2 = 4.d0 * t1
  ka = nint (t2)
  kb = nint (8.d0 * (t2 - ka))
else
  ka = 0
  kb = 0
endif
t1 = (8 * ka + kb) / 32.d0
call mpdmc (t1, 0, s(k1))
call mpsub (s(k3), s(k1), s(k2))
call mpmul (s(k0), s(k2), s(k1))

!   Compute cosine and sine of the reduced argument s using Taylor's series.

if (s(k1) .eq. 0.) then
  s(k0) = 0.
  s(k0+1) = 0.
  l1 = 0
  goto 110
endif
call mpeq (s(k1), s(k0))
call mpmul (s(k0), s(k0), s(k2))
l1 = 0

100  l1 = l1 + 1
if (l1 .eq. 10000) then
  if (mpker(29) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPCSSN: Iteration limit exceeded.')
    mpier = 29
    if (mpker(mpier) .eq. 2) call mpabrt
    mpnw = nws
    return
  endif
endif

t2 = - (2.d0 * l1) * (2.d0 * l1 + 1.d0)
call mpmul (s(k2), s(k1), s(k3))
call mpdivd (s(k3), t2, 0, s(k1))
call mpadd (s(k1), s(k0), s(k3))
call mpeq (s(k3), s(k0))

!   Check for convergence of the series.

if (s(k1) .ne. 0. .and. s(k1+1) .ge. s(k0+1) - mpnw) goto 100

!   Compute Cos (s) = Sqrt [1 - Sin^2 (s)].

110  call mpeq (s(k0), s(k1))
call mpmul (s(k0), s(k0), s(k2))
call mpsub (f, s(k2), s(k3))
call mpsqrt (s(k3), s(k0))

!   Compute cosine and sine of b * Pi / 16.

kc = abs (kb)
f(3) = 2.
if (kc .eq. 0) then
  s(k2) = 1.
  s(k2+1) = 0.
  s(k2+2) = 1.
  s(k3) = 0.
  s(k3+1) = 0.
else
  if (kc .eq. 1) then
    call mpsqrt (f, s(k4))
    call mpadd (f, s(k4), s(k5))
    call mpsqrt (s(k5), s(k4))
  elseif (kc .eq. 2) then
    call mpsqrt (f, s(k4))
  elseif (kc .eq. 3) then
    call mpsqrt (f, s(k4))
    call mpsub (f, s(k4), s(k5))
    call mpsqrt (s(k5), s(k4))
  elseif (kc .eq. 4) then
    s(k4) = 0.
    s(k4+1) = 0.
  endif
  call mpadd (f, s(k4), s(k5))
  call mpsqrt (s(k5), s(k3))
  call mpmuld (s(k3), 0.5d0, 0, s(k2))
  call mpsub (f, s(k4), s(k5))
  call mpsqrt (s(k5), s(k4))
  call mpmuld (s(k4), 0.5d0, 0, s(k3))
endif
if (kb .lt. 0) s(k3) = - s(k3)

!   Apply the trigonometric summation identities to compute cosine and sine of
!   s + b * Pi / 16.

call mpmul (s(k0), s(k2), s(k4))
call mpmul (s(k1), s(k3), s(k5))
call mpsub (s(k4), s(k5), s(k6))
call mpmul (s(k1), s(k2), s(k4))
call mpmul (s(k0), s(k3), s(k5))
call mpadd (s(k4), s(k5), s(k1))
call mpeq (s(k6), s(k0))


! GZ 
!y=0d0

!   This code in effect applies the trigonometric summation identities for
!   (s + b * Pi / 16) + a * Pi / 2.
!GZ 

!write(*,*) '==============================='
!write(*,*) 'Before calls to mpeq (in mpcssn)'
!write(*,*) 'mpnw',mpnw
!write(*,*) 'ka',ka       
!write(*,*) 's(k0),k0',s(k0),k0
!write(*,*) 's(k1),k1',s(k1),k1
!write(*,*) 'x',x
!write(*,*) 'y',y
!!write(*,*) 

if (ka .eq. 0) then
  call mpeq (s(k0), x)
  call mpeq (s(k1), y)
elseif (ka .eq. 1) then
  call mpeq (s(k1), x)
  x(1) = - x(1)
  call mpeq (s(k0), y)
elseif (ka .eq. -1) then
  call mpeq (s(k1), x)
  call mpeq (s(k0), y)
  y(1) = - y(1)
elseif (ka .eq. 2 .or. ka .eq. -2) then
  call mpeq (s(k0), x)
  x(1) = - x(1)
  call mpeq (s(k1), y)
  y(1) = - y(1)
endif
!write(*,*) 'after call to mpeq'
!write(*,*) 'x',x
!write(*,*) 'y',y
!write(*,*) 
!write(*,*) 

!   Restore original precision level.

mpnw = nws
call mproun (x)
!GZ this wil prob give wrong results, but at least it runs
call mproun (y)

120  if (mpidb .ge. 6) then
  write (mpldb, 3) l1
3 format ('Iteration count:',i5)
  call mpdeb ('MPCSSN O', x)
  call mpdeb ('MPCSSN O', y)
endif
return
end subroutine

subroutine mpexp (a, al2, b)

!   This computes the exponential function of the MP number A and returns the
!   MP result in B.  AL2 is the MP value of Log(2) produced by a prior call
!   to MPLOG.  For extra high levels of precision, use MPEXPX.  The last
!   word of the result is not reliable.  Debug output starts with MPIDB = 7.

!   Max SP space for B: MPNW + 4 cells.

!   This routine uses a modification of the Taylor's series for Exp (t):

!   Exp (t) =  (1 + r + r^2 / 2! + r^3 / 3! + r^4 / 4! ...) ^ q * 2 ^ n

!   where q = 256, r = t' / q, t' = t - n Log(2) and where n is chosen so
!   that -0.5 Log(2) < t' <= 0.5 Log(2).  Reducing t mod Log(2) and
!   dividing by 256 insures that -0.001 < r <= 0.001, which accelerates
!   convergence in the above series.

double precision alt, t1, t2
parameter (alt = 0.693147180559945309d0, nq = 8)
dimension a(mpnw+2), b(mpnw+4), al2(mpnw+2), f(8), s(4*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 7) call mpdeb ('MPEXP I', a)

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
call mpmdc (a, t1, n1)
t1 = t1 * 2.d0 ** n1

!   Unless the argument is near Log (2), Log(2) must be precomputed.  This
!   exception is necessary because MPLOG calls MPEXP to initialize Log (2).

if (abs (t1 - alt) .gt. mprdx) then
  call mpmdc (al2, t2, n2)
  if (n2 .ne. - mpnbt .or. abs (t2 * 0.5d0 ** mpnbt - alt) .gt. mprx2) then
    if (mpker(34) .ne. 0) then
      write (mpldb, 1)
1     format ('*** MPEXP: LOG (2) must be precomputed.')
      mpier = 34
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  endif
endif

!   Check for overflows and underflows.

if (abs (t1) .gt. 33271064.d0) then
  if (t1 .gt. 0.d0) then
    if (mpker(35) .ne. 0) then
      write (mpldb, 2) t1
2     format ('*** MPEXP: Argument is too large',1p,d25.16)
      mpier = 35
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  else
    b(1) = 0.
    b(2) = 0.
    l1 = 0
    goto 130
  endif
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
k3 = k2 + n5
nws = mpnw
mpnw = mpnw + 1
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.

!   Compute the reduced argument A' = A - Log(2) * Nint [A / Log(2)].  Save
!   NZ = Nint [A / Log(2)] for correcting the exponent of the final result.

if (abs (t1 - alt) .gt. mprdx) then
  call mpdiv (a, al2, s(k0))
  call mpnint (s(k0), s(k1))
  call mpmdc (s(k1), t1, n1)
  nz = t1 * 2.d0 ** n1 + sign (mprxx, t1)
  call mpmul (al2, s(k1), s(k2))
  call mpsub (a, s(k2), s(k0))
else
  call mpeq (a, s(k0))
  nz = 0
endif
tl = s(k0+1) - mpnw

!   Check if the reduced argument is zero.

if (s(k0) .eq. 0.d0) then
  s(k0) = 1.
  s(k0+1) = 0.
  s(k0+2) = 1.
  l1 = 0
  goto 120
endif

!   Divide the reduced argument by 2 ^ NQ.

call mpdivd (s(k0), 1.d0, nq, s(k1))

!   Compute Exp using the usual Taylor series.

call mpeq (f, s(k2))
call mpeq (f, s(k3))
l1 = 0

100  l1 = l1 + 1
if (l1 .eq. 10000) then
  if (mpker(36) .ne. 0) then
    write (mpldb, 3)
3   format ('*** MPEXP: Iteration limit exceeded.')
    mpier = 36
    if (mpker(mpier) .eq. 2) call mpabrt
    mpnw = nws
    return
  endif
endif

t2 = l1
call mpmul (s(k2), s(k1), s(k0))
call mpdivd (s(k0), t2, 0, s(k2))
call mpadd (s(k3), s(k2), s(k0))
call mpeq (s(k0), s(k3))

!   Check for convergence of the series.

if (s(k2) .ne. 0. .and. s(k2+1) .ge. tl) goto 100

!   Raise to the (2 ^ NQ)-th power.

do i = 1, nq
  call mpmul (s(k0), s(k0), s(k1))
  call mpeq (s(k1), s(k0))
enddo

!   Multiply by 2 ^ NZ.

120  call mpmuld (s(k0), 1.d0, nz, s(k1))
call mpeq (s(k1), b)

!   Restore original precision level.

mpnw = nws
call mproun (b)

130  if (mpidb .ge. 7) then
  write (mpldb, 4) l1
4 format ('Iteration count:',i5)
  call mpdeb ('MPEXP O', b)
endif
return
end subroutine

subroutine mpinrl (n, lx, x, rb, mt, lr, r, iq)

!   This routine searches for integer relations among the entries of the
!   N-long MP vector X.  An integer relation is an n-long vector r such that
!   r_1 x_1 + r_2 x_2 + ... + r_n x_n = 0.  The entries of x are assumed to
!   start at X(1), X(LX+1), X(2*LX+1), etc.  RB is the maximum Euclidean norm 
!   of an acceptable relation.  IQ is set to 1 if the routine succeeds in 
!   recovering a relation that (1) produces zero to within the input relative 
!   tolerance 10^MT and (2) has Euclidean norm less than RB.  If no relation 
!   is found that meets these standards, IQ is set to 0.  When a valid 
!   relation vector is recovered, it is placed in R, beginning at R(1), 
!   R(LR+1), R(2*LR+1), etc., where LR, like LX, is an input parameter.  Debug
!   output starts with MPIDB = 4.  When MPIDB = 5, norm bounds are output 
!   within which no relation can exist.

!   Max SP space for R: LR * N cells.

!   Be careful to set MT appropriately -- if it is set too highly negative,
!   the program may miss a relation.  MT is typically set to  15 - NDX, where
!   NDX is the number of accurate digits in the X vector.  Repeating a run
!   with somewhat higher precision is highly recommended to certify that bounds
!   results are valid. 

!   This routine uses the simplified PSLQ algorithm as described in a paper by
!   H.R.P. Ferguson, DHB and Arno, "Analysis of PSLQ, An Integer Relation
!   Finding Algorithm".  Contact DHB for a copy if desired.

!   The following parameters are set in this routine (and can be changed):
!     ipi = Iteration print interval when mpidb = 5.
!     itm = Maximum iteration count.  Run is aborted when this is exceeded.

double precision d1, d2, d3, d4, rb, rn
parameter (ipi = 100, itm = 1000000)
real a(mpnw+4,n,n), b(mpnw+4,n,n), h(mpnw+4,n,n-1), ss(mpnw+4,n,n), &
  r(lr,n), s(mpnw+4,n), x(lx,n), y(mpnw+4,n), t1(mpnw+4), &
  t2(mpnw+4), mpm1(mpnw+4), mpm2(mpnw+4), eps(mpnw+4)
double precision mpd1, mpd2, mpd3
logical mpl1

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  iq = 0
  return
endif
if (mpidb .ge. 4) then
  write (mpldb, 1) n, lx, rb, mt, lr
1 format ('MPINRL I',2i6,1pd15.6,2i6)
endif

!  Initialize.

iq = 0
it = 0
rn = 0.d0
call mpdmc (10.d0, 0, mpm1)
call mpnpwr (mpm1, mt, eps)
call mpinrp (n, a, b, h, s, lx, x, y)

!  Perform MP iterations.

100 it = it + 1
if (it .gt. itm) then
  if (mpidb .ge. 4) write (mpldb, 2) itm
2 format ('Iteration limit exceeded',i7)
  goto 120
endif
if (mpidb .ge. 6 .or. mpidb .eq. 5 .and. mod (it, ipi) .eq. 0) &
  write (mpldb, 3) it
3 format ('Iteration',i7)
  call mpinrq (eps, it, n, a, b, h, y, izm)
  if (izm .eq. 2) goto 120
  if (izm .gt. 0) goto 110
  if (mod (it, ipi) .eq. 0) then
  call mpdmc (1d100, 0, t1)
  call mpdmc (0.d0, 0, t2)

!  Compute min, max of y and norm bound.

  do i = 1, n
    call mpeq (y(1,i), mpm1)
    mpm1(1) = abs (mpm1(1))
    call mpcpr (t1, mpm1, mpi1)
    if (mpi1 .ge. 0) then
      call mpeq (mpm1, mpm2)
    else
      call mpeq (t1, mpm2)
    endif
    call mpeq (mpm2, t1)
    call mpeq (y(1,i), mpm1)
    mpm1(1) = abs (mpm1(1))
    call mpcpr (t2, mpm1, mpi1)
    if (mpi1 .ge. 0) then
      call mpeq (t2, mpm2)
    else
      call mpeq (mpm1, mpm2)
    endif
    call mpeq (mpm2, t2)
  enddo

  call mpmdc (t1, d3, n3)
  call dpdec (d3, n3, d1, n1)
  call mpmdc (t2, d3, n3)
  call dpdec (d3, n3, d2, n2)
  mpd3 = 0.d0

  do i = 1, n - 1
    call mpmdc (h(1,i,i), mpd2, mpi1)
    mpd1 = mpd2 * 2.d0 ** mpi1
    mpd2 = abs (mpd1)
    mpd3 = max (mpd3, mpd2)
  enddo

  d3 = 1.d0 / mpd3
  rn = max (rn, d3)
  if (mpidb .ge. 5) then
    write (mpldb, 4) d1, n1, d2, n2, d3, rn
4   format ('Min, max of y =',0p,f11.6,'D',i5,f11.6,'D',i5/ &
     'Norm bound =',1p,d15.6,4x,'Max. bound =',1p,d15.6)
  endif
  if (rn .gt. rb) then
    if (mpidb .ge. 4) write (mpldb, 5) rb
5   format ('Norm bound limit exceeded',1p,d15.6)
    goto 120
  endif
endif
goto 100

!  A relation has been detected.  Output the final norm bound and other info.

110   continue
if (mpidb .ge. 4) then
  call mpdmc (1d100, 0, t1)
  call mpdmc (0.d0, 0, t2)

  do i = 1, n
    call mpeq (y(1,i), mpm1)
    mpm1(1) = abs (mpm1(1))
    call mpcpr (t1, mpm1, mpi1)
    if (mpi1 .ge. 0) then
      call mpeq (mpm1, mpm2)
    else
      call mpeq (t1, mpm2)
    endif
    call mpeq (mpm2, t1)
    call mpeq (y(1,i), mpm1)
    mpm1(1) = abs (mpm1(1))
    call mpcpr (t2, mpm1, mpi1)
    if (mpi1 .ge. 0) then
      call mpeq (t2, mpm2)
    else
      call mpeq (mpm1, mpm2)
    endif
    call mpeq (mpm2, t2)
  enddo

  call mpmdc (t1, d3, n3)
  call dpdec (d3, n3, d1, n1)
  call mpmdc (t2, d3, n3)
  call dpdec (d3, n3, d2, n2)
  mpd3 = 0.d0

  do i = 1, n - 1
    call mpmdc (h(1,i,i), mpd2, mpi1)
    mpd1 = mpd2 * 2.d0 ** mpi1
    mpd2 = abs (mpd1)
    mpd3 = max (mpd3, mpd2)
  enddo

  d3 = 1.d0 / mpd3
  rn = max (rn, d3)
  write (mpldb, 6) it, d1, n1, d2, n2, rn
6 format ('Relation detected.',4x,'No. iterations =',i6/ &
  'Min, max of y =',0p,f11.6,'D',i5,f11.6,'D',i5/'Max. bound =',1p,d15.6)
endif
d3 = 1d100

!  If there is more than one relation, select the one with smallest norm.

do j = 1, n
  call mpeq (y(1,j), mpm1)
  mpm1(1) = abs (mpm1(1))
  call mpcpr (mpm1, eps, mpi1)
  mpl1 = mpi1 .le. 0
  if (mpl1) then
    d1 = 0.d0

    do i = 1, n
      call mpmdc (b(1,i,j), mpd2, mpi1)
      mpd1 = mpd2 * 2.d0 ** mpi1
      mpd2 = mpd1 ** 2
      mpd1 = d1 + mpd2
      d1 = mpd1
    enddo

    d1 = sqrt (d1)
    call mpmdc (y(1,j), d4, n4)
    call dpdec (d4, n4, d2, n2)
    if (mpidb .ge. 4) then
      write (mpldb, 7) j, d1, d2, n2
7     format ('Index of relation =',i4,3x,'Norm =',1p,d15.6,3x, &
        'Residual =',0p,f11.6,'D',i5)
    endif
    if (mpidb .ge. 5) then
      write (mpldb, 8)
8     format ('Relation:')
      call mpmout (1, n, mpnw + 4, b(1,1,j))
    endif
    if (d1 .lt. d3) then
      d3 = d1
      nws = mpnw
      mpnw = lr - 4

      do i = 1, n
        call mpeq (b(1,i,j), r(1,i))
      enddo
      mpnw = nws
    endif
  endif
enddo

if (d3 .le. rb) then
  iq = 1
else
  if (mpidb .ge. 5) write (mpldb, 9)
9 format ('Relation is too large.')
endif

120 return
end subroutine

subroutine mpinrp (n, a, b, h, s, lx, x, y)

!  Initializes MP arrays at the beginning.

real a(mpnw+4,n,n), b(mpnw+4,n,n), h(mpnw+4,n,n-1), s(mpnw+4,n), x(lx,n), &
  y(mpnw+4,n), t1(mpnw+4), mpm1(mpnw+4), mpm2(mpnw+4), mpm3(mpnw+4)
logical mpl1

if (mpidb .ge. 5) then
  write (mpldb, 1)
1 format ('Input x vector:')
  call mpmout (1, n, lx, x)
endif

!  Set a and b to the identity matrix.

do j = 1, n
  do i = 1, n
    call mpdmc (0.d0, 0, a(1,i,j))
    call mpdmc (0.d0, 0, b(1,i,j))
  enddo

  call mpdmc (1.d0, 0, a(1,j,j))
  call mpdmc (1.d0, 0, b(1,j,j))
enddo

call mpdmc (0.d0, 0, t1)

!  Compute the s vector, the square root of the partial sum of squares of x,
!  and compute the y vector, which is the normalized x vector.

do i = n, 1, -1
  call mpnpwr (x(1,i), 2, mpm1)
  call mpadd (t1, mpm1, mpm2)
  call mpeq (mpm2, t1)
  call mpsqrt (t1, mpm1)
  call mpeq (mpm1, s(1,i))
enddo

call mpdmc (1.d0, 0, mpm2)
call mpdiv (mpm2, s(1,1), mpm1)
call mpeq (mpm1, t1)

do i = 1, n
  call mpmul (t1, x(1,i), mpm1)
  call mpeq (mpm1, y(1,i))
  call mpmul (t1, s(1,i), mpm1)
  call mpeq (mpm1, s(1,i))
enddo

!  Compute the initial h matrix.

do i = 1, n
  do j = i + 1, n - 1
    call mpdmc (0.d0, 0, h(1,i,j))
  enddo

  if (i .le. n - 1) then
    call mpdiv (s(1,i+1), s(1,i), mpm1)
    call mpeq (mpm1, h(1,i,i))
  endif

  do j = 1, i - 1
    call mpmul (y(1,i), y(1,j), mpm1)
    call mpmul (s(1,j), s(1,j+1), mpm2)
    call mpdiv (mpm1, mpm2, mpm3)
    call mpeq (mpm3, mpm1)
    mpm1(1) = - mpm1(1)
    call mpeq (mpm1, h(1,i,j))
  enddo
enddo

!  Perform reduction on h, updating y, a, b and h.

do i = 2, n
  do j = i - 1, 1, -1
    call mpeq (h(1,i,j), mpm1)
    mpm1(1) = abs (mpm1(1))
    call mpeq (h(1,j,j), mpm2)
    mpm2(1) = abs (mpm2(1))
    call mpmuld (mpm2, 0.5d0, 0, mpm3)
    call mpcpr (mpm1, mpm3, mpi1)
    mpl1 = mpi1 .gt. 0
    if (mpl1) then
      call mpdiv (h(1,i,j), h(1,j,j), mpm1)
      call mpnint (mpm1, mpm2)
      call mpeq (mpm2, t1)
      call mpmul (t1, y(1,i), mpm1)
      call mpadd (y(1,j), mpm1, mpm2)
      call mpeq (mpm2, y(1,j))

      do k = i, n
        call mpmul (t1, a(1,j,k), mpm1)
        call mpsub (a(1,i,k), mpm1, mpm2)
        call mpeq (mpm2, a(1,i,k))
        call mpmul (t1, b(1,k,i), mpm1)
        call mpadd (b(1,k,j), mpm1, mpm2)
        call mpeq (mpm2, b(1,k,j))
      enddo

      do k = 1, j
        call mpmul (t1, h(1,j,k), mpm1)
        call mpsub (h(1,i,k), mpm1, mpm2)
        call mpeq (mpm2, h(1,i,k))
      enddo
    endif
  enddo
enddo

if (mpidb .ge. 6) then
  write (mpldb, 2)
2 format ('Initial y vector:')
  call mpmout (1, n, mpnw + 4, y)
  write (mpldb, 3)
3 format ('Initial a matrix:')
  call mpmout (n, n, mpnw + 4, a)
  write (mpldb, 4)
4 format ('Initial n matrix:')
  call mpmout (n, n, mpnw + 4, b)
  write (mpldb, 5)
5 format ('Initial h matrix:')
  call mpmout (n, n-1, mpnw + 4, h)
endif

return
end subroutine

subroutine mpinrq (eps, it, n, a, b, h, y, izm)

!  This performs one iteration of the PSLQ algorithm using MP arithmetic.

double precision d1, d2, d3, gam
parameter (gam = 1.154700538d0)
real a(mpnw+4,n,n), b(mpnw+4,n,n), h(mpnw+4,n,n-1), y(mpnw+4,n), &
  t1(mpnw+4), t2(mpnw+4), t3(mpnw+4), t4(mpnw+4), tl1(mpnw+4), &
  mpm1(mpnw+4), mpm2(mpnw+4), mpm3(mpnw+4), eps(mpnw+4)
double precision mpd1
logical mpl1

if (mpidb .ge. 6) write (mpldb, 1) it
1 format ('Iteration',i7,3x,'MP iteration')
 call mpdmc (1.d0, 0, mpm2)
 call mpdiv (mpm2, eps, mpm1)
 call mpeq (mpm1, tl1)

!  Select im = i such that gam^i * |h(i,i)| is maximal.

izm = 0
call mpdmc (0.d0, 0, t1)

do i = 1, n - 1
  mpd1 = gam ** i
  call mpeq (h(1,i,i), mpm1)
  mpm1(1) = abs (mpm1(1))
  call mpmuld (mpm1, mpd1, 0, mpm2)
  call mpeq (mpm2, t2)
  call mpcpr (t2, t1, mpi1)
  mpl1 = mpi1 .gt. 0
  if (mpl1) then
    im = i
    call mpeq (t2, t1)
  endif
enddo

if (mpidb .ge. 6) write (mpldb, 2) im
2 format ('im =',i4)
im1 = im + 1

!  Exchange the im and im+1 entries of y, rows of a and h, and columns of b.

call mpeq (y(1,im), t1)
call mpeq (y(1,im1), y(1,im))
call mpeq (t1, y(1,im1))

do i = 1, n
  call mpeq (a(1,im,i), t1)
  call mpeq (a(1,im1,i), a(1,im,i))
  call mpeq (t1, a(1,im1,i))
  call mpeq (b(1,i,im), t1)
  call mpeq (b(1,i,im1), b(1,i,im))
  call mpeq (t1, b(1,i,im1))
enddo

do i = 1, n - 1
  call mpeq (h(1,im,i), t1)
  call mpeq (h(1,im1,i), h(1,im,i))
  call mpeq (t1, h(1,im1,i))
enddo

!  Update h with permutation produced above.

if (im .le. n - 2) then
  call mpeq (h(1,im,im), t1)
  call mpeq (h(1,im,im1), t2)
  call mpnpwr (t1, 2, mpm1)
  call mpnpwr (t2, 2, mpm2)
  call mpadd (mpm1, mpm2, mpm3)
  call mpsqrt (mpm3, mpm1)
  call mpeq (mpm1, t3)
  call mpdiv (t1, t3, mpm1)
  call mpeq (mpm1, t1)
  call mpdiv (t2, t3, mpm1)
  call mpeq (mpm1, t2)

  do i = im, n
    call mpeq (h(1,i,im), t3)
    call mpeq (h(1,i,im1), t4)
    call mpmul (t1, t3, mpm1)
    call mpmul (t2, t4, mpm2)
    call mpadd (mpm1, mpm2, mpm3)
    call mpeq (mpm3, h(1,i,im))
    call mpmul (t2, t3, mpm1)
    call mpeq (mpm1, mpm2)
    mpm2(1) = - mpm2(1)
    call mpmul (t1, t4, mpm1)
    call mpadd (mpm2, mpm1, mpm3)
    call mpeq (mpm3, h(1,i,im1))
  enddo
endif

!  Perform reduction on h, updating y, a, b and h.

do i = im1, n
  jl = min (i - 1, im1)

  do j = jl, 1, -1
    call mpeq (h(1,i,j), mpm1)
    mpm1(1) = abs (mpm1(1))
    call mpeq (h(1,j,j), mpm2)
    mpm2(1) = abs (mpm2(1))
    call mpmuld (mpm2, 0.5d0, 0, mpm3)
    call mpcpr (mpm1, mpm3, mpi1)
    mpl1 = mpi1 .gt. 0
    if (mpl1) then
      call mpdiv (h(1,i,j), h(1,j,j), mpm1)
      call mpnint (mpm1, mpm2)
      call mpeq (mpm2, t1)
      call mpmul (t1, y(1,i), mpm1)
      call mpadd (y(1,j), mpm1, mpm2)
      call mpeq (mpm2, y(1,j))

      do k = 1, n
        call mpmul (t1, a(1,j,k), mpm1)
        call mpsub (a(1,i,k), mpm1, mpm2)
        call mpeq (mpm2, a(1,i,k))
        call mpmul (t1, b(1,k,i), mpm1)
        call mpadd (b(1,k,j), mpm1, mpm2)
        call mpeq (mpm2, b(1,k,j))
      enddo

      do k = 1, j
        call mpmul (t1, h(1,j,k), mpm1)
        call mpsub (h(1,i,k), mpm1, mpm2)
        call mpeq (mpm2, h(1,i,k))
      enddo
    endif
  enddo
enddo

!  Find the min and max of y and the max of a.

call mpdmc (1d100, 0, t1)
call mpdmc (0.d0, 0, t2)
call mpdmc (0.d0, 0, t3)

do i = 1, n
  call mpeq (y(1,i), mpm1)
  mpm1(1) = abs (mpm1(1))
  call mpcpr (t1, mpm1, mpi1)
  if (mpi1 .ge. 0) then
    call mpeq (mpm1, mpm2)
  else
    call mpeq (t1, mpm2)
  endif
  call mpeq (mpm2, t1)
  call mpeq (y(1,i), mpm1)
  mpm1(1) = abs (mpm1(1))
  call mpcpr (t2, mpm1, mpi1)
  if (mpi1 .ge. 0) then
    call mpeq (t2, mpm2)
  else
    call mpeq (mpm1, mpm2)
  endif
  call mpeq (mpm2, t2)
enddo

do i = im, n
  do j = 1, n
    call mpeq (a(1,i,j), mpm1)
    mpm1(1) = abs (mpm1(1))
    call mpcpr (t3, mpm1, mpi1)
    if (mpi1 .ge. 0) then
      call mpeq (t3, mpm2)
    else
      call mpeq (mpm1, mpm2)
    endif
    call mpeq (mpm2, t3)
  enddo
enddo

call mpcpr (t1, eps, mpi1)
mpl1 = mpi1 .le. 0
if (mpl1) then
  if (mpidb .ge. 5) then
    call mpmdc (t1, d3, n3)
    call dpdec (d3, n3, d1, n1)
    write (mpldb, 3) it, d1, n1
3   format ('Iteration',i7,3x,'itermp: Small value in y =',f11.6,'D',i5)
  endif
  izm = 1
endif

call mpcpr (t3, tl1, mpi1)
mpl1 = mpi1 .gt. 0
if (mpl1) then
  if (mpidb .ge. 4) then
    call mpmdc (t3, d3, n3)
    call dpdec (d3, n3, d1, n1)
    write (mpldb, 4) it, d1, n1
4   format ('Iteration',i7,3x,'itermp: Large value in a =',f11.6,'D',i5)
    write (mpldb, 5)
5   format ('Run aborted.')
  endif
  izm = 2
  return
endif

if (mpidb .ge. 6) then
  write (mpldb, 6)
6 format ('Updated y vector:')
  call mpmout (1, n, mpnw + 4, y)
  write (mpldb, 7)
7 format ('Updated a matrix:')
  call mpmout (n, n, mpnw + 4, a)
  write (mpldb, 8)
8 format ('Updated b matrix:')
  call mpmout (n, n, mpnw + 4, b)
  write (mpldb, 9)
9 format ('Updated h matrix:')
  call mpmout (n, n - 1, mpnw + 4, h)
endif

return
end subroutine

subroutine mplog (a, al2, b)

!   This computes the natural logarithm of the MP number A and returns the MP
!   result in B.  AL2 is the MP value of Log(2) produced by a prior call to
!   MPLOG.  For extra high levels of precision, use MPLOGX.  The last word of
!   the result is not reliable.  Debug output starts with MPIDB = 6.

!   Max SP space for B: MPNW + 4 cells.

!   The Taylor series for Log converges much more slowly than that of Exp.
!   Thus this routine does not employ Taylor series, but instead computes
!   logarithms by solving Exp (b) = a using the following Newton iteration,
!   which converges to b:

!     x_{k+1} = x_k + [a - Exp (x_k)] / Exp (x_k)

!   These iterations are performed with a maximum precision level MPNW that
!   is dynamically changed, approximately doubling with each iteration.
!   See the comment about the parameter NIT in MPDIVX.

double precision alt, cl2, t1, t2
parameter (alt = 0.693147180559945309d0, cl2 = 1.4426950408889633d0, nit = 3)
dimension a(mpnw+2), al2(mpnw+2), b(mpnw+4), s(3*mpnw+15)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 6) call mpdeb ('MPLOG I', a)

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)

if (ia .lt. 0 .or. na .eq. 0) then
  if (mpker(50) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPLOG: Argument is less than or equal to zero.')
    mpier = 50
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Unless the input is close to 2, Log (2) must have been precomputed.

call mpmdc (a, t1, n1)
if (abs (t1 - 2.d0) .gt. 1d-3 .or. n1 .ne. 0) then
  call mpmdc (al2, t2, n2)
  if (n2 .ne. - mpnbt .or. abs (t2 * 0.5d0 ** mpnbt - alt) .gt. mprx2) then
    if (mpker(51) .ne. 0) then
      write (mpldb, 2)
2     format ('*** MPLOG: LOG (2) must be precomputed.')
      mpier = 51
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  endif
endif

!   Check if input is exactly one.

if (a(1) .eq. 1. .and. a(2) .eq. 0. .and. a(3) .eq. 1.) then
  b(1) = 0.
  b(2) = 0.
  goto 120
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
nws = mpnw

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t2 = nws
mq = cl2 * log (t2) + 1.d0 - mprxx

!   Compute initial approximation of Log (A).

t1 = log (t1) + n1 * alt
call mpdmc (t1, 0, b)
mpnw = 3
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 2, mq
  mpnw = min (2 * mpnw - 2, nws) + 1
100  continue
  call mpexp (b, al2, s(k0))
  call mpsub (a, s(k0), s(k1))
  call mpdiv (s(k1), s(k0), s(k2))
  call mpadd (b, s(k2), s(k1))
  call mpeq (s(k1), b)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
enddo

!   Restore original precision level.

mpnw = nws
call mproun (b)

120  if (mpidb .ge. 6) call mpdeb ('MPLOG O', b)

return
end subroutine

subroutine mppi (pi)

!   This computes Pi to available precision (MPNW mantissa words).  For extra
!   high levels of precision, use MPPIX.  The last word of the result is not
!   reliable.  Debug output starts with MPIDB = 7.

!   Max SP space for PI: MPNW + 4 cells.

!   The algorithm that is used for computing Pi, which is due to Salamin
!   and Brent, is as follows:

!   Set  A_0 = 1,  B_0 = 1/Sqrt(2)  and  D_0 = Sqrt(2) - 1/2.

!   Then from k = 1 iterate the following operations:

!   A_k = 0.5 * (A_{k-1} + B_{k-1})
!   B_k = Sqrt (A_{k-1} * B_{k-1})
!   D_k = D_{k-1} - 2^k * (A_k - B_k) ^ 2

!   Then  P_k = (A_k + B_k) ^ 2 / D_k  converges quadratically to Pi.
!   In other words, each iteration approximately doubles the number of correct
!   digits, providing all iterations are done with the maximum precision.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0)
dimension f(8), pi(mpnw+4), s(5*mpnw+25)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  pi(1) = 0.
  pi(2) = 0.
  return
endif

!   Perform calculations to one extra word accuracy.

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
k3 = k2 + n5
k4 = k3 + n5
nws = mpnw
mpnw = mpnw + 1

!   Determine the number of iterations required for the given precision level.
!   This formula is good only for this Pi algorithm.

t1 = nws * log10 (mpbdx)
mq = cl2 * (log (t1) - 1.d0) + 1.d0

!   Initialize as above.

s(k0) = 1.
s(k0+1) = 0.
s(k0+2) = 1.
f(1) = 1.
f(2) = 0.
f(3) = 2.
f(4) = 0.
call mpsqrt (f, s(k2))
call mpmuld (s(k2), 0.5d0, 0, s(k1))
f(2) = -1.
f(3) = 0.5d0 * mpbdx
call mpsub (s(k2), f, s(k4))

!   Perform iterations as described above.

do k = 1, mq
  call mpadd (s(k0), s(k1), s(k2))
  call mpmul (s(k0), s(k1), s(k3))
  call mpsqrt (s(k3), s(k1))
  call mpmuld (s(k2), 0.5d0, 0, s(k0))
  call mpsub (s(k0), s(k1), s(k2))
  call mpmul (s(k2), s(k2), s(k3))
  t1 = 2.d0 ** k
  call mpmuld (s(k3), t1, 0, s(k2))
  call mpsub (s(k4), s(k2), s(k3))
  call mpeq (s(k3), s(k4))
enddo

!   Complete computation.

call mpadd (s(k0), s(k1), s(k2))
call mpmul (s(k2), s(k2), s(k3))
call mpdiv (s(k3), s(k4), s(k2))
call mpeq (s(k2), pi)

!   Restore original precision level.

mpnw = nws
call mproun (pi)

if (mpidb .ge. 7) call mpdeb ('MPPI O', pi)
return
end subroutine

subroutine mppol (n, l, a, x1, nx, x)

!   This finds a real root of the N-th degree polynomial whose MP coefficients
!   are in A by Newton-Raphson iterations, beginning at the DPE value (X1, NX)
!   and returns the MP root in X.  The N + 1 coefficients a_0, a_1, ..., a_N
!   are assumed to start in locations A(1), A(L+1), A(2*L+1), etc.  For extra
!   high levels of precision, use MPPOLX.  The last word of the result is not
!   reliable.  Debug output starts with MPIDB = 6.

!   Max SP space for X: MPNW + 4 cells.

!   One requirement for this routine to work is that the desired root is not
!   a repeated root.  If one wishes to apply this routine to find a repeated
!   root, it is first necessary to reduce the polynomial to one that has only
!   simple roots.  This can be done by performing the Euclidean algorithm in
!   the ring of polynomials to determine the greatest common divisor Q(t) of
!   P(t) and P'(t).  Here P(t) is the polynomial a_0 + a_1 t + a_2 t^2 +
!   ... + a_n t^n, and P'(t) is the derivative of P(t).  Then R(t) = P(t)/Q(t)
!   is a polynomial that has only simple roots.

!   This routine employs the standard form of the Newton-Raphson iteration:

!   X_{k+1} = X_k - P(X_k) / P'(X_k)

!   These iterations are performed with a maximum precision level MPNW that is
!   dynamically changed, approximately doubling with each iteration.

character*8 cx
double precision t1, x1
dimension a(l,n+1), x(mpnw+4), s(5*mpnw+25)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  x(1) = 0.
  x(2) = 0.
  return
endif
if (mpidb .ge. 6) then
  write (mpldb, 1) n
1 format ('MPPOL I',i4)

  do k = 0, n
    write (cx, '(I4)') k
    call mpdeb (cx, a(1,k+1))
  enddo

  write (mpldb, 2) x1, nx
2 format ('MPPOL I',f16.12,' x 10 ^',i6)
endif

!   Check if the polynomial is proper.

if (a(1,1) .eq. 0. .or. a(1,n+1) .eq. 0.) then
  if (mpker(63) .ne. 0) then
    write (mpldb, 3)
3   format ('*** MPPOL: Either the first or last input coefficient is zero.')
    mpier = 63
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
k3 = k2 + n5
k4 = k3 + n5
nws = mpnw
mpnw = mpnw + 1

!   Set the initial value.

call mpdmc (x1, nx, s(k0))
mpnw = 5
tl = -4.
l1 = 0
ls = -10

!   Perform MP Newton-Raphson iterations to solve P(x) = 0.

110  l1 = l1 + 1
if (l1 .eq. 50) then
  if (mpker(64) .ne. 0) then
    write (mpldb, 4)
4   format ('*** MPPOL: Iteration limit exceeded.')
    mpier = 64
    if (mpker(mpier) .eq. 2) call mpabrt
    mpnw = nws
    return
  endif
endif

!   Compute P(x).

call mpeq (a(1,n+1), s(k1))

do k = n - 1, 0, -1
  call mpmul (s(k0), s(k1), s(k2))
  call mpadd (s(k2), a(1,k+1), s(k1))
enddo

!   Compute P'(x).

t1 = n
call mpmuld (a(1,n+1), t1, 0, s(k2))

do k = n - 1, 1, -1
  call mpmul (s(k0), s(k2), s(k3))
  t1 = k
  call mpmuld (a(1,k+1), t1, 0, s(k4))
  call mpadd (s(k3), s(k4), s(k2))
enddo

!   Compute P(x) / P'(x) and update x.

call mpdiv (s(k1), s(k2), s(k3))
call mpsub (s(k0), s(k3), s(k4))

if (mpidb .ge. 7) then
  write (mpldb, 5) l1
5 format ('Iteration',i4)
  call mpdeb ('X', s(k0))
  call mpdeb ('P(X)', s(k1))
  call mpdeb ('P''(X)', s(k2))
  call mpdeb ('CORR', s(k3))
endif
call mpeq (s(k4), s(k0))

!   If this was the second iteration at full precision, there is no need to
!   continue (the adjusted value of x is correct); otherwise repeat.

if (l1 .eq. ls + 1) goto 140
if (s(k3) .ne. 0. .and. s(k3+1) .gt. tl) goto 110

!   Newton iterations have converged to current precision.  Increase precision
!   and continue.

if (mpnw .eq. nws + 1) goto 140
mpnw = min (2 * mpnw - 2, nws) + 1
if (mpnw .eq. nws + 1) ls = l1
tl = 1 - mpnw
if (mpidb .ge. 7) then
  write (mpldb, 6) mpnw
6 format (6x,'New MPNW =', i8)
endif
goto 110

140  call mpeq (s(k0), x)

!   Restore original precision level.

mpnw = nws
call mproun (x)

if (mpidb .ge. 6) then
  write (mpldb, 7) l1
7 format ('Iteration count:',i5)
  call mpdeb ('MPPOL O', x)
endif
return
end subroutine

end module

module mpfunf

!   This module defines complex arithmetic routines.

use mpfuna
use mpfunc
use mpfund
contains

subroutine mpcadd (l, a, b, c)

!   This computes the sum of the MPC numbers A and B and returns the MPC
!   result in C.  L is the offset between real and imaginary parts in A, B
!   and C.  L must be at least MPNW + 4.  Debug output starts with MPIDB = 9.

!   Max SP space for C: 2 * L cells.

dimension a(2*l), b(2*l), c(2*l)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  c(l+1) = 0.
  c(l+2) = 0.
  return
endif
if (mpidb .ge. 9) write (mpldb, 1)
1 format ('MPCADD')

if (l .lt. mpnw + 4) then
  if (mpker(11) .ne. 0) then
    write (mpldb, 2) l, mpnw + 4
2   format ('*** MPCADD: Offset parameter is too small',2i8)
    mpier = 11
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

l1 = l + 1
call mpadd (a, b, c)
call mpadd (a(l1), b(l1), c(l1))

return
end subroutine

subroutine mpcdiv (l, a, b, c)

!   This routine divides the MP complex numbers A and B to yield the MPC
!   quotient C.  L is the offset between real and imaginary parts in A, B
!   and the result C.  L must be at least MPNW + 4.  For extra high levels of
!   precision, use MPCDVX.  The last word of the result is not reliable.  
!   Debug output starts with MPIDB = 7

!   Max SP space for C: 2 * L cells.

!   This routine employs the formula described in MPCMUL to save multiprecision
!   multiplications.

dimension a(2*l), b(2*l), c(2*l), f(8), s(5*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  c(l+1) = 0.
  c(l+2) = 0.
  return
endif
l1 = l + 1
if (mpidb .ge. 7) then
  write (mpldb, 1) l
1 format ('MPCDIV I',i10)
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPCDIV I'/(6f12.0))
  no = min (int (abs (a(l1))), mpndb) + 2
  write (mpldb, 2) (a(l+i), i = 1, no)
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 2) (b(i), i = 1, no)
  no = min (int (abs (b(l1))), mpndb) + 2
  write (mpldb, 2) (b(l+i), i = 1, no)
endif

if (l .lt. mpnw + 4) then
  if (mpker(15) .ne. 0) then
    write (mpldb, 3) l, mpnw + 4
3   format ('*** MPCDIV: Offset parameter is too small',2i8)
    mpier = 15
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

if (b(1) .eq. 0. .and. b(l1) .eq. 0.) then
  if (mpker(16) .ne. 0) then
    write (mpldb, 4)
4   format ('*** MPCDIV: Divisor is zero.')
    mpier = 16
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
k3 = k2 + n4
k4 = k3 + n4
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.

call mpmul (a, b, s(k0))
call mpmul (a(l1), b(l1), s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpsub (s(k0), s(k1), s(k3))
call mpadd (a, a(l1), s(k0))
call mpsub (b, b(l1), s(k1))
call mpmul (s(k0), s(k1), s(k4))
call mpsub (s(k4), s(k3), s(k1))
call mpmul (b, b, s(k0))
call mpmul (b(l1), b(l1), s(k3))
call mpadd (s(k0), s(k3), s(k4))
call mpdiv (f, s(k4), s(k0))
call mpmul (s(k2), s(k0), c)
call mpmul (s(k1), s(k0), c(l1))

if (mpidb .ge. 7) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 5) (c(i), i = 1, no)
5 format ('MPCDIV O'/(6f12.0))
  no = min (int (abs (c(l1))), mpndb) + 2
  write (mpldb, 5) (c(l+i), i = 1, no)
endif
return
end subroutine

subroutine mpceq (l, a, b)

!   This sets the MPC number B equal to the MPC number A.  L is the offset
!   between real and imaginary parts in A and B.  Debug output starts with
!   MPIDB = 10.

!   Max SP space for B: 2 * L cells.

dimension a(2*l), b(2*l)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  b(l+1) = 0.
  b(l+2) = 0.
  return
endif
if (mpidb .ge. 10) write (mpldb, 1)
1 format ('MPCEQ')

i1 = sign (1., a(1))
n1 = min (int (abs (a(1))), mpnw, l - 2)
i2 = sign (1., a(l+1))
n2 = min (int (abs (a(l+1))), mpnw, l - 2)
b(1) = sign (n1, i1)
b(l+1) = sign (n2, i2)

do i = 2, n1 + 2
  b(i) = a(i)
enddo

do i = 2, n2 + 2
  b(l+i) = a(l+i)
enddo

return
end subroutine

subroutine mpcmul (l, a, b, c)

!   This routine multiplies the MP complex numbers A and B to yield the MPC
!   product C.  L is the offset between real and imaginary parts in A, B and
!   the result C.  L must be at least MPNW + 4.  For extra high levels of
!   precision, use MPCMLX.  The last word of the result is not reliable.  
!   Debug output starts with MPIDB = 7.

!   Max SP space for C: 2 * L cells.

!   This routine employs the formula

!   (a_1 + a_2 i) (b_1 + b_2 i)  =  [a_1 b_1 - a_2 b_2]  +
!    [(a_1 + a_2) (b_1 + b_2) - (a_1 b_1 + a_2 b_2)] i

!   Note that this formula can be implemented with only three multiplications
!   whereas the conventional formula requires four.

dimension a(2*l), b(2*l), c(2*l), s(4*mpnw+16)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  c(l+1) = 0.
  c(l+2) = 0.
  return
endif
l1 = l + 1
if (mpidb .ge.7) then
  write (mpldb, 1) l
1 format ('MPCMUL I',i10)
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPCMUL I'/(6f12.0))
  no = min (int (abs (a(l1))), mpndb) + 2
  write (mpldb, 2) (a(l+i), i = 1, no)
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 2) (b(i), i = 1, no)
  no = min (int (abs (b(l1))), mpndb) + 2
  write (mpldb, 2) (b(l+i), i = 1, no)
endif

if (l .lt. mpnw + 4) then
  if (mpker(20) .ne. 0) then
    write (mpldb, 3) l, mpnw + 4
3   format ('*** MPCMUL: Offset parameter is too small',2i8)
    mpier = 20
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
k3 = k2 + n4

call mpmul (a, b, s(k0))
call mpmul (a(l1), b(l1), s(k1))
call mpsub (s(k0), s(k1), c)
call mpadd (s(k0), s(k1), s(k2))
call mpadd (a, a(l1), s(k0))
call mpadd (b, b(l1), s(k1))
call mpmul (s(k0), s(k1), s(k3))
call mpsub (s(k3), s(k2), c(l1))

if (mpidb .ge. 7) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 4) (c(i), i = 1, no)
4 format ('MPCMUL O'/(6f12.0))
  no = min (int (abs (c(l1))), mpndb) + 2
  write (mpldb, 4) (c(l+i), i = 1, no)
endif
return
end subroutine

subroutine mpcpol (n, la, a, x1, nx, lx, x)

!   This routine finds a complex root of the N-th degree polynomial whose
!   MPC coefficients are in A by Newton-Raphson iterations, beginning
!   at the complex DPE value (X1(1), NX(1)) + i (X1(2), NX(2)), and returns
!   the MPC root in X.  The N + 1 coefficients a_0, a_1, ..., a_N are
!   assumed to start in locations A(1), A(2*LA+1), A(4*LA+1), etc.  LA is the
!   offset between the real and the imaginary parts of each input coefficient.
!   Typically LA = MPNW + 4.  LX, also an input parameter, is the offset
!   between the real and the imaginary parts of the result to be stored in X.
!   LX must be at least MPNW + 4.  For extra high levels of precision, use
!   MPCPLX.  Debug output starts with MPIDB = 5.

!   Max SP space for X: 2 * LX cells.

!   See the note about repeated roots in MPPOL.

!   This routine employs the complex form of the Newton-Raphson iteration:

!   X_{k+1} = X_k - P(X_k) / P'(X_k)

!   These iterations are performed with a maximum precision level MPNW that is
!   dynamically changed, approximately doubling with each iteration.

character*8 cx
double precision t1, x1
dimension a(2*la,n+1), nx(2), x(2*lx), x1(2), s(10*mpnw+50)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  x(1) = 0.
  x(2) = 0.
  x(lx+1) = 0.
  x(lx+2) = 0.
endif
if (mpidb .ge. 5) then
  write (mpldb, 1) n, lx
1 format ('MPCPOL I',2i6)

  do k = 0, n
    write (cx, '(I4)') k
    call mpdeb (cx, a(1,k+1))
    call mpdeb (cx, a(la+1,k+1))
  enddo

  write (mpldb, 2) x1(1), nx(1), x1(2), nx(2)
2 format ('MPCPOL I',f16.12,' x 10 ^',i6,f20.12,' x 10^',i6)
endif

!  Check if the polynomial is proper.

if (a(1,1) .eq. 0. .or. a(1,n+1) .eq. 0.) then
  if (mpker(23) .ne. 0) then
    write (mpldb, 3)
3   format ('*** MPCPOL: Either the first or last input coefficient is zero.')
    mpier = 23
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n5 = mpnw + 5
n10 = 2 * n5
k0 = 1
k1 = k0 + n10
k2 = k1 + n10
k3 = k2 + n10
k4 = k3 + n10
nws = mpnw
mpnw = mpnw + 1

!   Set the initial value.

call mpdmc (x1(1), nx(1), s(k0))
call mpdmc (x1(2), nx(2), s(k0+n5))
mpnw = 5
tl = -4.
l1 = 0
ls = -10

!   Perform MP Newton-Raphson iterations to solve P(x) = 0.

110  l1 = l1 + 1
if (l1 .eq. 50) then
  if (mpker(24) .ne. 0) then
    write (mpldb, 4)
4   format ('*** MPCPOL: Iteration limit exceeded.')
    mpier = 24
    if (mpker(mpier) .eq. 2) call mpabrt
    mpnw = nws
    return
  endif
endif

!   Compute P(x).

call mpmmpc (a(1,n+1), a(la+1,n+1), n5, s(k1))

do k = n - 1, 0, -1
  call mpcmul (n5, s(k0), s(k1), s(k2))
  call mpadd (s(k2), a(1,k+1), s(k1))
  call mpadd (s(k2+n5), a(la+1,k+1), s(k1+n5))
enddo

!   Compute P'(x).

t1 = n
call mpmuld (a(1,n+1), t1, 0, s(k2))
call mpmuld (a(la+1,n+1), t1, 0, s(k2+n5))

do k = n - 1, 1, -1
  call mpcmul (n5, s(k0), s(k2), s(k3))
  t1 = k
  call mpmuld (a(1,k+1), t1, 0, s(k4))
  call mpmuld (a(la+1,k+1), t1, 0, s(k4+n5))
  call mpcadd (n5, s(k3), s(k4), s(k2))
enddo

!   Compute P(x) / P'(x) and update x.

call mpcdiv (n5, s(k1), s(k2), s(k3))
call mpcsub (n5, s(k0), s(k3), s(k4))

if (mpidb .ge. 6) then
  write (mpldb, 5) l1
5 format ('Iteration',i4)
  call mpdeb ('X', s(k0))
  call mpdeb (' ', s(k0+n5))
  call mpdeb ('P(X)', s(k1))
  call mpdeb (' ', s(k1+n5))
  call mpdeb ('P''(X)', s(k2))
  call mpdeb (' ', s(k2+n5))
  call mpdeb ('CORR', s(k3))
  call mpdeb (' ', s(k3+n5))
endif
call mpceq (n5, s(k4), s(k0))

!   If this was the second iteration at full precision, there is no need to
!   continue (the adjusted value of x is correct); otherwise repeat.

if (l1 .eq. ls + 1) goto 140
if (s(k3) .ne. 0. .and. s(k3+1) .gt. tl .or. s(k3+n5) .ne. 0. &
  .and. s(k3+n5+1) .gt. tl) goto 110

!   Newton iterations have converged to current precision.  Increase precision
!   and continue.

if (mpnw .eq. nws + 1) goto 140
mpnw = min (2 * mpnw - 2, nws) + 1
if (mpnw .eq. nws + 1) ls = l1
tl = 1 - mpnw
if (mpidb .ge. 6) then
  write (mpldb, 6) mpnw
6 format (6x,'New MPNW =', i8)
endif
goto 110

140  call mpmmpc (s(k0), s(k0+n5), lx, x)

!   Restore original precision level.

mpnw = nws
call mproun (x)
call mproun (x(lx+1))

if (mpidb .ge. 5) then
  write (mpldb, 7) l1
7 format ('Iteration count:',i5)
  call mpdeb ('MPCPOL O', x)
  call mpdeb (' ', x(lx+1))
endif
return
end subroutine

subroutine mpcpwr (l, a, n, b)

!   This computes the N-th power of the MPC number A and returns the MPC
!   result C in B.  When N is zero, 1 is returned.  When N is negative, the
!   reciprocal of A ^ |N| is returned.  L is the offset between real and
!   imaginary parts in A and B.  L should be at least MPNW + 4.  For extra high
!   levels of precision, use MPCPWX.  Debug output starts with MPIDB = 7.

!   Max SP space for B: 2 * L cells.

!   This routine employs the binary method for exponentiation.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0)
dimension a(2*l), b(2*l), f1(8), f2(8), s(6*mpnw+30)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  b(l+1) = 0.
  b(l+2) = 0.
  return
endif
l1 = l + 1
if (mpidb .ge. 7) then
  write (mpldb, 1) l, n
1 format ('MPCPWR I',2i10)
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPCPWR I'/(6f12.0))
  no = min (int (abs (a(l1))), mpndb) + 2
  write (mpldb, 2) (a(l+i), i = 1, no)
endif

na1 = min (int (abs (a(1))), mpnw)
na2 = min (int (abs (a(l1))), mpnw)
if (na1 .eq. 0 .and. na2 .eq. 0) then
  if (n .ge. 0) then
    b(1) = 0.
    b(2) = 0.
    b(l1) = 0.
    b(l1+1) = 0.
    goto 120
  else
    if (mpker(25) .ne. 0) then
      write (mpldb, 3)
3     format ('*** MPCPWR: Argument is zero and N is negative or zero.')
      mpier = 25
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  endif
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + 2 * n5
k2 = k1 + 2 * n5
nws = mpnw
mpnw = mpnw + 1
nn = abs (n)
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
f2(1) = 0.
f2(2) = 0.
call mpmmpc (a, a(l1), n5, s(k0))
if (nn .eq. 0) then
  call mpmmpc (f1, f2, l, b)
  mpnw = nws
  goto 120
elseif (nn .eq. 1) then
  call mpceq (n5, s(k0), s(k2))
  goto 110
elseif (nn .eq. 2) then
  call mpcmul (n5, s(k0), s(k0), s(k2))
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN .GT. NN.

t1 = nn
mn = cl2 * log (t1) + 1.d0 + mprxx
call mpmmpc (f1, f2, n5, s(k2))
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn .ne. 2 * kk) then
    call mpcmul (n5, s(k2), s(k0), s(k1))
    call mpceq (n5, s(k1), s(k2))
  endif
  kn = kk
  if (j .lt. mn) then
    call mpcmul (n5, s(k0), s(k0), s(k1))
    call mpceq (n5, s(k1), s(k0))
  endif
enddo

!   Compute reciprocal if N is negative.

110  if (n .lt. 0) then
  call mpmmpc (f1, f2, n5, s(k1))
  call mpcdiv (n5, s(k1), s(k2), s(k0))
  call mpceq (n5, s(k0), s(k2))
endif
call mpmmpc (s(k2), s(n5+k2), l, b)

!   Restore original precision level.

mpnw = nws
call mproun (b)
call mproun (b(l1))

120  if (mpidb .ge. 7) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 4) (b(i), i = 1, no)
4 format ('MPCPWR O'/(6f12.0))
  no = min (int (abs (b(l1))), mpndb) + 2
  write (mpldb, 4) (b(l+i), i = 1, no)
endif
return
end subroutine

subroutine mpcsqt (l, a, b)

!   This routine computes the complex square root of the MPC number C.  L is
!   the offset between real and imaginary parts in A and B.  L must be at
!   least MPNW + 4.  For extra high levels of precision, use MPCSQX.  The last
!   word of the result is not reliable.  Debug output starts with MPIDB = 6.

!   Max SP space for B: 2 * L cells.

!   This routine uses the following formula, where A1 and A2 are the real and
!   imaginary parts of A, and where R = Sqrt [A1 ^ 2 + A2 ^2]:

!   B = Sqrt [(R + A1) / 2] + I Sqrt [(R - A1) / 2]

!   If the imaginary part of A is < 0, then the imaginary part of B is also
!   set to be < 0.

dimension a(2*l), b(2*l), s(3*mpnw+12)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  b(l+1) = 0.
  b(l+2) = 0.
  return
endif
l1 = l + 1
if (mpidb .ge. 6) then
  write (mpldb, 1) l
1 format ('MPCSQT I',i10)
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPCSQT I'/(6f12.0))
  no = min (int (abs (a(l1))), mpndb) + 2
  write (mpldb, 2) (a(l+i), i = 1, no)
endif

if (a(1) .eq. 0. .and. a(l+1) .eq. 0.) then
  b(1) = 0.
  b(2) = 0.
  b(l+1) = 0.
  b(l+2) = 0.
  goto 100
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4

call mpmul (a, a, s(k0))
call mpmul (a(l1), a(l1), s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpsqrt (s(k2), s(k0))
call mpeq (a, s(k1))
s(k1) = abs (s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpmuld (s(k2), 0.5d0, 0, s(k1))
call mpsqrt (s(k1), s(k0))
call mpmuld (s(k0), 2.d0, 0, s(k1))
if (a(1) .ge. 0.) then
  call mpeq (s(k0), b)
  call mpdiv (a(l1), s(k1), b(l1))
else
  call mpdiv (a(l1), s(k1), b)
  b(1) = abs (b(1))
  call mpeq (s(k0), b(l1))
  b(l1) = sign (b(l1), a(l1))
endif

100  if (mpidb .ge. 6) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPCSQT O'/(6f12.0))
  no = min (int (abs (b(l1))), mpndb) + 2
  write (mpldb, 3) (b(l+i), i = 1, no)
endif
return
end subroutine

subroutine mpcsub (l, a, b, c)

!   This subracts the MPC numbers A and B and returns the MPC difference in
!   C.  L is the offset between real and imaginary parts in A, B and C.  L
!   must be at least MPNW + 4.  Debug output starts with MPIDB = 9.

!   Max SP space for C: 2 * L cells.

dimension a(2*l), b(2*l), c(2*l)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  c(l+1) = 0.
  c(l+2) = 0.
  return
endif
if (mpidb .ge. 9) write (mpldb, 1)
1 format ('MPCSUB')

l1 = l + 1
call mpsub (a, b, c)
call mpsub (a(l1), b(l1), c(l1))

return
end subroutine

subroutine mpmmpc (a, b, l, c)

!   This converts MP numbers A and B to MPC form in C, i.e. C = A + B i.
!   L (an input parameter) is the offset between real and imaginary parts in
!   C.  Debug output starts with MPIDB = 10.

!   Max SP space for C: 2 * L cells.

dimension a(mpnw+2), b(mpnw+2), c(2*l)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  c(l+1) = 0.
  c(l+2) = 0.
  return
endif
if (mpidb .ge. 10) write (mpldb, 1)
1 format ('MPMMPC')

i1 = sign (1., a(1))
n1 = min (int (abs (a(1))), mpnw, l - 2)
i2 = sign (1., b(1))
n2 = min (int (abs (b(1))), mpnw, l - 2)
c(1) = sign (n1, i1)
c(l+1) = sign (n2, i2)

do i = 2, n1 + 2
  c(i) = a(i)
enddo

do i = 2, n2 + 2
  c(l+i) = b(i)
enddo

return
end subroutine

subroutine mpmpcm (l, a, b, c)

!   This converts the MPC number A to its MP real and imaginary parts, i.e.
!   B = Real (A) and C = Imag (A).  L is the offset between real and
!   imaginary parts in A.  Debug output starts with MPIDB = 10.

!   Max SP space for B and C: MPNW + 2 cells.

dimension a(2*l), b(mpnw+2), c(mpnw+2)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 10) write (mpldb, 1)
1 format ('MPMPCM')

i1 = sign (1., a(1))
n1 = min (int (abs (a(1))), mpnw, l - 2)
i2 = sign (1., a(l+1))
n2 = min (int (abs (a(l+1))), mpnw, l - 2)
b(1) = sign (n1, i1)
c(1) = sign (n2, i2)

do i = 2, n1 + 2
  b(i) = a(i)
enddo

do i = 2, n2 + 2
  c(i) = a(l+i)
enddo

return
end subroutine

end module

module mpfung

!   This module defines the core routines of the extra-high precision package,
!   including FFT and convolution routines.

use mpfuna
contains

subroutine mpinix (n)

!   This initializes the root of unity arrays UU1 and UU2, which are required
!   by the FFT routines that are called by MPMULX.  Before calling any of the
!   advanced MP routines (i.e. those whose names end in X), this routine must
!   be called with N set to the largest precision level MPNW that will be used
!   in the subsequent application.  It is not necessary for the user to call 
!   MPINIX if the advanced routines are not called.

double precision cl2, pi, t1, ti, tpn
parameter (cl2 = 1.4426950408889633d0, pi = 3.141592653589793238d0)

t1 = 0.75d0 * n
m = cl2 * log (t1) + 1.d0 - mprxx
mq = m + 2
nq = 2 ** mq
allocate (mpuu1(nq))
allocate (mpuu2(mq+nq))
mpuu1(1) = mq
ku = 2
ln = 1

do j = 1, mq
  t1 = pi / ln

  do i = 0, ln - 1
    ti = i * t1
    mpuu1(i+ku) = cmplx (cos (ti), sin (ti), mpkdp)
  enddo

  ku = ku + ln
  ln = 2 * ln
enddo

ku = mq + 1
mpuu2(1) = mq

do k = 2, mq - 1
  mpuu2(k) = ku
  mm = k
  nn = 2 ** mm
  mm1 = (mm + 1) / 2
  mm2 = mm - mm1
  nn1 = 2 ** mm1
  nn2 = 2 ** mm2
  tpn = 2.d0 * pi / nn

  do j = 0, nn2 - 1
    do i = 0, nn1 - 1
      iu = ku + i + j * nn1
      t1 = tpn * i * j
      mpuu2(iu) = cmplx (cos (t1), sin (t1), mpkdp)
    enddo
  enddo

  ku = ku + nn
enddo

return
end subroutine

subroutine mpfftcr (is, m, n, nsq, x, y)

!   This performs an N-point complex-to-real FFT, where N = 2^M.  X is the
!   double complex input array, and Y is the double precision output array.
!   The array X is used as a scratch array in MPFFT1, and so is overwritten.
!   X must be dimensioned with N/2+N1*NSP1+1 DC cells, and Y with N DP cells, 
!   where N = 2^M and N1 = 2^int(M/2).   This dimension requirement for X is 
!   somewhat greater than shown in the dimension statement below, because 
!   MPFFT1, which is called by this routine, requires more.  IS is the sign of
!   the transform.  Before calling MPFFTCR, the UU1 and UU2 arrays must be 
!   initialized by calling MPINIX.  This routine is not intended to be called 
!   directly by the user.

implicit double precision (a-h, o-z)
dimension y(n)
complex (mpkdp) dc1(n/2), x(n/2+nsq*mpnsp1+1), a1, a2, x1, x2

mx = mpuu1(1)

!   Check if input parameters are invalid.

if ((is .ne. 1 .and. is .ne. -1) .or. m .lt. 3 .or. m .gt. mx) then
  if (mpker(27) .ne. 0) then
    write (mpldb, 1)  is, m, mx
1   format ('*** MPFFTCR: Either the UU arrays have not been initialized'/ &
    'or else one of the input parameters is invalid', 3i5)
    mpier = 27
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif
n1 = 2 ** (m / 2)
n2 = n / 2
n21 = n2 + 1
n4 = n / 4

!   Construct the input to MPFFT1.

dc1(1) = 0.5d0 * cmplx (dble (x(1) + x(n2+1)), dble (x(1) - x(n2+1)), mpkdp)
if (is .eq. 1) then
  dc1(n4+1) = conjg (x(n4+1))
else
  dc1(n4+1) = x(n4+1)
endif
ku = n2

if (is .eq. 1) then
!dir$ ivdep
  do k = 2, n4
    x1 = x(k)
    x2 = conjg (x(n2+2-k))
    a1 = x1 + x2
    a2 = (0.d0, 1.d0) * mpuu1(k+ku) * (x1 - x2)
    dc1(k) = 0.5d0 * (a1 + a2)
    dc1(n2+2-k) = 0.5d0 * conjg (a1 - a2)
  enddo
else
!dir$ ivdep
  do k = 2, n4
    x1 = x(k)
    x2 = conjg (x(n2+2-k))
    a1 = x1 + x2
    a2 = (0.d0, 1.d0) * conjg (mpuu1(k+ku)) * (x1 - x2)
    dc1(k) = 0.5d0 * (a1 + a2)
    dc1(n2+2-k) = 0.5d0 * conjg (a1 - a2)
  enddo
endif

!   Perform a normal N/2-point FFT on DC1.

call mpfft1 (is, m - 1, n1, n2 / n1, dc1, x)

!   Copy DC1 to Y such that DC1(k) = Y(2k-1) + i Y(2k).

do k = 1, n / 2
  y(2*k-1) = dble (dc1(k))
  y(2*k) = aimag (dc1(k))
enddo

return
end subroutine

subroutine mpfftrc (is, m, n, nsq, x, y)

!   This performs an N-point real-to-complex FFT, where N = 2^M.  X is the
!   double precision input array, and Y is the double complex output array.
!   X must be dimensioned with N DP cells, and Y with N/2+N1*NSP1+1 DC cells,
!   where N = 2^M and N1 = 2^int(M/2).  This dimension requirement for Y is 
!   somewhat greater than that shown in the dimension statement below, because
!   MPFFT1, which is called by this routine, requires more.  IS is the sign of
!   the transform.  Before calling MPFFTRC, the UU1 and UU2 arrays must be 
!   initialized by calling MPINIX.  This routine is not intended to be called 
!   directly by the user.

implicit double precision (a-h, o-z)
dimension x(n)
complex (mpkdp) dc1(n/2), y(n/2+nsq*mpnsp1+1), a1, a2, z1, z2

mx = mpuu1(1)

!   Check if input parameters are invalid.

if ((is .ne. 1 .and. is .ne. -1) .or. m .lt. 3 .or. m .gt. mx) then
  if (mpker(67) .ne. 0) then
    write (mpldb, 1)  is, m, mx
1   format ('*** MPFFTRC: either the UU arrays have not been initialized'/ &
    'or else one of the input parameters is invalid',3i5)
    mpier = 67
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif
n1 = 2 ** (m / 2)
n2 = n / 2
n21 = n2 + 1
n4 = n / 4

!   Copy X to DC1 such that DC1(k) = X(2k-1) + i X(2k).

do k = 1, n2
  dc1(k) = cmplx (x(2*k-1), x(2*k), mpkdp)
enddo

!   Perform a normal N/2-point FFT on DC1.

call mpfft1 (is, m - 1, n1, n2 / n1, dc1, y)

!   Reconstruct the FFT of X.

y(1) = cmplx (2.d0 * (dble (dc1(1)) + aimag (dc1(1))), 0.d0, mpkdp)
if (is .eq. 1) then
  y(n4+1) = 2.d0 * dc1(n4+1)
else
  y(n4+1) = 2.d0 * conjg (dc1(n4+1))
endif
y(n2+1) = cmplx (2.d0 * (dble (dc1(1)) - aimag (dc1(1))), 0.d0, mpkdp)
ku = n2

if (is .eq. 1) then
!dir$ ivdep
  do k = 2, n4
    z1 = dc1(k)
    z2 = conjg (dc1(n2+2-k))
    a1 = z1 + z2
    a2 = (0.d0, -1.d0) * mpuu1(k+ku) * (z1 - z2)
    y(k) = a1 + a2
    y(n2+2-k) = conjg (a1 - a2)
  enddo
else
!dir$ ivdep
  do k = 2, n4
    z1 = dc1(k)
    z2 = conjg (dc1(n2+2-k))
    a1 = z1 + z2
    a2 = (0.d0, -1.d0) * conjg (mpuu1(k+ku)) * (z1 - z2)
    y(k) = a1 + a2
    y(n2+2-k) = conjg (a1 - a2)
  enddo
endif

return
end subroutine

subroutine mpfft1 (is, m, n1, n2, x, y)

!   This routine performs a complex-to-complex FFT.  IS is the sign of the
!   transform, N = 2^M is the size of the transform.  N1 = 2^M1 and N2 = 2^M2,
!   where M1 and M2 are defined as below.  X is the input and output array,
!   and Y is a scratch array.  X must have at N, and Y at least N + N1*MPNSP1,
!   double complex cells.  The arrays MPUU1 and MPUU2 must have been 
!   initialized by calling MPINIX.  This routine is not intended to be called 
!   directly by the user.

!   This employs the two-pass variant of the "four-step" FFT.  See the
!   article by David H. Bailey in J. of Supercomputing, March 1990, p. 23-35.

complex (mpkdp) x(n1,n2), y(n2+mpnsp1,n1), q1(mpnsp2), z1(mpnrow+mpnsp1,n1), &
  q2(mpnsp2), z2(mpnrow+mpnsp1,n1)

n = 2 ** m
m1 = (m + 1) / 2
m2 = m - m1
nr1 = min (n1, mpnrow)
nr2 = min (n2, mpnrow)
ku = mpuu2(m)

do i = 0, n1 - 1, nr1

!   Copy NR1 rows of X (treated as a N1 x N2 complex array) into Z1.

  do j = 1, n2
    do k = 1, nr1
      z1(k,j) = x(i+k,j)
    enddo
  enddo

!   Perform NR1 FFTs, each of length N2.

  call mpfft2 (is, nr1, m2, n2, z1, z2)

!   Multiply the resulting NR1 x N2 complex block by roots of unity and
!   store transposed into the appropriate section of Y.

  iu = i + ku - n1 - 1
  if (is .eq. 1) then
    do j = 1, n2
      do k = 1, nr1
        y(j,i+k) = mpuu2(iu+k+j*n1) * z1(k,j)
      enddo
    enddo
  else
    do j = 1, n2
      do k = 1, nr1
        y(j,i+k) = conjg (mpuu2(iu+k+j*n1)) * z1(k,j)
      enddo
    enddo
  endif
enddo

do i = 0, n2 - 1, nr2

!   Copy NR2 rows of the Y array into Z2.

  do j = 1, n1
    do k = 1, nr2
      z2(k,j) = y(i+k,j)
    enddo
  enddo

!   Perform NR2 FFTs, each of length N1.

  call mpfft2 (is, nr2, m1, n1, z2, z1)

!   Copy NR2 x N1 complex block back into X array.  It's a little more
!   complicated if M is odd.

  if (mod (m, 2) .eq. 0) then
    do j = 1, n1
      do k = 1, nr2
        x(i+k,j) = z2(k,j)
      enddo
    enddo
  else
    do j = 1, n1 / 2
      j2 = 2 * j - 1
!dir$ ivdep
      do k = 1, nr2
        x(i+k,j) = z2(k,j2)
        x(i+k+n2,j) = z2(k,j2+1)
      enddo
    enddo
  endif
enddo

return
end subroutine

subroutine mpfft2 (is, ns, m, n, x, y)

!   This performs NS simultaneous N-point complex-to-complex FFTs, where
!   N = 2^M.  X is the input and output array, UU1 is the root of unity array,
!   and Y is a scratch array.  X, Y and UU1 are double complex.   This routine
!   is not intended to be called directly by the user.

complex (mpkdp) x(mpnrow+mpnsp1,n), y(mpnrow+mpnsp1,n)

!   Perform the second variant of the Stockham FFT.

do l = 1, m, 2
  call mpfft3 (is, l, ns, m, n, x, y)
  if (l .eq. m) goto 100
  call mpfft3 (is, l + 1, ns, m, n, y, x)
enddo

goto 110

!   Copy Y to X.

100  do j = 1, n
  do i = 1, ns
    x(i,j) = y(i,j)
  enddo
enddo

110  continue

return
end subroutine

subroutine mpfft3 (is, l, ns, m, n, x, y)

!   This performs the L-th iteration of the second variant of the Stockham FFT
!   on the NS vectors in X.  Y is a scratch array, and UU1 is the root of
!   unity array.  X, Y and UU1 are double complex.  This routine is not
!   intended to be called directly by the user.

complex (mpkdp) x(mpnrow+mpnsp1,n), y(mpnrow+mpnsp1,n), u1, x1, x2

!   Set initial parameters.

n1 = n / 2
lk = 2 ** (l - 1)
li = 2 ** (m - l)
lj = 2 * lk
ku = li + 1

do i = 0, li - 1
  i11 = i * lk + 1
  i12 = i11 + n1
  i21 = i * lj + 1
  i22 = i21 + lk
  if (is .eq. 1) then
    u1 = mpuu1(i+ku)
  else
    u1 = conjg (mpuu1(i+ku))
  endif

  do k = 0, lk - 1
!dir$ ivdep
    do j = 1, ns
      x1 = x(j,i11+k)
      x2 = x(j,i12+k)
      y(j,i21+k) = x1 + x2
      y(j,i22+k) = u1 * (x1 - x2)
    enddo
  enddo
enddo

return
end subroutine

recursive subroutine mplconv (iq, n, nsq, a, b, c)

!  This computes the linear convolution of the N-long inputs A and B.  |IQ| is
!  the number of arguments (i.e., if IQ = 1, then B is ignored).  If IQ is
!  negative (and N < 64) then only the second half of the result vector is
!  required (i.e. this is a call by itself -- see below). NSQ = int(sqrt(3*N))
!  is an input required for the dimension of DC1 and DC2 (see below). 

!  This routine employs an advanced FFT-based scheme, except for small n.
!  This routine is not intended to be called directly by the user.
!>
!   Two machine-dependent parameters are set in this routine:
!     ERM = Maximum tolerated FFT roundoff error.  On IEEE systems ERM =
!     0.438D0.  It is not necessary to specify ERM for modest levels of
!     precision -- see comments below.
!     MBT = Number of mantissa bits in double precision data.  MBT = 53 on
!     IEEE systems, and MBT = 48 (i.e. single precision) on Crays.
!     It is not necessary to specify MBT for modest levels of precision.

double precision an, cl2, erm, t1, t2
double precision a(n), b(n), c(2*n), d1(3*n+2), d2(3*n+2), d3(3*n+2)
complex (mpkdp) dc1(3*n/2+nsq*mpnsp1+3), dc2(3*n/2+nsq*mpnsp1+3)
parameter (cl2 = 1.4426950408889633d0, erm = 0.438d0, mbt = 53)

!   Handle the case where N is less than NCR1 = 2 ** (MPMCR-1).  If IQ < 0, 
!   only the second half of the result vector is returned, since the first 
!   half won't be used.

ncr1 = 2 ** (mpmcr - 1)
if (n .lt. ncr1) then
  if (iq .eq. 1) then
    do k = 1, 2 * n
      t1 = 0.d0
      n1 = max (k - n + 1, 1)
      n2 = min (k, n)

      do j = n1, n2
        t1 = t1 + a(j) * a(k-j+1)
      enddo

      c(k) = t1
    enddo
  elseif (iq .eq. 2) then
    do k = 1, 2 * n
      t1 = 0.d0
      n1 = max (k - n + 1, 1)
      n2 = min (k, n)

      do j = n1, n2
        t1 = t1 + a(j) * b(k-j+1)
      enddo

      c(k) = t1
    enddo
  elseif (iq .eq. -1) then
    do k = 1, n - 1
      c(k) = 0.d0
    enddo

    do k = n, 2 * n
      t1 = 0.d0
      n1 = k - n + 1
      n2 = n

      do j = n1, n2
        t1 = t1 + a(j) * a(k-j+1)
      enddo

      c(k) = t1
    enddo
  elseif (iq .eq. -2) then
    do k = 1, n - 1
      c(k) = 0.d0
    enddo

    do k = n, 2 * n
      t1 = 0.d0
      n1 = k - n + 1
      n2 = n

      do j = n1, n2
        t1 = t1 + a(j) * b(k-j+1)
      enddo

      c(k) = t1
    enddo
  endif
  goto 100
endif

!   Determine M1 and N1.  Note that by this reckoning, N1 <= 1.5 N.  This is
!   the reason for the 3*n/2 dimensions above.

t1 = 0.75d0 * n
m1 = cl2 * log (t1) + 1.d0 - mprxx
n1 = 2 ** m1
m2 = m1 + 1
n2 = 2 * n1
n4 = 2 * n2
nm = min (2 * n, n2)

if (abs (iq) .eq. 1) then
  do i = 1, n
    d1(i) = a(i)
  enddo

  do i = n + 1, n2
    d1(i) = 0.d0
  enddo

!   Perform a forward real-to-complex FFT on the vector in a.

  call mpfftrc (1, m2, n2, nsq, d1, dc1)

!   Square the resulting complex vector.

  do i = 1, n1 + 1
    dc1(i) = dc1(i) ** 2
  enddo

else
  do i = 1, n
    d1(i) = a(i)
    d2(i) = b(i)
  enddo

  do i = n + 1, n2
    d1(i) = 0.d0
    d2(i) = 0.d0
  enddo

!   Perform forward real-to-complex FFTs on the vectors in a and b.

  call mpfftrc (1, m2, n2, nsq, d1, dc1)
  call mpfftrc (1, m2, n2, nsq, d2, dc2)

!   Multiply the resulting complex vectors.

  do i = 1, n1 + 1
    dc1(i) = dc1(i) * dc2(i)
  enddo
endif

!   Perform an inverse complex-to-real FFT on the resulting data.

call mpfftcr (-1, m2, n2, nsq, dc1, d3)

!   Divide by N4.

an = 1.d0 / n4

do i = 1, nm
  t1 = an * d3(i)
  t2 = anint (t1)
!  D1(I) = ABS (T2 - T1)
  c(i) = t2
enddo

!   Find the largest FFT roundoff error.  Roundoff error is minimal unless
!   exceedingly high precision (i.e. over one million digits) is used.  Thus
!   this test may be disabled in normal use.  To disable this test, uncomment
!   the next line of code and comment out the line of the previous loop
!   that begins D1(I) =.  To enable, do the reverse.

!   This code can be used as a rigorous system integrity test.  First set
!   MBT according to the system being used, and then set ERM to be fairly
!   small, say 0.001 or whatever is somewhat larger than the largest FFT
!   roundoff error typically encountered for a given precision level on the
!   computer being used.  Enable this test as explained in the previous
!   paragraph.  Then if an anomalously large roundoff error is detected, a
!   hardware or compiler error has likely occurred.

goto 180
t1 = 0.d0

do i = 1, nm
  if (d1(i) .gt. t1) then
    i1 = i
    t1 = d1(i)
  endif
enddo

!   Check if maximum roundoff error exceeds the limit ERM, which is set above.
!   Also determine the number of fractional bits and how large the error is in
!   terms of units in the last place (ulp).

if (t1 .gt. erm)  then
  if (mpker(55) .ne. 0) then
    t2 = an * d1(i1)
    i2 = cl2 * log (t1) + 1.d0 + mprxx
    i3 = cl2 * log (t2) + 1.d0 + mprxx
    i4 = mbt + i2 - i3
    i5 = t1 * 2 ** i4 + mprxx
    write (mpldb, 2) i1, t1, i4, i5
2   format ('*** MPLCONV: Excessive FFT roundoff error',i10,f10.6,2i6)
    mpier = 55
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
endif

180 continue

!   Handle case where n > n1.

if (n .gt. n1) then
  m = n - n1
  m2 = 2 * m
  m21 = 2 * m - 1
  ms = sqrt (3.d0 * m21) + mprxx
  k = n1 - m + 1

  if (abs (iq) .eq. 1) then
    do i = 1, m21
      d1(i) = a(k+i)
    enddo

    call mplconv (-1, m21, ms, d1, d2, d3)
  else
    do i = 1, m21
      d1(i) = a(k+i)
      d2(i) = b(k+i)
    enddo

    call mplconv (-2, m21, ms, d1, d2, d3)
  endif

  do i = 1, m2
    ii = i + m2 - 2
    c(i) = c(i) - d3(ii)
    c(i+n2) = d3(ii)
  enddo
endif

100 continue

return
end subroutine

end module

module mpfunh

!   This module defines the extra-high precision multiply and divide routines.

use mpfuna
use mpfunc
use mpfung
contains

subroutine mpcbrx (a, b)

!   This computes the cube root of the MP number A and returns the MP result
!   in B.  Before calling MPCBRX, the arrays UU1 and UU2 must be initialized by
!   calling MPINIX.  For modest levels of precision, use MPCBRT.  Debug output
!   starts with MPIDB = 6.

!   Max SP space for B: MPNW + 4 cells.

!   This routine uses basically the same Newton iteration algorithm as MPCBRT.
!   In fact, this routine calls MPCBRT to obtain an initial approximation.
!   See the comment about the parameter NIT in MPDIVX.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0, nit = 3)
dimension a(mpnw+2), b(mpnw+4), f(8), s(3*mpnw+15)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 6) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPCBRX I'/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ncr = 2 ** mpmcr

if (na .eq. 0) then
  b(1) = 0.
  b(2) = 0.
  goto 120
endif
if (ia .lt. 0.d0) then
  if (mpker(14) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPCBRX: Argument is negative.')
    mpier = 14
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if precision level is too low to justify the advanced routine.

if (mpnw .le. ncr) then
  call mpcbrt (a, b)
  goto 120
endif
n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
nws = mpnw

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Compute A^2 outside of the iteration loop.

call mpsqx (a, s(k0))

!   Compute the initial approximation of A ^ (-2/3).

mpnw = ncr + 1
call mpcbrt (a, s(k1))
call mpdiv (s(k1), a, b)
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (powers of two).

do k = mpmcr + 1, mq - 1
  nw1 = mpnw
  mpnw = min (2 * mpnw - 2, nws) + 1
  nw2 = mpnw
100  continue
  call mpsqx (b, s(k1))
  call mpmulx (b, s(k1), s(k2))
  call mpmulx (s(k0), s(k2), s(k1))
  call mpsub (f, s(k1), s(k2))
  mpnw = nw1
  call mpmulx (b, s(k2), s(k1))
  call mpdivd (s(k1), 3.d0, 0, s(k2))
  mpnw = nw2
  call mpadd (b, s(k2), s(k1))
  call mpeq (s(k1), b)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
 enddo

!   Perform last iteration using Karp's trick.

call mpmulx (a, b, s(k0))
nw1 = mpnw
mpnw = min (2 * mpnw - 2, nws) + 1
nw2 = mpnw
call mpsqx (s(k0), s(k1))
call mpmulx (s(k0), s(k1), s(k2))
call mpsub (a, s(k2), s(k1))
mpnw = nw1
call mpmulx (s(k1), b, s(k2))
call mpdivd (s(k2), 3.d0, 0, s(k1))
mpnw = nw2
call mpadd (s(k0), s(k1), s(k2))
call mpeq (s(k2), b)

!   Restore original precision level.

mpnw = nws
call mproun (b)

120  if (mpidb .ge. 6) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPCBRX O'/(6f12.0))
endif
return
end subroutine

subroutine mpdivx (a, b, c)

!   This divides the MP number A by the MP number B and returns the MP result
!   in C.  Before calling MPDIVX, the arrays UU1 and UU2 must be initialized by
!   calling MPINIX.  For modest levels of precision, use MPDIV.  Debug output 
!   starts with MPIDB = 7.

!   Max SP space for C: MPNW + 4 cells.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to 1 / B:

!    X_{k+1} = X_k + (1 - X_k * B) * X_k

!   where the muliplication () * X_k is performed with only half of the
!   normal level of precision.  These iterations are performed with a
!   maximum precision level MPNW that is dynamically changed, doubling with
!   each iteration.  The final iteration is performed as follows (this is
!   due to A. Karp):

!    A / B = (A * X_n) + [A - (A * X_n) * B] * X_n  (approx.)

!   where the multiplications A * X_n and [] * X_n are performed with only
!   half of the final level of precision.

!   One difficulty with this procedure is that errors often accumulate in the
!   trailing mantissa words.  This error can be controlled by repeating one of
!   the iterations.  The iteration that is repeated is controlled by setting
!   the parameter NIT below:  If NIT = 0, the last iteration is repeated (this
!   is most effective but most expensive).  If NIT = 1, then the next-to-last
!   iteration is repeated, etc.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0, nit = 3)
dimension a(mpnw+2), b(mpnw+2), c(mpnw+4), f(8), s(3*mpnw+15)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 7) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPDIVX I'/(6f12.0))
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 1) (b(i), i = 1, no)
endif

ia = sign (1., a(1))
ib = sign (1., b(1))
na = min (int (abs (a(1))), mpnw)
nb = min (int (abs (b(1))), mpnw)
ncr = 2 ** mpmcr

!   Check if dividend is zero.

if (na .eq. 0) then
  c(1) = 0.
  c(2) = 0.
  goto 120
endif

!   Check if divisor is zero.

if (nb .eq. 0)  then
  if (mpker(33) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPDIVX: Divisor is zero.')
    mpier = 33
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if precision level of divisor is too low to justify the advanced
!   routine.

if (nb .le. ncr) then
  call mpdiv (a, b, c)
  goto 120
endif
n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
nws = mpnw

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Compute the initial approximation of 1 / B to a precision of NCR words.

mpnw = ncr + 1
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.
call mpdiv (f, b, c)
iq = 0

!   Perform the Newton-Raphson iterations described above.

do k = mpmcr + 1, mq - 1
  nw1 = mpnw
  mpnw = min (2 * mpnw - 2, nws) + 1
  nw2 = mpnw
100  continue
  call mpmulx (b, c, s(k0))
  call mpsub (f, s(k0), s(k1))
  mpnw = nw1
  call mpmulx (c, s(k1), s(k0))
  mpnw = nw2
  call mpadd (c, s(k0), s(k1))
  call mpeq (s(k1), c)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
 enddo

!   Perform last iteration using Karp's trick.

call mpmulx (a, c, s(k0))
nw1 = mpnw
mpnw = min (2 * mpnw - 2, nws) + 1
nw2 = mpnw
call mpmulx (s(k0), b, s(k1))
call mpsub (a, s(k1), s(k2))
mpnw = nw1
call mpmulx (s(k2), c, s(k1))
mpnw = nw2
call mpadd (s(k0), s(k1), s(k2))
call mpeq (s(k2), c)

!   Restore original precision level.

mpnw = nws
call mproun (c)

120  if (mpidb .ge. 7) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 3) (c(i), i = 1, no)
3 format ('MPDIVX O'/(6f12.0))
endif
return
end subroutine

subroutine mpmulx (a, b, c)

!   This routine multiplies MP numbers A and B to yield the MP product C.
!   Before calling MPMULX, the arrays UU1 and UU2 must be initialized by
!   calling MPINIX.  For modest levels of precision, use MPMUL.  Debug output 
!   starts with MPIDB = 8.

!   Max SP space for C: MPNW + 4 cells.

!   This routine returns up to MPNW mantissa words of the product.  If the
!   complete double-long product of A and B is desired (for example in large
!   integer applications), then MPNW must be at least as large as the sum of
!   the mantissa lengths of A and B.  In other words, if the precision levels
!   of A and B are both 256 words, then MPNW must be at least 512 words to
!   obtain the complete double-long product in C.

!   This subroutine uses an advanced technique involving the fast Fourier
!   transform (FFT).  For high precision it is significantly faster than the
!   conventional scheme used in MPMUL.

double precision t1, t2, t3, t4
dimension a(mpnw+2), b(mpnw+2), c(mpnw+4)
double precision d1(2*mpnw+4), d2(2*mpnw+4), d3(4*mpnw+8)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  return
endif
if (mpidb .ge. 8)  then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPMULX I'/(6f12.0))
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 1) (b(i), i = 1, no)
endif

ia = sign (1., a(1))
ib = sign (1., b(1))
na = min (int (abs (a(1))), mpnw)
nb = min (int (abs (b(1))), mpnw)
ncr = 2 ** mpmcr
if (na .eq. 0 .or. nb .eq. 0) then

!   One of the inputs is zero -- result is zero.

  c(1) = 0.
  c(2) = 0.
  goto 190
endif

!   Check if precision level of one of the arguments is too low to justify the
!   advanced routine.

if (na .le. ncr .or. nb .le. ncr) then
  call mpmul (a, b, c)
  goto 190
endif

!   Place the input data in A and B into the scratch arrays DD1 and DD2.
!   This code also splits the input data into half-sized words.

!dir$ ivdep
do i = 1, na
  i2 = 2 * i - 1
  t1 = a(i+2)
  t2 = int (mprbx * t1)
  d1(i2) = t2
  d1(i2+1) = t1 - mpbbx * t2
enddo

do i = 2 * na + 1, 2 * nb
  d1(i) = 0.d0
enddo
    
!dir$ ivdep
do i = 1, nb
  i2 = 2 * i - 1
  t1 = b(i+2)
  t2 = int (mprbx * t1)
  d2(i2) = t2
  d2(i2+1) = t1 - mpbbx * t2
enddo

do i = 2 * nb + 1, 2 * na
  d2(i) = 0.d0
enddo

nn = 2 * max (na, nb)
nx = sqrt (3.d0 * nn) + mprxx
call mplconv (2, nn, nx, d1, d2, d3)

!   Recombine words and release carries.

nc = min (na + nb, mpnw)
nc1 = min (mpnw + 1, na + nb - 1)
d1(1) = sign (nc, ia * ib)
d1(2) = a(2) + b(2) + 1
d1(3) = d3(1)
d1(nc+3) = 0.d0
d1(nc+4) = 0.d0

!dir$ ivdep
do i = 1, nc1
  i2 = 2 * i
  t3 = d3(i2)
  t4 = d3(i2+1)
  t1 = int (mprdx * t3)
  t2 = t3 - mpbdx * t1
  t3 = int (mprdx * t4)
  t4 = t4 - mpbdx * t3
  d1(i+3) = mpbbx * t2 + t4
  d1(i+2) = d1(i+2) + mpbbx * t1 + t3
enddo

!   Fix up the result.

call mpnorm (d1, c)

190  if (mpidb .ge. 8) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 3) (c(i), i = 1, no)
3 format ('MPMULX O'/(6f12.0))
endif
return
end subroutine

subroutine mpnpwx (a, n, b)

!   This computes the N-th power of the MP number A and returns the MP result
!   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
!   of A ^ |N| is returned.  Before calling MPNPWX, the arrays UU1 and UU2
!   must be initialized by calling MPINIX.  For modest levels of precision, use
!   MPNPWR.  Debug output starts with MPIDB = 6.

!   Max SP space for B: MPNW + 4 cells.

!   This routine employs the binary method for exponentiation.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0)
dimension a(mpnw+2), b(mpnw+4), f1(8), s(2*mpnw+8)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 6) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) n, (a(i), i = 1, no)
1 format ('MPNPWX I',i5/(6f12.0))
endif

ncr = 2 ** mpmcr
na = min (int (abs (a(1))), mpnw)

!   Check if precision level of A is too low to justify the advanced routine.

if (na .le. ncr .and. n .ge. 0 .and. n .le. 4) then
  call mpnpwr (a, n, b)
  goto 120
endif
if (na .eq. 0) then
  if (n .ge. 0) then
    b(1) = 0.
    b(2) = 0.
    goto 120
  else
    if (mpker(58) .ne. 0) then
      write (mpldb, 2)
2     format ('*** MPNPWX: argument is zero and N is negative or zero.')
      mpier = 58
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  endif
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
nn = abs (n)
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
if (nn .eq. 0) then
  call mpeq (f1, b)
  goto 120
elseif (nn .eq. 1) then
  call mpeq (a, b)
  goto 110
elseif (nn .eq. 2) then
  call mpsqx (a, b)
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN .GT. NN.

t1 = nn
mn = cl2 * log (t1) + 1.d0 + mprxx
call mpeq (f1, b)
call mpeq (a, s(k0))
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn .ne. 2 * kk) then
    call mpmulx (b, s(k0), s(k1))
    call mpeq (s(k1), b)
  endif
  kn = kk
  if (j .lt. mn) then
    call mpsqx (s(k0), s(k1))
    call mpeq (s(k1), s(k0))
  endif
enddo

!   Compute reciprocal if N is negative.

110  if (n .lt. 0) then
  call mpdivx (f1, b, s(k0))
  call mpeq (s(k0), b)
endif

120  if (mpidb .ge. 6) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPNPWX O'/(6f12.0))
endif
return
end subroutine

subroutine mpnrtx (a, n, b)

!   This computes the N-th root of the MP number A and returns the MP result
!   in B.  N must be at least one and must not exceed 2 ^ 30.  Before calling
!   MPNRTX, the arrays UU1 and UU2 must be initialized by calling MPINIX.  For
!   modest levels of precision, use MPNRT.  Debug output starts with MPIDB = 6.

!   Max SP space for B: MPNW + 4 cells.

!   This routine uses basically the same Newton iteration algorithm as MPNRT.
!   In fact, this routine calls MPNRT to obtain an initial approximation.
!   See the comment about the parameter NIT in MPDIVX.

double precision cl2, t1, t2, tn
parameter (cl2 = 1.4426950408889633d0, nit = 3, n30 = 2 ** 30)
dimension a(mpnw+2), b(mpnw+4), f1(8), f2(8), s(4*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 6) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) n, (a(i), i = 1, no)
1 format ('MPNRTX I',i5/(6f12.0))
endif

ncr = 2 ** mpmcr
ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)

if (na .eq. 0) then
  b(1) = 0.
  b(2) = 0.
  goto 140
endif
if (ia .lt. 0) then
  if (mpker(61) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPNRTX: Argument is negative.')
    mpier = 61
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif
if (n .le. 0 .or. n .gt. n30) then
  if (mpker(62) .ne. 0) then
    write (mpldb, 3) n
3   format ('*** MPNRTX: Improper value of N',i10)
    mpier = 62
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if precision level is too low to justify the advanced routine.

if (mpnw .le. ncr) then
  call mpnrt (a, n, b)
  goto 140
endif

!   If N = 1, 2 or 3, call MPEQ, MPSQRX or MPCBRX.  These are faster.

if (n .eq. 1) then
  call mpeq (a, b)
  goto 140
elseif (n .eq. 2) then
  call mpsqrx (a, b)
  goto 140
elseif (n .eq. 3) then
  call mpcbrx (a, b)
  goto 140
endif

n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
k3 = k2 + n5
nws = mpnw
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Check how close A is to 1.

call mpsub (a, f1, s(k0))
if (s(k0) .eq. 0.) then
  call mpeq (f1, b)
  goto 140
endif
call mpmdc (s(k0), t1, n1)
n2 = cl2 * log (abs (t1))
t1 = t1 * 0.5d0 ** n2
n1 = n1 + n2
if (n1 .le. -30) then
  t2 = n
  n2 = cl2 * log (t2) + 1.d0 + mprxx
  n3 = - mpnbt * mpnw / n1
  if (n3 .lt. 1.25 * n2) then

!   A is so close to 1 that it is cheaper to use the binomial series.

    call mpdivd (s(k0), t2, 0, s(k1))
    call mpadd (f1, s(k1), s(k2))
    k = 0

100 k = k + 1
    t1 = 1 - k * n
    t2 = (k + 1) * n
    call mpmuld (s(k1), t1, 0, s(k3))
    call mpdivd (s(k3), t2, 0, s(k1))
    call mpmulx (s(k0), s(k1), s(k3))
    call mpeq (s(k3), s(k1))
    call mpadd (s(k1), s(k2), s(k3))
    call mpeq (s(k3), s(k2))
    if (s(k1) .ne. 0. .and. s(k1+1) .ge. - mpnw) goto 100

    call mpeq (s(k2), b)
    goto 130
  endif
endif

!   Compute the initial approximation of A ^ (-1/N).

mpnw = ncr + 1
call mpnrt (a, n, s(k0))
call mpdiv (f1, s(k0), b)
tn = n
call mpdmc (tn, 0, f2)
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (powers of two).

do k = mpmcr + 1, mq
  nw1 = mpnw
  mpnw = min (2 * mpnw - 1, nws) + 1
  nw2 = mpnw
110  continue
  call mpnpwx (b, n, s(k0))
  call mpmulx (a, s(k0), s(k1))
  call mpsub (f1, s(k1), s(k0))
  mpnw = nw1
  call mpmulx (b, s(k0), s(k1))
  call mpdivd (s(k1), tn, 0, s(k0))
  mpnw = nw2
  call mpadd (b, s(k0), s(k1))
  call mpeq (s(k1), b)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 110
  endif
enddo

!   Take the reciprocal to give final result.

call mpdivx (f1, b, s(k0))
call mpeq (s(k0), b)

!   Restore original precision level.

130  mpnw = nws
call mproun (b)

140  if (mpidb .ge. 6) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 4) (b(i), i = 1, no)
4 format ('MPNRTX O'/(6f12.0))
endif
return
end subroutine

subroutine mpsqrx (a, b)

!   This computes the square root of the MP number A and returns the MP result
!   in B.  Before calling MPSQRX, the arrays UU1 and UU2 must be initialized by
!   calling MPINIX.  For modest levels of precision, use MPSQRT.  Debug output
!   starts with MPIDB = 6.

!   Max SP space for B: MPNW + 4 cells.

!   This routine uses basically the same Newton iteration algorithm as MPSQRT.
!   In fact, this routine calls MPSQRT to obtain an initial approximation.
!   See the comment about the parameter NIT in MPDIVX.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0, nit = 3)
dimension a(mpnw+2), b(mpnw+4), f(8), s(3*mpnw+15)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 6) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPSQRX I'/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ncr = 2 ** mpmcr

if (na .eq. 0) then
  b(1) = 0.
  b(2) = 0.
  goto 120
endif
if (ia .lt. 0.d0) then
  if (mpker(71) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPSQRX: Argument is negative.')
    mpier = 71
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if precision level is too low to justify the advanced routine.

if (mpnw .le. ncr) then
  call mpsqrt (a, b)
  goto 120
endif
n5 = mpnw + 5
k0 = 1
k1 = k0 + n5
k2 = k1 + n5
nws = mpnw

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Compute the initial approximation of 1 / Sqrt(A).

mpnw = ncr + 1
call mpsqrt (a, s(k0))
call mpdiv (s(k0), a, b)
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = mpmcr + 1, mq - 1
  nw1 = mpnw
  mpnw = min (2 * mpnw - 2, nws) + 1
  nw2 = mpnw
100  continue
  call mpsqx (b, s(k0))
  call mpmulx (a, s(k0), s(k1))
  call mpsub (f, s(k1), s(k0))
  mpnw = nw1
  call mpmulx (b, s(k0), s(k1))
  call mpmuld (s(k1), 0.5d0, 0, s(k0))
  mpnw = nw2
  call mpadd (b, s(k0), s(k1))
  call mpeq (s(k1), b)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
enddo

!   Perform last iteration using Karp's trick.

call mpmulx (a, b, s(k0))
nw1 = mpnw
mpnw = min (2 * mpnw - 2, nws) + 1
nw2 = mpnw
call mpsqx (s(k0), s(k1))
call mpsub (a, s(k1), s(k2))
mpnw = nw1
call mpmulx (s(k2), b, s(k1))
call mpmuld (s(k1), 0.5d0, 0, s(k2))
mpnw = nw2
call mpadd (s(k0), s(k2), s(k1))
call mpeq (s(k1), b)

!   Restore original precision level.

mpnw = nws
call mproun (b)

120  if (mpidb .ge. 6) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPSQRX O'/(6f12.0))
endif

return
end subroutine

subroutine mpsqx (a, b)

!   This routine squares the MP number A to yield the MP product B.
!   Before calling MPSQX, the arrays UU1 and UU2 must be initialized by calling
!   MPINIX.  For modest levels of precision, use MPMUL.  MPNW should be a power
!   of two.  Debug output starts with MPIDB = 8.

!   Max SP space for B: MPNW + 4 cells.

!   This subroutine uses the same FFT technique as MPMULX.  It is faster
!   because only one forward FFT has to be computed.  See the comments in
!   MPMULX about obtaining the complete double-long result.

double precision t1, t2, t3, t4
dimension a(mpnw+2), b(mpnw+4)
double precision d1(2*mpnw+4), d2(4*mpnw+8)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 8)  then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPSQX I'/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ncr = 2 ** mpmcr

if (na .eq. 0) then
  b(1) = 0.
  b(2) = 0.
  goto 190
endif

!   Check if precision level of the argument is too low to justify the
!   advanced routine.

if (na .le. ncr) then
  call mpmul (a, a, b)
  goto 190
endif

!   Place the input data in A into the scratch array DD1.
!   This code also splits the input data into half-sized words.

!dir$ ivdep
do i = 1, na
  i2 = 2 * i - 1
  t1 = a(i+2)
  t2 = int (mprbx * t1)
  d1(i2) = t2
  d1(i2+1) = t1 - mpbbx * t2
enddo

nn = 2 * na
nx = sqrt (3.d0 * nn) + mprxx
call mplconv (1, nn, nx, d1, d1, d2)

!   Recombine words and release carries.

nc = min (2 * na, mpnw)
nc1 = min (mpnw + 1, 2 * na - 1)
d1(1) = nc
d1(2) = 2 * a(2) + 1
d1(3) = d2(1)
d1(nc+3) = 0.d0
d1(nc+4) = 0.d0

!dir$ ivdep
do i = 1, nc1
  i2 = 2 * i
  t3 = d2(i2)
  t4 = d2(i2+1)
  t1 = int (mprdx * t3)
  t2 = t3 - mpbdx * t1
  t3 = int (mprdx * t4)
  t4 = t4 - mpbdx * t3
  d1(i+3) = mpbbx * t2 + t4
  d1(i+2) = d1(i+2) + mpbbx * t1 + t3
enddo

!   Fix up the result.

call mpnorm (d1, b)

190  if (mpidb .ge. 8) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPSQX O'/(6f12.0))
endif
return
end subroutine

end module

module mpfuni

!   This module defines extra-high precision algebraic and transcendental
!   routines.

use mpfuna
use mpfunc
use mpfund
use mpfune
use mpfung
use mpfunh
contains

subroutine mpagmx (a, b)

!   This performs the arithmetic-geometric mean (AGM) iterations.  This routine
!   is called by MPLOGX.  It is not intended to be called directly by the user.

!   Max SP space for A and B: MPNW + 4 cells.

dimension a(mpnw+4), b(mpnw+4), s(2*mpnw+8)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  a(1) = 0.
  a(2) = 0.
  b(1) = 0.
  b(2) = 0.
  return
endif
n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
s(k0) = 0.
s(k0+1) = 0.
l1 = 0

100  l1 = l1 + 1
if (l1 .eq. 50) then
  if (mpker(5) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPAGMX: Iteration limit exceeded.')
    mpier = 5
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
endif

s1 = s(k0+1)
call mpadd (a, b, s(k0))
call mpmuld (s(k0), 0.5d0, 0, s(k1))
call mpmulx (a, b, s(k0))
call mpsqrx (s(k0), b)
call mpeq (s(k1), a)
call mpsub (a, b, s(k0))

!   Check for convergence.

if (s(k0) .ne. 0. .and. (s(k0+1) .lt. s1 .or. s(k0+1) .ge. -2)) goto 100

if (mpidb .ge. 6) write (mpldb, 2) l1, s(k0+1)
2 format ('MPAGMX: Iter., Tol. Achieved =',i5,f8.0)
return
end subroutine

subroutine mpcshx (a, pi, al2, x, y)

!   This computes the hyperbolic cosine and sine of the MP number A and
!   returns the two MP results in X and Y, respectively.  PI is the MP value
!   of Pi computed by a previous call to MPPI or MPPIX.  AL2 is the MP value
!   of Log (10) computed by a previous call to MPLOG or MPLOGX.  Before
!   calling MPCSHX, the arrays UU1 and Uu2 must be initialized by calling
!   MPINIX.  For modest levels of precision, use MPCSSH.  The last word 
!   of the result is not reliable.  Debug output starts with MPIDB = 5.

!   Max SP space for X and Y: MPNW + 4 cells.

dimension a(mpnw+2), f(8), al2(mpnw+2), pi(mpnw+2), x(mpnw+4), y(mpnw+4), &
  s(3*mpnw+12)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  x(1) = 0.
  x(2) = 0.
  y(1) = 0.
  y(2) = 0.
  return
endif
if (mpidb .ge. 5) call mpdeb ('MPCSHX I', a)

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.

call mpexpx (a, pi, al2, s(k0))
call mpdivx (f, s(k0), s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpmuld (s(k2), 0.5d0, 0, x)
call mpsub (s(k0), s(k1), s(k2))
call mpmuld (s(k2), 0.5d0, 0, y)

if (mpidb .ge. 5) then
  call mpdeb ('MPCSHX O', x)
  call mpdeb ('MPCSHX O', y)
endif
return
end subroutine

subroutine mpexpx (a, pi, al2, b)

!   This computes the exponential function of the MP number A and returns the
!   MP result in B.  PI is the MP value of Pi produced by a prior call to MPPI
!   or MPPIX.  AL2 is the MP value of Log(2) produced by a prior call to
!   MPLOG  or MPLOGX.  Before calling MPEXPX, the arrays UU1 and UU2 must be
!   initialized by calling MPINIX.  For modest levels of precision, use MPEXP.
!   The last word of the result is not reliable.  Debug output starts 
!   with MPIDB = 5.

!   Max SP space for B: MPNW + 4 cells.

!   This routine uses the Newton iteration

!     b_{k+1} = b_k [a + 1 - log b_k]

!   with a dynamically changing level of precision.  Logs are performed using
!   MPLOGX.  See the comment about the parameter NIT in MPDIVX.

double precision alt, cl2, cpi, t1, t2
parameter (alt = 0.693147180559945309d0, cl2 = 1.4426950408889633d0, &
  cpi = 3.141592653589793238d0, nit = 1)
dimension a(mpnw+2), al2(mpnw+2), b(mpnw+4), f1(8), pi(mpnw+2), s(3*mpnw+12)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 5) call mpdeb ('MPEXPX I', a)

ncr = 2 ** mpmcr
ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
call mpmdc (a, t1, n1)
t1 = t1 * 2.d0 ** n1

!   Check if precision level is too low to justify the advanced routine.

if (mpnw .le. ncr) then
  call mpexp (a, al2, b)
  goto 120
endif

!   Check if Log(2) has been precomputed.

call mpmdc (al2, t2, n2)
if (n2 .ne. - mpnbt .or. abs (t2 * 0.5d0 ** mpnbt - alt) .gt. mprx2) then
  if (mpker(37) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPEXPX: LOG (2) must be precomputed.')
    mpier = 37
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!     Check if Pi has been precomputed.

call mpmdc (pi, t2, n2)
if (n2 .ne. 0 .or. abs (t2 - cpi) .gt. mprx2) then
  if (mpker(38) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPEXPX: PI must be precomputed.')
    mpier = 38
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check for overflows and underflows.

if (abs (t1) .gt. 33271064.d0) then
  if (t1 .gt. 0.d0) then
    if (mpker(39) .ne. 0) then
      write (mpldb, 3) t1
3     format ('*** MPEXPX: Argument is too large',1p,d25.16)
      mpier = 39
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  else
    b(1) = 0.
    b(2) = 0.
    goto 120
  endif
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
nws = mpnw
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t2 = nws
mq = cl2 * log (t2) + 1.d0 - mprxx
call mpadd (a, f1, s(k0))

!   Compute initial approximation to Exp (A).

mpnw = ncr
call mpexp (a, al2, b)
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW.

do k = mpmcr + 1, mq
  mpnw = min (2 * mpnw, nws)
100  continue
  call mplogx (b, pi, al2, s(k1))
  call mpsub (s(k0), s(k1), s(k2))
  call mpmulx (b, s(k2), s(k1))
  call mpeq (s(k1), b)
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
enddo

120  if (mpidb .ge. 6) call mpdeb ('MPEXPX O', b)
return
end subroutine

subroutine mplogx (a, pi, al2, b)

!   This computes the natural logarithm of the MP number A and returns the MP
!   result in B.  PI is the MP value of Pi produced by a prior call to MPPI or
!   MPPIX.  AL2 is the MP value of Log(2) produced by a prior call to MPLOG
!   or MPLOGX.  Before calling MPLOGX, the arrays UU1 and UU2 must be
!   initialized by calling MPINIX.  For modest levels of precision, use MPLOG.
!   The last word of the result is not reliable.  Debug output starts 
!   with MPIDB = 6.

!   Max SP space for B: MPNW + 4 cells.

!   This uses the following algorithm, which is due to Salamin.  If a is
!   extremely close to 1, use a Taylor series.  Otherwise select n such that
!   z = x 2^n is at least 2^m, where m is the number of bits of desired
!   precision in the result.  Then

!   Log(x) = Pi / [2 AGM (1, 4/x)]

double precision alt, cpi, st, t1, t2, tn
parameter (mzl = -5, alt = 0.693147180559945309d0, cpi = 3.141592653589793d0)
dimension al2(mpnw+2), f1(8), f4(8), pi(mpnw+2), a(mpnw+4), b(mpnw+4), &
  s(4*mpnw+16)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  return
endif
if (mpidb .ge. 6) call mpdeb ('MPLOGX I', a)

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ncr = 2 ** mpmcr

!   Check if precision level is too low to justify the advanced routine.

if (mpnw .le. ncr) then
  call mplog (a, al2, b)
  goto 110
endif

if (ia .lt. 0 .or. na .eq. 0) then

!   Input is less than or equal to zero.

  if (mpker(52) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPLOGX: Argument is less than or equal to zero.')
    mpier = 52
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if Pi has been precomputed.

call mpmdc (pi, t1, n1)
if (n1 .ne. 0 .or. abs (t1 - cpi) .gt. mprx2) then
  if (mpker(53) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPLOGX: PI must be precomputed.')
    mpier = 53
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Unless the input is 2, Log (2) must have been precomputed.

if (a(1) .ne. 1. .or. a(2) .ne. 0. .or. a(3) .ne. 2.) then
  it2 = 0
  call mpmdc (al2, t2, n2)
  if (n2 .ne. - mpnbt .or. abs (t2 * 0.5d0 ** mpnbt - alt) .gt. mprx2) then
    if (mpker(54) .ne. 0) then
      write (mpldb, 3)
3     format ('*** MPLOGX: Log (2) must be precomputed.')
      mpier = 54
      if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  endif
else
  it2 = 1
endif

!   Define sections of the scratch array.

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
k3 = k2 + n4
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
f4(1) = 1.
f4(2) = 0.
f4(3) = 4.
f4(4) = 0.

!   If argument is 1, the result is zero.  If the argument is extremely close
!   to 1.  If so, employ a Taylor's series instead.

call mpsub (a, f1, s(k0))
if (s(k0) .eq. 0.) then
  b(1) = 0.
  b(2) = 0.
  goto 110
elseif (s(k0+1) .le. mzl) then
  call mpeq (s(k0), s(k1))
  call mpeq (s(k1), s(k2))
  i1 = 1
  is = 1
  tl = s(k0+1) - mpnw - 1

100  i1 = i1 + 1
  is = - is
  st = is * i1
  call mpmulx (s(k1), s(k2), s(k3))
  call mpdivd (s(k3), st, 0, s(k2))
  call mpadd (s(k0), s(k2), s(k3))
  call mpeq (s(k3), s(k0))
  if (s(k2+1) .ge. tl) goto 100

  call mpeq (s(k0), b)
  goto 110
endif

!   If input is exactly 2, set the exponent to a large value.  Otherwise
!   multiply the input by a large power of two.

call mpmdc (a, t1, n1)
n2 = mpnbt * (mpnw / 2 + 2) - n1
tn = n2
if (it2 .eq. 1) then
  call mpdmc (1.d0, n2, s(k0))
else
  call mpmuld (a, 1.d0, n2, s(k0))
endif

!   Perform AGM iterations.

call mpeq (f1, s(k1))
call mpdivx (f4, s(k0), s(k2))
call mpagmx (s(k1), s(k2))

!   Compute B = Pi / (2 * A), where A is the limit of the AGM iterations.

call mpmuld (s(k1), 2.d0, 0, s(k0))
call mpdivx (pi, s(k0), s(k1))

!   If the input was exactly 2, divide by TN.  Otherwise subtract TN * Log(2).

if (it2 .eq. 1) then
  call mpdivd (s(k1), tn, 0, s(k0))
else
  call mpmuld (al2, tn, 0, s(k2))
  call mpsub (s(k1), s(k2), s(k0))
endif
call mpeq (s(k0), b)

110  if (mpidb .ge. 6) call mpdeb ('MPLOGX O', b)
return
end subroutine

recursive subroutine mpoutx (a, b, n)

!   Converts the MP number A into character form in the CHARACTER*1 array B.
!   N (an output parameter) is the length of the output.  In other words, B is
!   contained in B(1), ..., B(N).  The format is analogous to the Fortran
!   exponential format (E format), except that the exponent is placed first.
!   Before calling MPOUTX, the arrays UU1 and UU2 must be initialized by
!   calling MPINIX.  For modest levels of precision, use MPOUTC.  Debug output
!   starts with MPIDB = 7.

!   Max CHARACTER*1 space for B: 7.225 * MPNW + 30 cells.

double precision al2, t1, t2, t3
character*1 b, b1, b2
character*10 dig
character*16 c1, c2
parameter (al2 = 0.301029995663981195d0, dig = '0123456789')
dimension a(mpnw+2), b(*), b1(8*mpnw+30), b2(8*mpnw+30), s(5*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = ' '
  n = 0
  return
endif
if (mpidb .ge. 7) then
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 1) (a(i), i = 1, no)
1 format ('MPOUTX I'/(6f12.0))
endif

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
k3 = k2 + n4
k4 = k3 + n4
ncr = 2 ** mpmcr

!   Check if actual precision level of argument is too low to justify the
!   advanced routine.

if (na .le. ncr) then
  call mpoutc (a, b, n)
  goto 110
endif

!   Normalize input to an integer by multiplying by a suitable power of 10.

t1 = a(3) + mprdx * a(4) + mprx2 * a(5)
t2 = log10 (t1)
m1 = max (al2 * mpnbt * (abs (a(1)) - a(2)) - t2, 0.d0)
call mpdmc (10.d0, 0, s(k0))
call mpnpwx (s(k0), m1, s(k2))
call mpmulx (a, s(k2), s(k1))
s(k1) = abs (s(k1))

!   Split large integer into two approximately equal decimal sections.

call mpmdc (s(k1), t1, i1)
call dpdec (t1, i1, t2, i2)
m2 = i2 / 2
call mpnpwx (s(k0), m2, s(k3))
call mpdivx (s(k1), s(k3), s(k0))
call mpinfr (s(k0), s(k2), s(k4))
call mpmulx (s(k2), s(k3), s(k0))
call mpsub (s(k1), s(k0), s(k3))

!   Recursively convert each section.

mpnws = mpnw
mpnw = s(k2) + 1
call mpoutx (s(k2), b1, nb1)
mpnw = s(k3) + 1
call mpoutx (s(k3), b2, nb2)
mpnw = mpnws

!   Obtain decimal exponents from each section.

c1 = ' '
c2 = ' '

do i = 1, 10
  c1(i:i) = b1(i+4)
  c2(i:i) = b2(i+4)
enddo

read (c1, '(I10)') ie1
read (c2, '(I10)') ie2

!   Set exponent of result.

ie = ie1 + m2 - m1
write (c1, '(I14)') ie

do i = 1, 4
  b(i) = b1(i)
enddo

do i = 5, 14
  b(i) = c1(i:i)
enddo

!   Copy mantissa of first section.

do i = 15, nb1
  b(i) = b1(i)
enddo

i1 = ie1 + m2 - ie2 + 19

!   If first section is too long, then round trailing digits (probably 9s).

if (nb1 .gt. i1) then
  i2 = index (dig, b(i1+1)) - 1
  if (i2 .ge. 5) then
    do i = i1, 21, -1
      if (b(i) .ne. '9') goto 100
      b(i) = '0'
    enddo

    write (mpldb, 2)
2   format ('*** MPOUTX: Exceptional case -- contact DHB.')
    stop

100 i2 = index (dig, b(i)) - 1
    write (c1, '(I1)') i2 + 1
    b(i) = c1(1:1)
  endif
elseif (nb1 .lt. i1) then

!   If first section is too short, then insert zeroes in gap.

  do i = nb1 + 1, i1
    b(i) = '0'
  enddo
endif

!   Copy mantissa of second section.

b(i1+1) = b2(19)
n = min (i1 + nb2 - 19, int (7.225 * mpnw + 30))

do i = i1 + 2, n
  b(i) = b2(i-i1+19)
enddo

!   Fix sign.

if (ia .eq. -1) b(18) = '-'

110 continue
if (mpidb .ge. 7) then
  no = min (n, 6 * mpndb + 20)
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPOUTX O'/(78a1))
endif

return
end subroutine

subroutine mppix (pi)

!   This computes Pi to available precision (MPNW mantissa words).  Before
!   calling MPPIX, the arrays UU1 and UU2 must be initialized by calling
!   MPINIX.  For modest levels of precision, use MPPI.  MPNW should be a power
!   of two.  The last word of the result is not reliable.  Debug
!   output starts with MPIDB = 7.

!   Max SP space for PI: MPNW + 4 cells.

!   This routine uses basically the same algorithm as MPPI.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0)
dimension f(8), pi(mpnw+4), s(5*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  pi(1) = 0.
  pi(2) = 0.
  return
endif
ncr = 2 ** mpmcr

!   Check if precision level is too low to justify the advanced routine.

if (mpnw .le. ncr) then
  call mppi (pi)
  goto 110
endif
n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
k3 = k2 + n4
k4 = k3 + n4

!   Determine the number of iterations required for the given precision level.
!   This formula is good only for this Pi algorithm.

t1 = mpnw * log10 (mpbdx)
mq = cl2 * (log (t1) - 1.d0) + 1.d0

!   Initialize as above.

s(k0) = 1.
s(k0+1) = 0.
s(k0+2) = 1.
f(1) = 1.
f(2) = 0.
f(3) = 2.
f(4) = 0.
call mpsqrx (f, s(k2))
call mpmuld (s(k2), 0.5d0, 0, s(k1))
f(2) = -1.
f(3) = 0.5d0 * mpbdx
call mpsub (s(k2), f, s(k4))

!   Perform iterations as described above.

do k = 1, mq
  call mpadd (s(k0), s(k1), s(k2))
  call mpmulx (s(k0), s(k1), s(k3))
  call mpsqrx (s(k3), s(k1))
  call mpmuld (s(k2), 0.5d0, 0, s(k0))
  call mpsub (s(k0), s(k1), s(k2))
  call mpsqx (s(k2), s(k3))
  t1 = 2.d0 ** k
  call mpmuld (s(k3), t1, 0, s(k2))
  call mpsub (s(k4), s(k2), s(k3))
  call mpeq (s(k3), s(k4))
enddo

!   Complete computation.

call mpadd (s(k0), s(k1), s(k2))
call mpsqx (s(k2), s(k3))
call mpdivx (s(k3), s(k4), s(k2))
call mpeq (s(k2), pi)

110  if (mpidb .ge. 7) call mpdeb ('MPPIX O', pi)

return
end subroutine

subroutine mppolx (n, l, a, x1, nx, x)

!   This finds a real root of the N-th degree polynomial whose MP coefficients
!   are in A by Newton-Raphson iterations, beginning at the DP value (X1, NX)
!   and returns the MP root in X.  The N + 1 coefficients a_0, a_1, ..., a_N
!   are assumed to start in locations A(1), A(L+1), A(2*L+1), etc.  Before
!   calling MPPOLX, the arrays UU1 and UU2 must be initialized by calling
!   MPINIX.  For modest levels of precision, use MPPOL.  The last word
!   of the result is not reliable.  Debug output starts with MPIDB = 5.

!   Max SP space for X: MPNW + 4 cells.

!   For a discussion of the algorithm and usage, see MPPOL.  This routine uses
!   basically the same Newton iteration algorithm as MPPOL.  In fact, this
!   routine calls MPPOL to obtain an initial approximation.

character*8 cx
double precision t1, x1
dimension a(l,n+1), x(mpnw+4), s(5*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  x(1) = 0.
  x(2) = 0.
  return
endif
if (mpidb .ge. 5) then
  write (mpldb, 1) n
1 format ('MPPOLX I',i4)

  do k = 0, n
    write (cx, '(I4)') k
    call mpdeb (cx, a(1,k+1))
  enddo

  write (mpldb, 2) x1, nx
2 format ('MPPOLX I',f16.12,' x 10 ^',i6)
endif

!   Check if precision level is too low to justify the advanced routine.

ncr = 2 ** mpmcr
if (mpnw .le. ncr) then
  call mppol (n, l, a, x1, nx, x)
  l1 = 0
  goto 150
endif

!   Check if the polynomial is proper.

if (a(1,1) .eq. 0. .or. a(1,n+1) .eq. 0.) then
  if (mpker(65) .ne. 0) then
    write (mpldb, 3)
3   format ('*** MPPOLX: Either the first or last input coefficient is zero.')
    mpier = 65
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
k3 = k2 + n4
k4 = k3 + n4
nws = mpnw

!   Compute the initial approximation.

mpnw = ncr
call mppol (n, l, a, x1, nx, x)
call mpeq (x, s(k0))
tl = 2 - mpnw
l1 = 0
ls = -10

!   Perform MP Newton-Raphson iterations to solve P(x) = 0.

110  l1 = l1 + 1
if (l1 .eq. 50) then
  if (mpker(66) .ne. 0) then
    write (mpldb, 4)
4   format ('*** MPPOLX: Iteration limit exceeded.')
    mpier = 66
    if (mpker(mpier) .eq. 2) call mpabrt
    mpnw = nws
    return
  endif
endif

!   Compute P(x).

call mpeq (a(1,n+1), s(k1))

do k = n - 1, 0, -1
  call mpmulx (s(k0), s(k1), s(k2))
  call mpadd (s(k2), a(1,k+1), s(k1))
enddo

!   Compute P'(x).

t1 = n
call mpmuld (a(1,n+1), t1, 0, s(k2))

do k = n - 1, 1, -1
  call mpmulx (s(k0), s(k2), s(k3))
  t1 = k
  call mpmuld (a(1,k+1), t1, 0, s(k4))
  call mpadd (s(k3), s(k4), s(k2))
enddo

!   Compute P(x) / P'(x) and update x.

call mpdivx (s(k1), s(k2), s(k3))
call mpsub (s(k0), s(k3), s(k4))

if (mpidb .ge. 6) then
  write (mpldb, 5) l1
5 format ('Iteration',i4)
  call mpdeb ('X', s(k0))
  call mpdeb ('P(X)', s(k1))
  call mpdeb ('P''(X)', s(k2))
  call mpdeb ('CORR', s(k3))
endif
call mpeq (s(k4), s(k0))

!   If this was the second iteration at full precision, there is no need to
!   continue (the adjusted value of x is correct); otherwise repeat.

if (l1 .eq. ls + 1) goto 140
if (s(k3) .ne. 0. .and. s(k3+1) .gt. tl) goto 110

!   Newton iterations have converged to current precision.  Increase precision
!   and continue.

if (mpnw .eq. nws) goto 140
mpnw = min (2 * mpnw, nws)
if (mpnw .eq. nws) ls = l1
if (mpnw .le. 32) then
  tl = 2 - mpnw
elseif (mpnw .le. 256) then
  tl = 3 - mpnw
else
  tl = 4 - mpnw
endif
if (mpidb .ge. 6) then
  write (mpldb, 6) mpnw
6 format (6x,'New MPNW =', i8)
endif
goto 110

140  call mpeq (s(k0), x)

150  if (mpidb .ge. 5) then
  write (mpldb, 7) l1
7 format ('Iteration count:',i5)
  call mpdeb ('MPPOLX O', x)
endif
return
end subroutine

end module

module mpfunj
use mpfuna
use mpfunc
use mpfund
use mpfune
use mpfunf
use mpfunh
use mpfuni
contains

subroutine mpangx (x, y, pi, a)

!   This computes the MP angle A subtended by the MP pair (X, Y) considered as
!   a point in the x-y plane.  This is more useful than an arctan or arcsin
!   routine, since it places the result correctly in the full circle, i.e.
!   -Pi < A <= Pi.  PI is the MP value of Pi computed by a previous call to
!   MPPI or MPPIX.  Before calling MPANGX, the arrays UU1 and UU2 must be
!   initialized by calling MPINIX.  For modest levels of precision, use MPANG.
!   The last word of the result is not reliable.  Debug output starts 
!   with MPIDB = 6.

!   Max SP space for A: MPNW + 4 cells.

!   This routine employs a complex arithmetic version of the MPLOGX alogirthm.

double precision cpi, t1
parameter (cpi = 3.141592653589793d0)
dimension a(mpnw+4), f0(8), f1(8), f4(8), pi(mpnw+2), x(mpnw+2), y(mpnw+2), &
  s(8*mpnw+32)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  a(1) = 0.
  a(2) = 0.
  return
endif
if (mpidb .ge. 6) then
  call mpdeb ('MPANGX I', x)
  call mpdeb ('MPANGX I', y)
endif

ix = sign (1., x(1))
nx = min (int (abs (x(1))), mpnw)
iy = sign (1., y(1))
ny = min (int (abs (y(1))), mpnw)
ncr = 2 ** mpmcr

!   Check if precision level is too low to justify the advanced routine.

if (mpnw .le. ncr) then
  call mpang (x, y, pi, a)
  goto 100
endif

!   Check if both X and Y are zero.

if (nx .eq. 0 .and. ny .eq. 0) then
  if (mpker(9) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPANGX: Both arguments are zero.')
    mpier = 9
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if Pi has been precomputed.

call mpmdc (pi, t1, n1)
if (n1 .ne. 0 .or. abs (t1 - cpi) .gt. mprx2) then
  if (mpker(10) .ne. 0) then
    write (mpldb, 2)
2   format ('*** MPANGX: PI must be precomputed.')
    mpier = 10
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

!   Check if one of X or Y is zero.

if (nx .eq. 0) then
  if (iy .gt. 0) then
    call mpmuld (pi, 0.5d0, 0, a)
  else
    call mpmuld (pi, -0.5d0, 0, a)
  endif
  goto 100
elseif (ny .eq. 0) then
  if (ix .gt. 0) then
    a(1) = 0.
    a(2) = 0.
  else
    call mpeq (pi, a)
  endif
  goto 100
endif

!   Define scratch space.

n4 = mpnw + 4
n42 = 2 * n4
k0 = 1
k1 = k0 + n42
k2 = k1 + n42
k3 = k2 + n42
f0(1) = 0.
f0(2) = 0.
f0(3) = 0.
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
f4(1) = 1.
f4(2) = 0.
f4(3) = 4.
f4(4) = 0.

!   Multiply the input by a large power of two.

call mpmdc (x, t1, n1)
n2 = mpnbt * (mpnw / 2 + 2) - n1
tn = n2
call mpmuld (x, 1.d0, n2, s(k1))
call mpmuld (y, 1.d0, n2, s(k2))
call mpmmpc (s(k1), s(k2), n4, s(k0))

!   Perform AGM iterations.

call mpmmpc (f1, f0, n4, s(k1))
call mpmmpc (f4, f0, n4, s(k3))
call mpcdvx (n4, s(k3), s(k0), s(k2))
call mpcagx (s(k1), s(k2))

!   Compute A = Imag (Pi / (2 * Z)), where Z is the limit of the complex AGM.

call mpmuld (s(k1), 2.d0, 0, s(k0))
call mpmuld (s(k1+n4), 2.d0, 0, s(k0+n4))
call mpmmpc (pi, f0, n4, s(k2))
call mpcdvx (n4, s(k2), s(k0), s(k1))
call mpeq (s(k1+n4), a)

100  if (mpidb .ge. 6) call mpdeb ('MPANGX O', a)

return
end subroutine

subroutine mpcagx (a, b)

!   This performs the arithmetic-geometric mean (AGM) iterations.  This routine
!   is called by MPANGX.  It is not intended to be called directly by the user.

!   Max SP space for A and B: 2*MPNW + 8 cells.

dimension a(2*mpnw+8), b(2*mpnw+8), s(4*mpnw+16)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  a(1) = 0.
  a(2) = 0.
  b(1) = 0.
  b(2) = 0.
  return
endif
n4 = mpnw + 4
k0 = 1
k1 = k0 + 2 * n4
s(k0) = 0.
s(k0+1) = 0.
l1 = 0

100  l1 = l1 + 1
if (l1 .eq. 50) then
  if (mpker(12) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPCAGX: Iteration limit exceeded.')
    mpier = 12
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
endif

s1 = s(k0+1)
call mpcadd (n4, a, b, s(k0))
call mpmuld (s(k0), 0.5d0, 0, s(k1))
call mpmuld (s(k0+n4), 0.5d0, 0, s(k1+n4))
call mpcmlx (n4, a, b, s(k0))
call mpcsqx (n4, s(k0), b)
call mpceq (n4, s(k1), a)
call mpsub (a, b, s(k0))

!   Check for convergence.

if (s(k0) .ne. 0. .and. (s(k0+1) .lt. s1 .or. s(k0+1) .ge. -2)) goto 100

if (mpidb .ge. 6) write (mpldb, 2) l1, s(k0+1)
2 format ('MPCAGX: Iter., Tol. Achieved =',i5,f8.0)
return
end subroutine

subroutine mpcdvx (l, a, b, c)

!   This routine divides the MP complex numbers A and B to yield the MPC
!   quotient C.  L is the offset between real and imaginary parts in A, B
!   the result C.  L must be at least MPNW + 4.  Before calling MPCDVX, the
!   arrays UU1 and UU2 must be initialized by calling MPINIX. For modest levels
!   of precision, use MPCDIV.  The last word of the result is not reliable.  
!   Debug output starts with MPIDB = 7

!   Max SP space for C: 2 * L cells.

!   This routine employs the same scheme as MPCDIV.

dimension a(2*l), b(2*l), c(2*l), f(8), s(5*mpnw+20)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  c(l+1) = 0.
  c(l+2) = 0.
  return
endif
l1 = l + 1
if (mpidb .ge. 7) then
  write (mpldb, 1) l
1 format ('MPCDVX I',i10)
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPCDVX I'/(6f12.0))
  no = min (int (abs (a(l1))), mpndb) + 2
  write (mpldb, 2) (a(l+i), i = 1, no)
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 2) (b(i), i = 1, no)
  no = min (int (abs (b(l1))), mpndb) + 2
  write (mpldb, 2) (b(l+i), i = 1, no)
endif

if (l .lt. mpnw + 4) then
  if (mpker(17) .ne. 0) then
    write (mpldb, 3) l, mpnw + 4
3   format ('*** MPCDVX: Offset parameter is too small',2i8)
    mpier = 17
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

if (b(1) .eq. 0. .and. b(l1) .eq. 0.) then
  if (mpker(18) .ne. 0) then
    write (mpldb, 4)
4   format ('*** MPCDVX: Divisor is zero.')
    mpier = 18
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
k3 = k2 + n4
k4 = k3 + n4
f(1) = 1.
f(2) = 0.
f(3) = 1.
f(4) = 0.

call mpmulx (a, b, s(k0))
call mpmulx (a(l1), b(l1), s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpsub (s(k0), s(k1), s(k3))
call mpadd (a, a(l1), s(k0))
call mpsub (b, b(l1), s(k1))
call mpmulx (s(k0), s(k1), s(k4))
call mpsub (s(k4), s(k3), s(k1))
call mpsqx (b, s(k0))
call mpsqx (b(l1), s(k3))
call mpadd (s(k0), s(k3), s(k4))
call mpdivx (f, s(k4), s(k0))
call mpmul (s(k2), s(k0), c)
call mpmul (s(k1), s(k0), c(l1))

if (mpidb .ge. 7) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 5) (c(i), i = 1, no)
5 format ('MPCDVX O'/(6f12.0))
  no = min (int (abs (c(l1))), mpndb) + 2
  write (mpldb, 5) (c(l+i), i = 1, no)
endif
return
end subroutine

subroutine mpcmlx (l, a, b, c)

!   This routine multiplies the MP complex numbers A and B to yield the MPC
!   product C.  L is the offset between real and imaginary parts in A, B and
!   the result C.  L must be at least MPNW + 4.  Before calling MPCMLX, the
!   arrays UU1 and UU2 must be initialized by calling MPINIX. For modest levels
!   of precision, use MPCMUL.  The last word of the result is not reliable.  
!   Debug output starts with MPIDB = 7.

!   Max SP space for C: 2 * L cells.

!   This routine employs the same scheme as MPCMUL.

dimension a(2*l), b(2*l), c(2*l), s(4*mpnw+16)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  c(1) = 0.
  c(2) = 0.
  c(l+1) = 0.
  c(l+2) = 0.
  return
endif
l1 = l + 1
if (mpidb .ge. 7) then
  write (mpldb, 1) l
1 format ('MPCMLX I',i10)
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPCMLX I'/(6f12.0))
  no = min (int (abs (a(l1))), mpndb) + 2
  write (mpldb, 2) (a(l+i), i = 1, no)
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 2) (b(i), i = 1, no)
  no = min (int (abs (b(l1))), mpndb) + 2
  write (mpldb, 2) (b(l+i), i = 1, no)
endif

if (l .lt. mpnw + 4) then
  if (mpker(19) .ne. 0) then
    write (mpldb, 3) l, mpnw + 4
3   format ('*** MPCMLX: Offset parameter is too small',2i8)
    mpier = 19
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4
k3 = k2 + n4

call mpmulx (a, b, s(k0))
call mpmulx (a(l1), b(l1), s(k1))
call mpsub (s(k0), s(k1), c)
call mpadd (s(k0), s(k1), s(k2))
call mpadd (a, a(l1), s(k0))
call mpadd (b, b(l1), s(k1))
call mpmulx (s(k0), s(k1), s(k3))
call mpsub (s(k3), s(k2), c(l1))

if (mpidb .ge. 7) then
  no = min (int (abs (c(1))), mpndb) + 2
  write (mpldb, 4) (c(i), i = 1, no)
4 format ('MPCMLX O'/(6f12.0))
  no = min (int (abs (c(l1))), mpndb) + 2
  write (mpldb, 4) (c(l+i), i = 1, no)
endif
return
end subroutine

subroutine mpcplx (n, la, a, x1, nx, lx, x)

!   This routine finds a complex root of the N-th degree polynomial whose
!   MPC coefficients are in A by Newton-Raphson iterations, beginning
!   at the complex DPE value (X1(1), NX(1)) + i (X1(2), NX(2)), and returns
!   the MPC root in X.  The N + 1 coefficients a_0, a_1, ..., a_N are
!   assumed to start in locations A(1), A(2*LA+1), A(4*LA+1), etc.  LA is the
!   offset between the real and the imaginary parts of each input coefficient.
!   Typically LA = MPNW + 4.  LX, also an input parameter, is the offset
!   between the real and the imaginary parts of the result to be stored in X.
!   LX should be at least MPNW + 4.  Before calling MPCPLX, the arrays UU1 and
!   UU2 must be initialized by calling MPINIX.  For modest levels of precision,
!   use MPCPOL.  The last word of the result is not reliable.  Debug 
!   output starts with MPIDB = 5.

!   Max SP space for X: 2 * LX cells.

!   See the note in MPPOL about repeated roots.

!   This routine employs the same scheme as MPCPOL.

character*8 cx
double precision t1, x1
dimension a(2*la,n+1), nx(2), x(2*lx), x1(2), s(10*mpnw+40)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  x(1) = 0.
  x(2) = 0.
  x(lx+1) = 0.
  x(lx+2) = 0.
endif
if (mpidb .ge. 5) then
  write (mpldb, 1) n, lx
1 format ('MPCPLX I',2i6)

  do k = 0, n
    write (cx, '(I4)') k
    call mpdeb (cx, a(1,k+1))
    call mpdeb (cx, a(la+1,k+1))
  enddo

  write (mpldb, 2) x1(2), nx(2)
2 format ('MPCPLX I',f16.12,' x 10 ^',i6,f20.12,' x 10^',i6)
endif

!   Check if precision level is too low to justify the advanced routine.

ncr = 2 ** mpmcr
if (mpnw .le. ncr) then
  call mpcpol (n, la, a, x1, nx, lx, x)
  l1 = 0
  goto 150
endif

!   Check if the polynomial is proper.

if (a(1,1) .eq. 0. .or. a(1,n+1) .eq. 0.) then
  if (mpker(21) .ne. 0) then
    write (mpldb, 3)
3   format ('*** MPCPLX: Either the first or last input coefficient is zero.')
    mpier = 21
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n4 = mpnw + 4
n8 = 2 * n4
k0 = 1
k1 = k0 + n8
k2 = k1 + n8
k3 = k2 + n8
k4 = k3 + n8
nws = mpnw

!   Set the initial value.

mpnw = ncr
call mpcpol (n, la, a, x1, nx, n4, s(k0))
tl = -4.
l1 = 0
ls = -10

!   Perform MP Newton-Raphson iterations to solve P(x) = 0.

110   l1 = l1 + 1
if (l1 .eq. 50) then
  if (mpker(22) .ne. 0) then
    write (mpldb, 4)
4   format ('*** MPCPLX: Iteration limit exceeded.')
    mpier = 22
    if (mpker(mpier) .eq. 2) call mpabrt
    mpnw = nws
    return
  endif
endif

!   Compute P(x).

call mpmmpc (a(1,n+1), a(la+1,n+1), n4, s(k1))

do k = n - 1, 0, -1
  call mpcmlx (n4, s(k0), s(k1), s(k2))
  call mpadd (s(k2), a(1,k+1), s(k1))
  call mpadd (s(k2+n4), a(la+1,k+1), s(k1+n4))
enddo

!   Compute P'(x).

t1 = n
call mpmuld (a(1,n+1), t1, 0, s(k2))
call mpmuld (a(la+1,n+1), t1, 0, s(k2+n4))

do k = n - 1, 1, -1
  call mpcmlx (n4, s(k0), s(k2), s(k3))
  t1 = k
  call mpmuld (a(1,k+1), t1, 0, s(k4))
  call mpmuld (a(la+1,k+1), t1, 0, s(k4+n4))
  call mpcadd (n4, s(k3), s(k4), s(k2))
enddo

!   Compute P(x) / P'(x) and update x.

call mpcdvx (n4, s(k1), s(k2), s(k3))
call mpcsub (n4, s(k0), s(k3), s(k4))

if (mpidb .ge. 6) then
  write (mpldb, 5) l1
5 format ('ITERATION',i4)
  call mpdeb ('X', s(k0))
  call mpdeb (' ', s(k0+n4))
  call mpdeb ('P(X)', s(k1))
  call mpdeb (' ', s(k1+n4))
  call mpdeb ('P''(X)', s(k2))
  call mpdeb (' ', s(k2+n4))
  call mpdeb ('CORR', s(k3))
  call mpdeb (' ', s(k3+n4))
endif
call mpceq (n4, s(k4), s(k0))

!   If this was the second iteration at full precision, there is no need to
!   continue (the adjusted value of x is correct); otherwise repeat.

if (l1 .eq. ls + 1) goto 140
if (s(k3) .ne. 0. .and. s(k3+1) .gt. tl .or. s(k3+n4) .ne. 0. &
  .and. s(k3+n4+1) .gt. tl) goto 110

!   Newton iterations have converged to current precision.  Increase precision
!   and continue.

if (mpnw .eq. nws) goto 140
mpnw = min (2 * mpnw, nws)
if (mpnw .eq. nws) ls = l1
if (mpnw .le. 32) then
  tl = 2 - mpnw
elseif (mpnw .le. 256) then
  tl = 3 - mpnw
else
  tl = 4 - mpnw
endif
if (mpidb .ge. 6) then
  write (mpldb, 6) mpnw
6 format (6x,'New MPNW =', i8)
endif
goto 110

140  call mpmmpc (s(k0), s(k0+n4), lx, x)

150  if (mpidb .ge. 5) then
  write (mpldb, 7) l1
7 format ('Iteration count:',i5)
  call mpdeb ('MPCPLX O', x)
  call mpdeb (' ', x(lx+1))
endif
return
end subroutine

subroutine mpcpwx (l, a, n, b)

!   This computes the N-th power of the MPC number A and returns the MPC
!   result C in B.  When N is zero, 1 is returned.  When N is negative, the
!   reciprocal of A ^ |N| is returned.  L is the offset between real and
!   imaginary parts in A and B.  L should be at least MPNW + 4.  Before calling
!   MPCPWX, the arrays UU1 and UU2 must be initialized by calling MPINIX.  For
!   modest levels of precision, use MPCPWR.  The last word of the result 
!   is not reliable.  Debug output starts with MPIDB = 6.

!   Max SP space for B: 2 * L cells.

!   This routine employs the binary method for exponentiation.

double precision cl2, t1
parameter (cl2 = 1.4426950408889633d0)
dimension a(2*l), b(2*l), f1(8), f2(8), s(6*mpnw+24)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  b(l+1) = 0.
  b(l+2) = 0.
  return
endif
l1 = l + 1
if (mpidb .ge. 6) then
  write (mpldb, 1) l, n
1 format ('MPCPWX I',2i10)
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPCPWX I'/(6f12.0))
  no = min (int (abs (a(l1))), mpndb) + 2
  write (mpldb, 2) (a(l+i), i = 1, no)
endif

na1 = min (int (abs (a(1))), mpnw)
na2 = min (int (abs (a(l1))), mpnw)
ncr = 2 ** mpmcr

!   Check if precision level of A is too low to justify advanced routine.

if (na1 .le. ncr .and. na2 .le. ncr) then
  call mpcpwr (l, a, n, b)
  goto 120
endif
if (na1 .eq. 0 .and. na2 .eq. 0) then
  if (n .ge. 0) then
    b(1) = 0.
    b(2) = 0.
    b(l1) = 0.
    b(l1+1) = 0.
    goto 120
  else
    if (mpker(26) .ne. 0) then
write (mpldb, 3)
3 format ('*** MPCPWX: Argument is zero and N is negative or zero.')
mpier = 26
if (mpker(mpier) .eq. 2) call mpabrt
    endif
    return
  endif
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + 2 * n4
k2 = k1 + 2 * n4
nn = abs (n)
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
f2(1) = 0.
f2(2) = 0.
call mpmmpc (a, a(l1), n4, s(k0))
if (nn .eq. 0) then
  call mpmmpc (f1, f2, l, b)
  goto 120
elseif (nn .eq. 1) then
  call mpceq (n4, s(k0), s(k2))
  goto 110
elseif (nn .eq. 2) then
  call mpcmlx (n4, s(k0), s(k0), s(k2))
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN .GT. NN.

t1 = nn
mn = cl2 * log (t1) + 1.d0 + mprxx
call mpmmpc (f1, f2, n4, s(k2))
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn .ne. 2 * kk) then
    call mpcmlx (n4, s(k2), s(k0), s(k1))
    call mpceq (n4, s(k1), s(k2))
  endif
  kn = kk
  if (j .lt. mn) then
    call mpcmlx (n4, s(k0), s(k0), s(k1))
    call mpceq (n4, s(k1), s(k0))
  endif
enddo

!   Compute reciprocal if N is negative.

110  if (n .lt. 0) then
  call mpmmpc (f1, f2, n4, s(k1))
  call mpcdvx (n4, s(k1), s(k2), s(k0))
  call mpceq (n4, s(k0), s(k2))
endif
call mpmmpc (s(k2), s(n4+k2), l, b)

120  if (mpidb .ge. 6) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 4) (b(i), i = 1, no)
4 format ('MPCPWX O'/(6f12.0))
  no = min (int (abs (b(l1))), mpndb) + 2
  write (mpldb, 4) (b(l+i), i = 1, no)
endif
return
end subroutine

subroutine mpcsqx (l, a, b)

!   This routine computes the complex square root of the MPC number C.  L is
!   the offset between real and imaginary parts in A and B.  L must be at
!   least MPNW + 4.  For modest levels of precision, use MPCSQT.  The last
!   word of the result is not reliable.  Debug output starts with MPIDB = 5.

!   Max SP space for B: 2 * L cells.

!   This routine uses the same algorithm as MPCSQT.

dimension a(2*l), b(2*l), s(3*mpnw+12)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  b(1) = 0.
  b(2) = 0.
  b(l+1) = 0.
  b(l+2) = 0.
  return
endif
l1 = l + 1
if (mpidb .ge. 5) then
  write (mpldb, 1) l
1 format ('MPCSQX I',i10)
  no = min (int (abs (a(1))), mpndb) + 2
  write (mpldb, 2) (a(i), i = 1, no)
2 format ('MPCSQX I'/(6f12.0))
  no = min (int (abs (a(l1))), mpndb) + 2
  write (mpldb, 2) (a(l+i), i = 1, no)
endif

if (a(1) .eq. 0. .and. a(l+1) .eq. 0.) then
  b(1) = 0.
  b(2) = 0.
  b(l+1) = 0.
  b(l+2) = 0.
  goto 100
endif

n4 = mpnw + 4
k0 = 1
k1 = k0 + n4
k2 = k1 + n4

call mpsqx (a, s(k0))
call mpsqx (a(l1), s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpsqrx (s(k2), s(k0))
call mpeq (a, s(k1))
s(k1) = abs (s(k1))
call mpadd (s(k0), s(k1), s(k2))
call mpmuld (s(k2), 0.5d0, 0, s(k1))
call mpsqrx (s(k1), s(k0))
call mpmuld (s(k0), 2.d0, 0, s(k1))
if (a(1) .ge. 0.) then
  call mpeq (s(k0), b)
  call mpdivx (a(l1), s(k1), b(l1))
else
  call mpdivx (a(l1), s(k1), b)
  b(1) = abs (b(1))
  call mpeq (s(k0), b(l1))
  b(l1) = sign (b(l1), a(l1))
endif

100  if (mpidb .ge. 5) then
  no = min (int (abs (b(1))), mpndb) + 2
  write (mpldb, 3) (b(i), i = 1, no)
3 format ('MPCSQX O'/(6f12.0))
  no = min (int (abs (b(l1))), mpndb) + 2
  write (mpldb, 3) (b(l+i), i = 1, no)
endif
return
end subroutine

subroutine mpcssx (a, pi, x, y)

!   This computes the cosine and sine of the MP number A and returns the two MP
!   results in X and Y, respectively.  PI is the MP value of Pi computed by a
!   previous call to MPPI or MPPIX.  Before calling MPCSSX, the arrays UU1 and
!   UU2 must be initialized by calling MPINIX.  For modest levels of
!   precision, use MPCSSN.  The last word of the result is not reliable.
!   Debug output starts with MPIDB = 5.

!   Max SP space for X and Y: MPNW + 4 cells.

!   This routine employs a complex arithmetic version of the scheme used in
!   MPEXPX.  See the comment about the parameter NIT in MPDIVX.

double precision cl2, cpi, t1, t2
parameter (cl2 = 1.4426950408889633d0, cpi = 3.141592653589793d0, nit = 1)
dimension a(mpnw+2), f1(8), pi(mpnw+2), x(mpnw+4), y(mpnw+4), s(8*mpnw+32)

if (mpier .ne. 0) then
  if (mpier .eq. 99) call mpabrt
  x(1) = 0.
  x(2) = 0.
  y(1) = 0.
  y(2) = 0.
  return
endif
if (mpidb .ge. 5) call mpdeb ('MPCSSX I', a)

ia = sign (1., a(1))
na = min (int (abs (a(1))), mpnw)
ncr = 2 ** mpmcr

!   Check if precision level is too low to justify advanced routine.

if (mpnw .le. ncr) then
  call mpcssn (a, pi, x, y)
  l1 = 0
  goto 120
endif

!   Check if input is zero.

if (na .eq. 0) then
  x(1) = 1.
  x(2) = 0.
  x(3) = 1.
  y(1) = 0.
  y(2) = 0.
  l1 = 0
  goto 120
endif

!   Check if Pi has been precomputed.

call mpmdc (pi, t1, n1)
if (n1 .ne. 0 .or. abs (t1 - cpi) .gt. mprx2) then
  if (mpker(30) .ne. 0) then
    write (mpldb, 1)
1   format ('*** MPCSSX: PI must be precomputed.')
    mpier = 30
    if (mpker(mpier) .eq. 2) call mpabrt
  endif
  return
endif

n4 = mpnw + 4
n42 = 2 * n4
k0 = 1
k1 = k0 + n42
k2 = k1 + n42
k3 = k2 + n42
f1(1) = 1.
f1(2) = 0.
f1(3) = 1.
f1(4) = 0.
nws = mpnw

!   Reduce argument to between - Pi and Pi.

call mpmuld (pi, 2.d0, 0, s(k0))
call mpdivx (a, s(k0), s(k1))
call mpnint (s(k1), s(k2))
call mpmul (s(k2), s(k0), s(k1))
call mpsub (a, s(k1), s(k0))

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t2 = nws
mq = cl2 * log (t2) + 1.d0 - mprxx
call mpeq (f1, s(k2))

!   Compute initial approximation to [Cos (A), Sin (A)].

mpnw = ncr
call mpcssn (s(k0), pi, s(k3), s(k3+n4))
iq = 0

!   Perform the Newton-Raphson iteration with a dynamically changing precision
!   level MPNW.

do k = mpmcr + 1, mq
  mpnw = min (2 * mpnw, nws)
100  continue
  call mpangx (s(k3), s(k3+n4), pi, s(k1))
  call mpsub (s(k0), s(k1), s(k2+n4))
  call mpcmlx (n4, s(k3), s(k2), s(k1))
  call mpceq (n4, s(k1), s(k3))
  if (k .eq. mq - nit .and. iq .eq. 0) then
    iq = 1
    goto 100
  endif
 enddo

!   The final (cos, sin) result must be normalized to have magnitude 1.

call mpsqx (s(k3), s(k0))
call mpsqx (s(k3+n4), s(k0+n4))
call mpadd (s(k0), s(k0+n4), s(k1))
call mpsqrx (s(k1), s(k2))
call mpdivx (s(k3), s(k2), s(k0))
call mpdivx (s(k3+n4), s(k2), s(k0+n4))
call mpmpcm (n4, s(k0), x, y)

120  if (mpidb .ge. 5) then
  call mpdeb ('MPCSSX O', x)
  call mpdeb ('MPCSSX O', y)
endif
return
end subroutine

end module

module mpfunmod
use mpfuna
use mpfunb
use mpfunc
use mpfund
use mpfune
use mpfunf
use mpfung
use mpfunh
use mpfuni
use mpfunj
end module

