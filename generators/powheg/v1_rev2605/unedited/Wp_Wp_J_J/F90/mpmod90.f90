module mpdefmod

!  These Fortran-90 modules perform translation of Fortran for multiprecision
!  execution.  Full details for usage are available in a separate document.
!  IEEE version.

!  David H. Bailey       1999-09-13

!  This code should work on any system employing the IEEE-754 arithmetic
!  standard, and on most others systems where the Fortran default real
!  datatype is a 32-bit word with 24 mantissa bits.

!  Code that needs to be changed for Cray vector systems is labeled with !>.
!  When compiling for Cray vector systems, disable double precision with -dp.

!  The following notational scheme is used to designate datatypes below:

!  A   Alphabetic [i.e. ASCII]
!  C   Complex
!  D   Double precision [i.e. REAL (KIND (0.D0))]
!  I   Integer
!  J   MP integer
!  Q   MP real
!  R   Real
!  X   Double complex  [i.e. COMPLEX (KIND (0.D0))]
!  Z   MP complex

!  The following parameters are all that need to be changed in normal usage:

!  MPIPL   Initial precision level, in digits.
!  MPIOU   Initial output precision level, in digits.
!  MPIEP   Log_10 of initial MP epsilon level.

use mpfunmod
! default values: 
parameter (mpipl = 2005, mpiou = 56, mpiep = 10 - mpipl)
!parameter (mpipl = 2005, mpiou = 15, mpiep = 10 - mpipl)

!----------------------------------------------------------------------------

!  *** No code below this point needs to be altered in normal usage.

!>  On Cray (non-IEEE) systems, change the constant 7.224719896 to 6.622659905.

parameter (mpwds = mpipl / 7.224719896d0 + 1)
private kdb, mp4, mp24, mp41
parameter (kdb = kind (0.d0), mp4 = mpwds + 4, mp24 = 2 * mp4, mp41 = mp4 + 1)
type mp_integer
  sequence
  real mpi(mp4)
end type
type mp_real
  sequence
  real mpr(mp4)
end type
type mp_complex
  sequence
  real mpc(mp24)
end type
type (mp_real), public:: mpl02, mpl10, mppic, mpeps, mplrg, mpsml
integer, public:: mpoud
real, private:: mpt0(mp24), mpt1(mp24), mpt2(mp24), mpt3(mp24), mpt4(mp24)

contains

  subroutine mpinit

!  MPINIT must be called at the start of execution in the user's main program.
!  It sets the numeric precision level, the MP epsilon level, the output
!  precision level, and computes the constants Pi, Log(2) and Log(10).  Note
!  that the numeric precision level MPNW is temporarily set to MPWDS + 1, in
!  order to permit extra-accurate calculations of the constants, and then is
!  reset to MPWDS.

    mpier = 0
    mpnw = mpwds + 1
    call mppi (mpt1)
    call mpdmc (2.d0, 0, mpt0)
    call mplog (mpt0, mpt2, mpt2)
    call mpdmc (10.d0, 0, mpt0)
    call mplog (mpt0, mpt2, mpt3)
    call mpnpwr (mpt0, mpiep, mpt4)
    mpnw = mpwds
    call mpeq (mpt1, mppic%mpr)
    call mpeq (mpt2, mpl02%mpr)
    call mpeq (mpt3, mpl10%mpr)
    call mpeq (mpt4, mpeps%mpr)
    mplrg%mpr(1) = 1.
    mplrg%mpr(2) = 2.e6
    mplrg%mpr(3) = 1.
    mplrg%mpr(4) = 0.
    mpsml%mpr(1) = 1.
    mpsml%mpr(2) = -2.e6
    mpsml%mpr(3) = 1.
    mpsml%mpr(4) = 0.
    mpoud = mpiou
    return
  end subroutine

  subroutine mpdexc (a, l, b)

!   This routine converts the character*1 string A, which
!   represents a multiprecision number in Fortran style, i.e.
!   '1234567890' or '1.23456789D-21', into standard MP binary format.
!   This routine is not intended to be called directly by the user.

    character*1 a(l), c(mpipl+100)
    dimension b(mpnw+4)

    do i = 1, l
      if (a(i) .eq. 'D' .or. a(i) .eq. 'E' .or. a(i) .eq. 'd' &
        .or. a(i) .eq. 'e') goto 100
    enddo

    call mpinpc (a, l, b)
    goto 110

100 i1 = i
    l1 = i - 1
    l2 = l - i
    c(1) = '1'
    c(2) = '0'
    c(3) = '^'

    do i = 1, l2
      c(i+3) = a(i+i1)
    enddo

    c(l2+4) = 'x'

    do i = 1, l1
      c(i+l2+4) = a(i)
    enddo

    call mpinpc (c, l1 + l2 + 4, b)
110 return
  end subroutine

  subroutine mpxzc (a, b)

!  This converts the DC variable A to the MPC variable B.
!  This routine is not intended to be called directly by the user.

    complex (kdb) a
    double precision da
    dimension b(mp24)
    da = a
    call mpdmc (da, 0, b)
    da = aimag (a)
    call mpdmc (da, 0, b(mp41))
    return
  end subroutine

  subroutine mpmzc (a, b)

!  This converts the MP real or MP integer variable A to the MPC variable B.
!  This routine is not intended to be called directly by the user.

    dimension a(mp4), b(mp24)
    call mpeq (a, b)
    b(mp41) = 0.
    b(mp4+2) = 0.
    return
  end subroutine

end module


module mpintmod

!  This Fortran-90 module defines operator extensions involving the
!  MP_INTEGER datatype.  For operations involving two MP data types,
!  those whose first argument is MP_INTEGER are included here.
!  Others are handled in other modules.

!  The subroutines and functions defined in this module are private
!  and not intended to be called directly by the user.

use mpfunmod
use mpdefmod
private kdb, mp4, mp24, mp41
parameter (kdb = kind (0.d0), mp4 = mpwds + 4, mp24 = 2 * mp4, mp41 = mp4 + 1)
real, private:: mpt0(mp24), mpt1(mp24), mpt2(mp24), mpt3(mp24), mpt4(mp24)
character*1, private:: az(mpipl+100)
private &
  mp_eqjj, mp_eqjq, mp_eqjz, mp_eqij, mp_eqji, mp_eqrj, mp_eqjr, &
  mp_eqcj, mp_eqjc, mp_eqdj, mp_eqjd, mp_eqxj, mp_eqjx, mp_eqja, &
  mp_addjj, mp_addjq, mp_addjz, mp_addij, mp_addji, mp_addrj, mp_addjr, &
  mp_addcj, mp_addjc, mp_adddj, mp_addjd, mp_addxj, mp_addjx, &
  mp_subjj, mp_subjq, mp_subjz, mp_subij, mp_subji, mp_subrj, mp_subjr, &
  mp_subcj, mp_subjc, mp_subdj, mp_subjd, mp_subxj, mp_subjx, mp_negj, &
  mp_muljj, mp_muljq, mp_muljz, mp_mulij, mp_mulji, mp_mulrj, mp_muljr, &
  mp_mulcj, mp_muljc, mp_muldj, mp_muljd, mp_mulxj, mp_muljx, &
  mp_divjj, mp_divjq, mp_divjz, mp_divij, mp_divji, mp_divrj, mp_divjr, &
  mp_divcj, mp_divjc, mp_divdj, mp_divjd, mp_divxj, mp_divjx, &
  mp_expjj, mp_expjq, mp_expij, mp_expji, mp_exprj, mp_expjr, &
  mp_expdj, mp_expjd, &
  mp_eqtjj, mp_eqtjq, mp_eqtjz, mp_eqtij, mp_eqtji, mp_eqtrj, mp_eqtjr, &
  mp_eqtcj, mp_eqtjc, mp_eqtdj, mp_eqtjd, mp_eqtxj, mp_eqtjx, &
  mp_netjj, mp_netjq, mp_netjz, mp_netij, mp_netji, mp_netrj, mp_netjr, &
  mp_netcj, mp_netjc, mp_netdj, mp_netjd, mp_netxj, mp_netjx, &
  mp_letjj, mp_letjq, mp_letij, mp_letji, mp_letrj, mp_letjr, &
  mp_letdj, mp_letjd, &
  mp_getjj, mp_getjq, mp_getij, mp_getji, mp_getrj, mp_getjr, &
  mp_getdj, mp_getjd, &
  mp_lttjj, mp_lttjq, mp_lttij, mp_lttji, mp_lttrj, mp_lttjr, &
  mp_lttdj, mp_lttjd, &
  mp_gttjj, mp_gttjq, mp_gttij, mp_gttji, mp_gttrj, mp_gttjr, &
  mp_gttdj, mp_gttjd

!  MPI operator extension interface blocks.

interface assignment (=)
  module procedure mp_eqjj
  module procedure mp_eqjq
  module procedure mp_eqjz
  module procedure mp_eqij
  module procedure mp_eqji
  module procedure mp_eqrj
  module procedure mp_eqjr
  module procedure mp_eqcj
  module procedure mp_eqjc
!>
  module procedure mp_eqdj
  module procedure mp_eqjd
  module procedure mp_eqxj
  module procedure mp_eqjx

  module procedure mp_eqja
end interface

interface operator (+)
  module procedure mp_addjj
  module procedure mp_addjq
  module procedure mp_addjz
  module procedure mp_addij
  module procedure mp_addji
  module procedure mp_addrj
  module procedure mp_addjr
  module procedure mp_addcj
  module procedure mp_addjc
!>
  module procedure mp_adddj
  module procedure mp_addjd
  module procedure mp_addxj
  module procedure mp_addjx
end interface

interface operator (-)
  module procedure mp_subjj
  module procedure mp_subjq
  module procedure mp_subjz
  module procedure mp_subij
  module procedure mp_subji
  module procedure mp_subrj
  module procedure mp_subjr
  module procedure mp_subcj
  module procedure mp_subjc
!>
  module procedure mp_subdj
  module procedure mp_subjd
  module procedure mp_subxj
  module procedure mp_subjx

  module procedure mp_negj
end interface

interface operator (*)
  module procedure mp_muljj
  module procedure mp_muljq
  module procedure mp_muljz
  module procedure mp_mulij
  module procedure mp_mulji
  module procedure mp_mulrj
  module procedure mp_muljr
  module procedure mp_mulcj
  module procedure mp_muljc
!>
  module procedure mp_muldj
  module procedure mp_muljd
  module procedure mp_mulxj
  module procedure mp_muljx
end interface

interface operator (/)
  module procedure mp_divjj
  module procedure mp_divjq
  module procedure mp_divjz
  module procedure mp_divij
  module procedure mp_divji
  module procedure mp_divrj
  module procedure mp_divjr
  module procedure mp_divcj
  module procedure mp_divjc
!>
  module procedure mp_divdj
  module procedure mp_divjd
  module procedure mp_divxj
  module procedure mp_divjx
end interface

interface operator (**)
  module procedure mp_expjj
  module procedure mp_expjq
  module procedure mp_expij
  module procedure mp_expji
  module procedure mp_exprj
  module procedure mp_expjr
!>
  module procedure mp_expdj
  module procedure mp_expjd
end interface

interface operator (.eq.)
  module procedure mp_eqtjj
  module procedure mp_eqtjq
  module procedure mp_eqtjz
  module procedure mp_eqtij
  module procedure mp_eqtji
  module procedure mp_eqtrj
  module procedure mp_eqtjr
  module procedure mp_eqtcj
  module procedure mp_eqtjc
!>
  module procedure mp_eqtdj
  module procedure mp_eqtjd
  module procedure mp_eqtxj
  module procedure mp_eqtjx
end interface

interface operator (.ne.)
  module procedure mp_netjj
  module procedure mp_netjq
  module procedure mp_netjz
  module procedure mp_netij
  module procedure mp_netji
  module procedure mp_netrj
  module procedure mp_netjr
  module procedure mp_netcj
  module procedure mp_netjc
!>
  module procedure mp_netdj
  module procedure mp_netjd
  module procedure mp_netxj
  module procedure mp_netjx
end interface

interface operator (.le.)
  module procedure mp_letjj
  module procedure mp_letjq
  module procedure mp_letij
  module procedure mp_letji
  module procedure mp_letrj
  module procedure mp_letjr
!>
  module procedure mp_letdj
  module procedure mp_letjd
end interface

interface operator (.ge.)
  module procedure mp_getjj
  module procedure mp_getjq
  module procedure mp_getij
  module procedure mp_getji
  module procedure mp_getrj
  module procedure mp_getjr
!>
  module procedure mp_getdj
  module procedure mp_getjd
end interface

interface operator (.lt.)
  module procedure mp_lttjj
  module procedure mp_lttjq
  module procedure mp_lttij
  module procedure mp_lttji
  module procedure mp_lttrj
  module procedure mp_lttjr
!>
  module procedure mp_lttdj
  module procedure mp_lttjd
end interface

interface operator (.gt.)
  module procedure mp_gttjj
  module procedure mp_gttjq
  module procedure mp_gttij
  module procedure mp_gttji
  module procedure mp_gttrj
  module procedure mp_gttjr
!>
  module procedure mp_gttdj
  module procedure mp_gttjd
end interface

contains

!  MPI assignment routines.

  subroutine mp_eqjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ja
    intent (in):: jb
    call mpeq (jb%mpi, ja%mpi)
    return
  end subroutine

  subroutine mp_eqjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ja
    intent (in):: qb
    call mpeq (qb%mpr, mpt1)
    call mpinfr (mpt1, ja%mpi, mpt2)
    return
  end subroutine

  subroutine mp_eqjz (ja, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ja
    intent (in):: zb
    call mpeq (zb%mpc, mpt1)
    call mpinfr (mpt1, ja%mpi, mpt2)
    return
  end subroutine

  subroutine mp_eqij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ia
    intent (in):: jb
    call mpmdc (jb%mpi, db, ib)
    ia = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ja
    intent (in):: ib
    db = ib
    call mpdmc (db, 0, ja%mpi)
    return
  end subroutine

  subroutine mp_eqrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ra
    intent (in):: jb
    call mpmdc (jb%mpi, db, ib)
    ra = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ja
    intent (in):: rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpinfr (mpt1, ja%mpi, mpt2)
    return
  end subroutine

  subroutine mp_eqcj (ca, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ca
    intent (in):: jb
    call mpmdc (jb%mpi, db, ib)
    ca = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqjc (ja, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ja
    intent (in):: cb
    db = cb
    call mpdmc (db, 0, mpt1)
    call mpinfr (mpt1, ja%mpi, mpt2)
    return
  end subroutine

  subroutine mp_eqdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: da
    intent (in):: jb
    call mpmdc (jb%mpi, db, ib)
    da = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ja
    intent (in):: db
    call mpdmc (db, 0, mpt1)
    call mpinfr (mpt1, ja%mpi, mpt2)
    return
  end subroutine

  subroutine mp_eqxj (xa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: xa
    intent (in):: jb
    call mpmdc (jb%mpi, db, ib)
    xa = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqjx (ja, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ja
    intent (in):: xb
    db = xb
    call mpdmc (db, 0, mpt1)
    call mpinfr (mpt1, ja%mpi, mpt2)
    return
  end subroutine

  subroutine mp_eqja (ja, ab)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    character*(*), intent (in):: ab
    intent (out):: ja
    l = len (ab)
    do i = 1, l
      az(i) = ab(i:i)
    enddo
    call mpdexc (az, l, mpt1)
    call mpinfr (mpt1, ja%mpi, mpt2)
    return
  end subroutine

!  MPI add routines.

  function mp_addjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_addjj
    intent (in):: ja, jb
    call mpadd (ja%mpi, jb%mpi, mpt1)
    call mpinfr (mpt1, mp_addjj%mpi, mpt2)
    return
  end function

  function mp_addjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addjq
    intent (in):: ja, qb
    call mpadd (ja%mpi, qb%mpr, mp_addjq%mpr)
    return
  end function

  function mp_addjz (ja, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addjz
    intent (in):: ja, zb
    call mpmzc (ja%mpi, mpt1)
    call mpcadd (mp4, mpt1, zb%mpc, mp_addjz%mpc)
    return
  end function

  function mp_addij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_addij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpadd (mpt1, jb%mpi, mpt2)
    call mpinfr (mpt2, mp_addij%mpi, mpt3)
    return
  end function

  function mp_addji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_addji
    intent (in):: ja, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpadd (ja%mpi, mpt1, mpt2)
    call mpinfr (mpt2, mp_addji%mpi, mpt3)
    return
  end function

  function mp_addrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpadd (mpt1, jb%mpi, mp_addrj%mpr)
    return
  end function

  function mp_addjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addjr
    intent (in):: ja, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpadd (ja%mpi, mpt1, mp_addjr%mpr)
    return
  end function

  function mp_addcj (ca, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addcj
    intent (in):: ca, jb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcadd (mp4, mpt1, mpt2, mp_addcj%mpc)
    return
  end function

  function mp_addjc (ja, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addjc
    intent (in):: ja, cb
    call mpmzc (ja%mpi, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcadd (mp4, mpt1, mpt2, mp_addjc%mpc)
    return
  end function

  function mp_adddj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_adddj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpadd (mpt1, jb%mpi, mp_adddj%mpr)
    return
  end function

  function mp_addjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addjd
    intent (in):: ja, db
    call mpdmc (db, 0, mpt1)
    call mpadd (ja%mpi, mpt1, mp_addjd%mpr)
    return
  end function

  function mp_addxj (xa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addxj
    intent (in):: xa, jb
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcadd (mp4, mpt1, mpt2, mp_addxj%mpc)
    return
  end function

  function mp_addjx (ja, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addjx
    intent (in):: ja, xb
    call mpmzc (ja%mpi, mpt1)
    call mpxzc (xb, mpt2)
    call mpcadd (mp4, mpt1, mpt2, mp_addjx%mpc)
    return
  end function

!  MPI subtract routines.

  function mp_subjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_subjj
    intent (in):: ja, jb
    call mpsub (ja%mpi, jb%mpi, mpt1)
    call mpinfr (mpt1, mp_subjj%mpi, mpt2)
    return
  end function

  function mp_subjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subjq
    intent (in):: ja, qb
    call mpsub (ja%mpi, qb%mpr, mp_subjq%mpr)
    return
  end function

  function mp_subjz (ja, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subjz
    intent (in):: ja, zb
    call mpmzc (ja%mpi, mpt1)
    call mpcsub (mp4, mpt1, zb%mpc, mp_subjz%mpc)
    return
  end function

  function mp_subij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_subij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpsub (mpt1, jb%mpi, mpt2)
    call mpinfr (mpt2, mp_subij%mpi, mpt3)
    return
  end function

  function mp_subji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_subji
    intent (in):: ja, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpsub (ja%mpi, mpt1, mpt2)
    call mpinfr (mpt2, mp_subji%mpi, mpt3)
    return
  end function

  function mp_subrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpsub (mpt1, jb%mpi, mp_subrj%mpr)
    return
  end function

  function mp_subjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subjr
    intent (in):: ja, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpsub (ja%mpi, mpt1, mp_subjr%mpr)
    return
  end function

  function mp_subcj (ca, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subcj
    intent (in):: ca, jb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcsub (mp4, mpt1, mpt2, mp_subcj%mpc)
    return
  end function

  function mp_subjc (ja, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subjc
    intent (in):: ja, cb
    call mpmzc (ja%mpi, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcsub (mp4, mpt1, mpt2, mp_subjc%mpc)
    return
  end function

  function mp_subdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpsub (mpt1, jb%mpi, mp_subdj%mpr)
    return
  end function

  function mp_subjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subjd
    intent (in):: ja, db
    call mpdmc (db, 0, mpt1)
    call mpsub (ja%mpi, mpt1, mp_subjd%mpr)
    return
  end function

  function mp_subxj (xa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subxj
    intent (in):: xa, jb
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcsub (mp4, mpt1, mpt2, mp_subxj%mpc)
    return
  end function

  function mp_subjx (ja, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subjx
    intent (in):: ja, xb
    call mpmzc (ja%mpi, mpt1)
    call mpxzc (xb, mpt2)
    call mpcsub (mp4, mpt1, mpt2, mp_subjx%mpc)
    return
  end function

!  MPI negation routine.

  function mp_negj (ja)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_negj
    intent (in):: ja
    call mpeq (ja%mpi, mp_negj%mpi)
    mp_negj%mpi(1) = - ja%mpi(1)
    return
  end function

!  MPI multiply routines.

  function mp_muljj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_muljj
    intent (in):: ja, jb
    call mpmul (ja%mpi, jb%mpi, mpt1)
    call mpinfr (mpt1, mp_muljj%mpi, mpt2)
    return
  end function

  function mp_muljq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_muljq
    intent (in):: ja, qb
    call mpmul (ja%mpi, qb%mpr, mp_muljq%mpr)
    return
  end function

  function mp_muljz (ja, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_muljz
    intent (in):: ja, zb
    call mpmzc (ja%mpi, mpt1)
    call mpcmul (mp4, mpt1, zb%mpc, mp_muljz%mpc)
    return
  end function

  function mp_mulij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_mulij
    intent (in):: ia, jb
    da = ia
    call mpmuld (jb%mpi, da, 0, mpt1)
    call mpinfr (mpt1, mp_mulij%mpi, mpt2)
    return
  end function

  function mp_mulji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_mulji
    intent (in):: ja, ib
    db = ib
    call mpmuld (ja%mpi, db, 0, mpt1)
    call mpinfr (mpt1, mp_mulji%mpi, mpt2)
    return
  end function

  function mp_mulrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_mulrj
    intent (in):: ra, jb
    da = ra
    call mpmuld (jb%mpi, da, 0, mp_mulrj%mpr)
    return
  end function

  function mp_muljr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_muljr
    intent (in):: ja, rb
    db = rb
    call mpmuld (ja%mpi, db, 0, mp_muljr%mpr)
    return
  end function

  function mp_mulcj (ca, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulcj
    intent (in):: ca, jb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcmul (mp4, mpt1, mpt2, mp_mulcj%mpc)
    return
  end function

  function mp_muljc (ja, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_muljc
    intent (in):: ja, cb
    call mpmzc (ja%mpi, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcmul (mp4, mpt1, mpt2, mp_muljc%mpc)
    return
  end function

  function mp_muldj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_muldj
    intent (in):: da, jb
    call mpmuld (jb%mpi, da, 0, mp_muldj%mpr)
    return
  end function

  function mp_muljd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_muljd
    intent (in):: ja, db
    call mpmuld (ja%mpi, db, 0, mp_muljd%mpr)
    return
  end function

  function mp_mulxj (xa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulxj
    intent (in):: xa, jb
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcmul (mp4, mpt1, mpt2, mp_mulxj%mpc)
    return
  end function

  function mp_muljx (ja, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_muljx
    intent (in):: ja, xb
    call mpmzc (ja%mpi, mpt1)
    call mpxzc (xb, mpt2)
    call mpcmul (mp4, mpt1, mpt2, mp_muljx%mpc)
    return
  end function

!  MPI divide routines.

  function mp_divjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_divjj
    intent (in):: ja, jb
    call mpdiv (ja%mpi, jb%mpi, mpt1)
    call mpinfr (mpt1, mp_divjj%mpi, mpt2)
    return
  end function

  function mp_divjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divjq
    intent (in):: ja, qb
    call mpdiv (ja%mpi, qb%mpr, mp_divjq%mpr)
    return
  end function

  function mp_divjz (ja, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divjz
    intent (in):: ja, zb
    call mpmzc (ja%mpi, mpt1)
    call mpcdiv (mp4, mpt1, zb%mpc, mp_divjz%mpc)
    return
  end function

  function mp_divij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_divij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpdiv (mpt1, jb%mpi, mpt2)
    call mpinfr (mpt2, mp_divij%mpi, mpt3)
    return
  end function

  function mp_divji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_divji
    intent (in):: ja, ib
    db = ib
    call mpdivd (ja%mpi, db, 0, mpt1)
    call mpinfr (mpt1, mp_divji%mpi, mpt2)
    return
  end function

  function mp_divrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpdiv (mpt1, jb%mpi, mp_divrj%mpr)
    return
  end function

  function mp_divjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divjr
    intent (in):: ja, rb
    db = rb
    call mpdivd (ja%mpi, db, 0, mp_divjr%mpr)
    return
  end function

  function mp_divcj (ca, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divcj
    intent (in):: ca, jb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcdiv (mp4, mpt1, mpt2, mp_divcj%mpc)
    return
  end function

  function mp_divjc (ja, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divjc
    intent (in):: ja, cb
    call mpmzc (ja%mpi, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcdiv (mp4, mpt1, mpt2, mp_divjc%mpc)
    return
  end function

  function mp_divdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpdiv (mpt1, jb%mpi, mp_divdj%mpr)
    return
  end function

  function mp_divjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divjd
    intent (in):: ja, db
    call mpdivd (ja%mpi, db, 0, mp_divjd%mpr)
    return
  end function

  function mp_divxj (xa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divxj
    intent (in):: xa, jb
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcdiv (mp4, mpt1, mpt2, mp_divxj%mpc)
    return
  end function

  function mp_divjx (ja, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divjx
    intent (in):: ja, xb
    call mpmzc (ja%mpi, mpt1)
    call mpxzc (xb, mpt2)
    call mpcdiv (mp4, mpt1, mpt2, mp_divjx%mpc)
    return
  end function

!  MPI exponentiation routines.

  function mp_expjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_expjj
    intent (in):: ja, jb
    call mplog (ja%mpi, mpl02%mpr, mpt1)
    call mpmul (mpt1, jb%mpi, mpt2)
    call mpexp (mpt2, mpl02%mpr, mpt1)
    call mpnint (mpt1, mp_expjj%mpi)
    return
  end function

  function mp_expjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expjq
    intent (in):: ja, qb
    call mplog (ja%mpi, mpl02%mpr, mpt1)
    call mpmul (mpt1, qb%mpr, mpt2)
    call mpexp (mpt2, mpl02%mpr, mp_expjq%mpr)
    return
  end function

  function mp_expij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_expij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mplog (mpt1, mpl02%mpr, mpt2)
    call mpmul (mpt2, jb%mpi, mpt3)
    call mpexp (mpt3, mpl02%mpr, mpt1)
    call mpnint (mpt1, mp_expij%mpi)
    return
  end function

  function mp_expji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_expji
    intent (in):: ja, ib
    call mpnpwr (ja%mpi, ib, mpt1)
    call mpnint (mpt1, mp_expji%mpi)
    return
  end function

  function mp_exprj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_exprj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mplog (mpt1, mpl02%mpr, mpt2)
    call mpmul (mpt2, jb%mpi, mpt3)
    call mpexp (mpt3, mpl02%mpr, mp_exprj%mpr)
    return
  end function

  function mp_expjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expjr
    intent (in):: ja, rb
    db = rb
    call mplog (ja%mpi, mpl02%mpr, mpt1)
    call mpmuld (mpt1, db, 0, mpt2)
    call mpexp (mpt2, mpl02%mpr, mp_expjr%mpr)
    return
  end function

  function mp_expdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mplog (mpt1, mpl02%mpr, mpt2)
    call mpmul (mpt2, jb%mpi, mpt3)
    call mpexp (mpt3, mpl02%mpr, mp_expdj%mpr)
    return
    end function

  function mp_expjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expjd
    intent (in):: ja, db
    call mplog (ja%mpi, mpl02%mpr, mpt1)
    call mpmuld (mpt1, db, 0, mpt2)
    call mpexp (mpt2, mpl02%mpr, mp_expjd%mpr)
    return
  end function

!  MPI .EQ. routines.

  function mp_eqtjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtjj
    intent (in):: ja, jb
    call mpcpr (ja%mpi, jb%mpi, ic)
    if (ic .eq. 0) then
      mp_eqtjj = .true.
    else
      mp_eqtjj = .false.
    endif
    return
  end function

  function mp_eqtjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtjq
    intent (in):: ja, qb
    call mpcpr (ja%mpi, qb%mpr, ic)
    if (ic .eq. 0) then
      mp_eqtjq = .true.
    else
      mp_eqtjq = .false.
    endif
    return
  end function

  function mp_eqtjz (ja, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtjz
    intent (in):: ja, zb
    call mpmzc (ja%mpi, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtjz = .true.
    else
      mp_eqtjz = .false.
    endif
    return
  end function

  function mp_eqtij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .eq. 0) then
      mp_eqtij = .true.
    else
      mp_eqtij = .false.
    endif
    return
  end function

  function mp_eqtji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtji
    intent (in):: ja, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .eq. 0) then
      mp_eqtji = .true.
    else
      mp_eqtji = .false.
    endif
    return
  end function

  function mp_eqtrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .eq. 0) then
      mp_eqtrj = .true.
    else
      mp_eqtrj = .false.
    endif
    return
  end function

  function mp_eqtjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtjr
    intent (in):: ja, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .eq. 0) then
      mp_eqtjr = .true.
    else
      mp_eqtjr = .false.
    endif
    return
  end function

  function mp_eqtcj (ca, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtcj
    intent (in):: ca, jb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtcj = .true.
    else
      mp_eqtcj = .false.
    endif
    return
  end function

  function mp_eqtjc (ja, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtjc
    intent (in):: ja, cb
    call mpmzc (ja%mpi, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtjc = .true.
    else
      mp_eqtjc = .false.
    endif
    return
  end function

  function mp_eqtdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .eq. 0) then
      mp_eqtdj = .true.
    else
      mp_eqtdj = .false.
    endif
    return
  end function

  function mp_eqtjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtjd
    intent (in):: ja, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .eq. 0) then
      mp_eqtjd = .true.
    else
      mp_eqtjd = .false.
    endif
    return
  end function

  function mp_eqtxj (xa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtxj
    intent (in):: xa, jb
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtxj = .true.
    else
      mp_eqtxj = .false.
    endif
    return
  end function

  function mp_eqtjx (ja, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtjx
    intent (in):: ja, xb
    call mpmzc (ja%mpi, mpt1)
    call mpxzc (xb, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtjx = .true.
    else
      mp_eqtjx = .false.
    endif
    return
  end function

!  MPI .NE. routines.

  function mp_netjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netjj
    intent (in):: ja, jb
    call mpcpr (ja%mpi, jb%mpi, ic)
    if (ic .ne. 0) then
      mp_netjj = .true.
    else
      mp_netjj = .false.
    endif
    return
  end function

  function mp_netjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netjq
    intent (in):: ja, qb
    call mpcpr (ja%mpi, qb%mpr, ic)
    if (ic .ne. 0) then
      mp_netjq = .true.
    else
      mp_netjq = .false.
    endif
    return
  end function

  function mp_netjz (ja, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netjz
    intent (in):: ja, zb
    call mpmzc (ja%mpi, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netjz = .true.
    else
      mp_netjz = .false.
    endif
    return
  end function

  function mp_netij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .ne. 0) then
      mp_netij = .true.
    else
      mp_netij = .false.
    endif
    return
  end function

  function mp_netji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netji
    intent (in):: ja, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .ne. 0) then
      mp_netji = .true.
    else
      mp_netji = .false.
    endif
    return
  end function

  function mp_netrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .ne. 0) then
      mp_netrj = .true.
    else
      mp_netrj = .false.
    endif
    return
  end function

  function mp_netjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netjr
    intent (in):: ja, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .ne. 0) then
      mp_netjr = .true.
    else
      mp_netjr = .false.
    endif
    return
  end function

  function mp_netcj (ca, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netcj
    intent (in):: ca, jb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netcj = .true.
    else
      mp_netcj = .false.
    endif
    return
  end function

  function mp_netjc (ja, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netjc
    intent (in):: ja, cb
    call mpmzc (ja%mpi, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netjc = .true.
    else
      mp_netjc = .false.
    endif
    return
  end function

  function mp_netdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .ne. 0) then
      mp_netdj = .true.
    else
      mp_netdj = .false.
    endif
    return
  end function

  function mp_netjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netjd
    intent (in):: ja, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .ne. 0) then
      mp_netjd = .true.
    else
      mp_netjd = .false.
    endif
    return
  end function

  function mp_netxj (xa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netxj
    intent (in):: xa, jb
    call mpxzc (xa, mpt1)
    call mpmzc (jb%mpi, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netxj = .true.
    else
      mp_netxj = .false.
    endif
    return
  end function

  function mp_netjx (ja, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netjx
    intent (in):: ja, xb
    call mpmzc (ja%mpi, mpt1)
    call mpxzc (xb, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netjx = .true.
    else
      mp_netjx = .false.
    endif
    return
  end function

!  MPI .LE. routines.

  function mp_letjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letjj
    intent (in):: ja, jb
    call mpcpr (ja%mpi, jb%mpi, ic)
    if (ic .le. 0) then
      mp_letjj = .true.
    else
      mp_letjj = .false.
    endif
    return
  end function

  function mp_letjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letjq
    intent (in):: ja, qb
    call mpcpr (ja%mpi, qb%mpr, ic)
    if (ic .le. 0) then
      mp_letjq = .true.
    else
      mp_letjq = .false.
    endif
    return
  end function

  function mp_letij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .le. 0) then
      mp_letij = .true.
    else
      mp_letij = .false.
    endif
    return
  end function

  function mp_letji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letji
    intent (in):: ja, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .le. 0) then
      mp_letji = .true.
    else
      mp_letji = .false.
    endif
    return
  end function

  function mp_letrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .le. 0) then
      mp_letrj = .true.
    else
      mp_letrj = .false.
    endif
    return
  end function

  function mp_letjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letjr
    intent (in):: ja, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .le. 0) then
      mp_letjr = .true.
    else
      mp_letjr = .false.
    endif
    return
  end function

  function mp_letdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .le. 0) then
      mp_letdj = .true.
    else
      mp_letdj = .false.
    endif
    return
  end function

  function mp_letjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letjd
    intent (in):: ja, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .le. 0) then
      mp_letjd = .true.
    else
      mp_letjd = .false.
    endif
    return
  end function

!  MPI .GE. routines.

  function mp_getjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getjj
    intent (in):: ja, jb
    call mpcpr (ja%mpi, jb%mpi, ic)
    if (ic .ge. 0) then
      mp_getjj = .true.
    else
      mp_getjj = .false.
    endif
    return
  end function

  function mp_getjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getjq
    intent (in):: ja, qb
    call mpcpr (ja%mpi, qb%mpr, ic)
    if (ic .ge. 0) then
      mp_getjq = .true.
    else
      mp_getjq = .false.
    endif
    return
  end function

  function mp_getij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .ge. 0) then
      mp_getij = .true.
    else
      mp_getij = .false.
    endif
    return
  end function

  function mp_getji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getji
    intent (in):: ja, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .ge. 0) then
      mp_getji = .true.
    else
      mp_getji = .false.
    endif
    return
  end function

  function mp_getrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .ge. 0) then
      mp_getrj = .true.
    else
      mp_getrj = .false.
    endif
    return
  end function

  function mp_getjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getjr
    intent (in):: ja, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .ge. 0) then
      mp_getjr = .true.
    else
      mp_getjr = .false.
    endif
    return
  end function

  function mp_getdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .ge. 0) then
      mp_getdj = .true.
    else
      mp_getdj = .false.
    endif
    return
  end function

  function mp_getjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getjd
    intent (in):: ja, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .ge. 0) then
      mp_getjd = .true.
    else
      mp_getjd = .false.
    endif
    return
  end function

!  MPI .LT. routines.

  function mp_lttjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttjj
    intent (in):: ja, jb
    call mpcpr (ja%mpi, jb%mpi, ic)
    if (ic .lt. 0) then
      mp_lttjj = .true.
    else
      mp_lttjj = .false.
    endif
    return
  end function

  function mp_lttjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttjq
    intent (in):: ja, qb
    call mpcpr (ja%mpi, qb%mpr, ic)
    if (ic .lt. 0) then
      mp_lttjq = .true.
    else
      mp_lttjq = .false.
    endif
    return
  end function

  function mp_lttij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .lt. 0) then
      mp_lttij = .true.
    else
      mp_lttij = .false.
    endif
    return
  end function

  function mp_lttji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttji
    intent (in):: ja, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .lt. 0) then
      mp_lttji = .true.
    else
      mp_lttji = .false.
    endif
    return
  end function

  function mp_lttrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .lt. 0) then
      mp_lttrj = .true.
    else
      mp_lttrj = .false.
    endif
    return
  end function

  function mp_lttjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttjr
    intent (in):: ja, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .lt. 0) then
      mp_lttjr = .true.
    else
      mp_lttjr = .false.
    endif
    return
  end function

  function mp_lttdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .lt. 0) then
      mp_lttdj = .true.
    else
      mp_lttdj = .false.
    endif
    return
  end function

  function mp_lttjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttjd
    intent (in):: ja, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .lt. 0) then
      mp_lttjd = .true.
    else
      mp_lttjd = .false.
    endif
    return
  end function

!  MPI .GT. routines.

  function mp_gttjj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttjj
    intent (in):: ja, jb
    call mpcpr (ja%mpi, jb%mpi, ic)
    if (ic .gt. 0) then
      mp_gttjj = .true.
    else
      mp_gttjj = .false.
    endif
    return
  end function

  function mp_gttjq (ja, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttjq
    intent (in):: ja, qb
    call mpcpr (ja%mpi, qb%mpr, ic)
    if (ic .gt. 0) then
      mp_gttjq = .true.
    else
      mp_gttjq = .false.
    endif
    return
  end function

  function mp_gttij (ia, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttij
    intent (in):: ia, jb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .gt. 0) then
      mp_gttij = .true.
    else
      mp_gttij = .false.
    endif
    return
  end function

  function mp_gttji (ja, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttji
    intent (in):: ja, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .gt. 0) then
      mp_gttji = .true.
    else
      mp_gttji = .false.
    endif
    return
  end function

  function mp_gttrj (ra, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttrj
    intent (in):: ra, jb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .gt. 0) then
      mp_gttrj = .true.
    else
      mp_gttrj = .false.
    endif
    return
  end function

  function mp_gttjr (ja, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttjr
    intent (in):: ja, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .gt. 0) then
      mp_gttjr = .true.
    else
      mp_gttjr = .false.
    endif
    return
  end function

  function mp_gttdj (da, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttdj
    intent (in):: da, jb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, jb%mpi, ic)
    if (ic .gt. 0) then
      mp_gttdj = .true.
    else
      mp_gttdj = .false.
    endif
    return
  end function

  function mp_gttjd (ja, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttjd
    intent (in):: ja, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (ja%mpi, mpt1, ic)
    if (ic .gt. 0) then
      mp_gttjd = .true.
    else
      mp_gttjd = .false.
    endif
    return
  end function

end module


module mprealmod

!  This Fortran-90 module defines operator extensions involving the
!  MP_REAL datatype.  For operations involving two MP data types,
!  those whose first argument is MP_REAL are included here.
!  Others are handled in other modules.

!  The subroutines and functions defined in this module are private
!  and not intended to be called directly by the user.

use mpfunmod
use mpdefmod
private kdb, mp4, mp24, mp41
parameter (kdb = kind (0.d0), mp4 = mpwds + 4, mp24 = 2 * mp4, mp41 = mp4 + 1)
real, private:: mpt0(mp24), mpt1(mp24), mpt2(mp24), mpt3(mp24), mpt4(mp24)
character*1, private:: az(mpipl+100)
private &
  mp_eqqj, mp_eqqq, mp_eqqz, mp_eqiq, mp_eqqi, mp_eqrq, mp_eqqr, &
  mp_eqcq, mp_eqqc, mp_eqdq, mp_eqqd, mp_eqxq, mp_eqqx, mp_eqqa, &
  mp_addqj, mp_addqq, mp_addqz, mp_addiq, mp_addqi, mp_addrq, mp_addqr, &
  mp_addcq, mp_addqc, mp_adddq, mp_addqd, mp_addxq, mp_addqx, &
  mp_subqj, mp_subqq, mp_subqz, mp_subiq, mp_subqi, mp_subrq, mp_subqr, &
  mp_subcq, mp_subqc, mp_subdq, mp_subqd, mp_subxq, mp_subqx, mp_negq, &
  mp_mulqj, mp_mulqq, mp_mulqz, mp_muliq, mp_mulqi, mp_mulrq, mp_mulqr, &
  mp_mulcq, mp_mulqc, mp_muldq, mp_mulqd, mp_mulxq, mp_mulqx, &
  mp_divqj, mp_divqq, mp_divqz, mp_diviq, mp_divqi, mp_divrq, mp_divqr, &
  mp_divcq, mp_divqc, mp_divdq, mp_divqd, mp_divxq, mp_divqx, &
  mp_expqj, mp_expqq, mp_expiq, mp_expqi, mp_exprq, mp_expqr, &
  mp_expdq, mp_expqd, &
  mp_eqtqj, mp_eqtqq, mp_eqtqz, mp_eqtiq, mp_eqtqi, mp_eqtrq, mp_eqtqr, &
  mp_eqtcq, mp_eqtqc, mp_eqtdq, mp_eqtqd, mp_eqtxq, mp_eqtqx, &
  mp_netqj, mp_netqq, mp_netqz, mp_netiq, mp_netqi, mp_netrq, mp_netqr, &
  mp_netcq, mp_netqc, mp_netdq, mp_netqd, mp_netxq, mp_netqx, &
  mp_letqj, mp_letqq, mp_letiq, mp_letqi, mp_letrq, mp_letqr, &
  mp_letdq, mp_letqd, &
  mp_getqj, mp_getqq, mp_getiq, mp_getqi, mp_getrq, mp_getqr, &
  mp_getdq, mp_getqd, &
  mp_lttqj, mp_lttqq, mp_lttiq, mp_lttqi, mp_lttrq, mp_lttqr, &
  mp_lttdq, mp_lttqd, &
  mp_gttqj, mp_gttqq, mp_gttiq, mp_gttqi, mp_gttrq, mp_gttqr, &
  mp_gttdq, mp_gttqd

!  MPR operator extension interface blocks.

interface assignment (=)
  module procedure mp_eqqj
  module procedure mp_eqqq
  module procedure mp_eqqz
  module procedure mp_eqiq
  module procedure mp_eqqi
  module procedure mp_eqrq
  module procedure mp_eqqr
  module procedure mp_eqcq
  module procedure mp_eqqc
!>
  module procedure mp_eqdq
  module procedure mp_eqqd
  module procedure mp_eqxq
  module procedure mp_eqqx

  module procedure mp_eqqa
end interface

interface operator (+)
  module procedure mp_addqj
  module procedure mp_addqq
  module procedure mp_addqz
  module procedure mp_addiq
  module procedure mp_addqi
  module procedure mp_addrq
  module procedure mp_addqr
  module procedure mp_addcq
  module procedure mp_addqc
!>
  module procedure mp_adddq
  module procedure mp_addqd
  module procedure mp_addxq
  module procedure mp_addqx
end interface

interface operator (-)
  module procedure mp_subqj
  module procedure mp_subqq
  module procedure mp_subqz
  module procedure mp_subiq
  module procedure mp_subqi
  module procedure mp_subrq
  module procedure mp_subqr
  module procedure mp_subcq
  module procedure mp_subqc
!>
  module procedure mp_subdq
  module procedure mp_subqd
  module procedure mp_subxq
  module procedure mp_subqx

  module procedure mp_negq
end interface

interface operator (*)
  module procedure mp_mulqj
  module procedure mp_mulqq
  module procedure mp_mulqz
  module procedure mp_muliq
  module procedure mp_mulqi
  module procedure mp_mulrq
  module procedure mp_mulqr
  module procedure mp_mulcq
  module procedure mp_mulqc
!>
  module procedure mp_muldq
  module procedure mp_mulqd
  module procedure mp_mulxq
  module procedure mp_mulqx
end interface

interface operator (/)
  module procedure mp_divqj
  module procedure mp_divqq
  module procedure mp_divqz
  module procedure mp_diviq
  module procedure mp_divqi
  module procedure mp_divrq
  module procedure mp_divqr
  module procedure mp_divcq
  module procedure mp_divqc
!>
  module procedure mp_divdq
  module procedure mp_divqd
  module procedure mp_divxq
  module procedure mp_divqx
end interface

interface operator (**)
  module procedure mp_expqj
  module procedure mp_expqq
  module procedure mp_expiq
  module procedure mp_expqi
  module procedure mp_exprq
  module procedure mp_expqr
!>
  module procedure mp_expdq
  module procedure mp_expqd
end interface

interface operator (.eq.)
  module procedure mp_eqtqj
  module procedure mp_eqtqq
  module procedure mp_eqtqz
  module procedure mp_eqtiq
  module procedure mp_eqtqi
  module procedure mp_eqtrq
  module procedure mp_eqtqr
  module procedure mp_eqtcq
  module procedure mp_eqtqc
!>
  module procedure mp_eqtdq
  module procedure mp_eqtqd
  module procedure mp_eqtxq
  module procedure mp_eqtqx
end interface

interface operator (.ne.)
  module procedure mp_netqj
  module procedure mp_netqq
  module procedure mp_netqz
  module procedure mp_netiq
  module procedure mp_netqi
  module procedure mp_netrq
  module procedure mp_netqr
  module procedure mp_netcq
  module procedure mp_netqc
!>
  module procedure mp_netdq
  module procedure mp_netqd
  module procedure mp_netxq
  module procedure mp_netqx
end interface

interface operator (.le.)
  module procedure mp_letqj
  module procedure mp_letqq
  module procedure mp_letiq
  module procedure mp_letqi
  module procedure mp_letrq
  module procedure mp_letqr
!>
  module procedure mp_letdq
  module procedure mp_letqd
end interface

interface operator (.ge.)
  module procedure mp_getqj
  module procedure mp_getqq
  module procedure mp_getiq
  module procedure mp_getqi
  module procedure mp_getrq
  module procedure mp_getqr
!>
  module procedure mp_getdq
  module procedure mp_getqd
end interface

interface operator (.lt.)
  module procedure mp_lttqj
  module procedure mp_lttqq
  module procedure mp_lttiq
  module procedure mp_lttqi
  module procedure mp_lttrq
  module procedure mp_lttqr
!>
  module procedure mp_lttdq
  module procedure mp_lttqd
end interface

interface operator (.gt.)
  module procedure mp_gttqj
  module procedure mp_gttqq
  module procedure mp_gttiq
  module procedure mp_gttqi
  module procedure mp_gttrq
  module procedure mp_gttqr
!>
  module procedure mp_gttdq
  module procedure mp_gttqd
end interface

contains

!  MPR assignment routines.

  subroutine mp_eqqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: qa
    intent (in):: jb
    call mpeq (jb%mpi, qa%mpr)
    return
  end subroutine

  subroutine mp_eqqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: qa
    intent (in):: qb
    call mpeq (qb%mpr, qa%mpr)
    return
  end subroutine

  subroutine mp_eqqz (qa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: qa
    intent (in):: zb
    call mpeq (zb%mpc, qa%mpr)
    return
  end subroutine

  subroutine mp_eqiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ia
    intent (in):: qb
    call mpmdc (qb%mpr, db, ib)
    ia = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: qa
    intent (in):: ib
    db = ib
    call mpdmc (db, 0, qa%mpr)
    return
  end subroutine

  subroutine mp_eqrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ra
    intent (in):: qb
    call mpmdc (qb%mpr, db, ib)
    ra = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: qa
    intent (in):: rb
    db = rb
    call mpdmc (db, 0, qa%mpr)
    return
  end subroutine

  subroutine mp_eqcq (ca, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ca
    intent (in):: qb
    call mpmdc (qb%mpr, db, ib)
    ca = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqqc (qa, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: qa
    intent (in):: cb
    db = cb
    call mpdmc (db, 0, qa%mpr)
    return
  end subroutine

  subroutine mp_eqdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: da
    intent (in):: qb
    call mpmdc (qb%mpr, db, ib)
    da = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: qa
    intent (in):: db
    call mpdmc (db, 0, qa%mpr)
    return
  end subroutine

  subroutine mp_eqxq (xa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: xa
    intent (in):: qb
    call mpmdc (qb%mpr, db, ib)
    xa = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqqx (qa, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: qa
    intent (in):: xb
    db = xb
    call mpdmc (db, 0, qa%mpr)
    return
  end subroutine

  subroutine mp_eqqa (qa, ab)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    character*(*), intent (in):: ab
    intent (out):: qa
    l = len (ab)
    do i = 1, l
      az(i) = ab(i:i)
    enddo
    call mpdexc (az, l, qa%mpr)
    return
  end subroutine

!  MPR add routines.

  function mp_addqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addqj
    intent (in):: qa, jb
    call mpadd (qa%mpr, jb%mpi, mp_addqj%mpr)
    return
  end function

  function mp_addqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addqq
    intent (in):: qa, qb
    call mpadd (qa%mpr, qb%mpr, mp_addqq%mpr)
    return
  end function

  function mp_addqz (qa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addqz
    intent (in):: qa, zb
    call mpmzc (qa%mpr, mpt1)
    call mpcadd (mp4, mpt1, zb%mpc, mp_addqz%mpc)
    return
  end function

  function mp_addiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpadd (mpt1, qb%mpr, mp_addiq%mpr)
    return
  end function

  function mp_addqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addqi
    intent (in):: qa, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpadd (qa%mpr, mpt1, mp_addqi%mpr)
    return
  end function

  function mp_addrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpadd (mpt1, qb%mpr, mp_addrq%mpr)
    return
  end function

  function mp_addqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addqr
    intent (in):: qa, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpadd (qa%mpr, mpt1, mp_addqr%mpr)
    return
  end function

  function mp_addcq (ca, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addcq
    intent (in):: ca, qb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcadd (mp4, mpt1, mpt2, mp_addcq%mpc)
    return
  end function

  function mp_addqc (qa, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addqc
    intent (in):: qa, cb
    call mpmzc (qa%mpr, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcadd (mp4, mpt1, mpt2, mp_addqc%mpc)
    return
  end function

  function mp_adddq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_adddq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpadd (mpt1, qb%mpr, mp_adddq%mpr)
    return
  end function

  function mp_addqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_addqd
    intent (in):: qa, db
    call mpdmc (db, 0, mpt1)
    call mpadd (qa%mpr, mpt1, mp_addqd%mpr)
    return
  end function

  function mp_addxq (xa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addxq
    intent (in):: xa, qb
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcadd (mp4, mpt1, mpt2, mp_addxq%mpc)
    return
  end function

  function mp_addqx (qa, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addqx
    intent (in):: qa, xb
    call mpmzc (qa%mpr, mpt1)
    call mpxzc (xb, mpt2)
    call mpcadd (mp4, mpt1, mpt2, mp_addqx%mpc)
    return
  end function

!  MPR subtract routines.

  function mp_subqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subqj
    intent (in):: qa, jb
    call mpsub (qa%mpr, jb%mpi, mp_subqj%mpr)
    return
  end function

  function mp_subqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subqq
    intent (in):: qa, qb
    call mpsub (qa%mpr, qb%mpr, mp_subqq%mpr)
    return
  end function

  function mp_subqz (qa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subqz
    intent (in):: qa, zb
    call mpmzc (qa%mpr, mpt1)
    call mpcsub (mp4, mpt1, zb%mpc, mp_subqz%mpc)
    return
  end function

  function mp_subiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpsub (mpt1, qb%mpr, mp_subiq%mpr)
    return
  end function

  function mp_subqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subqi
    intent (in):: qa, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpsub (qa%mpr, mpt1, mp_subqi%mpr)
    return
  end function

  function mp_subrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpsub (mpt1, qb%mpr, mp_subrq%mpr)
    return
  end function

  function mp_subqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subqr
    intent (in):: qa, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpsub (qa%mpr, mpt1, mp_subqr%mpr)
    return
  end function

  function mp_subcq (ca, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subcq
    intent (in):: ca, qb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcsub (mp4, mpt1, mpt2, mp_subcq%mpc)
    return
  end function

  function mp_subqc (qa, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subqc
    intent (in):: qa, cb
    call mpmzc (qa%mpr, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcsub (mp4, mpt1, mpt2, mp_subqc%mpc)
    return
  end function

  function mp_subdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpsub (mpt1, qb%mpr, mp_subdq%mpr)
    return
  end function

  function mp_subqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_subqd
    intent (in):: qa, db
    call mpdmc (db, 0, mpt1)
    call mpsub (qa%mpr, mpt1, mp_subqd%mpr)
    return
  end function

  function mp_subxq (xa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subxq
    intent (in):: xa, qb
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcsub (mp4, mpt1, mpt2, mp_subxq%mpc)
    return
  end function

  function mp_subqx (qa, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subqx
    intent (in):: qa, xb
    call mpmzc (qa%mpr, mpt1)
    call mpxzc (xb, mpt2)
    call mpcsub (mp4, mpt1, mpt2, mp_subqx%mpc)
    return
  end function

!  MPR negation routine.

  function mp_negq (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_negq
    intent (in):: qa
    call mpeq (qa%mpr, mp_negq%mpr)
    mp_negq%mpr(1) = - qa%mpr(1)
    return
  end function

!  MPR multiply routines.

  function mp_mulqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_mulqj
    intent (in):: qa, jb
    call mpmul (qa%mpr, jb%mpi, mp_mulqj%mpr)
    return
  end function

  function mp_mulqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_mulqq
    intent (in):: qa, qb
    call mpmul (qa%mpr, qb%mpr, mp_mulqq%mpr)
    return
  end function

  function mp_mulqz (qa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulqz
    intent (in):: qa, zb
    call mpmzc (qa%mpr, mpt1)
    call mpcmul (mp4, mpt1, zb%mpc, mp_mulqz%mpc)
    return
  end function

  function mp_muliq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_muliq
    intent (in):: ia, qb
    da = ia
    call mpmuld (qb%mpr, da, 0, mp_muliq%mpr)
    return
  end function

  function mp_mulqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_mulqi
    intent (in):: qa, ib
    db = ib
    call mpmuld (qa%mpr, db, 0, mp_mulqi%mpr)
    return
  end function

  function mp_mulrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_mulrq
    intent (in):: ra, qb
    da = ra
    call mpmuld (qb%mpr, da, 0, mp_mulrq%mpr)
    return
  end function

  function mp_mulqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_mulqr
    intent (in):: qa, rb
    db = rb
    call mpmuld (qa%mpr, db, 0, mp_mulqr%mpr)
    return
  end function

  function mp_mulcq (ca, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulcq
    intent (in):: ca, qb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcmul (mp4, mpt1, mpt2, mp_mulcq%mpc)
    return
  end function

  function mp_mulqc (qa, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulqc
    intent (in):: qa, cb
    call mpmzc (qa%mpr, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcmul (mp4, mpt1, mpt2, mp_mulqc%mpc)
    return
  end function

  function mp_muldq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_muldq
    intent (in):: da, qb
    call mpmuld (qb%mpr, da, 0, mp_muldq%mpr)
    return
  end function

  function mp_mulqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_mulqd
    intent (in):: qa, db
    call mpmuld (qa%mpr, db, 0, mp_mulqd%mpr)
    return
  end function

  function mp_mulxq (xa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulxq
    intent (in):: xa, qb
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcmul (mp4, mpt1, mpt2, mp_mulxq%mpc)
    return
  end function

  function mp_mulqx (qa, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulqx
    intent (in):: qa, xb
    call mpmzc (qa%mpr, mpt1)
    call mpxzc (xb, mpt2)
    call mpcmul (mp4, mpt1, mpt2, mp_mulqx%mpc)
    return
  end function

!  MPR divide routines.

  function mp_divqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divqj
    intent (in):: qa, jb
    call mpdiv (qa%mpr, jb%mpi, mp_divqj%mpr)
    return
  end function

  function mp_divqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divqq
    intent (in):: qa, qb
    call mpdiv (qa%mpr, qb%mpr, mp_divqq%mpr)
    return
  end function

  function mp_divqz (qa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divqz
    intent (in):: qa, zb
    call mpmzc (qa%mpr, mpt1)
    call mpcdiv (mp4, mpt1, zb%mpc, mp_divqz%mpc)
    return
  end function

  function mp_diviq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_diviq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpdiv (mpt1, qb%mpr, mp_diviq%mpr)
    return
  end function

  function mp_divqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divqi
    intent (in):: qa, ib
    db = ib
    call mpdivd (qa%mpr, db, 0, mp_divqi%mpr)
    return
  end function

  function mp_divrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpdiv (mpt1, qb%mpr, mp_divrq%mpr)
    return
  end function

  function mp_divqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divqr
    intent (in):: qa, rb
    db = rb
    call mpdivd (qa%mpr, db, 0, mp_divqr%mpr)
    return
  end function

  function mp_divcq (ca, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divcq
    intent (in):: ca, qb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcdiv (mp4, mpt1, mpt2, mp_divcq%mpc)
    return
  end function

  function mp_divqc (qa, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divqc
    intent (in):: qa, cb
    call mpmzc (qa%mpr, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcdiv (mp4, mpt1, mpt2, mp_divqc%mpc)
    return
  end function

  function mp_divdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpdiv (mpt1, qb%mpr, mp_divdq%mpr)
    return
  end function

  function mp_divqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_divqd
    intent (in):: qa, db
    call mpdivd (qa%mpr, db, 0, mp_divqd%mpr)
    return
  end function

  function mp_divxq (xa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divxq
    intent (in):: xa, qb
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcdiv (mp4, mpt1, mpt2, mp_divxq%mpc)
    return
  end function

  function mp_divqx (qa, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divqx
    intent (in):: qa, xb
    call mpmzc (qa%mpr, mpt1)
    call mpxzc (xb, mpt2)
    call mpcdiv (mp4, mpt1, mpt2, mp_divqx%mpc)
    return
  end function

!  MPR exponentiation routines.

  function mp_expqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expqj
    intent (in):: qa, jb
    call mplog (qa%mpr, mpl02%mpr, mpt1)
    call mpmul (mpt1, jb%mpi, mpt2)
    call mpexp (mpt2, mpl02%mpr, mp_expqj%mpr)
    return
  end function

  function mp_expqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expqq
    intent (in):: qa, qb
    call mplog (qa%mpr, mpl02%mpr, mpt1)
    call mpmul (mpt1, qb%mpr, mpt2)
    call mpexp (mpt2, mpl02%mpr, mp_expqq%mpr)
    return
  end function

  function mp_expiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mplog (mpt1, mpl02%mpr, mpt2)
    call mpmul (mpt2, qb%mpr, mpt3)
    call mpexp (mpt3, mpl02%mpr, mp_expiq%mpr)
    return
  end function

  function mp_expqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expqi
    intent (in):: qa, ib
    call mpnpwr (qa%mpr, ib, mp_expqi%mpr)
    return
  end function

  function mp_exprq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_exprq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mplog (mpt1, mpl02%mpr, mpt2)
    call mpmul (mpt2, qb%mpr, mpt3)
    call mpexp (mpt3, mpl02%mpr, mp_exprq%mpr)
    return
  end function

  function mp_expqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expqr
    intent (in):: qa, rb
    db = rb
    call mplog (qa%mpr, mpl02%mpr, mpt1)
    call mpmuld (mpt1, db, 0, mpt2)
    call mpexp (mpt2, mpl02%mpr, mp_expqr%mpr)
    return
  end function

  function mp_expdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mplog (mpt1, mpl02%mpr, mpt2)
    call mpmul (mpt2, qb%mpr, mpt3)
    call mpexp (mpt3, mpl02%mpr, mp_expdq%mpr)
    return
    end function

  function mp_expqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_expqd
    intent (in):: qa, db
    call mplog (qa%mpr, mpl02%mpr, mpt1)
    call mpmuld (mpt1, db, 0, mpt2)
    call mpexp (mpt2, mpl02%mpr, mp_expqd%mpr)
    return
  end function

!  MPR .EQ. routines.

  function mp_eqtqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtqj
    intent (in):: qa, jb
    call mpcpr (qa%mpr, jb%mpi, ic)
    if (ic .eq. 0) then
      mp_eqtqj = .true.
    else
      mp_eqtqj = .false.
    endif
    return
  end function

  function mp_eqtqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtqq
    intent (in):: qa, qb
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .eq. 0) then
      mp_eqtqq = .true.
    else
      mp_eqtqq = .false.
    endif
    return
  end function

  function mp_eqtqz (qa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtqz
    intent (in):: qa, zb
    call mpmzc (qa%mpr, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtqz = .true.
    else
      mp_eqtqz = .false.
    endif
    return
  end function

  function mp_eqtiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .eq. 0) then
      mp_eqtiq = .true.
    else
      mp_eqtiq = .false.
    endif
    return
  end function

  function mp_eqtqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtqi
    intent (in):: qa, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .eq. 0) then
      mp_eqtqi = .true.
    else
      mp_eqtqi = .false.
    endif
    return
  end function

  function mp_eqtrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .eq. 0) then
      mp_eqtrq = .true.
    else
      mp_eqtrq = .false.
    endif
    return
  end function

  function mp_eqtqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtqr
    intent (in):: qa, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .eq. 0) then
      mp_eqtqr = .true.
    else
      mp_eqtqr = .false.
    endif
    return
  end function

  function mp_eqtcq (ca, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtcq
    intent (in):: ca, qb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtcq = .true.
    else
      mp_eqtcq = .false.
    endif
    return
  end function

  function mp_eqtqc (qa, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtqc
    intent (in):: qa, cb
    call mpmzc (qa%mpr, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtqc = .true.
    else
      mp_eqtqc = .false.
    endif
    return
  end function

  function mp_eqtdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .eq. 0) then
      mp_eqtdq = .true.
    else
      mp_eqtdq = .false.
    endif
    return
  end function

  function mp_eqtqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtqd
    intent (in):: qa, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .eq. 0) then
      mp_eqtqd = .true.
    else
      mp_eqtqd = .false.
    endif
    return
  end function

  function mp_eqtxq (xa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtxq
    intent (in):: xa, qb
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtxq = .true.
    else
      mp_eqtxq = .false.
    endif
    return
  end function

  function mp_eqtqx (qa, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtqx
    intent (in):: qa, xb
    call mpmzc (qa%mpr, mpt1)
    call mpxzc (xb, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtqx = .true.
    else
      mp_eqtqx = .false.
    endif
    return
  end function

!  MPR .NE. routines.

  function mp_netqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netqj
    intent (in):: qa, jb
    call mpcpr (qa%mpr, jb%mpi, ic)
    if (ic .ne. 0) then
      mp_netqj = .true.
    else
      mp_netqj = .false.
    endif
    return
  end function

  function mp_netqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netqq
    intent (in):: qa, qb
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .ne. 0) then
      mp_netqq = .true.
    else
      mp_netqq = .false.
    endif
    return
  end function

  function mp_netqz (qa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netqz
    intent (in):: qa, zb
    call mpmzc (qa%mpr, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netqz = .true.
    else
      mp_netqz = .false.
    endif
    return
  end function

  function mp_netiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .ne. 0) then
      mp_netiq = .true.
    else
      mp_netiq = .false.
    endif
    return
  end function

  function mp_netqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netqi
    intent (in):: qa, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .ne. 0) then
      mp_netqi = .true.
    else
      mp_netqi = .false.
    endif
    return
  end function

  function mp_netrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .ne. 0) then
      mp_netrq = .true.
    else
      mp_netrq = .false.
    endif
    return
  end function

  function mp_netqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netqr
    intent (in):: qa, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .ne. 0) then
      mp_netqr = .true.
    else
      mp_netqr = .false.
    endif
    return
  end function

  function mp_netcq (ca, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netcq
    intent (in):: ca, qb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netcq = .true.
    else
      mp_netcq = .false.
    endif
    return
  end function

  function mp_netqc (qa, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netqc
    intent (in):: qa, cb
    call mpmzc (qa%mpr, mpt1)
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netqc = .true.
    else
      mp_netqc = .false.
    endif
    return
  end function

  function mp_netdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .ne. 0) then
      mp_netdq = .true.
    else
      mp_netdq = .false.
    endif
    return
  end function

  function mp_netqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netqd
    intent (in):: qa, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .ne. 0) then
      mp_netqd = .true.
    else
      mp_netqd = .false.
    endif
    return
  end function

  function mp_netxq (xa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netxq
    intent (in):: xa, qb
    call mpxzc (xa, mpt1)
    call mpmzc (qb%mpr, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netxq = .true.
    else
      mp_netxq = .false.
    endif
    return
  end function

  function mp_netqx (qa, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netqx
    intent (in):: qa, xb
    call mpmzc (qa%mpr, mpt1)
    call mpxzc (xb, mpt2)
    call mpcpr (mpt1, mpt2, ic1)
    call mpcpr (mpt1(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netqx = .true.
    else
      mp_netqx = .false.
    endif
    return
  end function

!  MPR .LE. routines.

  function mp_letqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letqj
    intent (in):: qa, jb
    call mpcpr (qa%mpr, jb%mpi, ic)
    if (ic .le. 0) then
      mp_letqj = .true.
    else
      mp_letqj = .false.
    endif
    return
  end function

  function mp_letqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letqq
    intent (in):: qa, qb
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .le. 0) then
      mp_letqq = .true.
    else
      mp_letqq = .false.
    endif
    return
  end function

  function mp_letiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .le. 0) then
      mp_letiq = .true.
    else
      mp_letiq = .false.
    endif
    return
  end function

  function mp_letqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letqi
    intent (in):: qa, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .le. 0) then
      mp_letqi = .true.
    else
      mp_letqi = .false.
    endif
    return
  end function

  function mp_letrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .le. 0) then
      mp_letrq = .true.
    else
      mp_letrq = .false.
    endif
    return
  end function

  function mp_letqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letqr
    intent (in):: qa, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .le. 0) then
      mp_letqr = .true.
    else
      mp_letqr = .false.
    endif
    return
  end function

  function mp_letdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .le. 0) then
      mp_letdq = .true.
    else
      mp_letdq = .false.
    endif
    return
  end function

  function mp_letqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_letqd
    intent (in):: qa, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .le. 0) then
      mp_letqd = .true.
    else
      mp_letqd = .false.
    endif
    return
  end function

!  MPR .GE. routines.

  function mp_getqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getqj
    intent (in):: qa, jb
    call mpcpr (qa%mpr, jb%mpi, ic)
    if (ic .ge. 0) then
      mp_getqj = .true.
    else
      mp_getqj = .false.
    endif
    return
  end function

  function mp_getqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getqq
    intent (in):: qa, qb
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .ge. 0) then
      mp_getqq = .true.
    else
      mp_getqq = .false.
    endif
    return
  end function

  function mp_getiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .ge. 0) then
      mp_getiq = .true.
    else
      mp_getiq = .false.
    endif
    return
  end function

  function mp_getqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getqi
    intent (in):: qa, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .ge. 0) then
      mp_getqi = .true.
    else
      mp_getqi = .false.
    endif
    return
  end function

  function mp_getrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .ge. 0) then
      mp_getrq = .true.
    else
      mp_getrq = .false.
    endif
    return
  end function

  function mp_getqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getqr
    intent (in):: qa, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .ge. 0) then
      mp_getqr = .true.
    else
      mp_getqr = .false.
    endif
    return
  end function

  function mp_getdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .ge. 0) then
      mp_getdq = .true.
    else
      mp_getdq = .false.
    endif
    return
  end function

  function mp_getqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_getqd
    intent (in):: qa, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .ge. 0) then
      mp_getqd = .true.
    else
      mp_getqd = .false.
    endif
    return
  end function

!  MPR .LT. routines.

  function mp_lttqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttqj
    intent (in):: qa, jb
    call mpcpr (qa%mpr, jb%mpi, ic)
    if (ic .lt. 0) then
      mp_lttqj = .true.
    else
      mp_lttqj = .false.
    endif
    return
  end function

  function mp_lttqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttqq
    intent (in):: qa, qb
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .lt. 0) then
      mp_lttqq = .true.
    else
      mp_lttqq = .false.
    endif
    return
  end function

  function mp_lttiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .lt. 0) then
      mp_lttiq = .true.
    else
      mp_lttiq = .false.
    endif
    return
  end function

  function mp_lttqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttqi
    intent (in):: qa, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .lt. 0) then
      mp_lttqi = .true.
    else
      mp_lttqi = .false.
    endif
    return
  end function

  function mp_lttrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .lt. 0) then
      mp_lttrq = .true.
    else
      mp_lttrq = .false.
    endif
    return
  end function

  function mp_lttqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttqr
    intent (in):: qa, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .lt. 0) then
      mp_lttqr = .true.
    else
      mp_lttqr = .false.
    endif
    return
  end function

  function mp_lttdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .lt. 0) then
      mp_lttdq = .true.
    else
      mp_lttdq = .false.
    endif
    return
  end function

  function mp_lttqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_lttqd
    intent (in):: qa, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .lt. 0) then
      mp_lttqd = .true.
    else
      mp_lttqd = .false.
    endif
    return
  end function

!  MPR .GT. routines.

  function mp_gttqj (qa, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttqj
    intent (in):: qa, jb
    call mpcpr (qa%mpr, jb%mpi, ic)
    if (ic .gt. 0) then
      mp_gttqj = .true.
    else
      mp_gttqj = .false.
    endif
    return
  end function

  function mp_gttqq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttqq
    intent (in):: qa, qb
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .gt. 0) then
      mp_gttqq = .true.
    else
      mp_gttqq = .false.
    endif
    return
  end function

  function mp_gttiq (ia, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttiq
    intent (in):: ia, qb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .gt. 0) then
      mp_gttiq = .true.
    else
      mp_gttiq = .false.
    endif
    return
  end function

  function mp_gttqi (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttqi
    intent (in):: qa, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .gt. 0) then
      mp_gttqi = .true.
    else
      mp_gttqi = .false.
    endif
    return
  end function

  function mp_gttrq (ra, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttrq
    intent (in):: ra, qb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .gt. 0) then
      mp_gttrq = .true.
    else
      mp_gttrq = .false.
    endif
    return
  end function

  function mp_gttqr (qa, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttqr
    intent (in):: qa, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .gt. 0) then
      mp_gttqr = .true.
    else
      mp_gttqr = .false.
    endif
    return
  end function

  function mp_gttdq (da, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttdq
    intent (in):: da, qb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, qb%mpr, ic)
    if (ic .gt. 0) then
      mp_gttdq = .true.
    else
      mp_gttdq = .false.
    endif
    return
  end function

  function mp_gttqd (qa, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_gttqd
    intent (in):: qa, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (qa%mpr, mpt1, ic)
    if (ic .gt. 0) then
      mp_gttqd = .true.
    else
      mp_gttqd = .false.
    endif
    return
  end function

end module


module mpcmpmod

!  This Fortran-90 module defines operator extensions involving the
!  MP_COMPLEX datatype.  For operations involving two MP data types,
!  those whose first argument is MP_COMPLEX are included here.
!  Others are handled in other modules.

!  The subroutines and functions defined in this module are private
!  and not intended to be called directly by the user.

use mpfunmod
use mpdefmod
private kdb, mp4, mp24, mp41
parameter (kdb = kind (0.d0), mp4 = mpwds + 4, mp24 = 2 * mp4, mp41 = mp4 + 1)
real, private:: mpt0(mp24), mpt1(mp24), mpt2(mp24), mpt3(mp24), mpt4(mp24)
character*1, private:: az(mpipl+100)
private &
  mp_eqzj, mp_eqzq, mp_eqzz, mp_eqiz, mp_eqzi, mp_eqrz, mp_eqzr, &
  mp_eqcz, mp_eqzc, mp_eqdz, mp_eqzd, mp_eqxz, mp_eqzx, &
  mp_addzj, mp_addzq, mp_addzz, mp_addiz, mp_addzi, mp_addrz, mp_addzr, &
  mp_addcz, mp_addzc, mp_adddz, mp_addzd, mp_addxz, mp_addzx, &
  mp_subzj, mp_subzq, mp_subzz, mp_subiz, mp_subzi, mp_subrz, mp_subzr, &
  mp_subcz, mp_subzc, mp_subdz, mp_subzd, mp_subxz, mp_subzx, mp_negz, &
  mp_mulzj, mp_mulzq, mp_mulzz, mp_muliz, mp_mulzi, mp_mulrz, mp_mulzr, &
  mp_mulcz, mp_mulzc, mp_muldz, mp_mulzd, mp_mulxz, mp_mulzx, &
  mp_divzj, mp_divzq, mp_divzz, mp_diviz, mp_divzi, mp_divrz, mp_divzr, &
  mp_divcz, mp_divzc, mp_divdz, mp_divzd, mp_divxz, mp_divzx, &
  mp_expzi, &
  mp_eqtzj, mp_eqtzq, mp_eqtzz, mp_eqtiz, mp_eqtzi, mp_eqtrz, mp_eqtzr, &
  mp_eqtcz, mp_eqtzc, mp_eqtdz, mp_eqtzd, mp_eqtxz, mp_eqtzx, &
  mp_netzj, mp_netzq, mp_netzz, mp_netiz, mp_netzi, mp_netrz, mp_netzr, &
  mp_netcz, mp_netzc, mp_netdz, mp_netzd, mp_netxz, mp_netzx

!  MPR operator extension interface blocks.

interface assignment (=)
  module procedure mp_eqzj
  module procedure mp_eqzq
  module procedure mp_eqzz
  module procedure mp_eqiz
  module procedure mp_eqzi
  module procedure mp_eqrz
  module procedure mp_eqzr
  module procedure mp_eqcz
  module procedure mp_eqzc
!>
  module procedure mp_eqdz
  module procedure mp_eqzd
  module procedure mp_eqxz
  module procedure mp_eqzx
end interface

interface operator (+)
  module procedure mp_addzj
  module procedure mp_addzq
  module procedure mp_addzz
  module procedure mp_addiz
  module procedure mp_addzi
  module procedure mp_addrz
  module procedure mp_addzr
  module procedure mp_addcz
  module procedure mp_addzc
!>
  module procedure mp_adddz
  module procedure mp_addzd
  module procedure mp_addxz
  module procedure mp_addzx
end interface

interface operator (-)
  module procedure mp_subzj
  module procedure mp_subzq
  module procedure mp_subzz
  module procedure mp_subiz
  module procedure mp_subzi
  module procedure mp_subrz
  module procedure mp_subzr
  module procedure mp_subcz
  module procedure mp_subzc
!>
  module procedure mp_subdz
  module procedure mp_subzd
  module procedure mp_subxz
  module procedure mp_subzx

  module procedure mp_negz
end interface

interface operator (*)
  module procedure mp_mulzj
  module procedure mp_mulzq
  module procedure mp_mulzz
  module procedure mp_muliz
  module procedure mp_mulzi
  module procedure mp_mulrz
  module procedure mp_mulzr
  module procedure mp_mulcz
  module procedure mp_mulzc
!>
  module procedure mp_muldz
  module procedure mp_mulzd
  module procedure mp_mulxz
  module procedure mp_mulzx
end interface

interface operator (/)
  module procedure mp_divzj
  module procedure mp_divzq
  module procedure mp_divzz
  module procedure mp_diviz
  module procedure mp_divzi
  module procedure mp_divrz
  module procedure mp_divzr
  module procedure mp_divcz
  module procedure mp_divzc
!>
  module procedure mp_divdz
  module procedure mp_divzd
  module procedure mp_divxz
  module procedure mp_divzx
end interface

interface operator (**)
  module procedure mp_expzi
end interface

interface operator (.eq.)
  module procedure mp_eqtzj
  module procedure mp_eqtzq
  module procedure mp_eqtzz
  module procedure mp_eqtiz
  module procedure mp_eqtzi
  module procedure mp_eqtrz
  module procedure mp_eqtzr
  module procedure mp_eqtcz
  module procedure mp_eqtzc
!>
  module procedure mp_eqtdz
  module procedure mp_eqtzd
  module procedure mp_eqtxz
  module procedure mp_eqtzx
end interface

interface operator (.ne.)
  module procedure mp_netzj
  module procedure mp_netzq
  module procedure mp_netzz
  module procedure mp_netiz
  module procedure mp_netzi
  module procedure mp_netrz
  module procedure mp_netzr
  module procedure mp_netcz
  module procedure mp_netzc
!>
  module procedure mp_netdz
  module procedure mp_netzd
  module procedure mp_netxz
  module procedure mp_netzx
end interface

contains

!  MPC assignment routines.

  subroutine mp_eqzj (za, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: za
    intent (in):: jb
    call mpmzc (jb%mpi, za%mpc)
    return
  end subroutine

  subroutine mp_eqzq (za, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: za
    intent (in):: qb
    call mpmzc (qb%mpr, za%mpc)
    return
  end subroutine

  subroutine mp_eqzz (za, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: za
    intent (in):: zb
    call mpceq (mp4, zb%mpc, za%mpc)
    return
  end subroutine

  subroutine mp_eqiz (ia, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ia
    intent (in):: zb
    call mpmdc (zb%mpc, db, ib)
    ia = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqzi (za, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: za
    intent (in):: ib
    xb = ib
    call mpxzc (xb, za%mpc)
    return
  end subroutine

  subroutine mp_eqrz (ra, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ra
    intent (in):: zb
    call mpmdc (zb%mpc, db, ib)
    ra = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqzr (za, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: za
    intent (in):: rb
    xb = rb
    call mpxzc (xb, za%mpc)
    return
  end subroutine

  subroutine mp_eqcz (ca, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: ca
    intent (in):: zb
    call mpmdc (zb%mpc, db, ib)
    call mpmdc (zb%mpc(mp41), dc, ic)
    ca = cmplx (db * 2.d0 ** ib, dc * 2.d0 ** ic)
    return
  end subroutine

  subroutine mp_eqzc (za, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: za
    intent (in):: cb
    xb = cb
    call mpxzc (xb, za%mpc)
    return
  end subroutine

  subroutine mp_eqdz (da, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: da
    intent (in):: zb
    call mpmdc (zb%mpc, db, ib)
    da = db * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqzd (za, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: za
    intent (in):: db
    xb = db
    call mpxzc (xb, za%mpc)
    return
  end subroutine

  subroutine mp_eqxz (xa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: xa
    intent (in):: zb
    call mpmdc (zb%mpc, db, ib)
    call mpmdc (zb%mpc(mp41), dc, ic)
    xa = cmplx (db * 2.d0 ** ib, dc * 2.d0 ** ic, kdb)
    return
  end subroutine

  subroutine mp_eqzx (za, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: za
    intent (in):: xb
    call mpxzc (xb, za%mpc)
    return
  end subroutine

!  MPC add routines.

  function mp_addzj (za, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addzj
    intent (in):: za, jb
    call mpmzc (jb%mpi, mpt1)
    call mpcadd (mp4, za%mpc, mpt1, mp_addzj%mpc)
    return
  end function

  function mp_addzq (za, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addzq
    intent (in):: za, qb
    call mpmzc (qb%mpr, mpt1)
    call mpcadd (mp4, za%mpc, mpt1, mp_addzq%mpc)
    return
  end function

  function mp_addzz (za, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addzz
    intent (in):: za, zb
    call mpcadd (mp4, za%mpc, zb%mpc, mp_addzz%mpc)
    return
  end function

  function mp_addiz (ia, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addiz
    intent (in):: ia, zb
    xa = ia
    call mpxzc (xa, mpt1)
    call mpcadd (mp4, mpt1, zb%mpc, mp_addiz%mpc)
    return
  end function

  function mp_addzi (za, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addzi
    intent (in):: za, ib
    xb = ib
    call mpxzc (xb, mpt1)
    call mpcadd (mp4, za%mpc, mpt1, mp_addzi%mpc)
    return
  end function

  function mp_addrz (ra, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addrz
    intent (in):: ra, zb
    xa = ra
    call mpxzc (xa, mpt1)
    call mpcadd (mp4, mpt1, zb%mpc, mp_addrz%mpc)
    return
  end function

  function mp_addzr (za, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addzr
    intent (in):: za, rb
    xb = rb
    call mpxzc (xb, mpt1)
    call mpcadd (mp4, za%mpc, mpt1, mp_addzr%mpc)
    return
  end function

  function mp_addcz (ca, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addcz
    intent (in):: ca, zb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpcadd (mp4, mpt1, zb%mpc, mp_addcz%mpc)
    return
  end function

  function mp_addzc (za, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addzc
    intent (in):: za, cb
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcadd (mp4, za%mpc, mpt2, mp_addzc%mpc)
    return
  end function

  function mp_adddz (da, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_adddz
    intent (in):: da, zb
    xa = da
    call mpxzc (xa, mpt1)
    call mpcadd (mp4, mpt1, zb%mpc, mp_adddz%mpc)
    return
  end function

  function mp_addzd (za, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addzd
    intent (in):: za, db
    xb = db
    call mpxzc (xb, mpt1)
    call mpcadd (mp4, za%mpc, mpt1, mp_addzd%mpc)
    return
  end function

  function mp_addxz (xa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addxz
    intent (in):: xa, zb
    call mpxzc (xa, mpt1)
    call mpcadd (mp4, mpt1, zb%mpc, mp_addxz%mpc)
    return
  end function

  function mp_addzx (za, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_addzx
    intent (in):: za, xb
    call mpxzc (xb, mpt2)
    call mpcadd (mp4, za%mpc, mpt2, mp_addzx%mpc)
    return
  end function

!  MPC subtract routines.

  function mp_subzj (za, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subzj
    intent (in):: za, jb
    call mpmzc (jb%mpi, mpt1)
    call mpcsub (mp4, za%mpc, mpt1, mp_subzj%mpc)
    return
  end function

  function mp_subzq (za, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subzq
    intent (in):: za, qb
    call mpmzc (qb%mpr, mpt1)
    call mpcsub (mp4, za%mpc, mpt1, mp_subzq%mpc)
    return
  end function

  function mp_subzz (za, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subzz
    intent (in):: za, zb
    call mpcsub (mp4, za%mpc, zb%mpc, mp_subzz%mpc)
    return
  end function

  function mp_subiz (ia, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subiz
    intent (in):: ia, zb
    xa = ia
    call mpxzc (xa, mpt1)
    call mpcsub (mp4, mpt1, zb%mpc, mp_subiz%mpc)
    return
  end function

  function mp_subzi (za, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subzi
    intent (in):: za, ib
    xb = ib
    call mpxzc (xb, mpt1)
    call mpcsub (mp4, za%mpc, mpt1, mp_subzi%mpc)
    return
  end function

  function mp_subrz (ra, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subrz
    intent (in):: ra, zb
    xa = ra
    call mpxzc (xa, mpt1)
    call mpcsub (mp4, mpt1, zb%mpc, mp_subrz%mpc)
    return
  end function

  function mp_subzr (za, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subzr
    intent (in):: za, rb
    xb = rb
    call mpxzc (xb, mpt1)
    call mpcsub (mp4, za%mpc, mpt1, mp_subzr%mpc)
    return
  end function

  function mp_subcz (ca, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subcz
    intent (in):: ca, zb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpcsub (mp4, mpt1, zb%mpc, mp_subcz%mpc)
    return
  end function

  function mp_subzc (za, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subzc
    intent (in):: za, cb
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcsub (mp4, za%mpc, mpt2, mp_subzc%mpc)
    return
  end function

  function mp_subdz (da, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subdz
    intent (in):: da, zb
    xa = da
    call mpxzc (xa, mpt1)
    call mpcsub (mp4, mpt1, zb%mpc, mp_subdz%mpc)
    return
  end function

  function mp_subzd (za, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subzd
    intent (in):: za, db
    xb = db
    call mpxzc (xb, mpt1)
    call mpcsub (mp4, za%mpc, mpt1, mp_subzd%mpc)
    return
  end function

  function mp_subxz (xa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subxz
    intent (in):: xa, zb
    call mpxzc (xa, mpt1)
    call mpcsub (mp4, mpt1, zb%mpc, mp_subxz%mpc)
    return
  end function

  function mp_subzx (za, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_subzx
    intent (in):: za, xb
    call mpxzc (xb, mpt2)
    call mpcsub (mp4, za%mpc, mpt2, mp_subzx%mpc)
    return
  end function

!  MPC negation routine.

  function mp_negz (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_negz
    intent (in):: za
    call mpceq (mp4, za%mpc, mp_negz%mpc)
    mp_negz%mpc(1) = - za%mpc(1)
    mp_negz%mpc(mp41) = - za%mpc(mp41)
    return
  end function

!  MPC multiply routines.

  function mp_mulzj (za, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulzj
    intent (in):: za, jb
    call mpmzc (jb%mpi, mpt1)
    call mpcmul (mp4, za%mpc, mpt1, mp_mulzj%mpc)
    return
  end function

  function mp_mulzq (za, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulzq
    intent (in):: za, qb
    call mpmzc (qb%mpr, mpt1)
    call mpcmul (mp4, za%mpc, mpt1, mp_mulzq%mpc)
    return
  end function

  function mp_mulzz (za, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulzz
    intent (in):: za, zb
    call mpcmul (mp4, za%mpc, zb%mpc, mp_mulzz%mpc)
    return
  end function

  function mp_muliz (ia, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_muliz
    intent (in):: ia, zb
    xa = ia
    call mpxzc (xa, mpt1)
    call mpcmul (mp4, mpt1, zb%mpc, mp_muliz%mpc)
    return
  end function

  function mp_mulzi (za, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulzi
    intent (in):: za, ib
    xb = ib
    call mpxzc (xb, mpt1)
    call mpcmul (mp4, za%mpc, mpt1, mp_mulzi%mpc)
    return
  end function

  function mp_mulrz (ra, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulrz
    intent (in):: ra, zb
    xa = ra
    call mpxzc (xa, mpt1)
    call mpcmul (mp4, mpt1, zb%mpc, mp_mulrz%mpc)
    return
  end function

  function mp_mulzr (za, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulzr
    intent (in):: za, rb
    xb = rb
    call mpxzc (xb, mpt1)
    call mpcmul (mp4, za%mpc, mpt1, mp_mulzr%mpc)
    return
  end function

  function mp_mulcz (ca, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulcz
    intent (in):: ca, zb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpcmul (mp4, mpt1, zb%mpc, mp_mulcz%mpc)
    return
  end function

  function mp_mulzc (za, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulzc
    intent (in):: za, cb
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcmul (mp4, za%mpc, mpt2, mp_mulzc%mpc)
    return
  end function

  function mp_muldz (da, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_muldz
    intent (in):: da, zb
    xa = da
    call mpxzc (xa, mpt1)
    call mpcmul (mp4, mpt1, zb%mpc, mp_muldz%mpc)
    return
  end function

  function mp_mulzd (za, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulzd
    intent (in):: za, db
    xb = db
    call mpxzc (xb, mpt1)
    call mpcmul (mp4, za%mpc, mpt1, mp_mulzd%mpc)
    return
  end function

  function mp_mulxz (xa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulxz
    intent (in):: xa, zb
    call mpxzc (xa, mpt1)
    call mpcmul (mp4, mpt1, zb%mpc, mp_mulxz%mpc)
    return
  end function

  function mp_mulzx (za, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_mulzx
    intent (in):: za, xb
    call mpxzc (xb, mpt2)
    call mpcmul (mp4, za%mpc, mpt2, mp_mulzx%mpc)
    return
  end function

!  MPC divide routines.

  function mp_divzj (za, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divzj
    intent (in):: za, jb
    call mpmzc (jb%mpi, mpt1)
    call mpcdiv (mp4, za%mpc, mpt1, mp_divzj%mpc)
    return
  end function

  function mp_divzq (za, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divzq
    intent (in):: za, qb
    call mpmzc (qb%mpr, mpt1)
    call mpcdiv (mp4, za%mpc, mpt1, mp_divzq%mpc)
    return
  end function

  function mp_divzz (za, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divzz
    intent (in):: za, zb
    call mpcdiv (mp4, za%mpc, zb%mpc, mp_divzz%mpc)
    return
  end function

  function mp_diviz (ia, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_diviz
    intent (in):: ia, zb
    xa = ia
    call mpxzc (xa, mpt1)
    call mpcdiv (mp4, mpt1, zb%mpc, mp_diviz%mpc)
    return
  end function

  function mp_divzi (za, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divzi
    intent (in):: za, ib
    xb = ib
    call mpxzc (xb, mpt1)
    call mpcdiv (mp4, za%mpc, mpt1, mp_divzi%mpc)
    return
  end function

  function mp_divrz (ra, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divrz
    intent (in):: ra, zb
    xa = ra
    call mpxzc (xa, mpt1)
    call mpcdiv (mp4, mpt1, zb%mpc, mp_divrz%mpc)
    return
  end function

  function mp_divzr (za, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divzr
    intent (in):: za, rb
    xb = rb
    call mpxzc (xb, mpt1)
    call mpcdiv (mp4, za%mpc, mpt1, mp_divzr%mpc)
    return
  end function

  function mp_divcz (ca, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divcz
    intent (in):: ca, zb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpcdiv (mp4, mpt1, zb%mpc, mp_divcz%mpc)
    return
  end function

  function mp_divzc (za, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divzc
    intent (in):: za, cb
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcdiv (mp4, za%mpc, mpt2, mp_divzc%mpc)
    return
  end function

  function mp_divdz (da, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divdz
    intent (in):: da, zb
    xa = da
    call mpxzc (xa, mpt1)
    call mpcdiv (mp4, mpt1, zb%mpc, mp_divdz%mpc)
    return
  end function

  function mp_divzd (za, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divzd
    intent (in):: za, db
    xb = db
    call mpxzc (xb, mpt1)
    call mpcdiv (mp4, za%mpc, mpt1, mp_divzd%mpc)
    return
  end function

  function mp_divxz (xa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divxz
    intent (in):: xa, zb
    call mpxzc (xa, mpt1)
    call mpcdiv (mp4, mpt1, zb%mpc, mp_divxz%mpc)
    return
  end function

  function mp_divzx (za, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_divzx
    intent (in):: za, xb
    call mpxzc (xb, mpt2)
    call mpcdiv (mp4, za%mpc, mpt2, mp_divzx%mpc)
    return
  end function

!  MPC exponentiation routines.

  function mp_expzi (za, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_expzi
    intent (in):: za, ib
    call mpcpwr (mp4, za%mpc, ib, mp_expzi%mpc)
    return
  end function

!  MPC .EQ. routines.

  function mp_eqtzj (za, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtzj
    intent (in):: za, jb
    call mpmzc (jb%mpi, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtzj = .true.
    else
      mp_eqtzj = .false.
    endif
    return
  end function

  function mp_eqtzq (za, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtzq
    intent (in):: za, qb
    call mpmzc (qb%mpr, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtzq = .true.
    else
      mp_eqtzq = .false.
    endif
    return
  end function

  function mp_eqtzz (za, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtzz
    intent (in):: za, zb
    call mpcpr (za%mpc, zb%mpc, ic1)
    call mpcpr (za%mpc(mp41), zb%mpc(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtzz = .true.
    else
      mp_eqtzz = .false.
    endif
    return
  end function

  function mp_eqtiz (ia, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtiz
    intent (in):: ia, zb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtiz = .true.
    else
      mp_eqtiz = .false.
    endif
    return
  end function

  function mp_eqtzi (za, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtzi
    intent (in):: za, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtzi = .true.
    else
      mp_eqtzi = .false.
    endif
    return
  end function

  function mp_eqtrz (ra, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtrz
    intent (in):: ra, zb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtrz = .true.
    else
      mp_eqtrz = .false.
    endif
    return
  end function

  function mp_eqtzr (za, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtzr
    intent (in):: za, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtzr = .true.
    else
      mp_eqtzr = .false.
    endif
    return
  end function

  function mp_eqtcz (ca, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtcz
    intent (in):: ca, zb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtcz = .true.
    else
      mp_eqtcz = .false.
    endif
    return
  end function

  function mp_eqtzc (za, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtzc
    intent (in):: za, cb
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcpr (za%mpc, mpt2, ic1)
    call mpcpr (za%mpc(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtzc = .true.
    else
      mp_eqtzc = .false.
    endif
    return
  end function

  function mp_eqtdz (da, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtdz
    intent (in):: da, zb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtdz = .true.
    else
      mp_eqtdz = .false.
    endif
    return
  end function

  function mp_eqtzd (za, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtzd
    intent (in):: za, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtzd = .true.
    else
      mp_eqtzd = .false.
    endif
    return
  end function

  function mp_eqtxz (xa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtxz
    intent (in):: xa, zb
    call mpxzc (xa, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtxz = .true.
    else
      mp_eqtxz = .false.
    endif
    return
  end function

  function mp_eqtzx (za, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_eqtzx
    intent (in):: za, xb
    call mpxzc (xb, mpt2)
    call mpcpr (za%mpc, mpt2, ic1)
    call mpcpr (za%mpc(mp41), mpt2(mp41), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      mp_eqtzx = .true.
    else
      mp_eqtzx = .false.
    endif
    return
  end function

!  MPC .NE. routines.

  function mp_netzj (za, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netzj
    intent (in):: za, jb
    call mpmzc (jb%mpi, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netzj = .true.
    else
      mp_netzj = .false.
    endif
    return
  end function

  function mp_netzq (za, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netzq
    intent (in):: za, qb
    call mpmzc (qb%mpr, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netzq = .true.
    else
      mp_netzq = .false.
    endif
    return
  end function

  function mp_netzz (za, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netzz
    intent (in):: za, zb
    call mpcpr (za%mpc, zb%mpc, ic1)
    call mpcpr (za%mpc(mp41), zb%mpc(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netzz = .true.
    else
      mp_netzz = .false.
    endif
    return
  end function

  function mp_netiz (ia, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netiz
    intent (in):: ia, zb
    da = ia
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netiz = .true.
    else
      mp_netiz = .false.
    endif
    return
  end function

  function mp_netzi (za, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netzi
    intent (in):: za, ib
    db = ib
    call mpdmc (db, 0, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netzi = .true.
    else
      mp_netzi = .false.
    endif
    return
  end function

  function mp_netrz (ra, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netrz
    intent (in):: ra, zb
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netrz = .true.
    else
      mp_netrz = .false.
    endif
    return
  end function

  function mp_netzr (za, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netzr
    intent (in):: za, rb
    db = rb
    call mpdmc (db, 0, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netzr = .true.
    else
      mp_netzr = .false.
    endif
    return
  end function

  function mp_netcz (ca, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netcz
    intent (in):: ca, zb
    xa = ca
    call mpxzc (xa, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netcz = .true.
    else
      mp_netcz = .false.
    endif
    return
  end function

  function mp_netzc (za, cb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netzc
    intent (in):: za, cb
    xb = cb
    call mpxzc (xb, mpt2)
    call mpcpr (za%mpc, mpt2, ic1)
    call mpcpr (za%mpc(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netzc = .true.
    else
      mp_netzc = .false.
    endif
    return
  end function

  function mp_netdz (da, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netdz
    intent (in):: da, zb
    call mpdmc (da, 0, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netdz = .true.
    else
      mp_netdz = .false.
    endif
    return
  end function

  function mp_netzd (za, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netzd
    intent (in):: za, db
    call mpdmc (db, 0, mpt1)
    call mpcpr (za%mpc, mpt1, ic1)
    call mpcpr (za%mpc(mp41), mpt1(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netzd = .true.
    else
      mp_netzd = .false.
    endif
    return
  end function

  function mp_netxz (xa, zb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netxz
    intent (in):: xa, zb
    call mpxzc (xa, mpt1)
    call mpcpr (mpt1, zb%mpc, ic1)
    call mpcpr (mpt1(mp41), zb%mpc(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netxz = .true.
    else
      mp_netxz = .false.
    endif
    return
  end function

  function mp_netzx (za, xb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    logical mp_netzx
    intent (in):: za, xb
    call mpxzc (xb, mpt2)
    call mpcpr (za%mpc, mpt2, ic1)
    call mpcpr (za%mpc(mp41), mpt2(mp41), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      mp_netzx = .true.
    else
      mp_netzx = .false.
    endif
    return
  end function

end module


module mpgenmod

!  This Fortran-90 module defines generic functions involving all
!  MP datatypes.

!  The subroutines and functions defined in this module are private
!  and not intended to be called directly by the user.  The generic
!  names (i.e. interface block names) are publicly accessible, though.

use mpfunmod
use mpdefmod
private kdb, mp4, mp24, mp41
parameter (kdb = kind (0.d0), mp4 = mpwds + 4, mp24 = 2 * mp4, mp41 = mp4 + 1)
real, private:: mpt0(mp24), mpt1(mp24), mpt2(mp24), mpt3(mp24), mpt4(mp24), &
  mpt5(mp24), mpt6(mp24)
character*1, private:: az(mpipl+100)
private &
  mp_absj, mp_absq, mp_absz, mp_acos, mp_imag, mp_aint, mp_anint, &
  mp_asin, mp_atan, mp_atan2, mp_jtoc, mp_qtoc, mp_ztoc, mp_conjg, &
  mp_cos, mp_cosz, mp_cosh, mp_jtod, mp_qtod, mp_ztod, mp_jtox, mp_qtox, &
  mp_ztox, mp_exp, mp_expz, mp_jtoi, mp_qtoi, mp_ztoi, mp_log, mp_logz, &
  mp_log10, mp_maxj, mp_maxq, mp_maxq3, mp_minj, mp_minq, mp_minq3, mp_modj, &
  mp_modq, mp_jtoz, mp_qtoz, mp_itoz, mp_rtoz, mp_ctoz, mp_dtoz, mp_xtoz, &
  mp_atoz, mp_jjtoz, mp_qqtoz, mp_iitoz, mp_rrtoz, mp_ddtoz, mp_aatoz, &
  mp_cssh, mp_cssn, mp_qtoj, mp_ztoj, mp_itoj, mp_rtoj, mp_ctoj, mp_dtoj, &
  mp_xtoj, mp_atoj, mp_nrt, mp_rand, mp_inpj, mp_inpq, mp_inpz, &
  mp_jtoq, mp_ztoq, mp_itoq, mp_rtoq, mp_ctoq, mp_dtoq, mp_xtoq, &
  mp_atoq, mp_outj, mp_outq, mp_outz, mp_nint, mp_jtor, mp_qtor, &
  mp_ztor, mp_signj, mp_signq, mp_sin, mp_sinz, mp_sinh, mp_sqrtq, &
  mp_sqrtz, mp_tan, mp_tanh

!  MP generic interface blocks.

interface abs
  module procedure mp_absj
  module procedure mp_absq
  module procedure mp_absz
end interface

interface acos
  module procedure mp_acos
end interface

interface aimag
  module procedure mp_imag
end interface

interface aint
  module procedure mp_aint
end interface

interface anint
  module procedure mp_anint
end interface

interface asin
  module procedure mp_asin
end interface

interface atan
  module procedure mp_atan
end interface

interface atan2
  module procedure mp_atan2
end interface

interface cmplx
  module procedure mp_jtoc
  module procedure mp_qtoc
  module procedure mp_ztoc
end interface

interface conjg
  module procedure mp_conjg
end interface

interface cos
  module procedure mp_cos
  module procedure mp_cosz
end interface

interface cosh
  module procedure mp_cosh
end interface

interface dble
  module procedure mp_jtod
  module procedure mp_qtod
  module procedure mp_ztod
end interface

interface dcmplx
  module procedure mp_jtox
  module procedure mp_qtox
  module procedure mp_ztox
end interface

interface exp
  module procedure mp_exp
  module procedure mp_expz
end interface

interface int
  module procedure mp_jtoi
  module procedure mp_qtoi
  module procedure mp_ztoi
end interface

interface log
  module procedure mp_log
  module procedure mp_logz
end interface

interface log10
  module procedure mp_log10
end interface

interface max
  module procedure mp_maxj
  module procedure mp_maxq
  module procedure mp_maxq3
end interface

interface min
  module procedure mp_minj
  module procedure mp_minq
  module procedure mp_minq3
end interface

interface mod
  module procedure mp_modj
  module procedure mp_modq
end interface

interface mpcmpl
  module procedure mp_jtoz
  module procedure mp_qtoz
  module procedure mp_itoz
  module procedure mp_rtoz
  module procedure mp_ctoz
!>
  module procedure mp_dtoz
  module procedure mp_xtoz

  module procedure mp_atoz

  module procedure mp_jjtoz
  module procedure mp_qqtoz
  module procedure mp_iitoz
  module procedure mp_rrtoz
!>
  module procedure mp_ddtoz

  module procedure mp_aatoz
end interface

interface mpcsshf
  module procedure mp_cssh
end interface

interface mpcssnf
  module procedure mp_cssn
end interface

interface mpint
  module procedure mp_qtoj
  module procedure mp_ztoj
  module procedure mp_itoj
  module procedure mp_rtoj
  module procedure mp_ctoj
!>
  module procedure mp_dtoj
  module procedure mp_xtoj

  module procedure mp_atoj
end interface

interface mpnrtf
  module procedure mp_nrt
end interface

interface mpranf
  module procedure mp_rand
end interface

interface mpread
  module procedure mp_inpj
  module procedure mp_inpq
  module procedure mp_inpz
end interface

interface mpreal
  module procedure mp_jtoq
  module procedure mp_ztoq
  module procedure mp_itoq
  module procedure mp_rtoq
  module procedure mp_ctoq
!>
  module procedure mp_dtoq
  module procedure mp_xtoq

  module procedure mp_atoq
end interface

interface mpwrite
  module procedure mp_outj
  module procedure mp_outq
  module procedure mp_outz
end interface

interface nint
  module procedure mp_nint
end interface

interface real
  module procedure mp_jtor
  module procedure mp_qtor
  module procedure mp_ztor
end interface

interface sign
  module procedure mp_signj
  module procedure mp_signq
end interface

interface sin
  module procedure mp_sin
  module procedure mp_sinz
end interface

interface sinh
  module procedure mp_sinh
end interface

interface sqrt
  module procedure mp_sqrtq
  module procedure mp_sqrtz
end interface

interface tan
  module procedure mp_tan
end interface

interface tanh
  module procedure mp_tanh
end interface

contains

  function mp_absj (ja)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_absj
    intent (in):: ja
    call mpeq (ja%mpi, mp_absj%mpi)
    mp_absj%mpi(1) = abs (ja%mpi(1))
    return
  end function

  function mp_absq (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_absq
    intent (in):: qa
    call mpeq (qa%mpr, mp_absq%mpr)
    mp_absq%mpr(1) = abs (qa%mpr(1))
    return
  end function

  function mp_absz (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_absz
    intent (in):: za
    call mpmul (za%mpc, za%mpc, mpt1)
    call mpmul (za%mpc(mp41), za%mpc(mp41), mpt2)
    call mpadd (mpt1, mpt2, mpt3)
    call mpsqrt (mpt3, mp_absz%mpr)
    return
  end function

  function mp_acos (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_acos
    intent (in):: qa
    call mpdmc (1.d0, 0, mpt1)
    call mpmul (qa%mpr, qa%mpr, mpt2)
    call mpsub (mpt1, mpt2, mpt3)
    call mpsqrt (mpt3, mpt1)
    call mpang (qa%mpr, mpt1, mppic%mpr, mp_acos%mpr)
    return
  end function

  function mp_aint (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_aint
    intent (in):: qa
    call mpinfr (qa%mpr, mp_aint%mpr, mpt1)
    return
  end function

  function mp_anint (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_anint
    intent (in):: qa
    call mpnint (qa%mpr, mp_anint%mpr)
    return
  end function

  function mp_asin (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_asin
    intent (in):: qa
    call mpdmc (1.d0, 0, mpt1)
    call mpmul (qa%mpr, qa%mpr, mpt2)
    call mpsub (mpt1, mpt2, mpt3)
    call mpsqrt (mpt3, mpt1)
    call mpang (mpt1, qa%mpr, mppic%mpr, mp_asin%mpr)
    return
  end function

  function mp_atan (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_atan
    intent (in):: qa
    call mpdmc (1.d0, 0, mpt1)
    call mpang (mpt1, qa%mpr, mppic%mpr, mp_atan%mpr)
    return
  end function

  function mp_atan2 (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_atan2
    intent (in):: qa, qb
    call mpang (qb%mpr, qa%mpr, mppic%mpr, mp_atan2%mpr)
    return
  end function

  function mp_jtoc (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    complex:: mp_jtoc
    intent (in):: ja, jb
    call mpmdc (ja%mpi, da, ia)
    call mpmdc (jb%mpi, db, ib)
    mp_jtoc = cmplx (da * 2.d0 ** ia, db * 2.d0 ** ib)
    return
  end function

  function mp_qtoc (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    complex:: mp_qtoc
    intent (in):: qa, qb
    call mpmdc (qa%mpr, da, ia)
    call mpmdc (qb%mpr, db, ib)
    mp_qtoc = cmplx (da * 2.d0 ** ia, db * 2.d0 ** ib)
    return
  end function

  function mp_ztoc (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    complex:: mp_ztoc
    intent (in):: za
    call mpmdc (za%mpc, da, ia)
    call mpmdc (za%mpc(mp41), db, ib)
    mp_ztoc = cmplx (da * 2.d0 ** ia, db * 2.d0 ** ib)
    return
  end function

  function mp_conjg (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_conjg
    intent (in):: za
    call mpceq (mp4, za%mpc, mp_conjg%mpc)
    mp_conjg%mpc(mp41) = - za%mpc(mp41)
    return
  end function

  function mp_cos (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_cos
    intent (in):: qa
    call mpcssn (qa%mpr, mppic%mpr, mp_cos%mpr, mpt1)
    return
  end function

  function mp_cosz (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_cosz
    intent (in):: za
    call mpeq (za%mpc(mp41), mpt2)
    mpt2(1) = - mpt2(1)
    call mpexp (mpt2, mpl02%mpr, mpt1)
    call mpdmc (1.d0, 0, mpt3)
    call mpdiv (mpt3, mpt1, mpt2)
    call mpcssn (za%mpc, mppic%mpr, mpt3, mpt4)
    call mpadd (mpt1, mpt2, mpt5)
    call mpmuld (mpt5, 0.5d0, 0, mpt6)
    call mpmul (mpt6, mpt3, mp_cosz%mpc)
    call mpsub (mpt1, mpt2, mpt5)
    call mpmuld (mpt5, 0.5d0, 0, mpt6)
    call mpmul (mpt6, mpt4, mp_cosz%mpc(mp41))
    return
  end function

  function mp_cosh (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_cosh
    intent (in):: qa
    call mpcssh (qa%mpr, mpl02%mpr, mp_cosh%mpr, mpt1)
    return
  end function

  function mp_jtod (ja)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: ja
    double precision mp_jtod
    call mpmdc (ja%mpi, da, ia)
    mp_jtod = da * 2.d0 ** ia
    return
  end function

  function mp_qtod (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: qa
    double precision:: mp_qtod, da
    call mpmdc (qa%mpr, da, ia)
    mp_qtod = da * 2.d0 ** ia
    return
  end function

  function mp_ztod (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: za
    double precision:: mp_ztod, da
    call mpmdc (za%mpc, da, ia)
    mp_ztod = da * 2.d0 ** ia
    return
  end function

  function mp_jtox (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    complex (kdb):: mp_jtox
    intent (in):: ja, jb
    call mpmdc (ja%mpi, da, ia)
    call mpmdc (jb%mpi, db, ib)
    mp_jtox = cmplx (da * 2.d0 ** ia, db * 2.d0 ** ib, kdb)
    return
  end function

  function mp_qtox (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    complex (kdb):: mp_qtox
    intent (in):: qa, qb
    call mpmdc (qa%mpr, da, ia)
    call mpmdc (qb%mpr, db, ib)
    mp_qtox = cmplx (da * 2.d0 ** ia, db * 2.d0 ** ib, kdb)
    return
  end function

  function mp_ztox (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    complex (kdb):: mp_ztox
    intent (in):: za
    call mpmdc (za%mpc, da, ia)
    call mpmdc (za%mpc(mp41), db, ib)
    mp_ztox = cmplx (da * 2.d0 ** ia, db * 2.d0 ** ib, kdb)
    return
  end function

  function mp_exp (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_exp
    intent (in):: qa
    call mpexp (qa%mpr, mpl02%mpr, mp_exp%mpr)
    return
  end function

  function mp_expz (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_expz
    intent (in):: za
    call mpexp (za%mpc, mpl02%mpr, mpt1)
    call mpcssn (za%mpc(mp41), mppic%mpr, mpt2, mpt3)
    call mpmul (mpt1, mpt2, mp_expz%mpc)
    call mpmul (mpt1, mpt3, mp_expz%mpc(mp41))
    return
  end function

  function mp_imag (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_imag
    intent (in):: za
    call mpeq (za%mpc(mp41), mp_imag%mpr)
    return
  end function

  function mp_jtoi (ja)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    integer:: mp_jtoi
    intent (in):: ja
    call mpmdc (ja%mpi, da, ia)
    mp_jtoi = da * 2.d0 ** ia
    return
  end function

  function mp_qtoi (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    integer:: mp_qtoi
    intent (in):: qa
    call mpmdc (qa%mpr, da, ia)
    mp_qtoi = da * 2.d0 ** ia
    return
  end function

  function mp_ztoi (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    integer:: mp_ztoi
    intent (in):: za
    call mpmdc (za%mpc, da, ia)
    mp_ztoi = da * 2.d0 ** ia
    return
  end function

  function mp_log (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_log
    intent (in):: qa
    call mplog (qa%mpr, mpl02%mpr, mp_log%mpr)
    return
  end function

  function mp_logz (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_logz
    intent (in):: za
    call mpmul (za%mpc, za%mpc, mpt1)
    call mpmul (za%mpc(mp41), za%mpc(mp41), mpt2)
    call mpadd (mpt1, mpt2, mpt3)
    call mplog (mpt3, mpl02%mpr, mpt4)
    call mpmuld (mpt4, 0.5d0, 0, mp_logz%mpc)
    call mpang (za%mpc, za%mpc(mp41), mppic%mpr, mp_logz%mpc(mp41))
    return
  end function

  function mp_log10 (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_log10
    intent (in):: qa
    call mplog (qa%mpr, mpl02%mpr, mpt1)
    call mpdiv (mpt1, mpl10%mpr, mp_log10%mpr)
    return
  end function

  function mp_maxj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_maxj
    intent (in):: ja, jb
    call mpcpr (ja%mpi, jb%mpi, ic)
    if (ic .ge. 0) then
      call mpeq (ja%mpi, mp_maxj%mpi)
    else
      call mpeq (jb%mpi, mp_maxj%mpi)
    endif
    return
  end function

  function mp_maxq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_maxq
    intent (in):: qa, qb
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .ge. 0) then
      call mpeq (qa%mpr, mp_maxq%mpr)
    else
      call mpeq (qb%mpr, mp_maxq%mpr)
    endif
    return
  end function

  function mp_maxq3 (qa, qb, qc)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_maxq3
    intent (in):: qa, qb, qc
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .ge. 0) then
      call mpeq (qa%mpr, mpt0)
    else
      call mpeq (qb%mpr, mpt0)
    endif
    call mpcpr (mpt0, qc%mpr, ic)
    if (ic .ge. 0) then
      call mpeq (mpt0, mp_maxq3%mpr)
    else
      call mpeq (qc%mpr, mp_maxq3%mpr)
    endif
    return
  end function

  function mp_minj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_minj
    intent (in):: ja, jb
    call mpcpr (ja%mpi, jb%mpi, ic)
    if (ic .lt. 0) then
      call mpeq (ja%mpi, mp_minj%mpi)
    else
      call mpeq (jb%mpi, mp_minj%mpi)
    endif
    return
  end function

  function mp_minq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_minq
    intent (in):: qa, qb
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .lt. 0) then
      call mpeq (qa%mpr, mp_minq%mpr)
    else
      call mpeq (qb%mpr, mp_minq%mpr)
    endif
    return
  end function

  function mp_minq3 (qa, qb, qc)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_minq3
    intent (in):: qa, qb, qc
    call mpcpr (qa%mpr, qb%mpr, ic)
    if (ic .lt. 0) then
      call mpeq (qa%mpr, mpt0)
    else
      call mpeq (qb%mpr, mpt0)
    endif
    call mpcpr (mpt0, qc%mpr, ic)
    if (ic .lt. 0) then
      call mpeq (mpt0, mp_minq3%mpr)
    else
      call mpeq (qc%mpr, mp_minq3%mpr)
    endif
    return
  end function

  function mp_modj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_modj
    intent (in):: ja, jb
    call mpdiv (ja%mpi, jb%mpi, mpt1)
    call mpinfr (mpt1, mpt2, mpt3)
    call mpmul (jb%mpi, mpt2, mpt1)
    call mpsub (ja%mpi, mpt1, mp_modj%mpi)
    return
  end function

  function mp_modq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_modq
    intent (in):: qa, qb
    call mpdiv (qa%mpr, qb%mpr, mpt1)
    call mpinfr (mpt1, mpt2, mpt3)
    call mpmul (qb%mpr, mpt2, mpt1)
    call mpsub (qa%mpr, mpt1, mp_modq%mpr)
    return
  end function

  function mp_jtoz (ja)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_jtoz
    intent (in):: ja
    call mpmzc (ja%mpi, mp_jtoz%mpc)
    return
  end function

  function mp_qtoz (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_qtoz
    intent (in):: qa
    call mpmzc (qa%mpr, mp_qtoz%mpc)
    return
  end function

  function mp_itoz (ia)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_itoz
    intent (in):: ia
    xa = ia
    call mpxzc (xa, mp_itoz%mpc)
    return
  end function

  function mp_rtoz (ra)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_rtoz
    intent (in):: ra
    xa = ra
    call mpxzc (xa, mp_rtoz%mpc)
    return
  end function

  function mp_ctoz (ca)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_ctoz
    intent (in):: ca
    xa = ca
    call mpxzc (xa, mp_ctoz%mpc)
    return
  end function

  function mp_dtoz (da)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_dtoz
    intent (in):: da
    xa = da
    call mpxzc (xa, mp_dtoz%mpc)
    return
  end function

  function mp_xtoz (xa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_xtoz
    intent (in):: xa
    call mpxzc (xa, mp_xtoz%mpc)
    return
  end function

  function mp_atoz (aa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    character*(*), intent (in):: aa
    type (mp_complex):: mp_atoz
    l = len (aa)
    do i = 1, l
      az(i) = aa(i:i)
    enddo
    call mpinpc (az, l, mpt1)
    call mpmzc (mpt1, mp_atoz%mpc)
    return
  end function

  function mp_jjtoz (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_jjtoz
    intent (in):: ja, jb
    call mpmmpc (ja%mpi, jb%mpi, mp4, mp_jjtoz%mpc)
    return
  end function

  function mp_qqtoz (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_qqtoz
    intent (in):: qa, qb
    call mpmmpc (qa%mpr, qb%mpr, mp4, mp_qqtoz%mpc)
    return
  end function

  function mp_iitoz (ia, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_iitoz
    intent (in):: ia, ib
    xa = cmplx (ia, ib, kdb)
    call mpxzc (xa, mp_iitoz%mpc)
    return
  end function

  function mp_rrtoz (ra, rb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_rrtoz
    intent (in):: ra, rb
    xa = cmplx (ra, rb, kdb)
    call mpxzc (xa, mp_rrtoz%mpc)
    return
  end function

  function mp_ddtoz (da, db)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_ddtoz
    intent (in):: da, db
    xa = cmplx (da, db, kdb)
    call mpxzc (xa, mp_ddtoz%mpc)
    return
  end function

  function mp_aatoz (aa, ab)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    character*(*), intent (in):: aa, ab
    type (mp_complex):: mp_aatoz
    l = len (aa)
    do i = 1, l
      az(i) = aa(i:i)
    enddo
    call mpinpc (az, l, mp_aatoz%mpc)
    l = len (ab)
    do i = 1, l
      az(i) = ab(i:i)
    enddo
    call mpinpc (az, l, mp_aatoz%mpc(mp41))
    return
  end function

  subroutine mp_cssh (qa, qb, qc)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: qa
    intent (out):: qb, qc
    call mpcssh (qa%mpr, mpl02%mpr, qb%mpr, qc%mpr)
    return
  end subroutine

  subroutine mp_cssn (qa, qb, qc)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: qa
    intent (out):: qb, qc
    call mpcssn (qa%mpr, mppic%mpr, qb%mpr, qc%mpr)
    return
  end subroutine

  function mp_qtoj (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_qtoj
    intent (in):: qa
    call mpeq (qa%mpr, mpt1)
    call mpinfr (mpt1, mp_qtoj%mpi, mpt2)
    return
  end function

  function mp_ztoj (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_ztoj
    intent (in):: za
    call mpeq (za%mpc, mpt1)
    call mpinfr (mpt1, mp_ztoj%mpi, mpt2)
    return
  end function

  function mp_itoj (ia)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_itoj
    intent (in):: ia
    da = ia
    call mpdmc (da, 0, mp_itoj%mpi)
    return
  end function

  function mp_rtoj (ra)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_rtoj
    intent (in):: ra
    da = ra
    call mpdmc (da, 0, mpt1)
    call mpinfr (mpt1, mp_rtoj%mpi, mpt2)
    return
  end function

  function mp_ctoj (ca)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_ctoj
    intent (in):: ca
    da = ca
    call mpdmc (da, 0, mpt1)
    call mpinfr (mpt1, mp_ctoj%mpi, mpt2)
    return
  end function

  function mp_dtoj (da)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_dtoj
    intent (in):: da
    call mpdmc (da, 0, mpt1)
    call mpinfr (mpt1, mp_dtoj%mpi, mpt2)
    return
  end function

  function mp_xtoj (xa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_xtoj
    intent (in):: xa
    da = xa
    call mpdmc (da, 0, mpt1)
    call mpinfr (mpt1, mp_xtoj%mpi, mpt2)
    return
  end function

  function mp_atoj (aa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    character*(*), intent (in):: aa
    type (mp_integer):: mp_atoj
    l = len (aa)
    do i = 1, l
      az(i) = aa(i:i)
    enddo
    call mpinpc (az, l, mpt1)
    call mpinfr (mpt1, mp_atoj%mpi, mpt2)
    return
  end function

  function mp_nrt (qa, ib)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_nrt
    intent (in):: qa, ib
    call mpnrt (qa%mpr, ib, mp_nrt%mpr)
    return
  end function

  function mp_rand ()
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_rand
    call mprand (mp_rand%mpr)
    return
  end function

  subroutine mp_inpj (iu, j1, j2, j3, j4, j5, j6, j7, j8, j9)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: j1, j2, j3, j4, j5, j6, j7, j8, j9
    optional:: j2, j3, j4, j5, j6, j7, j8, j9
    call mpinp (iu, j1%mpi, az)
    if (present (j2)) call mpinp (iu, j2%mpi, az)
    if (present (j3)) call mpinp (iu, j3%mpi, az)
    if (present (j4)) call mpinp (iu, j4%mpi, az)
    if (present (j5)) call mpinp (iu, j5%mpi, az)
    if (present (j6)) call mpinp (iu, j6%mpi, az)
    if (present (j7)) call mpinp (iu, j7%mpi, az)
    if (present (j8)) call mpinp (iu, j8%mpi, az)
    if (present (j9)) call mpinp (iu, j9%mpi, az)
    return
  end subroutine

  subroutine mp_inpq (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: q1, q2, q3, q4, q5, q6, q7, q8, q9
    optional:: q2, q3, q4, q5, q6, q7, q8, q9
    call mpinp (iu, q1%mpr, az)
    if (present (q2)) call mpinp (iu, q2%mpr, az)
    if (present (q3)) call mpinp (iu, q3%mpr, az)
    if (present (q4)) call mpinp (iu, q4%mpr, az)
    if (present (q5)) call mpinp (iu, q5%mpr, az)
    if (present (q6)) call mpinp (iu, q6%mpr, az)
    if (present (q7)) call mpinp (iu, q7%mpr, az)
    if (present (q8)) call mpinp (iu, q8%mpr, az)
    if (present (q9)) call mpinp (iu, q9%mpr, az)
    return
  end subroutine

  subroutine mp_inpz (iu, z1, z2, z3, z4, z5, z6, z7, z8, z9)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (out):: z1, z2, z3, z4, z5, z6, z7, z8, z9
    optional:: z2, z3, z4, z5, z6, z7, z8, z9
    call mpinp (iu, z1%mpc, az)
    call mpinp (iu, z1%mpc(mp41), az)
    if (present (z2)) call mpinp (iu, z2%mpc, az)
    if (present (z2)) call mpinp (iu, z2%mpc(mp41), az)
    if (present (z3)) call mpinp (iu, z3%mpc, az)
    if (present (z3)) call mpinp (iu, z3%mpc(mp41), az)
    if (present (z4)) call mpinp (iu, z4%mpc, az)
    if (present (z4)) call mpinp (iu, z4%mpc(mp41), az)
    if (present (z5)) call mpinp (iu, z5%mpc, az)
    if (present (z5)) call mpinp (iu, z5%mpc(mp41), az)
    if (present (z6)) call mpinp (iu, z6%mpc, az)
    if (present (z6)) call mpinp (iu, z6%mpc(mp41), az)
    if (present (z7)) call mpinp (iu, z7%mpc, az)
    if (present (z7)) call mpinp (iu, z7%mpc(mp41), az)
    if (present (z8)) call mpinp (iu, z8%mpc, az)
    if (present (z8)) call mpinp (iu, z8%mpc(mp41), az)
    if (present (z9)) call mpinp (iu, z9%mpc, az)
    if (present (z9)) call mpinp (iu, z9%mpc(mp41), az)
    return
  end subroutine

  function mp_jtoq (ja)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_jtoq
    intent (in):: ja
    call mpeq (ja%mpi, mp_jtoq%mpr)
    return
  end function

  function mp_ztoq (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_ztoq
    intent (in):: za
    call mpeq (za%mpc, mp_ztoq%mpr)
    return
  end function

  function mp_itoq (ia)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_itoq
    intent (in):: ia
    da = ia
    call mpdmc (da, 0, mp_itoq%mpr)
    return
  end function

  function mp_rtoq (ra)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_rtoq
    intent (in):: ra
    da = ra
    call mpdmc (da, 0, mp_rtoq%mpr)
    return
  end function

  function mp_ctoq (ca)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_ctoq
    intent (in):: ca
    da = ca
    call mpdmc (da, 0, mp_ctoq%mpr)
    return
  end function

  function mp_dtoq (da)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_dtoq
    intent (in):: da
    call mpdmc (da, 0, mp_dtoq%mpr)
    return
  end function

  function mp_xtoq (xa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_xtoq
    intent (in):: xa
    da = xa
    call mpdmc (da, 0, mp_xtoq%mpr)
    return
  end function

  function mp_atoq (aa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    character*(*), intent (in):: aa
    type (mp_real):: mp_atoq
    l = len (aa)
    do i = 1, l
      az(i) = aa(i:i)
    enddo
    call mpdexc (az, l, mp_atoq%mpr)
    return
  end function

  subroutine mp_outj (iu, j1, j2, j3, j4, j5, j6, j7, j8, j9)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: j1, j2, j3, j4, j5, j6, j7, j8, j9
    optional:: j2, j3, j4, j5, j6, j7, j8, j9
    call mpout (iu, j1%mpi, mpoud, az)
    if (present (j2)) call mpout (iu, j2%mpi, mpoud, az)
    if (present (j3)) call mpout (iu, j3%mpi, mpoud, az)
    if (present (j4)) call mpout (iu, j4%mpi, mpoud, az)
    if (present (j5)) call mpout (iu, j5%mpi, mpoud, az)
    if (present (j6)) call mpout (iu, j6%mpi, mpoud, az)
    if (present (j7)) call mpout (iu, j7%mpi, mpoud, az)
    if (present (j8)) call mpout (iu, j8%mpi, mpoud, az)
    if (present (j9)) call mpout (iu, j9%mpi, mpoud, az)
     return
  end subroutine

  subroutine mp_outq (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: q1, q2, q3, q4, q5, q6, q7, q8, q9
    optional:: q2, q3, q4, q5, q6, q7, q8, q9
    call mpout (iu, q1%mpr, mpoud, az)
    if (present (q2)) call mpout (iu, q2%mpr, mpoud, az)
    if (present (q3)) call mpout (iu, q3%mpr, mpoud, az)
    if (present (q4)) call mpout (iu, q4%mpr, mpoud, az)
    if (present (q5)) call mpout (iu, q5%mpr, mpoud, az)
    if (present (q6)) call mpout (iu, q6%mpr, mpoud, az)
    if (present (q7)) call mpout (iu, q7%mpr, mpoud, az)
    if (present (q8)) call mpout (iu, q8%mpr, mpoud, az)
    if (present (q9)) call mpout (iu, q9%mpr, mpoud, az)
     return
  end subroutine

  subroutine mp_outz (iu, z1, z2, z3, z4, z5, z6, z7, z8, z9)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: z1, z2, z3, z4, z5, z6, z7, z8, z9
    optional:: z2, z3, z4, z5, z6, z7, z8, z9
    call mpout (iu, z1%mpc, mpoud, az)
    call mpout (iu, z1%mpc(mp41), mpoud, az)
    if (present (z2)) call mpout (iu, z2%mpc, mpoud, az)
    if (present (z2)) call mpout (iu, z2%mpc(mp41), mpoud, az)
    if (present (z3)) call mpout (iu, z3%mpc, mpoud, az)
    if (present (z3)) call mpout (iu, z3%mpc(mp41), mpoud, az)
    if (present (z4)) call mpout (iu, z4%mpc, mpoud, az)
    if (present (z4)) call mpout (iu, z4%mpc(mp41), mpoud, az)
    if (present (z5)) call mpout (iu, z5%mpc, mpoud, az)
    if (present (z5)) call mpout (iu, z5%mpc(mp41), mpoud, az)
    if (present (z6)) call mpout (iu, z6%mpc, mpoud, az)
    if (present (z6)) call mpout (iu, z6%mpc(mp41), mpoud, az)
    if (present (z7)) call mpout (iu, z7%mpc, mpoud, az)
    if (present (z7)) call mpout (iu, z7%mpc(mp41), mpoud, az)
    if (present (z8)) call mpout (iu, z8%mpc, mpoud, az)
    if (present (z8)) call mpout (iu, z8%mpc(mp41), mpoud, az)
    if (present (z9)) call mpout (iu, z9%mpc, mpoud, az)
    if (present (z9)) call mpout (iu, z9%mpc(mp41), mpoud, az)
     return
  end subroutine

  function mp_nint (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_nint
    intent (in):: qa
    call mpnint (qa%mpr, mp_nint%mpi)
    return
  end function

  function mp_jtor (ja)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: ja
    real:: mp_jtor
    call mpmdc (ja%mpi, da, ia)
    mp_jtor = da * 2.d0 ** ia
    return
  end function

  function mp_qtor (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: qa
    real:: mp_qtor
    call mpmdc (qa%mpr, da, ia)
    mp_qtor = da * 2.d0 ** ia
    return
  end function

  function mp_ztor (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    intent (in):: za
    real:: mp_ztor
    call mpmdc (za%mpc, da, ia)
    mp_ztor = da * 2.d0 ** ia
    return
  end function

  function mp_signj (ja, jb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_integer):: mp_signj
    intent (in):: ja, jb
    call mpeq (ja%mpi, mp_signj%mpi)
    mp_signj%mpi(1) = sign (mp_signj%mpi(1), jb%mpi(1))
    return
  end function

  function mp_signq (qa, qb)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_signq
    intent (in):: qa, qb
    call mpeq (qa%mpr, mp_signq%mpr)
    mp_signq%mpr(1) = sign (mp_signq%mpr(1), qb%mpr(1))
    return
  end function

  function mp_sin (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_sin
    intent (in):: qa
    call mpcssn (qa%mpr, mppic%mpr, mpt1, mp_sin%mpr)
    return
  end function

  function mp_sinz (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_sinz
    intent (in):: za
    call mpeq (za%mpc(mp41), mpt2)
    mpt2(1) = - mpt2(1)
    call mpexp (mpt2, mpl02%mpr, mpt1)
    call mpdmc (1.d0, 0, mpt3)
    call mpdiv (mpt3, mpt1, mpt2)
    call mpcssn (za%mpc, mppic%mpr, mpt3, mpt4)
    call mpadd (mpt1, mpt2, mpt5)
    call mpmuld (mpt5, 0.5d0, 0, mpt6)
    call mpmul (mpt6, mpt4, mp_sinz%mpc)
    call mpsub (mpt1, mpt2, mpt5)
    call mpmuld (mpt5, -0.5d0, 0, mpt6)
    call mpmul (mpt6, mpt3, mp_sinz%mpc(mp41))
    return
  end function

  function mp_sinh (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_sinh
    intent (in):: qa
    call mpcssh (qa%mpr, mpl02%mpr, mpt1, mp_sinh%mpr)
    return
  end function

  function mp_sqrtq (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_sqrtq
    intent (in):: qa
    call mpsqrt (qa%mpr, mp_sqrtq%mpr)
    return
  end function

  function mp_sqrtz (za)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_complex):: mp_sqrtz
    intent (in):: za
    call mpcsqt (mp4, za%mpc, mp_sqrtz%mpc)
    return
  end function

  function mp_tan (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_tan
    intent (in):: qa
    call mpcssn (qa%mpr, mppic%mpr, mpt1, mpt2)
    call mpdiv (mpt2, mpt1, mp_tan%mpr)
    return
  end function

  function mp_tanh (qa)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
      type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
    type (mp_real):: mp_tanh
    intent (in):: qa
    call mpcssh (qa%mpr, mpl02%mpr, mpt1, mpt2)
    call mpdiv (mpt2, mpt1, mp_tanh%mpr)
    return
  end function

end module

module mpmodule
use mpfunmod
use mpintmod
use mprealmod
use mpcmpmod
use mpgenmod
end module
