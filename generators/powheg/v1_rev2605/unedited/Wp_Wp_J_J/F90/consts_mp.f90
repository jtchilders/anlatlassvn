!! For mp constants to be available one must first call fixmpconsts
module consts_mp
  use mpmodule
  implicit none
  private

  type(mp_real), public, save :: pi,twopi,zero,half,one,mone,two,three,four
  type(mp_real), public, save :: sqrt2,msqrt2,sqrt3,csqrt2,ln10

  type(mp_complex), public, save :: czero, cone, ctwo, chalf, ci 

  type(mp_real), public, save :: tol    ! tolerance related to precision 
  type(mp_real), public, save :: sq2tol ! sqrt(tolerance) related to precision 

  type(mp_real), public, save :: propcut  
  
  ! -- masses of various particles
  type(mp_real), public, save :: mt, mb, mw, mz  

  public :: fixmpconsts

contains

  subroutine fixmpconsts(init_prec)       
    integer, intent(in) :: init_prec 

    ! -- real consts
    half   = mpreal('0.5')
    two    = mpreal('2.0') 
    zero   = mpreal('0.0')
    one    = mpreal('1.0') 
    mone   = -one
    three  = mpreal('3.0')
    four   = mpreal('4.0')
    msqrt2 = -sqrt2
    sqrt3  = sqrt(three)
    pi     = mppic
    twopi  = two*pi
    sqrt2  = sqrt(two)

    ! -- complex consts
    czero  = mpcmpl('0.0','0.0')
    cone   = mpcmpl('1.0','0.0')
    ctwo   = mpcmpl('2.0','0.0')
    chalf  = mpcmpl('0.5','0.0')
    ci     = mpcmpl('0.0','1.0')
    csqrt2 = cone*sqrt2 
    ln10 = mpreal('2.3025850929940456840179914546843642')
    ! -- tolerances 
    tol = exp(mpreal('1.0')*(-(init_prec-2))*ln10)
    write(*,*) 'MPTOL IS',init_prec
    call mpwrite(0,tol)
    sq2tol = sqrt(tol)
    propcut = sq2tol 

    ! - masses 
    mt = mpreal('0.0')
    mb = mpreal('0.0') 
    mw = mpreal('0.8')
    

  end subroutine fixmpconsts

end module consts_mp
