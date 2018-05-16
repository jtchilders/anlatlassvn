C *** MODULE incnst ***
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      integer nmax                          !maximum number of particles
      parameter (nmax = 10)
!  integer*4, parameter :: nxop=10    !maximum number of particles, optimized version
      integer nlor         !maximum number of lorentz degrees of freedom
      parameter (nlor = 6)
      integer npmax     !maximum number of different particles "flavour"
      parameter (npmax = 43)
      integer nmrgmx
      parameter (nmrgmx = 5001)
C *** END MODULE incnst ***
