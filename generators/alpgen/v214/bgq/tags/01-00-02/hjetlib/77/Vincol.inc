C *** MODULE incol ***
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
c  ! TBMOM(I,J) gives the coulor of the combination of particles i and j
      integer tbmom(1058)
c        ! TBNEW(I,J) gives the number of different coulors from i and j
      integer tbnew(529)
!                                         combination
c                                  ! TBCOEFF(I,J) SU3 clebsh coefficient
      double precision tbcoeff(1058)
      common /incol/tbcoeff, tbmom, tbnew
C *** END MODULE incol ***
