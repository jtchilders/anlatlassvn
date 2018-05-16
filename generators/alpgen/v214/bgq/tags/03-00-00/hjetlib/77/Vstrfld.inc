C *** MODULE strfld ***
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
      integer momd(2047)                 !number of momenta contributing
!
! LORTYP(j) lables the lorentz representation of j; NLORDOF(j) returns the number of lorentz d.o.f. of j
! COLREP(j) lables the coulor representation of j
      integer lortyp(43), nlordof(43), colrep(43)
! FLDPTR is used to select the field amplitudes in the array OFFSH. FLDPTR(j1,j2,j3) returns :
! if j1=1 the number of off-shell fields of type J2 constructed out of J3 on shell external particles,
! if j1=2 the position in the array OFFSHCOL where relevant information for the fields are stored
! are stored (momenta, color ...)
      integer fldptr(2,43,5)
! NEXC,NEXCSM arrays for the iteration
      integer nexcmax
      parameter (nexcmax = 13)
      integer nexc(13,2:5), nexcsm(13,2:5)
      common /strfld/momd, lortyp, nlordof, colrep, fldptr, nexc, nexcsm
C *** END MODULE strfld ***
