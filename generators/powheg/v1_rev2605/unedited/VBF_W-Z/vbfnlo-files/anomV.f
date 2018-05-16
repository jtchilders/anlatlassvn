********************************************************************************
********************************************************************************
***                                                                          ***
*** anomV.F                                                                  ***
*** 1 June 2012                                                              ***
*** sophy@particle.uni-karlsruhe                                             ***
***                                                                          ***
*** Legacy NOTE: this is a consolidation of the 4 previous anomalous boson   ***
*** couplings read-in routines                                               ***
***                                                                          ***
********************************************************************************
********************************************************************************

      subroutine read_anomVcouplings

c reads in initial set of anomalous coupling constants 
c from file anomV.dat. The coupling constants
c are stored in common blocks in anomV_cplg.inc      

      implicit none

** Dummy parameter used to check consistency of Delta kappa_Z and _photon
      double precision consistency, massscale2

** Dummy parameters for filling HVV couplings
      double precision treefacW, treefacZ, loopfac, dum
      common /lhcoup/ treefacW, treefacZ, loopfac
      integer formfact1, formfact_ind1
      
      double precision anomVinput
      external anomVinput

       include "an_couplings.inc"   
       include "mssm.inc"
!        include "global.inc"
!        include "process.inc"


      print *," "
      print *,"  Information on anomalous coupling parameters  "
      print *,"------------------------------------------------"


** Read switches controlling whether we use a formfactor:
      formfact=.false.
      formfact1=anomVinput("FORMFAC")
      if(formfact1.eq.1d0) formfact=.true.

      formfac_ind=.false.
      formfact_ind1=anomVinput("FORMFAC_IND")
      if(formfact_ind1.eq.1d0) formfac_ind=.true.

** universal formfactors:
      if (formfact) then
         ffmassscale2=anomVinput("FFMASSSCALE")
         if (ffmassscale2 .lt. 1d-3) then
            write(*,*)'The input value of FFMASSSCALE is too small'
            write(*,*)"We will use FFMASSSCALE = 2000 GeV instead"
            ffmassscale2 = 2000d0
         end if            
         ffmassscale2 = ffmassscale2**2
         ffexponent=anomVinput("FFEXP")
         if (ffexponent .lt. 0) then
            write(*,*)"You cannot use a negative exponent for the"
            write(*,*)"formfactor.  We will instead use"
            write(*,*)"FFEXP = 2"
            ffexponent = 2
         end if
         if (.not. formfac_ind) then
            massscale2FB = ffmassscale2
            massscale2FW = ffmassscale2
            massscale2FWWW = ffmassscale2
            ffexpFWWW = ffexponent
            ffexpFW = ffexponent
            ffexpFB = ffexponent

            massscale2L = ffmassscale2
            massscale2KZ = ffmassscale2
            massscale2KA = ffmassscale2
            massscale2G = ffmassscale2
            ffexpL = ffexponent
            ffexpKZ = ffexponent
            ffexpKA = ffexponent
            ffexpG = ffexponent
         end if
      end if


** Now we move on to the dim-6 operators that parametrise the WWZ and WWgamma
** couplings:
      trianom=anomVinput("TRIANOM")

      if (trianom .eq. 1) then

         fwww_0=anomVinput("FWWW")
         fw_0=anomVinput("FW")
         fb_0=anomVinput("FB")

         if (formfac_ind .and. formfact) then
            massscale2FWWW=anomVinput("MASS_SCALE_FWWW")
            massscale2FWWW = massscale2FWWW**2
            if (massscale2FWWW .lt. 1d-3) then
               write(*,*)'The input value MASS_SCALE_FWWW is too small'
               write(*,*)"We will use the universal formfactor scale."
               massscale2FWWW = ffmassscale2
               write(*,*)"MASS_SCALE_FWWW =", sqrt(massscale2FWWW)
            end if            
            ffexpFWWW=anomVinput("FFEXP_FWWW")
            if (ffexpFWWW .lt. 0) then
               write(*,*)"You cannot use a negative exponent for the"
               write(*,*)"formfactor.  We will instead use value from"
               write(*,*)"the universal formfactor"
               ffexpFWWW = ffexponent 
               write(*,*)"FFEXP_FWWW =", ffexpFWWW
            end if

            massscale2FW=anomVinput("MASS_SCALE_FW")
            massscale2FW = massscale2FW**2
            if (massscale2FW .lt. 1d-3) then
               write(*,*)'The input value MASS_SCALE_FW is too small'
               write(*,*)"We will use the universal formfactor scale."
               massscale2FW = ffmassscale2
               write(*,*)"MASS_SCALE_FW =", sqrt(massscale2FW)
            end if            
            ffexpFW=anomVinput("FFEXP_FW") 
            if (ffexpFW .lt. 0) then
               write(*,*)"You cannot use a negative exponent for the"
               write(*,*)"formfactor.  We will instead use value from"
               write(*,*)"the universal formfactor"
               ffexpFW = ffexponent 
               write(*,*)"FFEXP_FW =", ffexpFW
            end if

            massscale2FB=anomVinput("MASS_SCALE_FB")
            massscale2FB = massscale2FB**2
            if (massscale2FB .lt. 1d-3) then
               write(*,*)'The input value MASS_SCALE_FB is too small'
               write(*,*)"We will use the universal formfactor scale."
               massscale2FB = ffmassscale2
               write(*,*)"MASS_SCALE_FB =", sqrt(massscale2FB)
            end if            
            ffexpFB= anomVinput("FFEXP_FB")
            if (ffexpFB .lt. 0) then
               write(*,*)"You cannot use a negative exponent for the"
               write(*,*)"formfactor.  We will instead use value from"
               write(*,*)"the universal formfactor"
               ffexpFB = ffexponent 
               write(*,*)"FFEXP_FB =", ffexpFB
            end if
         end if

** setting trianom=2 parametrisation
         aDkappa0_0 = 0.5d0*MW2*(fb_0 + fw_0)
         zDkappa0_0 = 0.5d0*MZ2*(CW2*fw_0 - SW2*fb_0)
         lambda0_0 = 3d0*MW2*EL*EL*fwww_0/(2d0*SW2)
         zDg0_0 = 0.5d0*MZ2*fw_0
** formfactors: note formfactor for kappa0(z/gamma) is not a simple conversion 
**              and is therefore done in the formfactor routine
         massscale2L = massscale2FWWW
         ffexpL = ffexpFWWW
         massscale2G = massscale2FW
         ffexpG = ffexpFW


      else if (trianom .eq. 2) then

         lambda0_0= anomVinput("LAMBDA0")
         zDg0_0 =anomVinput("ZDELTAG1")
         zDkappa0_0 =anomVinput("ZDELTAKAPPA0")
         aDkappa0_0 =anomVinput("ADELTAKAPPA0")

** If zdeltag, zdeltakappa0 or adeltakappa0 is zero, we set it to be consistent
         if ((zDkappa0_0.eq.0D0) .and. 
     -        ((zDg0_0.ne.0D0) .or. (aDkappa0_0.ne.0D0)))then
            zDkappa0_0 = zDg0_0-SW2/CW2*aDkappa0_0
            write(*,*)'Note!' 
            write(*,*)'ZDELTAKAPPA0 = ZDELTAG0-SW2/CW2*ADELTAKAPPA0' 
            write(*,*)'Thus setting ZDELTAKAPPA0 = ', zDkappa0_0
         else if ((zDg0_0.eq.0D0) .and. 
     -        ((zDkappa0_0.ne.0D0) .or. (aDkappa0_0.ne.0D0))) then
            zDg0_0 = zDkappa0_0 + SW2/CW2*aDkappa0_0
            write(*,*)'Note!' 
            write(*,*)'ZDELTAG0 = ZDELTAKAPPA0+SW2/CW2*ADELTAKAPPA0'
            write(*,*)'Thus setting ZDELTAG0 = ', zDg0_0
         else if ((aDkappa0_0.eq.0D0) .and. 
     -        ((zDkappa0_0.ne.0D0) .or. (zDg0_0.ne.0D0))) then
            aDkappa0_0 = CW2/SW2*(zDg0_0-zDkappa0_0)
            write(*,*)'Note!' 
            write(*,*)'ADELTAKAPPA0 = CW2/SW2*(ZDELTAG0-ZDELTAKAPPA0)'
            write(*,*)'Thus setting ADELTAKAPPA0 = ', aDkappa0_0
         end if

** Consistency check on kappa0, lambda0, g
         consistency = (zDkappa0_0/(SW2*MZ2)) + (aDkappa0_0/MW2)
         consistency = consistency - (zDg0_0/(SW2*MZ2))
         consistency = sqrt(consistency**2)
         if (consistency .gt. 1D-7) then
            aDkappa0_0 = CW2/SW2*(zDg0_0-zDkappa0_0)
            write(*,*)'WARNING! The values for ZDELTAKAPPA0,'
            write(*,*)'ADELTAKAPPA0 AND ZDELTAG0 are not consistent.'
            write(*,*)'We will use ZDELTAKAPPA0 AND ZDELTAG0 as input.'
            write(*,*)'This gives ADELTAKAPPA0 =', aDkappa0_0
            write(*,*)'If you want to use ADELTAKAPPA0 as input'
            write(*,*)'set either ZDELTAKAPPA0 or ZDELTAG0 to 0.'
         end if

** Setting individual formfactors
         if (formfact .and. formfac_ind) then
            massscale2L =anomVinput("MASS_SCALE_LAMBDA")
            if(massscale2L.lt.0d0) massscale2L=0d0
            ffexpL =anomVinput("FFEXP_LAMBDA")
            massscale2L = massscale2L**2
            massscale2G =anomVinput("MASS_SCALE_G")
            if(massscale2G.lt.0d0) massscale2G=0d0
            ffexpG=anomVinput("FFEXP_G")
            massscale2G = massscale2G**2
            massscale2KZ =anomVinput("MASS_SCALE_ZKAPPA")
            if(massscale2KZ.lt.0d0) massscale2KZ=0d0
            ffexpKZ =anomVinput("FFEXP_ZKAPPA")
            massscale2KZ = massscale2KZ**2
            massscale2KA =anomVinput("MASS_SCALE_AKAPPA")
            if(massscale2KA.lt.0d0) massscale2KA=0d0
            ffexpKA =anomVinput("FFEXP_AKAPPA")
            massscale2KA = massscale2KA**2

* Check for formfactor mass scales set to zero
            if (massscale2L .eq. 0d0) then
               write(*,*)"You have input 0 GeV mass scales for lambda "
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", ffmassscale2
               massscale2L = ffmassscale2
            end if
               
            if ((massscale2KA .eq. 0d0) .and. (massscale2KZ .eq. 0d0)
     &           .and. (massscale2G .eq. 0d0)) then
               write(*,*)"You have input 0 GeV mass scales for all "
               write(*,*)"individual formfactors."
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", sqrt(ffmassscale2)
               write(*,*)"for MASS_SCALE_ZKAPPA and MASS_SCALE_G"
               massscale2KZ = ffmassscale2
               massscale2G = ffmassscale2
            else if ((massscale2KA .eq. 0d0) .and. 
     &              (massscale2KZ .eq. 0d0)) then
               write(*,*)"You have input 0 GeV mass scales for"
               write(*,*)"both MASS_SCALE_ZKAPPA and MASS_SCALE_AKAPPA"
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", sqrt(ffmassscale2)
               write(*,*)"for MASS_SCALE_ZKAPPA"
               massscale2KZ = ffmassscale2
            else if ((massscale2KA .eq. 0d0) .and.
     &           (massscale2G .eq. 0d0)) then
               write(*,*)"You have input 0 GeV mass scales for"
               write(*,*)"both MASS_SCALE_AKAPPA and MASS_SCALE_G"
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", sqrt(ffmassscale2)
               write(*,*)"for MASS_SCALE_G"
               massscale2G = ffmassscale2
            else if ((massscale2KZ .eq. 0d0) .and.
     &           (massscale2G .eq. 0d0)) then
               write(*,*)"You have input 0 GeV mass scales for"
               write(*,*)"both MASS_SCALE_ZKAPPA and MASS_SCALE_G"
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", sqrt(ffmassscale2)
               write(*,*)"for MASS_SCALE_G"
               massscale2G = ffmassscale2
            end if

** Checking for negative exponents in formfactor
            if (ffexpKA .lt. 0) then
               write(*,*)"You cannot use negative exponents in the"
               write(*,*)"formfactors.  We will instead use"
               write(*,*)"FFEXPKA = 2"
               ffexpKA = 2
            end if
            if (ffexpKZ .lt. 0) then
               write(*,*)"You cannot use negative exponents in the"
               write(*,*)"formfactors.  We will instead use"
               write(*,*)"FFEXPKZ = 2"
               ffexpKZ = 2
            end if
            if (ffexpG .lt. 0) then
               write(*,*)"You cannot use negative exponents in the"
               write(*,*)"formfactors.  We will instead use"
               write(*,*)"FFEXPKZ = 2"
               ffexpG = 2
            end if
            if (ffexpL .lt. 0) then
               write(*,*)"You cannot use negative exponents in the"
               write(*,*)"formfactors.  We will instead use"
               write(*,*)"FFEXPKZ = 2"
               ffexpL = 2
            end if

** Consistency checks and settings:
            if ((massscale2KA .eq. 0d0) .and. 
     &           (aDkappa0_0 .ne. 0d0)) then 
               massscale2 = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &              dble(ffexpG)) - zDkappa0_0/((1d0 + 1000d0/
     &              massscale2KZ)**dble(ffexpKZ)))*CW2/(SW2*aDkappa0_0)
               massscale2 = 1000d0/(massscale2**
     &              dble(-1d0/ffexpKA) - 1d0)
               if (massscale2 .gt. 1d0) then
                  massscale2KA = massscale2
                  write(*,*)'MASS_SCALE_AKAPPA is set to', 
     &                 sqrt(massscale2KA)
                  write(*,*)'for consistency'
               else
                  write(*,*)"Sorry! We cannot make the input scales"
                  write(*,*)"MASS_SCALE_G and MASS_SCALE_ZKAPPA"
                  write(*,*)"consistent.  Please change a value or"
                  write(*,*)"use a universal formfactor."
                  stop
               end if
            else if ((massscale2KZ .eq. 0d0) .and.
     &              (zDkappa0_0 .ne. 0d0)) then
               massscale2 = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &              dble(ffexpG)) - (SW2*aDkappa0_0/CW2)/((1d0 + 1000d0/
     &              massscale2KA)**dble(ffexpKA)))/zDkappa0_0
               massscale2 = 1000d0/(massscale2**dble(-1d0/ffexpKZ) - 
     &              1d0)
               if (massscale2 .gt. 1d0) then
                  massscale2KZ = massscale2
                  write(*,*)'MASS_SCALE_ZKAPPA is set to', 
     &                 sqrt(abs(massscale2KZ))
                  write(*,*)'for consistency'
               else
                  write(*,*)"Sorry! We cannot make the input scales"
                  write(*,*)"MASS_SCALE_G and MASS_SCALE_AKAPPA"
                  write(*,*)"consistent.  Please change a value or"
                  write(*,*)"use a universal formfactor."
                  stop
               end if
            else if ((massscale2G .eq. 0d0) .and. 
     &              (zDg0_0 .ne. 0d0)) then
               massscale2 = ((SW2*aDkappa0_0/CW2)/((1d0 + 1000d0/
     &              massscale2KA)**dble(ffexpKA)) + zDkappa0_0/((1d0 + 
     &              1000d0/massscale2KZ)**dble(ffexpKZ)))/zDg0_0
               massscale2 = 1000d0/(massscale2**dble(-1d0/ffexpG) - 
     &              1d0)  
               if (massscale2G .gt. 1d0) then
                  massscale2G = massscale2
                  write(*,*)'MASS_SCALE_G is set to', 
     &                 sqrt(abs(massscale2G))
                  write(*,*)'for consistency'
               else
                  write(*,*)"Sorry! We cannot make the input scales"
                  write(*,*)"MASS_SCALE_ZKAPPA and MASS_SCALE_AKAPPA"
                  write(*,*)"consistent.  Please change a value or"
                  write(*,*)"use a universal formfactor."
                  stop
               end if
            end if

            consistency = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &              dble(ffexpG)) - zDkappa0_0/((1d0 + 1000d0/
     &              massscale2KZ)**dble(ffexpKZ)))*CW2/SW2
            consistency = consistency - aDkappa0_0/((1d0 + 1000d0/
     &              massscale2KA)**dble(ffexpKA))
            if (abs(consistency) .gt. 1d-7) then
               massscale2 = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &              dble(ffexpG)) - zDkappa0_0/((1d0 + 1000d0/
     &              massscale2KZ)**dble(ffexpKZ)))*CW2/(SW2*aDkappa0_0)
               massscale2 = 1000d0/(massscale2**
     &              dble(-1d0/ffexpKA) - 1d0)
               if (massscale2 .gt. 1d0) then
                  massscale2KA = massscale2
                  write(*,*)"The mass scales MASS_SCALE_G,MASS_SCALE_KZ"
                  write(*,*)"and MASS_SCALE_KA are inconsistent."
                  write(*,*)"We will use MASS_SCALE_G AND MASS_SCALE_KZ"
                  write(*,*)"as inputs, which gives"
                  write(*,*)'MASS_SCALE_AKAPPA =', sqrt(massscale2KA)
                  write(*,*)"If you want to use MASS_SCALE_AKAPPA as an"
                  write(*,*)"input, set either MASS_SCALE_G or "
                  write(*,*)"MASS_SCALE_ZKAPPA to zero."
               else
                  massscale2 = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &                 dble(ffexpG)) - (SW2*aDkappa0_0/CW2)/((1d0 + 
     &                 1000d0/massscale2KA)**dble(ffexpKA)))/zDkappa0_0
                  massscale2 = 1000d0/(massscale2**
     &                 dble(-1d0/ffexpKZ) - 1d0)
                  if (massscale2 .gt. 1d0) then
                     massscale2KZ = massscale2
                     write(*,*)"The mass scales MASS_SCALE_G,"
                     write(*,*)"MASS_SCALE_KZ and MASS_SCALE_KA are" 
                     write(*,*)" inconsistent.  We will use "
                     write(*,*) "MASS_SCALE_G AND MASS_SCALE_KA as"
                     write(*,*)"inputs, which gives"
                     write(*,*)'MASS_SCALE_ZKAPPA =', sqrt(massscale2KZ)
                     write(*,*)"If you want to use MASS_SCALE_AKAPPA "
                     write(*,*)"as an input, set either MASS_SCALE_G"
                     write(*,*)"or MASS_SCALE_ZKAPPA to zero."
                  else
                     massscale2 = ((SW2*aDkappa0_0/CW2)/((1d0 + 1000d0/
     &                    massscale2KA)**dble(ffexpKA)) + zDkappa0_0/ 
     &                    ((1d0 + 1000d0/massscale2KZ)**dble(ffexpKZ)))/
     &                    zDg0_0
                     massscale2 = 1000d0/(massscale2**
     &                    dble(-1d0/ffexpG) - 1d0)  
                     if (massscale2 .gt. 1d0) then
                        massscale2G = massscale2
                        write(*,*)"The mass scales MASS_SCALE_G,"
                        write(*,*)"MASS_SCALE_KZ and MASS_SCALE_KA are" 
                        write(*,*)"inconsistent.  We will use "
                        write(*,*) "MASS_SCALE_KZ AND MASS_SCALE_KA as"
                        write(*,*)"inputs, which gives"
                        write(*,*)'MASS_SCALE_G =', sqrt(massscale2G)
                        write(*,*)"If you want to use MASS_SCALE_G "
                        write(*,*)"as an input, set either "
                        write(*,*)"MASS_SCALE_AKAPPA or "
                        write(*,*)"MASS_SCALE_ZKAPPA to zero."
                     else
                        write(*,*)"Sorry! We cannot make the input"
                        write(*,*)"formfactor mass scales MASS_SCALE_G,"
                        write(*,*)"MASS_SCALE_AKAPPA and"
                        write(*,*)"MASS_SCALE_ZKAPPA consistent. Please"
                        write(*,*)"change one or more values or use a"
                        write(*,*)"universal formfactor."
                        stop
                     end if
                  end if
               end if
                  
            end if               

         end if


** Setting trianom=1 parameterisation
         fwww_0 = 2D0*SW2*lambda0_0/(3D0*MW2*EL*EL)
         fw_0 = 2D0*zDg0_0/MZ2
         fb_0 = (CW2*zDg0_0 - zDkappa0_0)*2D0/(SW2*MZ2)
* Formfactors: note formfactor for FB is not a simple conversion and is thus
*              set in the formfactor subroutine
         massscale2FWWW = massscale2L
         ffexpFWWW = ffexpL
         massscale2FW = massscale2G
         ffexpFW = ffexpG

      else 
         write(*,*)'Invalid entry for TRIANOM:'
         write(*,*)'gauge boson anomalous coupling parametrisation'
         write(*,*)'Please set TRIANOM = 1 or 2'
         stop
      end if




** converting to HVV coupling notation, for Higgs width calculation
      treefacZ = 1d0
      treefacW = 1d0
      loopfac = 1d0
      dum = 0d0
      call anomH_convert(4,4,(trianom .eq. 2),(trianom .eq. 1),
     &     dum, dum, dum, dum, zDg0_0, aDkappa0_0, dum,
     &     dum, dum, dum, dum, fw_0, fb_0, dum)
      


** Setting formfactors
      fwww=fwww_0
      fw=fw_0
      fb=fb_0

      lambda0 = lambda0_0
      zDg0 = zDg0_0
      zDkappa0 = zDkappa0_0
      aDkappa0 = aDkappa0_0



      end



********************************************************************************
********************************************************************************

