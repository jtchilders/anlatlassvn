      subroutine hdecayparser(in_higgs_mass,higgs_width)
      implicit none
      include 'PhysPars.h'
      character* 50 hdecayfname
      character* 1 dummy
      real* 8 in_higgs_mass,higgs_width
      real* 8 higgsmass,gammagammawidth,zgammawidth
      real* 8 hphmwidth,hhwidth
      real* 8 wwwidth,zzwidth,totalwidth
      real* 8 tmp_mass,tmp_width
      integer i
      logical ifirst

c     h
      if (ih.eq.1) then
         hdecayfname = 'br.l3_2HDM'
c     H
      else if (ih.eq.2) then
         hdecayfname = 'br.h3_2HDM'
c     A
      else if (ih.eq.3) then
         hdecayfname = 'br.a3_2HDM'
      else
         write(*,*) 'Higgs type not recognized'
         stop
      endif

      write(*,*) 'Opening hdecay file'
      open(21, FILE=hdecayfname, STATUS="old", ACTION="read")

c     Skip the first lines
      read(21,*) dummy
      read(21,*) dummy

      ifirst = .true.
c     cycle on all the lines 
      do
         if (ih.eq.1) then
            read(21,*,END=400) higgsmass,gammagammawidth,zgammawidth,
     $           wwwidth,hphmwidth,totalwidth
         else if (ih.eq.2) then
            read(21,*,END=400) higgsmass,hhwidth,gammagammawidth,
     $            zgammawidth,wwwidth,hphmwidth,totalwidth
         else if (ih.eq.3) then
            read(21,*,END=400) higgsmass,wwwidth,totalwidth
         endif
c      write(*,*) higgsmass,gammagammawidth,zgammawidth
c     $            ,wwwidth,zzwidth,totalwidth
c     Find the best possible approximation for the total width

      if (higgsmass.ge.in_higgs_mass) then
         if (ifirst.eqv..false.) then
            if (abs(higgsmass-in_higgs_mass).lt.
     $           abs(tmp_mass-in_higgs_mass)) then
               higgs_width = totalwidth
            else
               higgs_width = tmp_width
            end if
         else
            higgs_width = totalwidth
         endif

         goto 401
      end if

      tmp_mass  = higgsmass
      tmp_width = totalwidth
      end do
c     We exited without finding an upper bound. We take the last value
 400  higgs_width = totalwidth
 401  write(*,*) 'Finished reading hdecay file'
      write(*,*) 'Best total width for Higgs mass equal to ',
     $               in_higgs_mass, 'GeV is ', higgs_width, 'GeV',
     $ ' ( computed at ', higgsmass, ')'
      close(21)
c      stop
      end subroutine
