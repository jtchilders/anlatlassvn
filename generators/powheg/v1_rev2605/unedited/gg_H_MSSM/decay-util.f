      subroutine hdecayparser(in_higgs_mass,higgs_width)
      implicit none
      character* 50 hdecayfname
      character* 1 dummy
      real* 8 in_higgs_mass,higgs_width
      real* 8 higgsmass,ggwidth,gammagammawidth,zgammawidth
      real* 8 wwwidth,zzwidth,totalwidth
      real* 8 tmp_mass,tmp_width
      
      hdecayfname = 'br.sm2'

      write(*,*) 'Opening hdecay file'
      open(21, FILE=hdecayfname, STATUS="old", ACTION="read")

c     Skip the first lines
      read(21,*) dummy
      read(21,*) dummy
      
 300  read(21,*,END=400) higgsmass,ggwidth,gammagammawidth,zgammawidth
     $            ,wwwidth,zzwidth,totalwidth
c      write(*,*) higgsmass,ggwidth,gammagammawidth,zgammawidth
c     $            ,wwwidth,zzwidth,totalwidth
c     Find the best possible approximation for the total width
      if (higgsmass.ge.in_higgs_mass) then
         if (abs(higgsmass-in_higgs_mass).lt.
     $       abs(tmp_mass-in_higgs_mass)) then
            higgs_width = totalwidth
         else
            higgs_width = tmp_width
         end if

         write(*,*) 'Best total width for Higgs mass equal to ',
     $               in_higgs_mass, 'GeV is ', higgs_width, 'GeV'
         goto 400
      end if
      tmp_mass  = higgsmass
      tmp_width = totalwidth
      goto 300
 400  write(*,*) 'Finished reading hdecay file'
      end subroutine
