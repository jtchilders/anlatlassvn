      subroutine reweightifneeded(dsig0,dsig)
      implicit none
      real * 8 dsig0,dsig
      real * 8 powheginput
      logical ini
      include 'pwhg_weights.h'
      integer j
      real * 8 facscfact,renscfact
      data ini/.true./
      save ini,facscfact,renscfact
      if(ini) then
         facscfact = powheginput("#facscfact")
         renscfact = powheginput("#renscfact")
         if(facscfact.lt.0) facscfact=1
         if(renscfact.lt.0) renscfact=1
         write(*,*) 'Found ',weights_num,' reweighted cross sections'
         ini=.false.
      endif
      do j=1,weights_num
         if(abs(facscfact-weights_facfac(j)).lt.1d-5
     1        .and.
     2      abs(renscfact-weights_renfac(j)).lt.1d-5) then
            dsig=weights_val(j)
            return
         endif
      enddo
      dsig=dsig0
      end
