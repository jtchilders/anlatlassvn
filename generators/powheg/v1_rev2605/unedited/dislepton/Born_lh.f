c     Sets up the colour for the given flavour configuration
c     already filled in the Les Houches interface.
c     In case there are several colour structures, one
c     should pick one with a probability proportional to
c     the value of the corresponding cross section, for the
c     kinematics defined in the Les Houches interface
      subroutine borncolour_lh
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      ! neutral particles
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0
      ! colored particles
      if((idup(1).gt.0).and.(idup(2).lt.0)) then
         icolup(1,1)=501
         icolup(2,1)=0
         icolup(1,2)=0
         icolup(2,2)=501
      elseif((idup(1).lt.0).and.(idup(2).gt.0)) then
         icolup(1,1)=0
         icolup(2,1)=501
         icolup(1,2)=501
         icolup(2,2)=0
      else
         write(*,*) ' invalid flavour'
         stop
      endif
      end



c     Sets up the resonances whose mass must be preserved
c     on the Les Houches interface; not needed for slepton production.
      subroutine finalize_lh
      implicit none
      end
