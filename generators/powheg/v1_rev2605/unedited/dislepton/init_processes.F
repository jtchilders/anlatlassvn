c     setup Born and real subprocesses
c     set process id using slepton types
      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'LesHouches.h'
      integer i1,i2,i3,i4,i5,k,ii(5)
      equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),(i4,ii(4)),
     &              (i5,ii(5))
      logical debug
      parameter (debug=.true.)
      integer j
      real *8 powheginput
      external powheginput
      logical condition
      integer max_flav,encode_sfermion_pair

      ! set number of active flavors for incoming partons
      max_flav = int(powheginput('incomingflavors'))
      if (max_flav.eq.5)  write(*,*) "WARNING: ignoring sbottom mixing"
      if (max_flav.gt.5)  stop "can't handle more than 5 flavours"

      ! choose process ------------------------------------------------

      i3 = int(powheginput('sleptona'))
      i4 = int(powheginput('sleptonb'))

      if (i3.lt.0) then
        write(*,*) 'WARNING in init_processes:'
     &    //' swapping order of sleptons A and B'
        call swapi(i3,i4)
      endif

      if (mod(i3,100).ne.mod(-i4,100)) stop
     &    'sleptona and sleptonb must be '
     &    //'opposite charge and same generation'

      write(*,*) '****************************************************'
      write(*,*) ' POWHEG: slepton production'
      write(*,*) ' final state: ', i3, ', ', i4
      write(*,*) '****************************************************'
      write(*,*) 
      write(*,*) ' number of incoming flavors: ', max_flav

      ! set process id: encode types of slepton A and B
      lprup(1) = encode_sfermion_pair(i3,i4)
      write(*,*) ' encoded final state to process ID: ', lprup(1)

      ! index of the first coloured particle in the final state
      ! (all subsequent particles are coloured)
      flst_lightpart=5


      ! Born subprocesses ---------------------------------------------

      flst_nborn=0
      do i1=-max_flav,max_flav
         do i2=-max_flav,max_flav
            if(i1.ne.0.and.i1+i2.eq.0) then  ! q qbar
               flst_nborn=flst_nborn+1
               if(flst_nborn.gt.maxprocborn) goto 999
               do k=1,nlegborn
                  flst_born(k,flst_nborn)=ii(k)
               enddo
            endif
         enddo
      enddo
      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif

      ! real subprocesses ---------------------------------------------

      flst_nreal=0
      do i1=-max_flav,max_flav
         do i2=-max_flav,max_flav
            do i5=-max_flav,max_flav
               condition=.false.
               if(.not.(i1.eq.0.and.i2.eq.0)) then ! exclude gg
                  if((i1.ne.0).and.(i1+i2.eq.0).and.(i5.eq.0)) then ! q qbar -> g
                     condition=.true.
                  elseif((i1.eq.0).and.(i2.eq.i5)) then ! g q
                     condition=.true.
                  elseif((i2.eq.0).and.(i1.eq.i5)) then ! q g
                     condition=.true.
                  endif
               endif
               if(condition) then
                  flst_nreal=flst_nreal+1
                  if(flst_nreal.gt.maxprocreal) goto 998
                  do k=1,nlegreal
                     flst_real(k,flst_nreal)=ii(k)
                  enddo
               endif
            enddo
         enddo
      enddo
      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
         enddo
      endif

      return
 998  write(*,*) 'init_processes: increase maxprocreal'
      stop
 999  write(*,*) 'init_processes: increase maxprocborn'
      stop
      end
