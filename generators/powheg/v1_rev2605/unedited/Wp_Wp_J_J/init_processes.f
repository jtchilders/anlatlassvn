      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      include 'pwhg_flg-add.h'
      include 'pwhg_par.h'
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,k,ii(nlegreal)
      equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),
     #  (i4,ii(4)),(i5,ii(5)),(i6,ii(6)),(i7,ii(7)),
     #  (i8,ii(8)),(i9,ii(9))
      logical debug
!      parameter (debug=.false.)
      parameter (debug=.true.)
      integer j
      integer charge3(-6:6)
      data charge3 /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      integer fam(-6:6)
      data fam /-3,-3,-2,-2,-1,-1,0,1,1,2,2,3,3/
      logical condition
      real * 8 powheginput
      external powheginput
c     vector boson id and decay
      include 'cvecbos.h'
c     lepton masses
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      flg_fastbtlbound=.true.

      par_isrtinycsi = 1d-6
      par_isrtinyy = 1d-6
      par_fsrtinycsi = 1d-5
      par_fsrtinyy = 1d-5


c******************************************************
c     Choose the process to be implemented
c******************************************************
c    ID of vector boson produced
c   decay products of the vector boson
      vdecaymodeW1=powheginput('vdecaymodeW1')
      vdecaymodeW2=powheginput('vdecaymodeW2')
      if(vdecaymodeW1*vdecaymodeW2.lt.0) then
         write(*,*) ' incompatible decay modes for W1 and W2;'
         write(*,*) ' Must have the same sign! Aborting'
         call exit(-1)
      endif
      if(vdecaymodeW1.gt.0) then
         idvecbos=-24
      else
         idvecbos=24
      endif

      if (lepmass(1).ne.0.51099891d-3) then
         write(*,*) 'block data lepmass not loaded. stop running' 
         stop
      endif
      
      if(idvecbos.eq.24) then
         if ((vdecaymodeW1.ne.-11).and.(vdecaymodeW1.ne.-13)
     $        .and.(vdecaymodeW1.ne.-15)) then
            write(*,*) 'ERROR: The decay mode for W1 you selected' /
     $           /' is not allowed '
            stop 
         endif 
         if ((vdecaymodeW2.ne.-11).and.(vdecaymodeW2.ne.-13)
     $        .and.(vdecaymodeW2.ne.-15)) then
            write(*,*) 'ERROR: The decay mode for W2 you selected' /
     $           /' is not allowed '
            
            stop
         endif
         write(*,*) 
         write(*,*) ' POWHEG: W+ W+ + 2j production and decay ' 
         if (vdecaymodeW1.eq.-11) write(*,*) '         to e+ ve '
         if (vdecaymodeW1.eq.-13) write(*,*) '         to mu+ vmu'
         if (vdecaymodeW1.eq.-15) write(*,*) '         to tau+ vtau'
         write(*,*)                          '            and'
         if (vdecaymodeW2.eq.-11) write(*,*) '         to e+ ve '
         if (vdecaymodeW2.eq.-13) write(*,*) '         to mu+ vmu'
         if (vdecaymodeW2.eq.-15) write(*,*) '         to tau+ vtau'
         write(*,*) 
      elseif(idvecbos.eq.-24) then
         if ((vdecaymodeW1.ne.11).and.(vdecaymodeW1.ne.13)
     $        .and.(vdecaymodeW1.ne.15)) then
            write(*,*) 'ERROR: The decay mode for W1 you selected' /
     $           /' is not allowed '
            stop
         endif
         if ((vdecaymodeW2.ne.11).and.(vdecaymodeW2.ne.13)
     $        .and.(vdecaymodeW2.ne.15)) then
            write(*,*) 'ERROR: The decay mode for W2 you selected' /
     $           /' is not allowed '
            stop
         endif

         write(*,*) 
         write(*,*) ' POWHEG: W- W- + 2j production and decay '
         if (vdecaymodeW1.eq.11) write(*,*) '         to e- ve~ '
         if (vdecaymodeW1.eq.13) write(*,*) '         to mu- vmu~'
         if (vdecaymodeW1.eq.15) write(*,*) '         to tau- vtau~'
         write(*,*)                          '            and'
         if (vdecaymodeW2.eq.11) write(*,*) '         to e- ve~ '
         if (vdecaymodeW2.eq.13) write(*,*) '         to mu- vmu~'
         if (vdecaymodeW2.eq.15) write(*,*) '         to tau- vtau~'
         write(*,*)    
      else
         write(*,*) 'ERROR: The ID of vector boson you selected' 
     $        //' is not allowed (24: W+ -24: W-)'
         stop
      endif

c     change the LHUPI id of the process according to vector boson id
c     and decay
      lprup(1)=10000+vdecaymodeW1 ! 10000+idup of charged decay product of the W
      lprup(2)=10000+vdecaymodeW2 ! 10000+idup of charged decay product of the W
      
      if(lprup(1).eq.10011) then
         decmass=lepmass(1)
         
      elseif(lprup(1).eq.(10000-11)) then
         decmass=lepmass(1)
        
      elseif(lprup(1).eq.10013) then
         decmass=lepmass(2)
         
      elseif(lprup(1).eq.(10000-13)) then
         decmass=lepmass(2)

      elseif(lprup(1).eq.10015) then
         decmass=lepmass(3)
         
      elseif(lprup(1).eq.(10000-15)) then
         decmass=lepmass(3) 
  
      else
c     not yet implemented
         write(*,*) 'non leptonic W decays '//
     #        'not yet implemented'
         stop
      endif   

      if(lprup(2).eq.10011) then
         decmass=lepmass(1)
         
      elseif(lprup(2).eq.(10000-11)) then
         decmass=lepmass(1)
        
      elseif(lprup(2).eq.10013) then
         decmass=lepmass(2)
         
      elseif(lprup(2).eq.(10000-13)) then
         decmass=lepmass(2)

      elseif(lprup(2).eq.10015) then
         decmass=lepmass(3)
         
      elseif(lprup(2).eq.(10000-15)) then
         decmass=lepmass(3) 
  
      else
c     not yet implemented
         write(*,*) 'non leptonic W decays '//
     #        'not yet implemented'
         stop
      endif   
c*********************************************************     
c
c     index of the first coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=7
      i4=vdecaymodeW1
      if ((idvecbos.eq.24).and.(vdecaymodeW1.lt.0)) then
         i3=-vdecaymodeW1+1
      elseif ((idvecbos.eq.-24).and.(vdecaymodeW1.gt.0)) then
         i3=-(vdecaymodeW1+1)
      endif

      i6=vdecaymodeW2
      if ((idvecbos.eq.24).and.(vdecaymodeW2.lt.0)) then
         i5=-vdecaymodeW2+1
      elseif ((idvecbos.eq.-24).and.(vdecaymodeW2.gt.0)) then
         i5=-(vdecaymodeW2+1)
      endif

!      i5 = i3 
!      i6 = i4 


      if(vdecaymodeW1.eq.vdecaymodeW2) then
         vsymfact=0.5d0
      else
         vsymfact=1d0
      endif

c     Born graphs
      flst_nborn=0
      condition=.false.

      do i1=-4,4 ! no b in initial state 
      do i2=-4,4
      do i7=-4,4
      do i8=i7,4

C     charge conservation 
C     construct the Born flavour list for Wp Wp (then reverse signs for Wm Wm) 
                  condition=(charge3(i1)+charge3(i2))
C     .                 .eq.(2*sign(3,idvecbos)+charge3(i7)+charge3(i8))
     .             .eq.(2*sign(3,abs(idvecbos))+charge3(i7)+charge3(i8))
C     assume diagonal CKM 
                  condition = condition .and. (
     .                 (-fam(i1)== fam(i2)  .and. fam(i7)==-fam(i8)).or. 
     .                 (-fam(i1)==-fam(i7) .and. -fam(i2)==-fam(i8)).or. 
     .                 (-fam(i1)==-fam(i8) .and. -fam(i2)==-fam(i7)))
                  if(condition) then
c     q qbar'
                     flst_nborn=flst_nborn+1
                     if(flst_nborn.gt.maxprocborn) goto 999
                     do k=1,nlegborn
                        flst_born(k,flst_nborn)=ii(k)
                     enddo
                  endif
               enddo
            enddo
         enddo
      enddo
C     reverse signs for Wm Wm
      if (idvecbos == -24 ) then 
         flst_born(1:2,:) = -flst_born(1:2,:)
         flst_born(7:8,:) = -flst_born(7:8,:)
      endif

      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif

c     Real graphs
      flst_nreal=0
      condition=.false.
      do i1=-4,4 ! no b in initial state 
      do i2=-4,4
      do i7=-4,4
      do i8=i7,4
      do i9=i8,4


C     charge conservation 
C     construct the real flavour list for Wp Wp (then reverse signs for Wm Wm) 
         condition=(charge3(i1)+charge3(i2))
     .        .eq.(2*sign(3,abs(idvecbos))+charge3(i7)+charge3(i8)
     .        +charge3(i9))
C     assume diagonal CKM 
         condition = condition .and. (
     .    (-fam(i1)== fam(i2).and. fam(i7)==-fam(i8).and.fam(i9)==0).or. 
     .    (-fam(i1)==-fam(i7).and.-fam(i2)==-fam(i8).and.fam(i9)==0).or. 
     .    (-fam(i1)==-fam(i8).and.-fam(i2)==-fam(i7).and.fam(i9)==0).or. 

     .    (fam(i9)== fam(i2).and. fam(i7)==-fam(i8).and.-fam(i1)==0).or. 
     .    (fam(i9)==-fam(i7).and.-fam(i2)==-fam(i8).and.-fam(i1)==0).or. 
     .    (fam(i9)==-fam(i8).and.-fam(i2)==-fam(i7).and.-fam(i1)==0).or. 

     .    (-fam(i1)==-fam(i9).and.fam(i7)==-fam(i8).and.-fam(i2)==0).or. 
     .    (-fam(i1)==-fam(i7).and.fam(i9)==-fam(i8).and.-fam(i2)==0).or. 
     .    (-fam(i1)==-fam(i8).and.fam(i9)==-fam(i7).and.-fam(i2)==0).or. 

     .    (-fam(i1)== fam(i2).and. fam(i9)==-fam(i8).and.fam(i7)==0).or. 
     .    (-fam(i1)==-fam(i9).and.-fam(i2)==-fam(i8).and.fam(i7)==0).or. 
     .    (-fam(i1)==-fam(i8).and.-fam(i2)==-fam(i9).and.fam(i7)==0).or. 

     .    (-fam(i1)== fam(i2).and. fam(i7)==-fam(i9).and.fam(i8)==0).or. 
     .    (-fam(i1)==-fam(i7).and.-fam(i2)==-fam(i9).and.fam(i8)==0).or. 
     .    (-fam(i1)==-fam(i9).and.-fam(i2)==-fam(i7).and.fam(i8)==0)


     .)    
         if(condition) then
c     q qbar'
            flst_nreal=flst_nreal+1
            if(flst_nreal.gt.maxprocreal) goto 998
            do k=1,nlegreal
               flst_real(k,flst_nreal)=ii(k)
            enddo
         endif
      enddo
      enddo
      enddo
      enddo
      enddo

C     reverse signs for Wm Wm
      if (idvecbos == -24 ) then 
         flst_real(1:2,:) = -flst_real(1:2,:)
         flst_real(7:9,:) = -flst_real(7:9,:)
      endif

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
      end
      
      block data lepmass_data 
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      data lepmass /0.51099891d-3,0.1056583668d0,1.77684d0/
      end 
