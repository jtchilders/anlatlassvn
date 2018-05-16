      subroutine init_processes
      implicit none
      include 'pwhg_math.h'
      include 'process.inc'
      include 'vbfnlo-files/global.inc'
      real* 8 vbfnloinput
      external vbfnloinput
      
      decaymode = vbfnloinput("DECAYMODE") 
      procID = vbfnloinput("PROC_ID")      
      
      select case(procID)
      
      case(Zjj_l) 
          Write(*,*) "Process ", procID
          write(*,*) "Z production in VBF with leptonic decay"
          write(*,*) "Decaymode: " , decaymode
          if(decaymode.eq.11) then
             write(*,*) "into e+ e-"
          elseif(decaymode.eq.13) then
             write(*,*) "into mu+ mu-"
          else
             write(*,*) "DECAYMODE should be 11 or 13, STOP"
             stop
          endif
          call init_process_Z
      case(Wpjj) 
          Write(*,*) "Process ", procID
          write(*,*) "W+ production in VBF with leptonic decay"
          write(*,*) "Decaymode: " , decaymode
          if(decaymode.eq.11) then
             write(*,*) "into e+ nu_e"
          elseif(decaymode.eq.13) then
             write(*,*) "into mu+ nu_mu"
          else
             write(*,*) "DECAYMODE should be 11 or 13, STOP"
             stop
          endif          
          call init_process_Wp
      case(Wmjj) 
          Write(*,*) "Process ", procID
          write(*,*) "W- production in VBF with leptonic decay"
          write(*,*) "Decaymode: " , decaymode
          call init_process_Wm          
          if(decaymode.eq.11) then
             write(*,*) "into e- nu_ebar"
          elseif(decaymode.eq.13) then
             write(*,*) "into mu- nu_mubar"
          else
             write(*,*) "DECAYMODE should be 11 or 13, STOP"
             stop
          endif            
      end select
      
      end
      
      subroutine init_process_Z
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'LesHouches.h'
      include 'process.inc'
      logical debug
      parameter (debug=.false.)
      integer j,i,ii,jj,k
      integer charge3(-6:6)
      data charge3 /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      logical condition
      integer ferm_charge(5)
      character *3 flav(-5:5)
      data (flav(i),i=-5,5) 
     #     /'b~','c~','s~','u~','d~','g','d','u','s','c','b'/
      integer max_flav
      logical emit_Wp_upper,emit_Wm_upper,emit_Wp_lower,emit_Wm_lower
      integer flst_nreal_WW, flst_nborn_WW
      integer ZCC,ZNC
      logical CKM_diag
      integer flst_real_tmp(nlegreal,maxprocreal)
      logical flavequiv
      external flavequiv
      logical tag,newtag
      real * 8 powheginput
      external powheginput



      tag = .true.
      newtag = .true.

      if (.not.tag) then
         do i=1,nlegborn
            do j=1,maxprocborn
               flst_borntags(i,j)=0
            enddo
         enddo
         do i=1,nlegreal
            do j=1,maxprocreal
               flst_realtags(i,j)=0
            enddo
         enddo
      endif


c     index of the first light coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=5
      max_flav = 4
      flst_nborn=0
      flst_nreal=0
      ZCC=21001
      ZNC=21002

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC           BORN GRAPHS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      CKM_diag = .true.

c     WW -> Z case
      condition=.false.

      emit_Wp_upper = .true.
      emit_Wm_upper = .true.
      emit_Wp_lower = .true.
      emit_Wm_lower = .true.
 
      do i=-max_flav,max_flav
         do j=-max_flav,max_flav
            do ii=-max_flav,max_flav
               do jj=-max_flav,max_flav               
                  if (.not.((i.eq.0).or.(j.eq.0).or.
     #                    (ii.eq.0).or.(jj.eq.0))) then ! NOT a gluon
                     ferm_charge(1) = charge3(i)
                     ferm_charge(2) = charge3(j)
                     ferm_charge(3) = charge3(ii)
                     ferm_charge(4) = charge3(jj)

                     if (CKM_diag) then
                        emit_Wp_upper = (ii.eq.i-1)
                        emit_Wm_upper = (ii.eq.i+1)
                        emit_Wp_lower = (jj.eq.j-1)
                        emit_Wm_lower = (jj.eq.j+1)                 
                     endif

                     condition = 
c     W+ emission from upper leg                        
     #                    (((ferm_charge(1)-(ferm_charge(3)+3)
     #                    .eq.0).and.(emit_Wp_upper)) .and.
c     W- emission from lower leg                        
     #                    ((ferm_charge(2)-(ferm_charge(4)-3)
     #                    .eq.0).and.(emit_Wm_lower)))
     #                    .or.
c     W- emission from upper leg                        
     #                    (((ferm_charge(1)-(ferm_charge(3)-3)
     #                    .eq.0).and.(emit_Wm_upper)) .and.
c     W+ emission from lower leg                        
     #                    ((ferm_charge(2)-(ferm_charge(4)+3)
     #                    .eq.0).and.(emit_Wp_lower)))
                     
                     if (condition) then
                        flst_nborn=flst_nborn+1
                        if(flst_nborn.gt.maxprocborn) goto 999
                        flst_born(1,flst_nborn)=i
                        flst_born(2,flst_nborn)=j
                        flst_born(3,flst_nborn)=ZCC
                        flst_born(4,flst_nborn)=ii
                        flst_born(5,flst_nborn)=jj
                        if (tag) then
                           flst_borntags(1,flst_nborn)=1
                           flst_borntags(2,flst_nborn)=2
                           flst_borntags(3,flst_nborn)=0
                           flst_borntags(4,flst_nborn)=4
                           flst_borntags(5,flst_nborn)=5
                           if (newtag) then
                              flst_borntags(1,flst_nborn)=1
                              flst_borntags(2,flst_nborn)=2
                              flst_borntags(3,flst_nborn)=0
                              flst_borntags(4,flst_nborn)=1 !4
                              flst_borntags(5,flst_nborn)=2 !5
                           endif
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
      
      if (debug) then
         write(*,*) ' born processes: CC ',flst_nborn
         do j=1,flst_nborn
            write(*,*) 'proc ',j,' ', (flst_born(k,j),k=1,nlegborn)
         enddo
      endif
      flst_nborn_WW = flst_nborn
      
c     NC case
      do i=-max_flav,max_flav
         do j=-max_flav,max_flav
            if (.not.((i.eq.0).or.(j.eq.0))) then
               ferm_charge(1) = charge3(i)
               ferm_charge(2) = charge3(j)
               ferm_charge(3) = charge3(i)
               ferm_charge(4) = charge3(j)
                     
               flst_nborn=flst_nborn+1
               if(flst_nborn.gt.maxprocborn) goto 999
               flst_born(1,flst_nborn)=i
               flst_born(2,flst_nborn)=j
               flst_born(3,flst_nborn)=ZNC
               flst_born(4,flst_nborn)=i
               flst_born(5,flst_nborn)=j
               if (tag) then
                  flst_borntags(1,flst_nborn)=1
                  flst_borntags(2,flst_nborn)=2
                  flst_borntags(3,flst_nborn)=0
                  flst_borntags(4,flst_nborn)=4
                  flst_borntags(5,flst_nborn)=5
                  if (newtag) then
                     flst_borntags(1,flst_nborn)=1
                     flst_borntags(2,flst_nborn)=2
                     flst_borntags(3,flst_nborn)=0
                     flst_borntags(4,flst_nborn)=1 
                     flst_borntags(5,flst_nborn)=2 
                  endif
               endif
            endif
         enddo
      enddo
      
      if (debug) then
         write(*,*) ' born processes: NC ',flst_nborn
         do j=flst_nborn_WW+1,flst_nborn
            write(*,*) 'proc ',j,' ', (flst_born(k,j),k=1,nlegborn)
         enddo
      endif
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC                REAL GRAPHS    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     WW -> Z case
c     q q -> Z q q g
      do i=-max_flav,max_flav
         do j=-max_flav,max_flav
            do ii=-max_flav,max_flav
               do jj=-max_flav,max_flav               
                  if (.not.((i.eq.0).or.(j.eq.0).or.
     #                    (ii.eq.0).or.(jj.eq.0))) then ! NOT a gluon
                     ferm_charge(1) = charge3(i)
                     ferm_charge(2) = charge3(j)
                     ferm_charge(3) = charge3(ii)
                     ferm_charge(4) = charge3(jj)

                     if (CKM_diag) then
                        emit_Wp_upper = (ii.eq.i-1)
                        emit_Wm_upper = (ii.eq.i+1)
                        emit_Wp_lower = (jj.eq.j-1)
                        emit_Wm_lower = (jj.eq.j+1)                 
                     endif
                     
                     condition = 
c     W+ emission from upper leg                        
     #                    (((ferm_charge(1)-(ferm_charge(3)+3)
     #                    .eq.0).and.(emit_Wp_upper)) .and.
c     W- emission from lower leg                        
     #                    ((ferm_charge(2)-(ferm_charge(4)-3)
     #                    .eq.0).and.(emit_Wm_lower)))
     #                    .or.
c     W- emission from upper leg                        
     #                    (((ferm_charge(1)-(ferm_charge(3)-3)
     #                    .eq.0).and.(emit_Wm_upper)) .and.
c     W+ emission from lower leg                        
     #                    ((ferm_charge(2)-(ferm_charge(4)+3)
     #                    .eq.0).and.(emit_Wp_lower)))
                     
                     if (condition) then
                        flst_nreal=flst_nreal+1
                        if(flst_nreal.gt.maxprocreal) goto 998
                        flst_real(1,flst_nreal)=i
                        flst_real(2,flst_nreal)=j
                        flst_real(3,flst_nreal)=ZCC
                        flst_real(4,flst_nreal)=ii
                        flst_real(5,flst_nreal)=jj
                        flst_real(6,flst_nreal)=0 ! gluon
                        if (tag) then
                           flst_realtags(1,flst_nreal)=1
                           flst_realtags(2,flst_nreal)=2
                           flst_realtags(3,flst_nreal)=0
                           flst_realtags(4,flst_nreal)=4
                           flst_realtags(5,flst_nreal)=5
                           flst_realtags(6,flst_nreal)=0
                           if (newtag) then
                              flst_realtags(1,flst_nreal)=1
                              flst_realtags(2,flst_nreal)=2
                              flst_realtags(3,flst_nreal)=0
                              flst_realtags(4,flst_nreal)=1 !4
                              flst_realtags(5,flst_nreal)=2 !5
                              flst_realtags(6,flst_nreal)=0
                           endif
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
c     g q -> Z q q q
c     loop on only HALF of the incoming upper-line quark, not to double count!
c     In fact, the real-radiation term contains TWO Feynman diagrams.
      do i=1,max_flav
         do j=-max_flav,max_flav
            do ii=-max_flav,max_flav
               do jj=-max_flav,max_flav               
                  if (.not.((i.eq.0).or.(j.eq.0).or.
     #                    (ii.eq.0).or.(jj.eq.0))) then ! NOT a gluon
                     ferm_charge(1) = 0
                     ferm_charge(2) = charge3(j)
                     ferm_charge(3) = charge3(ii)
                     ferm_charge(4) = charge3(jj)
                     ferm_charge(5) = charge3(-i)

                     if (CKM_diag) then
                        emit_Wp_upper = (ii.eq.i-1)
                        emit_Wm_upper = (ii.eq.i+1)
                        emit_Wp_lower = (jj.eq.j-1)
                        emit_Wm_lower = (jj.eq.j+1)                 
                     endif

                     condition = 
c     W+ emission from upper leg                        
     #                    (((-ferm_charge(5)-(ferm_charge(3)+3)
     #                    .eq.0).and.(emit_Wp_upper)) .and.
c     W- emission from lower leg                        
     #                    ((ferm_charge(2)-(ferm_charge(4)-3)
     #                    .eq.0).and.(emit_Wm_lower)))
     #                    .or.
c     W- emission from upper leg                        
     #                    (((-ferm_charge(5)-(ferm_charge(3)-3)
     #                    .eq.0).and.(emit_Wm_upper)) .and.
c     W+ emission from lower leg                        
     #                    ((ferm_charge(2)-(ferm_charge(4)+3)
     #                    .eq.0).and.(emit_Wp_lower)))
                     
                     if (condition) then
                        flst_nreal=flst_nreal+1
                        if(flst_nreal.gt.maxprocreal) goto 998
                        flst_real(1,flst_nreal)=0 ! gluon
                        flst_real(2,flst_nreal)=j
                        flst_real(3,flst_nreal)=ZCC
                        flst_real(4,flst_nreal)=ii
                        flst_real(5,flst_nreal)=jj
                        flst_real(6,flst_nreal)=-i
                        if (tag) then
                           flst_realtags(1,flst_nreal)=0
                           flst_realtags(2,flst_nreal)=2
                           flst_realtags(3,flst_nreal)=0
                           flst_realtags(4,flst_nreal)=4
                           flst_realtags(5,flst_nreal)=5
                           flst_realtags(6,flst_nreal)=1
                           if (newtag) then
                              flst_realtags(1,flst_nreal)=0
                              flst_realtags(2,flst_nreal)=2
                              flst_realtags(3,flst_nreal)=0
                              flst_realtags(4,flst_nreal)=1 !4
                              flst_realtags(5,flst_nreal)=2 !5
                              flst_realtags(6,flst_nreal)=1
                           endif
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo

c     q g -> Z q q q
c     loop on only HALF of the incoming lower-line quark, not to double count!
c     In fact, the real-radiation term contains TWO Feynman diagrams.
      do i=-max_flav,max_flav
         do j=1,max_flav             
            do ii=-max_flav,max_flav
               do jj=-max_flav,max_flav               
                  if (.not.((i.eq.0).or.(j.eq.0).or.
     #                    (ii.eq.0).or.(jj.eq.0))) then ! NOT a gluon
                     ferm_charge(1) = charge3(i)
                     ferm_charge(2) = 0
                     ferm_charge(3) = charge3(ii)
                     ferm_charge(4) = charge3(jj)
                     ferm_charge(5) = charge3(-j)

                     if (CKM_diag) then
                        emit_Wp_upper = (ii.eq.i-1)
                        emit_Wm_upper = (ii.eq.i+1)
                        emit_Wp_lower = (jj.eq.j-1)
                        emit_Wm_lower = (jj.eq.j+1)                 
                     endif
                     
                     condition = 
c     W+ emission from upper leg                        
     #                    (((ferm_charge(1)-(ferm_charge(3)+3)
     #                    .eq.0).and.(emit_Wp_upper)) .and.
c     W- emission from lower leg                        
     #                    ((-ferm_charge(5)-(ferm_charge(4)-3)
     #                    .eq.0).and.(emit_Wm_lower)))
     #                    .or.
c     W- emission from upper leg                        
     #                    (((ferm_charge(1)-(ferm_charge(3)-3)
     #                    .eq.0).and.(emit_Wm_upper)) .and.
c     W+ emission from lower leg                        
     #                    ((-ferm_charge(5)-(ferm_charge(4)+3)
     #                    .eq.0).and.(emit_Wp_lower)))
                     if (condition) then
                        flst_nreal=flst_nreal+1
                        if(flst_nreal.gt.maxprocreal) goto 998
                        flst_real(1,flst_nreal)=i
                        flst_real(2,flst_nreal)=0 ! gluon
                        flst_real(3,flst_nreal)=ZCC
                        flst_real(4,flst_nreal)=ii
                        flst_real(5,flst_nreal)=jj
                        flst_real(6,flst_nreal)=-j
                        if (tag) then
                           flst_realtags(1,flst_nreal)=1
                           flst_realtags(2,flst_nreal)=0
                           flst_realtags(3,flst_nreal)=0
                           flst_realtags(4,flst_nreal)=4
                           flst_realtags(5,flst_nreal)=5
                           flst_realtags(6,flst_nreal)=2
                           if (newtag) then
                              flst_realtags(1,flst_nreal)=1
                              flst_realtags(2,flst_nreal)=0
                              flst_realtags(3,flst_nreal)=0
                              flst_realtags(4,flst_nreal)=1 !4
                              flst_realtags(5,flst_nreal)=2 !5
                              flst_realtags(6,flst_nreal)=2
                           endif                           
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo              


      if (debug) then
         write(*,*) ' real processes: CC ',flst_nreal
         do j=1,flst_nreal
            write(*,*) 'proc ',j,' ',(flst_real(k,j),k=1,nlegreal)
         enddo
      endif
      flst_nreal_WW = flst_nreal


c     NC case
c     q q -> Z q q g
      do i=-max_flav,max_flav
         do j=-max_flav,max_flav
            if (.not.((i.eq.0).or.(j.eq.0))) then
               ferm_charge(1) = charge3(i)
               ferm_charge(2) = charge3(j)
               ferm_charge(3) = charge3(i)
               ferm_charge(4) = charge3(j)
                     
               flst_nreal=flst_nreal+1
               
               if(flst_nreal.gt.maxprocreal) goto 998
               flst_real(1,flst_nreal)=i
               flst_real(2,flst_nreal)=j
               flst_real(3,flst_nreal)=ZNC
               flst_real(4,flst_nreal)=i
               flst_real(5,flst_nreal)=j
               flst_real(6,flst_nreal)=0 ! gluon
               if (tag) then
                  flst_realtags(1,flst_nreal)=1
                  flst_realtags(2,flst_nreal)=2
                  flst_realtags(3,flst_nreal)=0
                  flst_realtags(4,flst_nreal)=4
                  flst_realtags(5,flst_nreal)=5
                  flst_realtags(6,flst_nreal)=0
                  if (newtag) then
                     flst_realtags(1,flst_nreal)=1
                     flst_realtags(2,flst_nreal)=2
                     flst_realtags(3,flst_nreal)=0
                     flst_realtags(4,flst_nreal)=1 !4
                     flst_realtags(5,flst_nreal)=2 !5
                     flst_realtags(6,flst_nreal)=0
                  endif
               endif
            endif
         enddo
      enddo
c     g q -> Z q q q
c     loop on only HALF of the incoming upper-line quark, not to double count!
c     In fact, the real-radiation term contains TWO Feynman diagrams.
      do i=1,max_flav
         do j=-max_flav,max_flav
            if (.not.((i.eq.0).or.(j.eq.0))) then
               ferm_charge(1) = 0
               ferm_charge(2) = charge3(j)
               ferm_charge(3) = charge3(i)
               ferm_charge(4) = charge3(j)
               ferm_charge(5) = -charge3(i)
               
               flst_nreal=flst_nreal+1

               if(flst_nreal.gt.maxprocreal) goto 998
               flst_real(1,flst_nreal)=0 ! gluon
               flst_real(2,flst_nreal)=j
               flst_real(3,flst_nreal)=ZNC
               flst_real(4,flst_nreal)=i
               flst_real(5,flst_nreal)=j
               flst_real(6,flst_nreal)=-i
               if (tag) then
                  flst_realtags(1,flst_nreal)=0
                  flst_realtags(2,flst_nreal)=2
                  flst_realtags(3,flst_nreal)=0
                  flst_realtags(4,flst_nreal)=4
                  flst_realtags(5,flst_nreal)=5
                  flst_realtags(6,flst_nreal)=1
                  if (newtag) then
                     flst_realtags(1,flst_nreal)=0
                     flst_realtags(2,flst_nreal)=2
                     flst_realtags(3,flst_nreal)=0
                     flst_realtags(4,flst_nreal)=1 !4
                     flst_realtags(5,flst_nreal)=2 !5
                     flst_realtags(6,flst_nreal)=1
                  endif
               endif
            endif
         enddo
      enddo

c     q g -> Z q q q
c     loop on only HALF of the incoming lower-line quark, not to double count!
c     In fact, the real-radiation term contains TWO Feynman diagrams.
      do i=-max_flav,max_flav
         do j=1,max_flav  
            if (.not.((i.eq.0).or.(j.eq.0))) then
               ferm_charge(1) = charge3(i)
               ferm_charge(2) = 0
               ferm_charge(3) = charge3(i)
               ferm_charge(4) = charge3(j)
               ferm_charge(5) = -charge3(j)
               
               flst_nreal=flst_nreal+1

               if(flst_nreal.gt.maxprocreal) goto 998
               flst_real(1,flst_nreal)=i
               flst_real(2,flst_nreal)=0 ! gluon
               flst_real(3,flst_nreal)=ZNC
               flst_real(4,flst_nreal)=i
               flst_real(5,flst_nreal)=j
               flst_real(6,flst_nreal)=-j
               if (tag) then
                  flst_realtags(1,flst_nreal)=1
                  flst_realtags(2,flst_nreal)=0
                  flst_realtags(3,flst_nreal)=0
                  flst_realtags(4,flst_nreal)=4
                  flst_realtags(5,flst_nreal)=5
                  flst_realtags(6,flst_nreal)=2
                  if (newtag) then
                     flst_realtags(1,flst_nreal)=1
                     flst_realtags(2,flst_nreal)=0
                     flst_realtags(3,flst_nreal)=0
                     flst_realtags(4,flst_nreal)=1 !4
                     flst_realtags(5,flst_nreal)=2 !5
                     flst_realtags(6,flst_nreal)=2
                  endif
               endif
            endif
         enddo
      enddo              

      if (debug) then
         write(*,*) ' real processes: NC ',flst_nreal-flst_nreal_WW
         do j=flst_nreal_WW+1,flst_nreal
            write(*,*) 'proc ',j-flst_nreal_WW,' ',
     #           (flst_real(k,j),k=1,nlegreal)
         enddo
      endif

c      stop
      
      do i=1,flst_nborn
        flst_born(6,i)=flst_born(5,i)
        flst_born(5,i)=flst_born(4,i)
        flst_born(4,i)=decaymode
        flst_born(3,i)=-decaymode
        
        flst_borntags(6,i)=flst_borntags(5,i)
        flst_borntags(5,i)=flst_borntags(4,i)
        flst_borntags(4,i)=0
        flst_borntags(3,i)=0        
      enddo

      do i=1,flst_nreal
        flst_real(7,i)=flst_real(6,i)
        flst_real(6,i)=flst_real(5,i)
        flst_real(5,i)=flst_real(4,i)
        flst_real(4,i)=decaymode
        flst_real(3,i)=-decaymode
        
        flst_realtags(7,i)=flst_realtags(6,i)
        flst_realtags(6,i)=flst_realtags(5,i)
        flst_realtags(5,i)=flst_realtags(4,i)
        flst_realtags(4,i)=0
        flst_realtags(3,i)=0      
      enddo
      
         
         
      return
 998  write(*,*) 'init_processes: increase maxprocreal:',maxprocreal
      stop
 999  write(*,*) 'init_processes: increase maxprocborn:',maxprocborn
      stop
      end
      

      subroutine init_process_Wp
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'LesHouches.h'
      include 'process.inc'
      logical debug
      parameter (debug=.false.)
      integer j,i,ii,jj,k,l
      integer charge3(-6:6)
      data charge3 /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      logical condition
      integer ferm_charge(5)
      character *3 flav(-5:5)
      data (flav(i),i=-5,5) 
     #     /'b~','c~','s~','u~','d~','g','d','u','s','c','b'/
      integer max_flav
      logical emit_Wp_upper,emit_Wm_upper,emit_Wp_lower,emit_Wm_lower
      integer flst_nreal_WW, flst_nborn_WW
      logical CKM_diag
      integer flst_real_tmp(nlegreal,maxprocreal), 
     # flst_born_tmp(nlegborn,maxprocborn)
      logical flavequiv
      external flavequiv
      logical tag,newtag
      integer hdecaymode
      real * 8 powheginput
      external powheginput
      integer lepton, neutrino
      
      if(decaymode.eq.11) then 
         lepton = -11
         neutrino = 12
      
      elseif(decaymode.eq.13) then
          lepton = -13
         neutrino = 14     
      else
        write(*,*) "Error in init_process_Wp, DECAYMODE"
        stop
      endif


      tag = .true.
      newtag = .true.

      if (.not.tag) then
         do i=1,nlegborn
            do j=1,maxprocborn
               flst_borntags(i,j)=0
            enddo
         enddo
         do i=1,nlegreal
            do j=1,maxprocreal
               flst_realtags(i,j)=0
            enddo
         enddo
      endif


c     index of the first light coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=5
      max_flav = 4
      flst_nborn=0
      flst_nreal=0
      j=0

ccccccupper quarkline unchanged
      do i=-4,4
        do k=-3,-1,2
        if((i.ne.0).and.(k.ne.0)) then
        
               j=j+1
               flst_born(1,j)=i
               flst_born(2,j)=k
               flst_born(3,j)=neutrino
               flst_born(4,j)=lepton
               flst_born(5,j)=i  
               flst_born(6,j)=k-1

               flst_borntags(1,j)=1
               flst_borntags(2,j)=2
               flst_borntags(3,j)=0
               flst_borntags(4,j)=0
               flst_borntags(5,j)=1  
               flst_borntags(6,j)=2
               
        endif       
        enddo
        do k=2,4,2     
           if((i.ne.0).and.(k.ne.0)) then

               j=j+1
               flst_born(1,j)=i
               flst_born(2,j)=k
               flst_born(3,j)=neutrino
               flst_born(4,j)=lepton
               flst_born(5,j)=i  
               flst_born(6,j)=k-1  
               
               flst_borntags(1,j)=1
               flst_borntags(2,j)=2
               flst_borntags(3,j)=0
               flst_borntags(4,j)=0
               flst_borntags(5,j)=1  
               flst_borntags(6,j)=2               
               

           endif                 
        enddo  
      enddo  
      
cccccclower quarkline unchanged
      do i=-4,4
        do k=-3,-1,2
        if((i.ne.0).and.(k.ne.0)) then

               j=j+1
               flst_born(1,j)=k
               flst_born(2,j)=i
               flst_born(3,j)=neutrino
               flst_born(4,j)=lepton
               flst_born(5,j)=k-1  
               flst_born(6,j)=i


               flst_borntags(1,j)=1
               flst_borntags(2,j)=2
               flst_borntags(3,j)=0
               flst_borntags(4,j)=0
               flst_borntags(5,j)=1  
               flst_borntags(6,j)=2
        endif       
        enddo
        do k=2,4,2     
           if((i.ne.0).and.(k.ne.0)) then
               j=j+1
               flst_born(1,j)=k
               flst_born(2,j)=i
               flst_born(3,j)=neutrino
               flst_born(4,j)=lepton
               flst_born(5,j)=k-1  
               flst_born(6,j)=i    
               
               flst_borntags(1,j)=1
               flst_borntags(2,j)=2
               flst_borntags(3,j)=0
               flst_borntags(4,j)=0
               flst_borntags(5,j)=1  
               flst_borntags(6,j)=2               
               
           endif                 
        enddo  
      enddo        
      flst_nborn=j
      
      
      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif


      k=0

      do j=1, flst_nborn

               k=k+1
               flst_real(1,k)=flst_born(1,j)
               flst_real(2,k)=flst_born(2,j)
               flst_real(3,k)=flst_born(3,j)
               flst_real(4,k)=flst_born(4,j)
               flst_real(5,k)=flst_born(5,j) 
               flst_real(6,k)=flst_born(6,j)         
               flst_real(7,k)=0
               
               flst_realtags(1,k)=flst_borntags(1,j)
               flst_realtags(2,k)=flst_borntags(2,j)
               flst_realtags(3,k)=flst_borntags(3,j)
               flst_realtags(4,k)=flst_borntags(4,j)
               flst_realtags(5,k)=flst_borntags(5,j) 
               flst_realtags(6,k)=flst_borntags(6,j)         
               flst_realtags(7,k)=0               
               
      enddo


      do j=-4,-1
        do l=-3,-1,2
               k=k+1        ! upper quarkline q q -> j only (anti)quark
               flst_real(1,k)=0   
               flst_real(2,k)=l
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=j
               flst_real(6,k)=l-1
               flst_real(7,k)=-j

               flst_realtags(1,k)=0
               flst_realtags(2,k)=2
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1 
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=1                   
               

         enddo
         
        do l=2,4,2
                       k=k+1
               flst_real(1,k)=0   
               flst_real(2,k)=l
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=j
               flst_real(6,k)=l-1
               flst_real(7,k)=-j

               flst_realtags(1,k)=0
               flst_realtags(2,k)=2
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1 
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=1   

         enddo   
         
        do l=-3,-1,2
               k=k+1        ! lower quarkline q q -> j only (anti)quark
               flst_real(1,k)=l               
               flst_real(2,k)=0   
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=l-1               
               flst_real(6,k)=j
               flst_real(7,k)=-j

               flst_realtags(1,k)=1
               flst_realtags(2,k)=0
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0               
               flst_realtags(5,k)=1        
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=2                   
               

         enddo
         
        do l=2,4,2
                       k=k+1
               flst_real(1,k)=l
               flst_real(2,k)=0   
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=l-1               
               flst_real(6,k)=j
               flst_real(7,k)=-j

               flst_realtags(1,k)=1
               flst_realtags(2,k)=0
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1 
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=2   

         enddo           
         
         
       enddo



        do j=-4,4
            if (j.ne.0) then
            do l=-3,-1,2        !upper quarkline qq-> only 1 -2
                           k=k+1
               flst_real(1,k)=j                           
               flst_real(2,k)=0

               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               
               flst_real(5,k)=j 
               flst_real(6,k)=l-1
               flst_real(7,k)=-l


               flst_realtags(1,k)=1
               flst_realtags(2,k)=0
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1 
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=2                 

            enddo
            
            
            do l=-3,-1,2        !lower quarkline qq-> only 1 -2
                           k=k+1
               flst_real(1,k)=0
               flst_real(2,k)=j                           
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=l-1               
               flst_real(6,k)=j 
               flst_real(7,k)=-l


               flst_realtags(1,k)=0
               flst_realtags(2,k)=2
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1                
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=1                 

            enddo            
            

            
            endif
         enddo   


      flst_nreal=k


      ii=1

      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
         enddo  
!          stop
      endif
      
      return
 998  write(*,*) 'init_processes: increase maxprocreal:',maxprocreal
      stop
 999  write(*,*) 'init_processes: increase maxprocborn:',maxprocborn
      stop
      end
      
      subroutine init_process_Wm
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'LesHouches.h'
      include 'process.inc'
      logical debug
      parameter (debug=.false.)
      integer j,i,ii,jj,k,l
      integer charge3(-6:6)
      data charge3 /-2,1,-2,1,-2,1,0,-1,2,-1,2,-1,2/
      logical condition
      integer ferm_charge(5)
      character *3 flav(-5:5)
      data (flav(i),i=-5,5) 
     #     /'b~','c~','s~','u~','d~','g','d','u','s','c','b'/
      integer max_flav
      logical emit_Wp_upper,emit_Wm_upper,emit_Wp_lower,emit_Wm_lower
      integer flst_nreal_WW, flst_nborn_WW
      logical CKM_diag
      integer flst_real_tmp(nlegreal,maxprocreal), 
     # flst_born_tmp(nlegborn,maxprocborn)
      logical flavequiv
      external flavequiv
      logical tag,newtag
      integer hdecaymode
      real * 8 powheginput
      external powheginput
      integer lepton,neutrino


      if(decaymode.eq.11) then 
         lepton = 11
         neutrino = -12
      
      elseif(decaymode.eq.13) then
          lepton = 13
         neutrino = -14     
      else
        write(*,*) "Error in init_process_Wp, DECAYMODE"
        stop
      endif
      tag = .true.
      newtag = .true.

      if (.not.tag) then
         do i=1,nlegborn
            do j=1,maxprocborn
               flst_borntags(i,j)=0
            enddo
         enddo
         do i=1,nlegreal
            do j=1,maxprocreal
               flst_realtags(i,j)=0
            enddo
         enddo
      endif


c     index of the first light coloured particle in the final state
c     (all subsequent particles are coloured)
      flst_lightpart=5
      max_flav = 4
      flst_nborn=0
      flst_nreal=0
      j=0

ccccccupper quarkline unchanged
      do i=-4,4
        do k=-3,-1,2
        if((i.ne.0).and.(k.ne.0)) then
        

               j=j+1
               flst_born(1,j)=i
               flst_born(2,j)=k-1
               flst_born(3,j)=neutrino
               flst_born(4,j)=lepton
               flst_born(5,j)=i  
               flst_born(6,j)=k

               flst_borntags(1,j)=1
               flst_borntags(2,j)=2
               flst_borntags(3,j)=0
               flst_borntags(4,j)=0
               flst_borntags(5,j)=1  
               flst_borntags(6,j)=2
               

        endif       
        enddo
        do k=2,4,2     
           if((i.ne.0).and.(k.ne.0)) then

               j=j+1
               flst_born(1,j)=i
               flst_born(2,j)=k-1
               flst_born(3,j)=neutrino
               flst_born(4,j)=lepton
               flst_born(5,j)=i  
               flst_born(6,j)=k  
               
               flst_borntags(1,j)=1
               flst_borntags(2,j)=2
               flst_borntags(3,j)=0
               flst_borntags(4,j)=0
               flst_borntags(5,j)=1  
               flst_borntags(6,j)=2               
               

           endif                 
        enddo  
      enddo  
      
cccccclower quarkline unchanged
      do i=-4,4
        do k=-3,-1,2
        if((i.ne.0).and.(k.ne.0)) then

               j=j+1
               flst_born(1,j)=k-1
               flst_born(2,j)=i
               flst_born(3,j)=neutrino
               flst_born(4,j)=lepton
               flst_born(5,j)=k  
               flst_born(6,j)=i

               flst_borntags(1,j)=1
               flst_borntags(2,j)=2
               flst_borntags(3,j)=0
               flst_borntags(4,j)=0
               flst_borntags(5,j)=1  
               flst_borntags(6,j)=2
        endif       
        enddo
        do k=2,4,2     
           if((i.ne.0).and.(k.ne.0)) then
               j=j+1
               flst_born(1,j)=k-1
               flst_born(2,j)=i
               flst_born(3,j)=neutrino
               flst_born(4,j)=lepton
               flst_born(5,j)=k  
               flst_born(6,j)=i    
               
               flst_borntags(1,j)=1
               flst_borntags(2,j)=2
               flst_borntags(3,j)=0
               flst_borntags(4,j)=0
               flst_borntags(5,j)=1  
               flst_borntags(6,j)=2               
               
           endif                 
        enddo  
      enddo        
      flst_nborn=j
      
      
      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif


      k=0

      do j=1, flst_nborn

               k=k+1
               flst_real(1,k)=flst_born(1,j)
               flst_real(2,k)=flst_born(2,j)
               flst_real(3,k)=flst_born(3,j)
               flst_real(4,k)=flst_born(4,j)
               flst_real(5,k)=flst_born(5,j) 
               flst_real(6,k)=flst_born(6,j)         
               flst_real(7,k)=0
               
               flst_realtags(1,k)=flst_borntags(1,j)
               flst_realtags(2,k)=flst_borntags(2,j)
               flst_realtags(3,k)=flst_borntags(3,j)
               flst_realtags(4,k)=flst_borntags(4,j)
               flst_realtags(5,k)=flst_borntags(5,j) 
               flst_realtags(6,k)=flst_borntags(6,j)         
               flst_realtags(7,k)=0               
               
 
      enddo


      do j=-4,-1
        do l=-3,-1,2
               k=k+1        ! upper quarkline q q -> j only (anti)quark
               flst_real(1,k)=0   
               flst_real(2,k)=l-1
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=j
               flst_real(6,k)=l
               flst_real(7,k)=-j

               flst_realtags(1,k)=0
               flst_realtags(2,k)=2
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1 
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=1                   
               

         enddo
         
        do l=2,4,2
                       k=k+1
               flst_real(1,k)=0   
               flst_real(2,k)=l-1
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=j
               flst_real(6,k)=l
               flst_real(7,k)=-j

               flst_realtags(1,k)=0
               flst_realtags(2,k)=2
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1 
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=1   

         enddo   
         
        do l=-3,-1,2
               k=k+1        ! lower quarkline q q -> j only (anti)quark
               flst_real(1,k)=l-1              
               flst_real(2,k)=0   
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=l               
               flst_real(6,k)=j
               flst_real(7,k)=-j

               flst_realtags(1,k)=1
               flst_realtags(2,k)=0
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0               
               flst_realtags(5,k)=1        
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=2                   
               

         enddo
         
        do l=2,4,2
                       k=k+1
               flst_real(1,k)=l-1
               flst_real(2,k)=0   
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=l               
               flst_real(6,k)=j
               flst_real(7,k)=-j

               flst_realtags(1,k)=1
               flst_realtags(2,k)=0
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1 
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=2   

         enddo           
         
         
       enddo



        do j=-4,4
            if (j.ne.0) then
            do l=-3,-1,2        !upper quarkline qq-> only 1 -2
                           k=k+1
               flst_real(1,k)=j                           
               flst_real(2,k)=0

               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=j 
               flst_real(6,k)=l
               flst_real(7,k)=-(l-1)


               flst_realtags(1,k)=1
               flst_realtags(2,k)=0
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1 
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=2                 

            enddo
            
            
            do l=-3,-1,2        !lower quarkline qq-> only 1 -2
                           k=k+1
               flst_real(1,k)=0
               flst_real(2,k)=j                           
               flst_real(3,k)=neutrino
               flst_real(4,k)=lepton
               flst_real(5,k)=l               
               flst_real(6,k)=j 
               flst_real(7,k)=-(l-1)


               flst_realtags(1,k)=0
               flst_realtags(2,k)=2
               flst_realtags(3,k)=0
               flst_realtags(4,k)=0
               flst_realtags(5,k)=1                
               flst_realtags(6,k)=2        
               flst_realtags(7,k)=1                 

            enddo            
            

            
            endif
         enddo   


      flst_nreal=k


      ii=1

      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
         enddo  
      endif
      
      return
 998  write(*,*) 'init_processes: increase maxprocreal:',maxprocreal
      stop
 999  write(*,*) 'init_processes: increase maxprocborn:',maxprocborn
      stop
      end
