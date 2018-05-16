!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECopp.f90.

!         Procedures for residues & cuts

module qpopp
  use types; use consts_qp; use define_ampl
  use qppol_int; use qpaux_functions; use qpvvn
  use match1; use qpamplitude; use qpglobal 
  use qpcut_utils 
  implicit none 
  private 

  public :: resid5, resid4, resid3, resid2
  public :: pentcut, quadcut, tripcut, doubcut


  logical :: use_mynivecs = .true. 

contains

  ! Routine computing the pentagon residue. Arguments are       
  ! Ds is the dimensionlity of the loop momentum (only 5 allowed now!)
  ! lv is the 5-dimensional loop momentum (solution of unitarity eqs.)
  ! k1,k2,k3,k4 : 4 independent combinations of external momenta entering the cut
  ! lc5: cut labeling 
  ! fc5: specifies the flavour 
  ! ij : specifies the external momentum ordering  
  ! res: total residue summed over internal states 
  subroutine resid5(Ds,lv,k1,k2,k3,k4,l5c,f5c,ij,res)
    integer, intent(in)        :: Ds
    integer, intent(in)        :: l5c(5)
    character, intent(in)      :: f5c(5)*3
    integer, intent(in)        :: ij
    complex(qp), intent(in)  :: k1(4),k2(4),k3(4),k4(4)
    complex(qp), intent(in)  :: lv(5)
    complex(qp), intent(out) :: res
    ! ------------------------------------------------------------       
    complex(qp) :: q1(5),q2(5),q3(5),q4(5),q5(5)
    complex(qp) :: BPOL1(8,16),POL1(8,16)
    complex(qp) :: BPOL2(8,16),POL2(8,16)
    complex(qp) :: BPOL3(8,16),POL3(8,16)
    complex(qp) :: BPOL4(8,16),POL4(8,16)
    complex(qp) :: BPOL5(8,16),POL5(8,16)
    complex(qp) :: mur1(10,10),mur2(10,10)
    complex(qp) :: mur3(10,10),mur4(10,10)
    complex(qp) :: mur5(10,10)
    complex(qp) :: res6,res8
    integer       :: i,tag_pol, tag_f
    integer       :: j,j1,j2
    integer       :: i1,i2,i3,i4,ia,ib,j3,j4,j5
    integer       :: Nj1,Nj2,Nj3,Nj4,Nj5
    character     :: lab*3, lab1*3, lab2*3

    tag_pol = 0

    do i=1,4
       q1(i)=lv(i)
       q2(i)=lv(i)+k1(i)
       q3(i)=lv(i)+k2(i)
       q4(i)=lv(i)+k3(i)
       q5(i)=lv(i)+k4(i)
    enddo

    !     I assume here that Ds = 5, always, to have a pentuple cut
    q1(5)=lv(5)
    q2(5)=lv(5)
    q3(5)=lv(5)
    q4(5)=lv(5)
    q5(5)=lv(5)

    ! tag_f means that we have a fermion loop 
    if (any(f5c(1:5) .eq. 'glu')) then 
       tag_f = 0
    else
       tag_f = 1
    endif

    if (ferm_loops .or. ferm_loops_Z .or. ferm_loops_Z_sbs) tag_f = 2
    if (extra_ferm_pair_nf) tag_f = 2
    if (gluons_ferm_loops) tag_f = 2
          

    if (Ds == 5) then 

       ! compute polarization for the 5 cut lines
       ! givepol returns Nj, BPOL and POL 

       lab=f5c(1)
       call givepol(lab,q1,6,8,tag_pol,Nj1,BPOL1,POL1)

       lab = f5c(2)
       call givepol(lab,q2,6,8,tag_pol,Nj2,BPOL2,POL2)

       lab= f5c(3)
       call givepol(lab,q3,6,8,tag_pol,Nj3,BPOL3,POL3)

       lab=f5c(4)
       call givepol(lab,q4,6,8,tag_pol,Nj4,BPOL4,POL4)

       lab= f5c(5)
       call givepol(lab,q5,6,8,tag_pol,Nj5,BPOL5,POL5)


       !-------- calculate product of amplitudes
       !-------- first  organize 4 lists that will be arguments 
       !-------- for amplitude procedure 
       ! --------start with the first label on the cut list ll
       !-------- Lab_ex and Lab_in

       res=czero

       i1=l5c(2)-l5c(1)
       i2=l5c(3)-l5c(2)
       i3=l5c(4)-l5c(3)
       i4=l5c(5)-l5c(4)

       lab1=f5c(1)
       lab2=f5c(2)

       ! q1 and q2 are the momenta depending on the loop-momentum (i.e. the not external)
       call ampl(6,8,Nj1,Nj2,POL1,q1,lab1,l5c(1),1,i1,BPOL2,q2,lab2,ij,mur1)



       !----------next loop
       ia = i1+1
       ib = i1+i2

       lab1=f5c(2)
       lab2=f5c(3)

       call ampl(6,8,Nj2,Nj3,POL2,q2,lab1,l5c(2),ia,ib,BPOL3,q3,lab2,ij,mur2)


       !----------next loop
       ia=i1+i2+1
       ib=i1+i2+i3

       lab1=f5c(3)
       lab2=f5c(4)


       call ampl(6,8,Nj3,Nj4,POL3,q3,lab1,l5c(3),ia,ib,BPOL4,q4,lab2,ij,mur3)


       !----------next loop

       ia=i1+i2+i3+1
       ib=i1+i2+i3+i4

       lab1=f5c(4)
       lab2=f5c(5)


       call ampl(6,8,Nj4,Nj5,POL4,q4,lab1,l5c(4),ia,ib,BPOL5,q5,lab2,ij,mur4)



       !----------  next loop
       lab1=f5c(5)
       lab2=f5c(1)

       ia=i1+i2+i3+i4+1
       ib=Npoint


       call ampl(6,8,Nj5,Nj1,POL5,q5,lab1,l5c(5),ia,ib,BPOL1,q1,lab2,ij,mur5)



       !--------sum
       do j1=1,Nj1
          do j2=1,Nj2
             do j3=1,Nj3
                do j4=1,Nj4
                   do j5=1,Nj5
                      res=res + mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j5)*mur5(j5,j1)
                   enddo
                enddo
             enddo
          enddo
       enddo

       !           to multiply by proper power of I:

       do j=1,5
          if (quark_flavour(f5c(j))) then 
             res = ci*res
          else
             res = -ci*res
          endif
       enddo




       ! -- intermediate 6-dim result

       res6 = res

       ! -- 8-dim calculation

       lab=f5c(1)

       call givepol(lab,q1,8,16,tag_pol,Nj1,BPOL1,POL1)

       lab=f5c(2)
       call givepol(lab,q2,8,16,tag_pol,Nj2,BPOL2,POL2)

       lab=f5c(3)
       call givepol(lab,q3,8,16,tag_pol,Nj3,BPOL3,POL3)

       lab=f5c(4)
       call givepol(lab,q4,8,16,tag_pol,Nj4,BPOL4,POL4)

       lab=f5c(5)
       call givepol(lab,q5,8,16,tag_pol,Nj5,BPOL5,POL5)

       res=czero


       !----------
       lab1=f5c(1)
       lab2=f5c(2)

       call ampl(8,16,Nj1,Nj2,POL1,q1,lab1,l5c(1),1,i1,BPOL2,q2,lab2,ij,mur1)


       !------------
       ia =i1+1
       ib = i1+i2

       lab1=f5c(2)
       lab2=f5c(3)


       call ampl(8,16,Nj2,Nj3,POL2,q2,lab1,l5c(2),ia,ib,BPOL3,q3,lab2,ij,mur2)



       !-----------next loop
       ia=i1+i2+1
       ib=i1+i2+i3

       lab1=f5c(3)
       lab2=f5c(4)


       call ampl(8,16,Nj3,Nj4,POL3,q3,lab1,l5c(3),ia,ib,BPOL4,q4,lab2,ij,mur3)


       !---------
       ia=i1+i2+i3+1
       ib=i1+i2+i3+i4

       lab1=f5c(4)
       lab2=f5c(5)

       call ampl(8,16,Nj4,Nj5,POL4,q4,lab1,l5c(4),ia,ib,BPOL5,q5,lab2,ij,mur4)


       !-------
       ia=i1+i2+i3+i4+1
       ib=Npoint

       lab1=f5c(5)
       lab2=f5c(1)


       call ampl(8,16,Nj5,Nj1,POL5,q5,lab1,l5c(5),ia,ib,BPOL1,q1,lab2,ij,mur5)



       do j1=1,Nj1
          do j2=1,Nj2
             do j3=1,Nj3
                do j4=1,Nj4
                   do j5=1,Nj5
                      res=res + mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j5)*mur5(j5,j1)
                   enddo
                enddo
             enddo
          enddo
       enddo

       ! --  to multiply by proper power of I

       do j=1,5
          if (quark_flavour(f5c(j))) then 
             res = ci*res
          else
             res =-ci*res
          endif
       enddo


       ! --  intermediate 8-dim result           

       res8 = res


       !---------- final result 

       if (tag_f.eq.1) then 
          res=two*res6 -res8 
       elseif (tag_f.eq.2) then 
          res=res6 -res8/four 
       else
          res = res6 - res8
       endif

    else
       stop 'resid5: Ds/=5'
    end if ! Ds  

  end subroutine resid5


  subroutine resid4(Ds,lv,k1,k2,k3,l4c,f4c,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res)
    integer, intent(in)        :: N5,N4,N3,N2
    integer, intent(in)        :: Ds, l4c(4)
    character, intent(in)      :: f4c(4)*3
    integer, intent(in)        :: ij
    complex(qp), intent(in)  :: k1(4),k2(4),k3(4)
    complex(qp), intent(in)  :: lv(Ds)
    integer, intent(in)        :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    character, intent(in)      :: F5(:,:)*3,F4(:,:)*3,F3(:,:)*3,F2(:,:)*3
    complex(qp), intent(out) :: res
    ! --------------------------------------------------------------------       
    complex(qp) ::  q1(5),q2(5),q3(5),q4(5)
    complex(qp) ::  BPOL1(8,16),POL1(8,16)
    complex(qp) ::  BPOL2(8,16),POL2(8,16)
    complex(qp) ::  BPOL3(8,16),POL3(8,16)
    complex(qp) ::  BPOL4(8,16),POL4(8,16)
    complex(qp) ::  mur1(10,10),mur2(10,10)
    complex(qp) ::  mur3(10,10),mur4(10,10)
    complex(qp) ::  res6,res8
    complex(qp) ::  krefa(4),lvt(5)
    complex(qp) ::  e(1),vprop(5),prop
    character     ::  lab*3, lab1*3, lab2*3
    integer       ::  Nj1,Nj2,Nj3,Nj4
    integer       ::  n45,lmatch45(N5),im,pos,pos1
    integer       ::  j,j1,j2, tag_f,i,tag_pol
    integer       ::  i1,i2,i3,ia,ib,j3,j4,lpos0(1)

    tag_pol = 0

    do i=1,4
       q1(i)=lv(i)
       q2(i)=lv(i)+k1(i)
       q3(i)=lv(i)+k2(i)
       q4(i)=lv(i)+k3(i)
    enddo

    if (Ds.eq.4) then 
       q1(5)=czero
       q2(5)=czero
       q3(5)=czero
       q4(5)=czero
    else
       q1(5)=lv(5)
       q2(5)=lv(5)
       q3(5)=lv(5)
       q4(5)=lv(5)
    endif

    ! tag_f means that we have a fermion loop 
    if (any(f4c .eq. 'glu')) then 
       tag_f = 0
    else
       tag_f = 1
    endif

    if (ferm_loops .or. ferm_loops_Z .or. ferm_loops_Z_sbs) tag_f = 2
    if (extra_ferm_pair_nf) tag_f = 2
    if (gluons_ferm_loops) tag_f = 2

    !    residues for different space dimensionalities


    if (Ds.eq.4) then 

       lab=f4c(1)
       call givepol(lab,q1,4,4,tag_pol,Nj1,BPOL1,POL1)

       lab=f4c(2)
       call givepol(lab,q2,4,4,tag_pol,Nj2,BPOL2,POL2)

       lab=f4c(3)
       call givepol(lab,q3,4,4,tag_pol,Nj3,BPOL3,POL3)

       lab=f4c(4)
       call givepol(lab,q4,4,4,tag_pol,Nj4,BPOL4,POL4)

       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero


       i1=l4c(2)-l4c(1) 
       i2=l4c(3)-l4c(2)
       i3=l4c(4)-l4c(3)

       !----------

       lab1=f4c(1)
       lab2=f4c(2)



       call ampl(4,4,Nj1,Nj2,POL1,q1,lab1,l4c(1),1,i1,BPOL2,q2,lab2,ij,mur1)

       !-----------

       ia =i1+1
       ib = i1+i2

       lab1=f4c(2)
       lab2=f4c(3)




       call ampl(4,4,Nj2,Nj3,POL2,q2,lab1,l4c(2),ia,ib,BPOL3,q3,lab2,ij,mur2)



       !----------
       ia=i1+i2+1
       ib=i1+i2+i3

       lab1=f4c(3)
       lab2=f4c(4)


       call ampl(4,4,Nj3,Nj4,POL3,q3,lab1,l4c(3),ia,ib,BPOL4,q4,lab2,ij,mur3)





       !-----------new loop
       ia=i1+i2+i3+1
       ib=Npoint

       lab1=f4c(4)
       lab2=f4c(1)




       call ampl(4,4,Nj4,Nj1,POL4,q4,lab1,l4c(4),ia,ib,BPOL1,q1,lab2,ij,mur4)





       do j1=1,Nj1
          do j2=1,Nj2
             do j3=1,Nj3
                do j4=1,Nj4
                   res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j1)
                enddo
             enddo
          enddo
       enddo

       !           to multiply by proper power of I:

       do j=1,4
          if (quark_flavour(f4c(j))) then 
             res = ci*res
          else
             res = (-ci)*res
          endif
       enddo



       !ccccc!ccccc!ccccccc endif for Ds = 4 
    endif



    if (Ds.eq.5) then 

       lab=f4c(1)

       call givepol(lab,q1,6,8,tag_pol,Nj1,BPOL1,POL1)

       lab=f4c(2)

       call givepol(lab,q2,6,8,tag_pol,Nj2,BPOL2,POL2)

       lab=f4c(3)
       call givepol(lab,q3,6,8,tag_pol,Nj3,BPOL3,POL3)

       lab=f4c(4)
       call givepol(lab,q4,6,8,tag_pol,Nj4,BPOL4,POL4)


       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero

       i1=l4c(2)-l4c(1) 
       i2=l4c(3)-l4c(2)
       i3=l4c(4)-l4c(3)

       !---------

       lab1=f4c(1)
       lab2=f4c(2)

       call ampl(6,8,Nj1,Nj2,POL1,q1,lab1,l4c(1),1,i1,BPOL2,q2,lab2,ij,mur1)

       !-----------
       ia =i1+1
       ib = i1+i2

       lab1=f4c(2)
       lab2=f4c(3)


       call ampl(6,8,Nj2,Nj3,POL2,q2,lab1,l4c(2),ia,ib,BPOL3,q3,lab2,ij,mur2)


       !----------
       ia=i1+i2+1
       ib=i1+i2+i3

       lab1=f4c(3)
       lab2=f4c(4)




       call ampl(6,8,Nj3,Nj4,POL3,q3,lab1,l4c(3),ia,ib,BPOL4,q4,lab2,ij,mur3)



       !----------

       ia=i1+i2+i3+1
       ib=Npoint

       lab1=f4c(4)
       lab2=f4c(1)


       call ampl(6,8,Nj4,Nj1,POL4,q4,lab1,l4c(4),ia,ib,BPOL1,q1,lab2,ij,mur4)

       do j1=1,Nj1
          do j2=1,Nj2
             do j3=1,Nj3
                do j4=1,Nj4
                   res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j1)
                enddo
             enddo
          enddo
       enddo

       !           to multiply by proper power of I:

       do j=1,4
          if (quark_flavour(f4c(j))) then 
             res = ci*res
          else
             res = (-ci)*res
          endif
       enddo


       !           intermediate 6-dim result

       res6 = res




       !           8-dim calculation

       lab=f4c(1)
       call givepol(lab,q1,8,16,tag_pol,Nj1,BPOL1,POL1)

       lab=f4c(2)
       call givepol(lab,q2,8,16,tag_pol,Nj2,BPOL2,POL2)

       lab=f4c(3)
       call givepol(lab,q3,8,16,tag_pol,Nj3,BPOL3,POL3)

       lab=f4c(4)
       call givepol(lab,q4,8,16,tag_pol,Nj4,BPOL4,POL4)


       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero


       lab1=f4c(1)
       lab2=f4c(2)

       call ampl(8,16,Nj1,Nj2,POL1,q1,lab1,l4c(1),1,i1,BPOL2,q2,lab2,ij,mur1)


       !---------
       ia =i1+1
       ib = i1+i2

       lab1=f4c(2)
       lab2=f4c(3)

       call ampl(8,16,Nj2,Nj3,POL2,q2,lab1,l4c(2),ia,ib,BPOL3,q3,lab2,ij,mur2)



       !----------
       ia=i1+i2+1
       ib=i1+i2+i3

       lab1=f4c(3)
       lab2=f4c(4)


       call ampl(8,16,Nj3,Nj4,POL3,q3,lab1,l4c(3),ia,ib,BPOL4,q4,lab2,ij,mur3)

       !---------
       ia=i1+i2+i3+1
       ib=Npoint

       lab1=f4c(4)
       lab2=f4c(1)

       call ampl(8,16,Nj4,Nj1,POL4,q4,lab1,l4c(4),ia,ib,BPOL1,q1,lab2,ij,mur4)


       do j1=1,Nj1
          do j2=1,Nj2
             do j3=1,Nj3
                do j4=1,Nj4
                   res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j4)*mur4(j4,j1)
                enddo
             enddo
          enddo
       enddo

       !           to multiply by proper power of I:

       do j=1,4
          if (quark_flavour(f4c(j))) then 
             res = ci*res
          else
             res = (-ci)*res
          endif
       enddo


       res8 = res

       !---------- final result 

       if (tag_f.eq.1) then 
          res=two*res6 -res8 
       elseif (tag_f.eq.2) then  ! this is for internal fermion loop
          res=res6 -res8/four 
       else
          res = res6 - res8
       endif



       !ccccc!ccccc!ccccccc endif for Ds = 5 
    endif



    !     subtracting  the 5-cut contribution

    call match(4,5,l4c,f4c,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,n45,lmatch45)



    do i=1,n45

       im = lmatch45(i)



       j=1

49     if (l4c(1).ne.Lc5(im,j)) then 
          j=j+1
          go to 49
       endif

       pos=j


       call mismatch(4,5,l4c,im,N5,Lc5,N4,Lc4,N3,Lc3,N2,Lc2,lpos0) 

       pos1 = lpos0(1)

       do j=1,1
          e(j)=coeff5(im,j)
       enddo


       do j=1,4
          krefa(j)=propv5(im,4*(pos-1)+j)
       enddo


       do j=1,4
          lvt(j)=q1(j)-krefa(j)
       enddo
       lvt(5)=q1(5)

       do j=1,4 
          vprop(j)=lvt(j)+propv5(im,j+4*(pos1-1))
       enddo
       vprop(5)=lvt(5)

       prop =sc(vprop,vprop)

       prop = prop-mass5(im,pos1)**2
       
       res=res - e(1)*lvt(5)**2/prop 


    enddo

  end subroutine resid4


  subroutine resid3(Ds,lv,k1,k2,l3c,f3c,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res)
    integer, intent(in)        :: Ds
    integer, intent(in)        :: l3c(3), ij
    character, intent(in)      :: f3c(3)*3
    complex(qp), intent(in)  :: k1(4),k2(4),lv(Ds)
    integer, intent(in)        :: N5,N4,N3,N2
    integer, intent(in)        :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    character, intent(in)      :: F5(:,:)*3,F4(:,:)*3,F3(:,:)*3,F2(:,:)*3
    complex(qp), intent(out) :: res
    ! ----------------------------------------------------------------
    complex(qp) :: lvt(5)
    complex(qp) :: v45(4),vprop(5),prop,r1,r2,krefa(4),d(5)
    complex(qp) :: q1(5),q2(5),q3(5)
    complex(qp) :: vne(5)
    complex(qp) :: BPOL1(8,16),POL1(8,16)
    complex(qp) :: BPOL2(8,16),POL2(8,16)
    complex(qp) :: BPOL3(8,16),POL3(8,16)
    complex(qp) :: mur1(10,10),mur2(10,10)
    complex(qp) :: mur3(10,10)
    complex(qp) :: e(1)
    complex(qp) :: res6,res8
    complex(qp) :: prop1,prop2,r22
    integer       :: j,j1,j2,tag_f, lpos0(1)
    integer       :: i1,i2,ia,ib,j3,n34,lmatch34(N4)
    integer       :: lmatch35(N5),i,tag_pol,pos,pos1,pos2
    integer       :: Nj1,Nj2,Nj3,n35,im,lpos1(2)
    character     :: lab*3, lab1*3, lab2*3

    tag_pol=0

    do i=1,4
       q1(i)=lv(i)
       q2(i)=lv(i)+k1(i)
       q3(i)=lv(i)+k2(i)
    enddo

    if (Ds.eq.4) then 
       q1(5)=czero
       q2(5)=czero
       q3(5)=czero
    else
       q1(5)=lv(5)
       q2(5)=lv(5)
       q3(5)=lv(5)
    endif

    !---------required to properly calculate diagrams 
    !---------without virtual gluons

    ! tag_f means that we have a fermion loop 
    if (any(f3c .eq. 'glu')) then 
       tag_f = 0
    else
       tag_f = 1
    endif

    if (ferm_loops .or. ferm_loops_Z .or. ferm_loops_Z_sbs) tag_f = 2
    if (extra_ferm_pair_nf) tag_f = 2
    if (gluons_ferm_loops) tag_f = 2

    !     amplitudes for different space dimensions

    if (Ds.eq.4) then 

       lab=f3c(1)

       call givepol(lab,q1,4,4,tag_pol,Nj1,BPOL1,POL1)

       lab=f3c(2)
       call givepol(lab,q2,4,4,tag_pol,Nj2,BPOL2,POL2)

       lab=f3c(3)
       call givepol(lab,q3,4,4,tag_pol,Nj3,BPOL3,POL3)

       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero

       i1=l3c(2)-l3c(1)
       i2=l3c(3)-l3c(2)

       lab1=f3c(1)
       lab2=f3c(2)

       call ampl(4,4,Nj1,Nj2,POL1,q1,lab1,l3c(1),1,i1,BPOL2,q2,lab2,ij,mur1)


       !----------
       lab1=f3c(2)
       lab2=f3c(3)


       ia =i1+1
       ib = i1+i2


       call ampl(4,4,Nj2,Nj3,POL2,q2,lab1,l3c(2),ia,ib,BPOL3,q3,lab2,ij,mur2)

       !------------
       lab1=f3c(3)
       lab2=f3c(1)


       ia=i1+i2+1
       ib=Npoint

       call ampl(4,4,Nj3,Nj1,POL3,q3,lab1,l3c(3),ia,ib,BPOL1,q1,lab2,ij,mur3)


       do j1=1,Nj1
          do j2=1,Nj2
             do j3=1,Nj3
                res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j1)
             enddo
          enddo
       enddo

       !           to multiply by proper power of I:

       do j=1,3
          if (quark_flavour(f3c(j))) then 
             res = ci*res
          else 
             res = (-ci)*res
          endif
       enddo




       !ccccc!ccccc!ccccccc endif for Ds = 4 
    endif



    if (Ds.eq.5) then 

       lab=f3c(1)

       call givepol(lab,q1,6,8,tag_pol,Nj1,BPOL1,POL1)

       lab=f3c(2)
       call givepol(lab,q2,6,8,tag_pol,Nj2,BPOL2,POL2)

       lab=f3c(3)
       call givepol(lab,q3,6,8,tag_pol,Nj3,BPOL3,POL3)

       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero

       i1=l3c(2)-l3c(1) 
       i2=l3c(3)-l3c(2)


       !-----------
       lab1=f3c(1)
       lab2=f3c(2)


       call ampl(6,8,Nj1,Nj2,POL1,q1,lab1,l3c(1),1,i1,BPOL2,q2,lab2,ij,mur1)



       !-----------

       ia =i1+1
       ib = i1+i2

       lab1=f3c(2)
       lab2=f3c(3)

       call ampl(6,8,Nj2,Nj3,POL2,q2,lab1,l3c(2),ia,ib,BPOL3,q3,lab2,ij,mur2)


       !----------
       lab1=f3c(3)
       lab2=f3c(1)

       ia=i1+i2+1
       ib=Npoint


       call ampl(6,8,Nj3,Nj1,POL3,q3,lab1,l3c(3),ia,ib,BPOL1,q1,lab2,ij,mur3)



       do j1=1,Nj1
          do j2=1,Nj2
             do j3=1,Nj3
                res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j1)
             enddo
          enddo
       enddo


       !           to multiply by proper power of I:

       do j=1,3
          if (quark_flavour(f3c(j))) then 
             res = ci*res
          else
             res = (-ci)*res
          endif
       enddo


       !           intermediate 6-dim result

       res6 = res


       !           8-dim calculation

       lab=f3c(1)

       call givepol(lab,q1,8,16,tag_pol,Nj1,BPOL1,POL1)

       lab=f3c(2)
       call givepol(lab,q2,8,16,tag_pol,Nj2,BPOL2,POL2)

       lab=f3c(3)
       call givepol(lab,q3,8,16,tag_pol,Nj3,BPOL3,POL3)


       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero

       i1=l3c(2)-l3c(1) 
       i2=l3c(3)-l3c(2)

       lab1=f3c(1)
       lab2=f3c(2)


       call ampl(8,16,Nj1,Nj2,POL1,q1,lab1,l3c(1),1,i1,BPOL2,q2,lab2,ij,mur1)


       !------------new loop

       ia =i1+1
       ib = i1+i2


       lab1=f3c(2)
       lab2=f3c(3)


       call ampl(8,16,Nj2,Nj3,POL2,q2,lab1,l3c(2),ia,ib,BPOL3,q3,lab2,ij,mur2)




       !----------

       ia=i1+i2+1
       ib=Npoint

       lab1=f3c(3)
       lab2=f3c(1)


       call ampl(8,16,Nj3,Nj1,POL3,q3,lab1,l3c(3),ia,ib,BPOL1,q1,lab2,ij,mur3)


       do j1=1,Nj1
          do j2=1,Nj2
             do j3=1,Nj3
                res=res+ mur1(j1,j2)*mur2(j2,j3)*mur3(j3,j1)
             enddo
          enddo
       enddo


       !           to multiply by proper power of I:

       do j=1,3
          if (quark_flavour(f3c(j))) then 
             res = ci*res
          else
             res = (-ci)*res
          endif
       enddo


       res8 = res



       !           final result 



       if (tag_f.eq.1) then 
          res=two*res6-res8
       else
          if (tag_f.eq.2) then    ! this is for internal fermion loop
             res=res6 - res8/4._qp
          else
             res=res6 -res8 
          endif
       endif


       !ccccc!ccccc!ccccccc endif for Ds = 5 
    endif


    !     subtracting the 5-cut 

    call match(3,5,l3c,f3c,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,n35,lmatch35)

    do i=1,n35

       im = lmatch35(i)

       j=1
61     if (l3c(1).ne.Lc5(im,j)) then 
          j=j+1
          go to 61
       endif

       pos=j


       call mismatch(3,5,l3c,im,N5,Lc5,N4,Lc4,N3,Lc3,N2,Lc2,lpos1) 

       pos1=lpos1(1)
       pos2=lpos1(2)

       do j=1,1
          e(j)=coeff5(im,j)
       enddo

       do j=1,4
          krefa(j)=propv5(im,4*(pos-1)+j)
       enddo


       do j=1,4
          lvt(j)=q1(j)-krefa(j)
       enddo
       lvt(5)=q1(5)

       do j=1,4 
          vprop(j)=lvt(j)+propv5(im,j+4*(pos1-1))
       enddo
       vprop(5)=lvt(5)

       prop = sc(vprop,vprop)

       prop1 = prop-mass5(im,pos1)**2


       do j=1,4 
          vprop(j)=lvt(j)+propv5(im,j+4*(pos2-1))
       enddo
       vprop(5)=lvt(5)

       prop =  sc(vprop,vprop)

       prop2 = prop-mass5(im,pos2)**2

       res=res - e(1)*lvt(5)**2/prop1/prop2 


    enddo


    !     subtracting the 4-cut

    call match(3,4,l3c,f3c,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,n34,lmatch34)


    do i=1,n34

       im=lmatch34(i)

       j=1
50     if (l3c(1).ne.Lc4(im,j)) then 
          j=j+1
          go to 50
       endif

       pos=j


       call mismatch(3,4,l3c,im,N5,Lc5,N4,Lc4,N3,Lc3,N2,Lc2,lpos0) 

       pos1 = lpos0(1)


       do j=1,5
          d(j)=coeff4(im,j)
       enddo

       do j=1,4
          krefa(j)=propv4(im,4*(pos-1)+j)
          v45(j)=refvect4(im,j)   
       enddo

       do j=1,4
          lvt(j)=q1(j)-krefa(j)
       enddo
       lvt(5)=q1(5)

       do j=1,4 
          vprop(j)=lvt(j)+propv4(im,j+4*(pos1-1))
       enddo
       vprop(5)=lvt(5)

       prop=sc(vprop,vprop)

       prop = prop-mass4(im,pos1)**2


       do j=1,5
          vne(j)=czero
       enddo
       vne(5)=ci

       r1=sc(v45,lvt)
       r2=sc(vne,lvt)

       r22=r2**2

       res=res-(d(1)+d(2)*r1+r22*(d(3)+d(4)*r1+d(5)*r22))/prop 

    enddo

  end subroutine resid3

  !------------------------------------------------------------------------
  subroutine resid2(Ds,lv,k1,l2c,f2c,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res)
    integer, intent(in)        :: Ds, l2c(2), ij
    character, intent(in)      :: f2c(2)*3
    complex(qp), intent(in)  :: lv(Ds)
    complex(qp), intent(in)  :: k1(4)
    integer, intent(in)        :: N5,N4,N3,N2
    integer, intent(in)        :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    character, intent(in)      :: F5(:,:)*3,F4(:,:)*3,F3(:,:)*3,F2(:,:)*3
    complex(qp), intent(out) :: res
    ! ------------------------------------------------------------       
    complex(qp) :: lvt(5),e(1)
    complex(qp) :: v45(4),vprop(5),prop,r1,r2,krefa(4)
    complex(qp) :: d(5),c(10),re
    complex(qp) :: prop1,prop2
    complex(qp) :: q1(5),q2(5)
    complex(qp) :: r3,r4
    complex(qp) :: vne(5),trikoeff
    complex(qp) :: v3(4),v4(4)
    complex(qp) :: BPOL1(8,16),POL1(8,16)
    complex(qp) :: BPOL2(8,16),POL2(8,16)
    complex(qp) :: mur1(10,10),mur2(10,10)
    complex(qp) :: res6,res8
    complex(qp) :: q12(5),prop3,r22,r32,r42
    character     :: lab*3, lab1*3, lab2*3
    integer       :: i,tag_pol,pos,pos1,pos2,n25
    integer       :: j,j1,j2,lpos(2),lpos2(3),pos3
    integer       :: i1,ia,ib,n24,n23, lpos0(1)
    integer       :: lmatch25(N5),lmatch24(N4)
    integer       :: lmatch23(N3)
    integer       :: tag_f
    integer       :: Nj1,Nj2
    integer       :: im



    do i=1,4
       q1(i)=lv(i)
       q2(i)=lv(i)+k1(i)
    enddo

    if (Ds.eq.4) then 
       q1(5)=czero
       q2(5)=czero
    else
       q1(5)=lv(5)
       q2(5)=lv(5)
    endif

    do i=1,5
       q12(i)=q1(i)+q2(i)
    enddo

    r1 =  sc(q12,q12)
    if(abs(r1-cone*mt**2).lt.propcut) then 
       tag_pol=1
    else
       tag_pol=0
    endif

    ! tag_f means that we have a fermion loop 
    if (any(f2c .eq. 'glu')) then 
       tag_f = 0
    else
       tag_f = 1
    endif



    if (ferm_loops .or. ferm_loops_Z .or. ferm_loops_Z_sbs) tag_f = 2
    if (extra_ferm_pair_nf) tag_f = 2
    if (gluons_ferm_loops) tag_f = 2


    !         residues for various values of Ds


    if (Ds.eq.4) then 

       lab=f2c(1)

       call givepol(lab,q1,4,4,tag_pol,Nj1,BPOL1,POL1)

       lab=f2c(2)
       call givepol(lab,q2,4,4,tag_pol,Nj2,BPOL2,POL2)

       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero

       i1=l2c(2)-l2c(1)

       !----------new loop

       lab1=f2c(1)
       lab2=f2c(2)



       call ampl(4,4,Nj1,Nj2,POL1,q1,lab1,l2c(1),1,i1,BPOL2,q2,lab2,ij,mur1)




       !-----------new loop
       ia =i1+1
       ib=Npoint

       lab1=f2c(2)
       lab2=f2c(1)

       call ampl(4,4,Nj2,Nj1,POL2,q2,lab1,l2c(2),ia,ib,BPOL1,q1,lab2,ij,mur2)


       do j1=1,Nj1
          do j2=1,Nj2
             res=res+ mur1(j1,j2)*mur2(j2,j1)
          enddo
       enddo
       !--------- to multiply by proper power of I:

       do j=1,2
          if (quark_flavour(f2c(j))) then 
             res = ci*res
          else
             res = (-ci)*res
          endif
       enddo

       !ccccc!ccccc!ccccccc endif for Ds = 4 
    endif




    if (Ds.eq.5) then 

       lab=f2c(1)

       call givepol(lab,q1,6,8,tag_pol,Nj1,BPOL1,POL1)

       lab=f2c(2)
       call givepol(lab,q2,6,8,tag_pol,Nj2,BPOL2,POL2)


       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero

       i1=l2c(2)-l2c(1)

       !------------new loop

       lab1=f2c(1)
       lab2=f2c(2)


       call ampl(6,8,Nj1,Nj2,POL1,q1,lab1,l2c(1),1,i1,BPOL2,q2,lab2,ij,mur1)



       !------------new loop
       ia =i1+1
       ib = Npoint

       lab1=f2c(2)
       lab2=f2c(1)


       call ampl(6,8,Nj2,Nj1,POL2,q2,lab1,l2c(2),ia,ib,BPOL1,q1,lab2,ij,mur2)


       do j1=1,Nj1
          do j2=1,Nj2
             res=res+ mur1(j1,j2)*mur2(j2,j1)
          enddo
       enddo



       !           to multiply by proper power of I:

       do j=1,2
          if (quark_flavour(f2c(j))) then 
             res = ci*res
          else
             res = (-ci)*res
          endif
       enddo


       !----------- intermediate 6-dim result

       res6 = res

       !----------8-dim calculation

       lab=f2c(1)

       call givepol(lab,q1,8,16,tag_pol,Nj1,BPOL1,POL1)

       lab=f2c(2)
       call givepol(lab,q2,8,16,tag_pol,Nj2,BPOL2,POL2)



       !         calculate product of amplitudes
       !         first  organize 4 lists that will be arguments 
       !         for amplitude procedure 
       !         start with the first label on the cut list ll
       !         Lab_ex and Lab_in

       res=czero

       !----------new loop

       lab1=f2c(1)
       lab2=f2c(2)


       call ampl(8,16,Nj1,Nj2,POL1,q1,lab1,l2c(1),1,i1,BPOL2,q2,lab2,ij,mur1)


       !------------new loop
       ia =i1+1
       ib = Npoint


       lab1=f2c(2)
       lab2=f2c(1)

       call ampl(8,16,Nj2,Nj1,POL2,q2,lab1,l2c(2),ia,ib,BPOL1,q1,lab2,ij,mur2)




       do j1=1,Nj1
          do j2=1,Nj2
             res=res+ mur1(j1,j2)*mur2(j2,j1)
          enddo
       enddo


       !---------- multiply by proper power of I:
       do j=1,2
          if (quark_flavour(f2c(j))) then 
             res = ci*res
          else
             res = (-ci)*res
          endif
       enddo


       res8 = res


       !           final result 

       if (tag_f.eq.1) then 
          res=2*res6 -res8 
       else
          if (tag_f.eq.2) then 
             res = res6 - res8/4._qp
          else
             res = res6 - res8
          endif
       endif



       !ccccc!ccccc!ccccccc endif for Ds = 5 
    endif




    !     subtracting the 5-cut 

    call match(2,5,l2c,f2c,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,n25,lmatch25)


    do i=1,n25

       im = lmatch25(i)



       j=1
62     if (l2c(1).ne.Lc5(im,j)) then 
          j=j+1
          go to 62
       endif

       pos=j

       call mismatch(2,5,l2c,im,N5,Lc5,N4,Lc4,N3,Lc3,N2,Lc2,lpos2) 

       pos1=lpos2(1)
       pos2=lpos2(2)
       pos3=lpos2(3)


       do j=1,1
          e(j)=coeff5(im,j)
       enddo

       do j=1,4
          krefa(j)=propv5(im,4*(pos-1)+j)
       enddo


       do j=1,4
          lvt(j)=q1(j)-krefa(j)
       enddo
       lvt(5)=q1(5)

       do j=1,4 
          vprop(j)=lvt(j)+propv5(im,j+4*(pos1-1))
       enddo
       vprop(5)=lvt(5)

       prop=sc(vprop,vprop)

       prop1 = prop-mass5(im,pos1)**2



       do j=1,4 
          vprop(j)=lvt(j)+propv5(im,j+4*(pos2-1))
       enddo
       vprop(5)=lvt(5)

       prop=sc(vprop,vprop)

       prop2 = prop-mass5(im,pos2)**2

       do j=1,4 
          vprop(j)=lvt(j)+propv5(im,j+4*(pos3-1))
       enddo
       vprop(5)=lvt(5)

       prop=sc(vprop,vprop)

       prop3 = prop-mass5(im,pos3)**2

       res=res - e(1)*lvt(5)**2/prop1/prop2/prop3 

    enddo



    !        subtracting the 4-cut

    call match(2,4,l2c,f2c,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,n24,lmatch24)



    do i=1,n24

       im=lmatch24(i)

       j=1
51     if (l2c(1).ne.Lc4(im,j)) then 
          j=j+1
          go to 51
       endif

       pos=j


       call mismatch(2,4,l2c,im,N5,Lc5,N4,Lc4,N3,Lc3,N2,Lc2,lpos) 

       pos1=lpos(1)
       pos2=lpos(2)

       do j=1,5
          d(j)=coeff4(im,j)
       enddo

       do j=1,4
          krefa(j)=propv4(im,4*(pos-1)+j)
          v45(j)=refvect4(im,j)   
       enddo

       do j=1,4
          lvt(j)=q1(j)-krefa(j)
       enddo
       lvt(5)=q1(5)

       do j=1,4 
          vprop(j)=lvt(j)+propv4(im,j+4*(pos1-1))
       enddo
       vprop(5)=lvt(5)

       prop=sc(vprop,vprop)

       prop1 = prop-mass4(im,pos1)**2

       do j=1,4 
          vprop(j)=lvt(j)+propv4(im,j+4*(pos2-1))
       enddo
       vprop(5)=lvt(5)

       prop=sc(vprop,vprop)

       prop2 = prop-mass4(im,pos2)**2

       do j=1,5
          vne(j)=czero
       enddo
       vne(5)=ci

       r1 = sc(v45,lvt)
       r2 = sc(vne,lvt)

       r22=r2**2

       res=res-(d(1)+d(2)*r1+r22*(d(3) +d(4)*r1+d(5)*r22))/prop1/prop2 

    enddo

    !        subtracting the triple -cut contribution

    call match(2,3,l2c,f2c,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,n23,lmatch23)
    
    do i=1,n23

       im = lmatch23(i)

       j=1
52     if (l2c(1).ne.Lc3(im,j)) then 
          j=j+1
          go to 52
       endif

       pos=j

       call mismatch(2,3,l2c,im,N5,Lc5,N4,Lc4,N3,Lc3,N2,Lc2,lpos0) 
       pos1=lpos0(1)

       do j=1,10
          c(j)=coeff3(im,j)
       enddo


       do j=1,4
          krefa(j)=propv3(im,4*(pos-1)+j)
          v3(j)=refvect3(im,j)   
          v4(j)=refvect3(im,4+j)   
       enddo


       do j=1,4
          lvt(j)=q1(j)-krefa(j)
       enddo
       lvt(5)=q1(5)

       do j=1,4 
          vprop(j)=lvt(j)+propv3(im,j+4*(pos1-1))
       enddo
       vprop(5)=lvt(5)

       prop=sc(vprop,vprop)

       prop = prop-mass3(im,pos1)**2


       do j=1,5
          vne(j)=czero
       enddo
       vne(5)=ci

       r3 = sc(v3,lvt)
       r4 = sc(v4,lvt)
       re = sc(vne,lvt)

       r42=r4**2
       r32=r3**2

       trikoeff=c(1)+c(2)*r3+c(3)*r4+c(4)*r3*r4 +c(5)*(r32-r42)+&
            &c(6)*r32*r4 +c(7)*r3*r42+re**2*(c(8) +c(9)*r3+c(10)*r4)


       res=res-trikoeff/prop

    enddo

  end subroutine resid2



  !=========================================================================
  !      procedures for computing cuts 
  !=========================================================================
  subroutine pentcut(Lc5,F5,Yc5)
    integer, intent(in)   :: Lc5(:,:)
    character, intent(in) :: F5(:,:)*3
    integer, intent(in)   :: Yc5(:,:)
    ! -----------------------------------------------------------------
    complex(qp) :: k1(4),k2(4),k3(4),k4(4)
    complex(qp) :: v(4,4), v1(4),v2(4), v3(4), v4(4)
    complex(qp) :: lv5(5)
    complex(qp) :: V1234(4)
    complex(qp) :: x1,x2,x3,r1,r2,r3,r4,ltr,x4
    complex(qp) :: summ,sump
    complex(qp) :: e(1)
    real(qp)    :: m(5)
    integer       :: i,j,ij      
    integer       :: l5cut(5)
    integer       :: N5
    character     :: f5cut(5)*3
    ! ------------------------------------
    complex(qp) :: pred(4,4), vi(4,4)  


    N5 = size(Lc5,dim=1)


    do i=1,N5

       do j=1,5
          l5cut(j) = Lc5(i,j)
          f5cut(j) = F5(i,j)
       enddo


       ij = Yc5(i,1) !--- which momentum ordering to choose

       do j=1,4
          k1(j)=momline(ij,l5cut(2),j)-momline(ij,l5cut(1),j)
          k2(j)=momline(ij,l5cut(3),j)-momline(ij,l5cut(1),j)
          k3(j)=momline(ij,l5cut(4),j)-momline(ij,l5cut(1),j)
          k4(j)=momline(ij,l5cut(5),j)-momline(ij,l5cut(1),j)
       enddo

       do j=1,5

          if (quark_flavour(f5cut(j))) m(j)=zero
          if (f5cut(j).eq.'glu') then 
             m(j)=zero
          elseif (f5cut(j).eq.'top') then 
             m(j)=mt
          elseif (f5cut(j).eq.'bot') then 
             m(j)=mb
          endif
       enddo


    if (use_mynivecs) then 
          pred(:,1) = k1
          pred(:,2) = k2
          pred(:,3) = k3
          pred(:,4) = k4
          call compute_vi(pred(:,:4),vi(:,:4))

          v1 = vi(:,1)
          v2 = vi(:,2)
          v3 = vi(:,3)
          v4 = vi(:,4)
       else

          call give4to4vect(k1,k2,k3,k4,v)
          
          v1=v(1,:)
          v2=v(2,:)
          v3=v(3,:)
          v4=v(4,:)
       endif

       r1= sc(k1,k1)
       r2= sc(k2,k2)
       r3= sc(k3,k3)
       r4= sc(k4,k4)

       x1=-half*(r1+m(1)**2-m(2)**2)
       x2=-half*(r2+m(1)**2-m(3)**2)
       x3=-half*(r3+m(1)**2-m(4)**2)
       x4=-half*(r4+m(1)**2-m(5)**2)


       V1234=x1*v1+x2*v2+x3*v3+x4*v4

       r1 = sc(V1234,V1234)

       ltr=sqrt(m(1)**2-r1)

       do j=1,4
          lv5(j)=V1234(j)
       enddo

       lv5(5)=ltr*ci

       

       call resid5(5,lv5,k1,k2,k3,k4,l5cut,f5cut,ij,sump)

       lv5(5)=-ltr*ci
       call resid5(5,lv5,k1,k2,k3,k4,l5cut,f5cut,ij,summ)

!------changed, KM, divided by ltr**2

       e(1)=half*(sump + summ)/lv5(5)**2

       do j=1,1
          coeff5(i,j)=e(j)
       enddo

       do j=1,5
          mass5(i,j) = m(j)
       enddo

       do j=1,4
          propv5(i,j)  = czero
          propv5(i,j+4)= k1(j)
          propv5(i,j+8)= k2(j)
          propv5(i,j+12)=k3(j)
          propv5(i,j+16)=k4(j)
       enddo

    enddo

  end subroutine pentcut


  subroutine quadcut(Lc5,F5,Yc5,Lc4,F4,Yc4,Lc3,F3,Yc3,Lc2,F2,Yc2)  ! quadrupole cut
    integer, intent(in)   :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    character, intent(in) :: F5(:,:)*3,F4(:,:)*3,F3(:,:)*3, F2(:,:)*3
    integer, intent(in)   :: Yc5(:,:),Yc4(:,:),Yc3(:,:),Yc2(:,:)
    ! ------------------------------------------------------------      
    complex(qp) ::  k1(4),k2(4),k3(4)
    complex(qp) ::  v(4,4), v1(4),v2(4), v3(4), v4(4)
    complex(qp) ::  lv4(4),lv5(5)
    complex(qp) ::  V123(4)
    complex(qp) ::  x1,x2,x3,r1,r2,r3,ltr,xtr,xe,x4
    complex(qp) ::  koeffp, koeffm, koeff1, koeff2
    complex(qp) ::  summ,sump
    complex(qp) ::  d(5)
    real(qp)    :: m(4)
    integer       ::  i,j,ij      
    integer       ::  l4cut(4)
    character     :: f4cut(4)*3
    integer       :: N5,N4,N3,N2
    ! ------------------------------------
    complex(qp) :: ni(4,1), pred(4,3), vi(4,3)  


    N5 = size(Lc5,dim=1)
    N4 = size(Lc4,dim=1)
    N3 = size(Lc3,dim=1)
    N2 = size(Lc2,dim=1)

    do i=1,N4

       do j=1,4
          l4cut(j)=Lc4(i,j)
          f4cut(j)=F4(i,j)
       enddo

       ij = Yc4(i,1)

       do j=1,4
          k1(j)=momline(ij,l4cut(2),j)-momline(ij,l4cut(1),j)
          k2(j)=momline(ij,l4cut(3),j)-momline(ij,l4cut(1),j)
          k3(j)=momline(ij,l4cut(4),j)-momline(ij,l4cut(1),j)
       enddo




       do j=1,4
          if (quark_flavour(f4cut(j))) m(j)=zero
          if (f4cut(j).eq.'glu') m(j)=zero
          if (f4cut(j).eq.'bot') m(j)=mb
          if (f4cut(j).eq.'top') m(j)=mt
       enddo


       if (use_mynivecs) then 
          pred(:,1) = k1
          pred(:,2) = k2
          pred(:,3) = k3

          call compute_ni(pred(:,:3),bvec(:,:1),ni(:,:1)) 
          call compute_vi(pred(:,:3),vi(:,:3))

          v1 = vi(:,1)
          v2 = vi(:,2)
          v3 = vi(:,3)
          v4 = ni(:,1)

       else

          call give3to4vect(k1,k2,k3,v)
          
          v1=v(1,:)
          v2=v(2,:)
          v3=v(3,:)
          v4=v(4,:)
       endif

       r1=sc(k1,k1)
       r2=sc(k2,k2)
       r3=sc(k3,k3)


       x1=-half*(r1+m(1)**2-m(2)**2)
       x2=-half*(r2+m(1)**2-m(3)**2)
       x3=-half*(r3+m(1)**2-m(4)**2)


       V123=x1*v1+x2*v2+x3*v3 


       r1 = sc(V123,V123)
       ltr=sqrt(m(1)**2-r1)

       !----------------------------------------     4-dim coefficients

       do j=1,4
          lv4(j)=ltr*v4(j) + V123(j)
       enddo



       call resid4(4,lv4,k1,k2,k3,l4cut,f4cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)


       do j=1,4
          lv4(j)= -ltr*v4(j) + V123(j)
       enddo

       call resid4(4,lv4,k1,k2,k3,l4cut,f4cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,summ)


       d(1)=half*(sump +summ )
       d(2)=half/ltr*(sump  - summ)


       !-------------------------------------- d-dim coefficients

       xtr=ltr/sqrt2
       xe=sqrt(ltr**2-xtr**2)

       do j=1,4
          lv5(j)=xtr*v4(j)+V123(j)
       enddo

       lv5(5) = xe*ci


       call resid4(5,lv5,k1,k2,k3,l4cut,f4cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)


       r1=sc(lv5,v4)

       do j=1,4
          lv5(j)= -xtr*v4(j)+V123(j)
       enddo

       lv5(5) = xe*ci

       call resid4(5,lv5,k1,k2,k3,l4cut,f4cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,summ)

       r2=sc(lv5,v4)

       koeffp=sump-d(1)-d(2)*r1
       koeffm=summ-d(1)-d(2)*r2

       d(4)=half/r1/xe**2*(koeffp - koeffm)


       x4=ltr/sqrt3
       xe=sqrt(ltr**2-x4**2)

       do j=1,4
          lv5(j)=x4*v4(j)+V123(j)
       enddo

       lv5(5)=xe*ci

       r1 = sc(lv5,v4)

       call resid4(5,lv5,k1,k2,k3,l4cut,f4cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)


       koeff1=sump-d(1)-d(2)*r1-d(4)*r1*xe**2


       x4=ltr/sqrt2
       xe=sqrt(ltr**2-x4**2)

       do j=1,4
          lv5(j)=x4*v4(j)+V123(j)
       enddo

       lv5(5)=xe*ci

       r1 =sc(lv5,v4)

       call resid4(5,lv5,k1,k2,k3,l4cut,f4cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)

       koeff2=sump-d(1)-d(2)*r1- d(4)*r1*xe**2

       d(3)=-half*(9.0_qp*koeff1-16.0_qp*koeff2)/ltr**2
       d(5)= 3.0_qp*(3.0_qp*koeff1-4.0_qp*koeff2)/ltr**4


       do j=1,5
          coeff4(i,j)=d(j)
       enddo

       do j=1,4
          refvect4(i,j)=v4(j)
       enddo

       do j=1,4
          mass4(i,j) = m(j)
       enddo

       do j=1,4
          propv4(i,j)=czero
          propv4(i,j+4)=k1(j)
          propv4(i,j+8)=k2(j)
          propv4(i,j+12)=k3(j)
       enddo

    enddo

  end subroutine quadcut



  subroutine tripcut(Lc5,F5,Yc5,Lc4,F4,Yc4,Lc3,F3,Yc3,Lc2,F2,Yc2)
    integer, intent(in)   :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    character, intent(in) :: F5(:,:)*3,F4(:,:)*3,F3(:,:)*3, F2(:,:)*3
    integer, intent(in)   :: Yc5(:,:),Yc4(:,:),Yc3(:,:),Yc2(:,:)
    ! ----------------------------------------------------------       
    complex(qp) :: k1(4),k2(4)
    complex(qp) :: v(4,4), v1(4),v2(4), v3(4), v4(4)
    complex(qp) :: lv4(4),lva(5),lvb(5),lv4a(4),lv4b(4)
    complex(qp) :: V12(4)
    complex(qp) :: x1,x2,r1,r2,ltr,xtr,xe,x4
    complex(qp) :: summ,sump
    complex(qp) :: c(10)
    complex(qp) :: lhs_eqs(7), lhseq_new(4),lhseq(4)
    complex(qp) :: lhseqa,lhseqb
    complex(qp) :: phi1,r3a,r4a,r3b,r4b
    complex(qp) :: bcoeff(-3:3)
    complex(qp) :: dcp,dsp
    real(qp)    :: m(3)
    real(qp)    :: buff, phi
    integer       :: i,j,kk,ss,i1,ij      
    integer       :: l3cut(3)
    character     :: f3cut(3)*3
    integer       :: N5,N4,N3,N2
    ! ------------------------------------
    complex(qp) :: ni(4,2), pred(4,2),vi(4,2)
    ! -- variables needed for stability test 
    complex(qp) :: rnx,lv_rn(5),res_v1,res_v2,res_err
    complex(qp) :: ldv3,ldv4,lvsq
    real(qp)    :: rn1,rn2,rn3


    N5 = size(Lc5,dim=1)
    N4 = size(Lc4,dim=1)
    N3 = size(Lc3,dim=1)
    N2 = size(Lc2,dim=1)

    buff=1d-1

    do i=1,N3

       do j=1,3
          l3cut(j)=Lc3(i,j)
          f3cut(j)=F3(i,j)
       enddo

       ij = Yc3(i,1)
       do j=1,4
          k1(j)=momline(ij,l3cut(2),j)-momline(ij,l3cut(1),j)
          k2(j)=momline(ij,l3cut(3),j)-momline(ij,l3cut(1),j)
       enddo

       do j=1,3
          if (quark_flavour(f3cut(j))) m(j)=zero
          if (f3cut(j).eq.'glu') m(j)=zero
          if (f3cut(j).eq.'top') m(j)=mt
          if (f3cut(j).eq.'bot') m(j)=mb
       enddo
       
       if (use_mynivecs) then 
          pred(:,1) = k1
          pred(:,2) = k2
          call compute_ni(pred(:,:2),bvec(:,:2),ni(:,:2)) 
          call compute_vi(pred(:,:2),vi)

          v1 = vi(:,1)
          v2 = vi(:,2)
          v3 = ni(:,1)
          v4 = ni(:,2)
          
       else
          call give2to4vect(k1,k2,v)
          
          v1=v(1,:)
          v2=v(2,:)
          v3=v(3,:)
          v4=v(4,:)
       endif
       r1 = sc(k1,k1)
       r2 = sc(k2,k2)

       x1=-half*(r1+m(1)**2-m(2)**2)
       x2=-half*(r2+m(1)**2-m(3)**2)


       V12=x1*v1+x2*v2

       r1 = sc(V12,V12)

       ltr=sqrt(m(1)**2-r1)


       !------ two cases -- large and small ltr separated by buff
       if (abs(ltr).ge.buff) then 

          do kk=1,7
             phi=2.0_qp*pi*(kk-1)/7.0_qp
             dcp=ltr*cos(phi)
             dsp=ltr*sin(phi)

             lv4=V12+dcp*v3+dsp*v4

             call resid3(4,lv4,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)
             lhs_eqs(kk)=sump
          enddo


          do j=1,7
             ss=j-4
             bcoeff(ss)=czero
             do kk=1,7
                phi1=-2.0_qp*pi*ci/7.0_qp*(kk-1)*ss
                bcoeff(ss)= bcoeff(ss)+1.0_qp/7.0_qp*lhs_eqs(kk)*exp(phi1)
             enddo
          enddo


          c(1)=bcoeff(0)
          c(4)=-2.0_qp*ci*(bcoeff(-2)-bcoeff(2))/ltr**2
          c(5)=(bcoeff(-2)+bcoeff(2))/ltr**2
          c(6)=-4.0_qp*ci*(bcoeff(-3)-bcoeff(3))/ltr**3
          c(7)=-4.0_qp*(bcoeff(-3)+bcoeff(3))/ltr**3
          c(2)= 0.25_qp*(-c(7)*ltr**3+4*bcoeff(-1)+4.0_qp*bcoeff(1))/ltr
          c(3)= 0.25_qp*ci*(c(6)*ltr**3*ci   &
               &-4.0_qp*bcoeff(-1)+4.0_qp*bcoeff(1))/ltr

          ! -- stability check in 4 dimensions 
          !rnx = 17._qp/19._qp
          !lv_rn(:4) = V12+ltr*(rnx*v3+sqrt(one-rnx**2)*v4)
          !call resid3(4,lv_rn(:4),k1,k2,l3cut,f3cut,ij,&
          !     &N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res_v1)
          !
          !ldv3 = dot(lv_rn(:4),v3)
          !ldv4 = dot(lv_rn(:4),v4)
          !res_v2 = c(1)+c(2)*ldv3+c(3)*ldv4+c(5)*(ldv3**2-ldv4**2)+&
          !     &ldv3*ldv4*(c(4)+c(6)*ldv3+c(7)*ldv4)
          !
          !if (abs(res_v1-res_v2) > sq2tol) then 
          !   write(*,*) 'not ZERO?',(res_v1-res_v2)
          !endif
          !------now determine  ep-dimensional structures

          xtr=ltr/sqrt2
          xe=sqrt(ltr**2 - xtr**2)

          do j=1,4
             lva(j)=V12(j)+xtr*v3(j)
          enddo

          lva(5)=xe*ci 

          call resid3(5,lva,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

          xtr=-xtr
          xe=sqrt(ltr**2 - xtr**2)

          do j=1,4
             lvb(j)=V12(j)+xtr*v3(j)
          enddo

          lvb(5)=xe*ci 

          call resid3(5,lvb,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqb)

          r3a= sc(lva,v3)
          r4a= sc(lva,v4)

          lhseqa=lhseqa - (c(1) + c(2)*r3a+c(3)*r4a+c(4)*r3a*r4a   &
               &+c(5)*(r3a**2-r4a**2)+c(6)*r3a**2*r4a+c(7)*r3a*r4a**2)

          r3b = sc(lvb,v3)
          r4b = sc(lvb,v4)


          lhseqb=lhseqb - (c(1) + c(2)*r3b+c(3)*r4b+c(4)*r3b*r4b &
               &+c(5)*(r3b**2-r4b**2)+c(6)*r3b**2*r4b+c(7)*r3b*r4b**2)

          c(9)=(lhseqb - lhseqa)/2.0_qp/r3b/xe**2
          c(8)=(lhseqb - c(9)*r3b*xe**2)/xe**2

          xtr=ltr/sqrt2
          xe=sqrt(ltr**2-xtr**2)

          do j=1,4
             lva(j)=V12(j) + xtr*v4(j)
          enddo

          lva(5)=xe*ci

          call resid3(5,lva,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)


          r3a =  sc(lva,v3)
          r4a =  sc(lva,v4)


          lhseqa=lhseqa - (c(1) + c(2)*r3a+c(3)*r4a+c(4)*r3a*r4a &
               &+c(5)*(r3a**2-r4a**2)+c(6)*r3a**2*r4a+c(7)*r3a*r4a**2)

          c(10)=(lhseqa - c(8)*xe**2 -c(9)*r3a*xe**2)/r4a/xe**2


          ! -- stability check in D dimension 
          rn1=17._qp/19._qp
          rn2=-5._qp/7._qp
          rn3=23._qp/29._qp
          xtr = ltr/sqrt(rn1**2+rn2**2+rn3**2)
          do j=1,4
             lv_rn(j)=V12(j)+xtr*(rn1*v3(j)+rn2*v4(j))
          enddo
          lv_rn(5)=xtr*rn3*ci

          ldv3 = dot(lv_rn(:4),v3)
          ldv4 = dot(lv_rn(:4),v4)
          lvsq = (xtr*rn3)**2
          
          call resid3(5,lv_rn,k1,k2,l3cut,f3cut,ij,&
               &N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res_v1)
          res_v2 = c(1)+c(2)*ldv3+c(3)*ldv4+c(5)*(ldv3**2-ldv4**2)+&
               &ldv3*ldv4*(c(4)+c(6)*ldv3+c(7)*ldv4)+&
               &lvsq*(c(8)+c(9)*ldv3+c(10)*ldv4)

       elseif (abs(ltr).lt.buff) then 

          do kk=1,4
             phi=2.0_qp*pi/four*(kk-1)
             xtr=exp(ci*phi)+ci*10._qp*tol

             x4=sqrt(ltr**2-xtr**2)

             lv4a=V12 + xtr*v3+x4*v4
             lv4b=V12 + xtr*v3-x4*v4

             call resid3(4,lv4a,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)

             call resid3(4,lv4b,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,summ)

             lhseq(kk)=half*(sump+summ) 

             lhseq_new(kk)=half/x4*(sump-summ)

          enddo


          do i1=1,4
             ss=i1-1
             bcoeff(ss) = czero

             do kk=1,4  
                phi1=-2.0_qp*pi*ci/four*(kk-1)*ss
                bcoeff(ss)=bcoeff(ss)+0.25_qp*lhseq(kk)*exp(phi1)
             enddo
          enddo

          c(7)=-bcoeff(3)
          c(5)=half*bcoeff(2)
          c(2)=bcoeff(1) - c(7)*ltr**2
          c(1)=bcoeff(0) + c(5)*ltr**2

          c(4)=(lhseq_new(1)-lhseq_new(3))/2.0_qp
          c(3)=(lhseq_new(1)+lhseq_new(2)-(cone+ci)*c(4))/2.0_qp
          c(6)=(lhseq_new(1)-lhseq_new(2)-(cone-ci)*c(4))/2.0_qp


          !! -- stability check in 4 dimension 
          !rnx = exp(ci*(17._qp/19._qp))
          !x4 = sqrt(ltr**2-rnx**2)
          !lv_rn(:4) = V12+rnx*v3+x4*v4
          !call resid3(4,lv_rn(:4),k1,k2,l3cut,f3cut,ij,&
          !     &N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res_v1)
          !
          !ldv3 = dot(lv_rn(:4),v3)
          !ldv4 = dot(lv_rn(:4),v4)
          !res_v2 = c(1)+c(2)*ldv3+c(3)*ldv4+c(5)*(ldv3**2-ldv4**2)+&
          !     &ldv3*ldv4*(c(4)+c(6)*ldv3+c(7)*ldv4)
          !
          !if (abs(res_v1-res_v2) > sq2tol) then 
          !   write(*,*) 'small buff: not ZERO?',(res_v1-res_v2),&
          !        &(res_v1),(res_v2)
          !endif

          !------now determine  ep-dimensional structures

          xtr=cone-ci*10._qp*tol

          xe=sqrt(ltr**2 - xtr**2)

          do j=1,4
             lva(j)=V12(j) + xtr*v3(j)
          enddo

          lva(5)=xe*ci

          call resid3(5,lva,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

          xtr=-cone

          do j=1,4
             lvb(j)=V12(j) + xtr*v3(j)
          enddo

          lvb(5)=xe*ci

          call resid3(5,lvb,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqb)

          r3a = sc(lva,v3)
          r4a = sc(lva,v4)

          lhseqa=lhseqa - (c(1) + c(2)*r3a+c(3)*r4a+c(4)*r3a*r4a   &
               &+c(5)*(r3a**2-r4a**2)+c(6)*r3a**2*r4a+c(7)*r3a*r4a**2)

          r3b =  sc(lvb,v3)
          r4b =  sc(lvb,v4)

          lhseqb=lhseqb - (c(1) + c(2)*r3b+c(3)*r4b+c(4)*r3b*r4b &
               &+c(5)*(r3b**2-r4b**2)+c(6)*r3b**2*r4b+c(7)*r3b*r4b**2)

          c(9)=(lhseqb - lhseqa)/two/r3b/xe**2
          c(8)=(lhseqb - c(9)*r3b*xe**2)/xe**2


          xtr=cone +ci*10._qp*tol
          
          xe=sqrt(ltr**2-xtr**2)

          do j=1,4
             lva(j)=V12(j) + xtr*v4(j)
          enddo

          lva(5)=xe*ci

          call resid3(5,lva,k1,k2,l3cut,f3cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

          r3a =  sc(lva,v3)
          r4a =  sc(lva,v4)

          lhseqa=lhseqa - (c(1) + c(2)*r3a+c(3)*r4a+c(4)*r3a*r4a &
               &+c(5)*(r3a**2-r4a**2)+c(6)*r3a**2*r4a+c(7)*r3a*r4a**2)

          c(10)=(lhseqa - c(8)*xe**2 -c(9)*r3a*xe**2)/r4a/xe**2

          ! -- stability check in D dimensions 
          rn1=17._qp/19._qp
          rn2=-5._qp/7._qp
          xe = sqrt(ltr**2-rn1**2-rn2**2)
          do j=1,4
             lv_rn(j)=V12(j)+(rn1*v3(j)+rn2*v4(j))
          enddo
          lv_rn(5)=xe*ci
          ldv3 = dot(lv_rn(:4),v3)
          ldv4 = dot(lv_rn(:4),v4)
          lvsq = (xe)**2

          call resid3(5,lv_rn,k1,k2,l3cut,f3cut,ij,&
               &N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res_v1)

          res_v2 = c(1)+c(2)*ldv3+c(3)*ldv4+c(5)*(ldv3**2-ldv4**2)+&
               &ldv3*ldv4*(c(4)+c(6)*ldv3+c(7)*ldv4)+&
               &lvsq*(c(8)+c(9)*ldv3+c(10)*ldv4)

       endif

       !if (abs(res_v2) > tol) then 
       if (abs(res_v2) > sq2tol) then ! changed 4/12/2008 
          res_err = (res_v1-res_v2)/res_v2
       else
          res_err = (res_v1-res_v2)
       endif
       
       do j=1,10
          coeff3(i,j)=c(j)
          dcoeff3(i,j)=coeff3(i,j)*res_err 
       enddo

       do j=1,4
          refvect3(i,j)=v3(j)
          refvect3(i,j+4)=v4(j)
       enddo

       do j=1,3
          mass3(i,j) = m(j)
       enddo

       do j=1,4
          propv3(i,j)=czero
          propv3(i,j+4)=k1(j)
          propv3(i,j+8)=k2(j)
       enddo

    enddo
  end subroutine tripcut


  !---- the double cut
  subroutine doubcut(Lc5,F5,Yc5,Lc4,F4,Yc4,Lc3,F3,Yc3,Lc2,F2,Yc2)
    integer, intent(in)   :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    character, intent(in) :: F5(:,:)*3,F4(:,:)*3,F3(:,:)*3, F2(:,:)*3
    integer, intent(in)   :: Yc5(:,:),Yc4(:,:),Yc3(:,:),Yc2(:,:)
    ! --------------------------------------------------------------       
    complex(qp) :: k1(4),k1d(4)
    complex(qp) :: d1k1_test
    complex(qp) :: v1(4),v2(4), v3(4), v4(4)
    complex(qp) :: V1234(4,4),lv(5),V12(4)
    complex(qp) :: x1,x2,x3,r1,r2,r3,r4,ltr,xtr,xe,x4
    complex(qp) :: ltr2a,ltr2b
    complex(qp) :: sump,lhseqc
    complex(qp) :: b(10),phi1
    complex(qp) :: lhs_eqs(7)
    complex(qp) :: lhseqa,lhseqb,lhseqd,lhs_eqa(3), lhs_eqb(3)
    complex(qp) :: lhseq1,lhseq47a,lhseq47b,lhseq2new
    complex(qp) :: lhseq4new,lhseq3,lhseq36a,lhseq36b
    complex(qp) :: lhseq1new,lhseq2
    complex(qp) :: lhs_eqm(3),lhs_eqp(3),b1eff
    complex(qp) :: bcoeff(-3:3)
    complex(qp) :: dcp,dsp
    real(qp)    :: m(2)
    real(qp)    :: buff, phi
    integer       :: i,j,kk,ss,i1,ij,k      
    integer       :: l2cut(2)
    character     :: f2cut(2)*3
    integer       :: N5,N4,N3,N2
    ! ------------------------------------
    complex(qp) :: ni(4,3), pred(4,1) 
    ! -- variables needed to test stability 
    complex(qp) :: rnx,lv_rn(5),res_v1,res_v2,res_err
    complex(qp) :: ldv2,ldv3,ldv4,lvsq
    real(qp)    :: rn1,rn2,rn3,rn4


    N5 = size(Lc5,dim=1)
    N4 = size(Lc4,dim=1)
    N3 = size(Lc3,dim=1)
    N2 = size(Lc2,dim=1)

    buff=0.1_qp

    do i=1,N2
       do j=1,2
          l2cut(j)=Lc2(i,j)
          f2cut(j)=F2(i,j)
       enddo

       ij = Yc2(i,1)

       do j=1,4
          k1(j)=momline(ij,l2cut(2),j)-momline(ij,l2cut(1),j)
       enddo

       do j=1,2
          if (quark_flavour(f2cut(j))) m(j)=zero
          if (f2cut(j).eq.'glu') m(j)=zero
          if (f2cut(j).eq.'bot') m(j)=mb
          if (f2cut(j).eq.'top') m(j)=mt
       enddo



       !----- a check if a k1 vector is light-like or not: 
       d1k1_test = sc(k1,k1)



       if (abs(d1k1_test).lt.1d-15) then 
          d1k1_test=czero
       endif




       !------ the light like external vector
       if (abs(d1k1_test).eq.zero) then 
          call give1to4vect_light(k1,v1234)
          !----- v1234(1) is just k1, so we have it already


          k1d=v1234(2,:)
          v3=v1234(3,:)
          v4=v1234(4,:)

          !----------------------------------------------------------------------
          ! The loop  momentum is parameterized as
          ! q = x1*k1 + x2*k1d + qtr*(v3*cos(phi)+v4*sin(phi)):
          !
          ! If we are on the cut, we can determine x2 but
          ! x1 and phi are arbitrary; there is a relation between the two, though:
          !       qtr**2 = aD2[2]**2*(1+x1) - x1*aD1[2]**2:
          ! note that, if the two masses are zero or equal, then qtr is 
          ! independent of x1:
          !
          ! The tensor parameterization for the ligh-like vector : 
          !  b1 + b2*sc(k1d,q) + b3*sc(v3,q)+b4*sc(v4,q)+b5*sc(k1d,q)*sc(k1d,q)
          ! + b6*sc(k1d,q)*sc(v3,q)+b7*sc(k1d,q)*sc(v4,q) 
          ! + b8*(sc(v3,q)**2 - sc(v4,q)**2)+b9*sc(v3,q)*sc(v4,q)
          !------------------------------------------------------------------

          r1 = sc(k1,k1d)
          x2=half*(m(2)**2 - m(1)**2)/r1

          x1=chalf
          ltr=sqrt(m(1)**2*(cone+x1)-x1*m(2)**2)

          x3=cone
          x4=sqrt(ltr**2 - x3**2)

          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo



          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

          x3=-x3
          x4=x4
          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo
          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqb)

          x3=-x3
          x4=-x4
          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo
          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqc)

          x3=-x3
          x4=x4
          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo
          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqd)

          lhseq1=0.25_qp*(lhseqa+lhseqb+lhseqc+lhseqd)

          b(9)=0.25_qp/(-x4)*(lhseqa-lhseqc-lhseqb+lhseqd)

          lhseq47a=0.25_qp/(-x4)*(lhseqa-lhseqc+lhseqb-lhseqd)
          lhseq36a=(lhseqa-lhseqb)/2.0_qp - b(9)*(-x4)


          x3=chalf
          x4=sqrt(ltr**2 - x3**2)
          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo
          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseq2new)

          x1=-chalf
          ltr=sqrt(m(1)**2*(cone+x1) - x1*m(2)**2)

          x3=cone
          x4=sqrt(ltr**2 - x3**2)
          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo

          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

          x3=-x3
          x4=x4
          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo

          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqb)

          x3=-x3
          x4=-x4
          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo

          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqc)

          x3=-x3
          x4=x4

          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo

          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqd)

          lhseq3=0.25_qp*(lhseqa+lhseqb+lhseqc+lhseqd)

          lhseq47b=0.25_qp/(-x4)*(lhseqa-lhseqc+lhseqb-lhseqd)
          lhseq36b=(lhseqa-lhseqb)/2.0_qp-b(9)*(-x4)


          b(4)=half*(lhseq47a+lhseq47b)
          b(7)=(lhseq47a-lhseq47b)/r1
          b(3)=half*(lhseq36a+lhseq36b)
          b(6)=(lhseq36a-lhseq36b)/r1

          x1=chalf
          ltr2a=m(1)**2*(cone+x1) - x1*m(2)**2 

          x3=chalf
          x4=sqrt(ltr2a - x3**2)

          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo

          r2= sc(lv,k1d)
          r3 =sc(lv,v3)
          r4 = sc(lv,v4)

          lhseq2new=lhseq2new - b(3)*r3-b(4)*r4  -b(6)*r2*r3 -b(7)*r2*r4 -b(9)*r3*r4


          x1=-chalf
          ltr2b=m(1)**2*(cone+x1) - x1*m(2)**2 

          b(8)=2.0_qp/3.0_qp*(lhseq1-lhseq2new)

          b(2)= ( (lhseq1-b(8)*(ctwo-ltr2a)) &
               &   - (lhseq3-b(8)*(ctwo-ltr2b)) )/r1


          !----------now we repeat this again for a different choice of x1
          x1=chalf*chalf
          ltr=sqrt(m(1)**2*(cone+x1) - x1*m(2)**2)

          x3=cone
          x4=sqrt(ltr**2 - x3**2)

          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo

          call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseq4new)

          r2 =  sc(lv,k1d)
          r3=   sc(lv,v3)
          r4 =  sc(lv,v4)

          lhseq4new=lhseq4new-b(2)*r2-b(3)*r3 -b(4)*r4-b(6)*r2*r3 &
               &-b(7)*r2*r4 -b(8)*(r3**2-r4**2)  -b(9)*r3*r4

          lhseq1new=lhseq1-b(8)*(2.0_qp-ltr2a) -b(2)/two*r1

          b(1)=(4.0_qp*lhseq4new - lhseq1new)/3.0_qp
          b(5)=16.0_qp/r1*(lhseq4new - b(1))


          !-------- the ep-dependent coefficient
          x1=-chalf
          ltr=sqrt(m(1)**2*(cone+x1) - x1*m(2)**2)

          x3=czero

          xe=czero+0.5768943_qp 

          x4=sqrt(ltr**2 - xe**2 - x3**2)

          do j=1,4
             lv(j) = x1*k1(j)+x2*k1d(j)+x3*v3(j)+x4*v4(j)
          enddo

          lv(5)=xe*ci

          call resid2(5,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

          b(10) = one/xe**2*(lhseqa - b(1) - b(2)*r1*x1  - b(3)*x3-b(4)*x4  &
               &-b(5)*r1**2*x1**2 -b(6)*r1*x1*x3 - b(7)*r1*x1*x4 &
               &-b(8)*(x3**2 - x4**2) - b(9)*x3*x4 )



          tagdcut(i,1)=999

          do j=1,10
             coeff2(i,j)=b(j)
          enddo

          do j=1,4
             refvect2(i,j)=k1d(j)
             refvect2(i,j+4)=v3(j)
             refvect2(i,j+8)=v4(j)
          enddo

          do j=1,2
             mass2(i,j) = m(j)
          enddo

          do j=1,4
             propv2(i,j)=czero
             propv2(i,j+4)=k1(j)
          enddo

          !---------------"else is for k1 NOT being a light-like vector
       else 
          if (use_mynivecs) then 
             pred(:,1) = k1
             call compute_ni(pred,bvec,ni) 
             v1= k1/sc(k1,k1)
             v2=ni(:,1)
             v3=ni(:,2)
             v4=ni(:,3)
          else
             
             call give1to4vect(k1,v1234)
             v1=v1234(1,:)
             v2=v1234(2,:)
             v3=v1234(3,:)
             v4=v1234(4,:)
          endif

          r1 =  sc(k1,k1)

          x1=-half*(r1+m(1)**2-m(2)**2)

          V12=x1*v1

          r2 =  sc(V12,V12)

          ltr=sqrt(m(1)**2 - r2)




          !---------we now distinguish two cases -- ``large'' and ``small'' ltr 

          if (abs(ltr) > buff ) then 

             do kk=1,6

                phi=2.0_qp*pi/6.0_qp*(kk-1)
                dcp=ltr*cos(phi)
                dsp=ltr*sin(phi)

                do j=1,4
                   lv(j)=V12(j)+ dcp*v2(j)+ dsp*v3(j)
                enddo

                call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)


                lhs_eqs(kk)=sump
             enddo


             do i1=1,6
                ss=i1-3

                bcoeff(ss) = czero

                do kk=1,6  
                   phi1=-2.0_qp*pi*ci/6.0_qp*(kk-1)*ss
                   bcoeff(ss)=bcoeff(ss)+one/6._qp*lhs_eqs(kk)*exp(phi1)
                enddo
             enddo

             !-------------------------------------b1eff is b1 + b6*ltr**2
             b1eff=bcoeff(0)
             b(2)=1.0_qp/ltr*(bcoeff(1)+bcoeff(-1))
             b(3)=ci/ltr*(bcoeff(1)-bcoeff(-1))
             b(5)=1.0_qp/ltr**2*(bcoeff(2)+bcoeff(-2))
             b(7)=2.0_qp*ci/ltr**2*(bcoeff(2)-bcoeff(-2))


             !-------- take projection on v3 to be 0 and compute three things:

             phi=pi/four
             dcp=ltr*cos(phi)
             dsp=ltr*sin(phi)

             do j=1,4         
                lv(j)= V12(j)+dcp*v2(j)+dsp*v4(j)
             enddo

             call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseq1)

             r2 = sc(lv,v2)
             r3 = sc(lv,v3)

             lhseq1=lhseq1-b(2)*r2-b(3)*r3-b(5)*(r2**2-r3**2)-b(7)*r2*r3


             do j=1,4         
                lv(j)= V12(j)-dcp*v2(j)+dsp*v4(j)
             enddo

             call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseq2)

             r2 = sc(lv,v2)
             r3 = sc(lv,v3)

             lhseq2=lhseq2-b(2)*r2-b(3)*r3-b(5)*(r2**2-r3**2)-b(7)*r2*r3

             do j=1,4         
                lv(j)= V12(j)+dcp*v2(j)-dsp*v4(j)
             enddo

             call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseq3)

             r2 =  sc(lv,v2)
             r3 =  sc(lv,v3)

             lhseq3=lhseq3-b(2)*r2-b(3)*r3-b(5)*(r2**2-r3**2)-b(7)*r2*r3

             b(8)=(lhseq1 - lhseq2)/ltr**2
             b(4)=sqrt2/ltr*half*(lhseq2-lhseq3)
             b(6)=2.0_qp/3.0_qp/ltr**2*(b1eff-(lhseq2+lhseq3)/2.0_qp -b(8)*ltr**2/2.0_qp)

             b(1)=b1eff - b(6)*ltr**2

             do j=1,4         
                lv(j)= V12(j)+dcp*v3(j)+dsp*v4(j)
             enddo

             call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseq1)

             r2 =  sc(lv,v2)
             r3 =  sc(lv,v3)
             r4 =  sc(lv,v4)


             b(9)=2.0_qp/ltr**2*(lhseq1-b(1)-b(2)*r2-b(3)*r3-b(4)*r4-b(5)&
                  &*(r2**2-r3**2)-b(6)*(r2**2 + r3**2 - 2*r4**2)-b(7)*r2*r3 -b(8)*r2*r4)


             !! -- stability check in 4 dimension 
             !rn1=7._qp/19._qp
             !rn2=-5._qp/17._qp
             !rn3=sqrt(one-rn1**2-rn2**2)
             !do j=1,4
             !   lv_rn(j)=V12(j)+ltr*(rn1*v2(j)+rn2*v3(j)+rn3*v4(j))
             !enddo
             !ldv2 = dot(lv_rn(:4),v2)
             !ldv3 = dot(lv_rn(:4),v3)
             !ldv4 = dot(lv_rn(:4),v4)
             !
             !call resid2(4,lv_rn(:4),k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,&
             !     &Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res_v1)
             !
             !!       ----------> GOT HERE 
             !res_v2 = b(1)+b(2)*ldv2+b(3)*ldv3+b(4)*ldv4+b(5)*(ldv2**2-ldv3**2)&
             !     &+b(6)*(ldv2**2+ldv3**2-two*ldv4**2)+b(7)*ldv2*ldv3&
             !     &+b(8)*ldv2*ldv4+b(9)*ldv3*ldv4
             !
             !if (abs(res_v1-res_v2) > sq2tol) then 
             !   write(*,*) 'bub: big buff ZERO2?',(abs(res_v1-res_v2))&
             !        &,(abs(res_v1)),(abs(res_v2))
             !endif
             

             !  this `else' is for  ltr <> 0 condition
          else

             !   now we focus on the case of small |ltr|  case

             r1 =  sc(V12,V12)
             ltr=sqrt(m(1)**2-r1)

             do kk=1,3
                phi=2.0_qp*pi/three*(kk-1)
                xtr=exp(ci*phi)

                x3=sqrt(ltr**2 - xtr**2-ci*10._qp*tol)

                do j=1,4
                   lv(j) = V12(j) + xtr*v2(j) + x3*v3(j) 
                enddo


                call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)

                lhs_eqa(kk)=sump

                do j=1,4
                   lv(j) = V12(j) + xtr*v2(j) - x3*v3(j) 
                enddo

                call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,sump)
                lhs_eqb(kk)=sump

                lhs_eqp(kk)=(lhs_eqa(kk)+lhs_eqb(kk))/two

                lhs_eqm(kk)=(lhs_eqa(kk)-lhs_eqb(kk))/two/x3

             enddo

             do i1=1,3
                ss=i1-1

                bcoeff(ss)=czero

                do kk=1,3
                   bcoeff(ss)=bcoeff(ss)+1.0_qp/three*lhs_eqp(kk)&
                        &*exp(-2.0_qp*pi*ci/three*(kk-1)*ss)
                enddo
             enddo

             b(5)=bcoeff(2)/2.0_qp
             b(2)=bcoeff(1)

             b1eff=(bcoeff(0)+b(5)*ltr**2)


             b(7)=(lhs_eqm(1)-lhs_eqm(2))/(cone -exp(2.0_qp*pi/three*ci))
             b(3)=lhs_eqm(1) - b(7)


             !------------------now we get the rest, considering v2-v4 plane

             x4=cone
             x2=sqrt(ltr**2-x4**2)
             do j=1,4
                lv(j)=V12(j)+x2*v2(j)+x4*v4(j)
             enddo

             call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

             r2 = sc(lv,v2)
             r3 = sc(lv,v3)

             lhseqa=lhseqa-b(2)*r2-b(3)*r3-b(5)*(r2**2-r3**2)-b(7)*r2*r3

             x4=x4
             x2=-x2
             do j=1,4
                lv(j)=V12(j)+x2*v2(j)+x4*v4(j)
             enddo

             call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqb)


             r2 = sc(lv,v2)
             r3 = sc(lv,v3)

             lhseqb=lhseqb-b(2)*r2-b(3)*r3-b(5)*(r2**2-r3**2)-b(7)*r2*r3

             x4=-x4
             x2=x2
             do j=1,4
                lv(j)=V12(j)+x2*v2(j)+x4*v4(j)
             enddo

             call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqc)

             r2 =  sc(lv,v2)
             r3 =  sc(lv,v3)

             lhseqc=lhseqc-b(2)*r2-b(3)*r3-b(5)*(r2**2-r3**2)-b(7)*r2*r3

             b(8)=one/(-x2)*(lhseqa - lhseqb)/two
             b(4)=half*(lhseqa-lhseqc)
             b(6)=(lhseqa-b1eff-b(4)-b(8)*(-x2))*(-one/three)
             b(1)=b1eff - ltr**2*b(6)

             !          calculation of b(9)

             x4=sqrt(ltr**2-cone)
             do j=1,4
                lv(j)=V12(j)+v3(j)+x4*v4(j)
             enddo

             call resid2(4,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

             r2 = sc(lv,v2)
             r3 = sc(lv,v3)
             r4 = sc(lv,v4)

             b(9)= (lhseqa - (b(1) + b(2)*r2 + b(3)*r3+b(4)*r4 &
                  &+ b(5)*(r2**2 - r3**2) + b(6)*(r2**2+r3**2-2*r4**2) &
                  &+ b(7)*r2*r3 + b(8)*r2*r4 ) )/x4 

             !----------------------- endif for ltr <> 0 condition
          endif


          !----------- now  ep-dependent coefficient for both ltr <> 0 cases

          xe= cone
          x4=sqrt(ltr**2-xe**2)
          do j=1,4
             lv(j)=V12(j)+ x4*v4(j)
          enddo
          lv(5)=ci*xe

          call resid2(5,lv,k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,lhseqa)

          r2 = sc(lv,v2)
          r3 = sc(lv,v3)
          r4 = sc(lv,v4)


          b(10)= (lhseqa - (b(1) + b(2)*r2 + b(3)*r3+b(4)*r4  &
               & + b(5)*(r2**2 - r3**2) + b(6)*(r2**2+r3**2-two*r4**2) &
               & + b(7)*r2*r3+ b(8)*r2*r4+ b(9)*r3*r4) )/xe**2


          ! -- stability check in D dimension 
          if (abs(ltr) > buff ) then 
             
             rn1=17._qp/19._qp
             rn2=-5._qp/7._qp
             rn3=23._qp/29._qp
             rn4=-31._qp/39._qp
             xtr = ltr/sqrt(rn1**2+rn2**2+rn3**2+rn4**2)
             do j=1,4
                lv_rn(j)=V12(j)+xtr*(rn1*v2(j)+rn2*v3(j)+rn3*v4(j))
             enddo
             lv_rn(5)=xtr*rn4*ci
             
             ldv2 = dot(lv_rn(:4),v2)
             ldv3 = dot(lv_rn(:4),v3)
             ldv4 = dot(lv_rn(:4),v4)
             lvsq = (xtr*rn4)**2
             
             call resid2(5,lv_rn(:5),k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,&
                  &Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res_v1)
             
             res_v2 = b(1)+b(2)*ldv2+b(3)*ldv3+b(4)*ldv4+b(5)*(ldv2**2-ldv3**2)&
                  &+b(6)*(ldv2**2+ldv3**2-two*ldv4**2)+b(7)*ldv2*ldv3&
                  &+b(8)*ldv2*ldv4+b(9)*ldv3*ldv4+b(10)*lvsq
             

          else
             rn1=11._qp/19._qp
             rn2=-1._qp/7._qp
             rn3= 3._qp/29._qp
             xe =sqrt(ltr**2-rn1**2-rn2**2-rn3**2) 
             do j=1,4
                lv_rn(j)=V12(j)+(rn1*v2(j)+rn2*v3(j)+rn3*v4(j))
             enddo
             lv_rn(5)=xe*ci
             
             ldv2 = dot(lv_rn(:4),v2)
             ldv3 = dot(lv_rn(:4),v3)
             ldv4 = dot(lv_rn(:4),v4)
             lvsq = xe**2
             
             call resid2(5,lv_rn(:5),k1,l2cut,f2cut,ij,N5,Lc5,F5,N4,&
                  &Lc4,F4,N3,Lc3,F3,N2,Lc2,F2,res_v1)
             
             res_v2 = b(1)+b(2)*ldv2+b(3)*ldv3+b(4)*ldv4+b(5)*(ldv2**2-ldv3**2)&
                  &+b(6)*(ldv2**2+ldv3**2-two*ldv4**2)+b(7)*ldv2*ldv3&
                  &+b(8)*ldv2*ldv4+b(9)*ldv3*ldv4+b(10)*lvsq
          endif
          
       if (abs(res_v2) > sq2tol) then 
             res_err = (res_v1-res_v2)/res_v2
          else
             res_err = (res_v1-res_v2)
          endif
          
          tagdcut(i,1)=666

          do j=1,10
             coeff2(i,j)=b(j)
             dcoeff2(i,j)=coeff2(i,j)*res_err
          enddo

          do j=1,4
             refvect2(i,j)=v2(j)
             refvect2(i,j+4)=v3(j)
             refvect2(i,j+8)=v4(j)
          enddo

          do j=1,2
             mass2(i,j) = m(j)
          enddo

          do j=1,4
             propv2(i,j)=czero
             propv2(i,j+4)=k1(j)
          enddo


       endif ! ends light-like. or not light-like

    enddo ! ends loop over double cuts

  end subroutine doubcut



end module qpopp
