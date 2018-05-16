      subroutine gosam_momenta(p,pgosam)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 p(0:3,nlegborn)
c     real * 8 pgosam(5*nlegborn)
c     In GoSam the array of the momenta has dimension 50.
c     It accounts for 10 particles at most
      integer dim_mom_array
      parameter (dim_mom_array=50)
      real * 8 pgosam(dim_mom_array)
      integer i
      
      if (nlegborn*5 .gt. dim_mom_array) then
         write(*,*) 'The dimension of the pgosam array in the '//
     $        'pwhg_gosam.f file NEEDS to be increased'
         write(*,*) 'PROGRAM ABORT'
         call exit(1)
      endif

      do i=1,nlegborn
         pgosam(1+5*(i-1))=p(0,i)
         pgosam(2+5*(i-1))=p(1,i)
         pgosam(3+5*(i-1))=p(2,i)
         pgosam(4+5*(i-1))=p(3,i)
c         write(*,*) i,p(0,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
         pgosam(5+5*(i-1))=kn_masses(i)
      enddo
      end

      
      subroutine gosam_flst_born(AlphasPower,AlphaPower)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'LesHouches.h'
      include 'pwhg_st.h'
      include 'pwhg_par.h'      
      integer AlphasPower,AlphaPower
      integer flst_gosam(nlegborn),j,k
      integer iun,error
      integer date(3), time(3)
      logical file_exists
      character * 20 stringa
      character* 1 stringa1
      if (AlphasPower.eq.0.and.AlphaPower.eq.0) then
         write(*,*) 'AlphasPower and AlphaPower are zero'
         write(*,*) 'Please fix their value in the '//
     $        'init_processes.f file'
         write(*,*) 'CANNOT PROCEED. POWHEG abort'
         call exit(1)
      endif
      call idate(date)
      call itime(time)
c     check if orderfile.lh exists. If not, create it      
      inquire(file='orderfile.lh', exist=file_exists)
      if (file_exists .eqv. .false.) then
C         call newunit(iun)
         iun=11
         open(unit=iun,file='orderfile.lh',err=200)     
         write(*,*) '************************************************'//
     $        '******'
         write(*,*) 'Writing the orderfile...'
         
         write(iun,*) '# orderfile.lh'
         write(iun,*) '# Created by POWHEG-BOX'
         write(iun,1000) date(1), date(2), date(3), time 
         write(iun,*)
         write(iun,*) 'MatrixElementSquareType CHaveraged'
         write(iun,*) 'CorrectionType          QCD'
         write(iun,*) 'IRregularisation        CDR'
         write(iun,'(1x,a,i2)') 'AlphasPower            ',AlphasPower
         write(iun,'(1x,a,i2)') 'AlphaPower             ',AlphaPower
         write(iun,*) 'SubdivideSubprocess     no'
         write(iun,*)
         write(iun,*) '# processes  list'
         do j=1,flst_nborn
            do k=1,nlegborn
               flst_gosam(k)=flst_born(k,j)
               if (flst_gosam(k).eq.0) flst_gosam(k)=21
            enddo
c            stringa='(2i4,a3,3i4)'
            write(stringa1, '(i1)' )  nlegborn-2  
            stringa="(2i4,a3,"//stringa1//"i4)"
            write(iun,stringa) (flst_gosam(k),k=1,2),' ->',
     $           (flst_gosam(k),k=3,nlegborn)
         enddo
         
         write(*,*) 'GoSam orderfile written in the GoSam_POWHEG'//
     $        ' directory'
C         write(*,*) 'Enter GoSam_POWHEG and build virtual'      
         write(*,*) '************************************************'//
     $        '******'
         stop
 200     write(*,*) '******************************************'
         write(*,*) 'Problem in writing the order file'
         write(*,*) 'CANNOT PROCEED. POWHEG execution abort'
         write(*,*) '******************************************'
         stop
      else
         write(*,*) '************************************************'//
     $        '************************'
         write(*,*) 'GoSam_POWHEG/orderfile.lh exists.'
         write(*,*) 'The POWHEG BOX assumes that the virtual.f file '//
     $        'has been correctly generated.'
         write(*,*) 'Execution continues'
         write(*,*) '************************************************'//
     $        '************************'
         
      endif
 1000 format ( ' # Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &     i2.2, ':', i2.2, ':', i2.2 )
      
      end

