*--   Author :   Andrea Messina
C-----------------------------------------------------------------------------
      subroutine inithbook
C-----------------------------------------------------------------------------

      
      implicit none
      include 'hnt.inc'
C...EVENT FILE DATA 
      INTEGER NUNIT,NUNITOUT 
C...PROCESS CODE 
      INTEGER IHRD 
C...PDF SET TYPE 
      INTEGER NDNS 
      CHARACTER PDFTYP*25 
C...TOTAL NUMBER OF PARTONS 
      INTEGER NPART 
      CHARACTER*100 FILENAME 
      COMMON /EVFLUP/ NUNIT,NUNITOUT,IHRD,NPART,NDNS,PDFTYP,FILENAME 
      character*120 name
      integer IST


      call hlimit(nwpawc)
      CALL STRCATH(FILENAME,'_py.ntpl',name)
      call hropen(200,'HNTPL',name,'N',1024,IST)
      if (IST .NE. 0) then
         print*,'Error opening ntuple'
         stop
      endif
      

      call hbnt(200,'HNTPL',' ')
       
      call hbname(200,'HEPG',ngen,
     +     'ngen[0,900],E_gn(ngen),Px_gn(ngen),Py_gn(ngen),Pz_gn(ngen),
     +      Vt_gn(ngen),Vx_gn(ngen),Vy_gn(ngen),Vz_gn(ngen),id_gn(ngen),
     +      istd_gn(ngen),jda1_gn(ngen),jda2_gn(ngen),jmo1_gn(ngen),
     +      jmo2_gn(ngen),m_gn(ngen):r')

      end


C-----------------------------------------------------------------------------
      subroutine fillntuple
C-----------------------------------------------------------------------------

      include 'hnt.inc'
      include 'hepevt.inc'
C      integer bevt
C
      integer i

      ngen=0
C      bevt=0
C      do i=1,50
C         if(abs(idhep(i)).eq.5) bevt=1
C      enddo
      do i=1,nhep
c         if (ISTHEP(i).eq.1.and.bevt.eq.1) then
         ngen=ngen+1
         E_gn(ngen) =real(phep(4,i))
         Px_gn(ngen)=real(phep(1,i))
         Py_gn(ngen)=real(phep(2,i))
         Pz_gn(ngen)=real(phep(3,i))
         Vt_gn(ngen)=real(vhep(4,i))
         Vx_gn(ngen)=real(vhep(1,i))
         Vy_gn(ngen)=real(vhep(2,i))
         Vz_gn(ngen)=real(vhep(3,i))
         id_gn(ngen)=idhep(i)
         istd_gn(ngen)=isthep(i)
         jda1_gn(ngen)=jdahep(1,i)
         jda2_gn(ngen)=jdahep(2,i)
         jmo1_gn(ngen)=jmohep(1,i)
         jmo2_gn(ngen)=jmohep(2,i)
         m_gn(ngen)=real(phep(5,i))
c     endif
      enddo

      call hfnt(200)

      end

C-----------------------------------------------------------------------------
      subroutine closentuple
C-----------------------------------------------------------------------------
      implicit none

      integer ICYCLE
      call hrout(0,ICYCLE,' ')
      call hrend('HNTPL')

      return
      end
