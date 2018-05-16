      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'     
      double precision p(0:3,nlegreal)
      integer rflav(nlegreal)
      double precision amp2
      integer mu
      logical ini
      data ini/.true./
      save ini
cccccccccccccccccccc
c$$$      logical debug
c$$$      parameter(debug=.false.)
      logical compare_with_DUW_paper
      parameter(compare_with_DUW_paper=.false.)
      double precision kr_mad(0:3,nlegreal)
c$$$      double precision amp2uwer,tmp
      integer ileg
      double precision tiny
      parameter (tiny=1d-4)
c colored output
      character*1 red(5)
      character*1 reset(4) 
      data red /' ','[','3','1','m'/ 
      data reset /' ','[','0','m'/ 
      red(1) = char(27)     ! Escape character (ASCII 27).
      reset(1) = char(27)
c
      if(ini) then
c     Set MADGRAPH parameters at first call
         call set_madgraph_parameters
c$$$          if (debug) then
c$$$c     compare with madgraph from P.Uwer
c$$$             write(*,*) "REAL: PERFORMING CHECK AGAINST P.UWER MAdGRAPH"
c$$$             write(*,*) "      AND USING THEM FOR REALS "
c$$$c     setting madgraph parameters (needed for madgraph subroutines)
c$$$            call set_madgraph_parameters_uwer
c$$$         endif
      endif
      ini=.false.     
      amp2 = 0d0
      if (abs(rflav(3)).le.5.or.abs(rflav(4)).le.5) then
         write(*,*) 'real_ampsq: ERROR in flavor assignement'
         write(*,*) rflav(1),' ',rflav(2),'->',rflav(3),' ',rflav(4),' '
     $        ,rflav(5),' ',rflav(6)
         call exit(1)
      endif
ccccccccccccccccccc
      if(compare_with_DUW_paper) then
         call check_DUW_paper_real
      endif
cccccccccccccccccccccc

c     Reset MADGRAPH parameters (those 
c     that depends on the kinematics) at each call
      call mad_setparam
      do ileg=1,nlegreal
         do mu=0,3
            kr_mad(mu,ileg)=p(mu,ileg)
         enddo
      enddo

c     to avoid bugs in HELAS, restore exact masslessness of incoming partons 
      kr_mad(0,1)=dabs(kr_mad(3,1))
      kr_mad(0,2)=dabs(kr_mad(3,2))
      call compreal(kr_mad,rflav,amp2)
      
      if(amp2.eq.0d0) then
        write(*,*) 'WARNING real_ampsq: returning 0 amplitude!'
        write(*,*) rflav(1),' ',rflav(2),'->',rflav(3),' ',rflav(4),' '
     $        ,rflav(5),' ',rflav(6)
      endif


c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$$$c compare with madgraph from P.Uwer library
c$$$      if (debug) then
c$$$c     set madgraph parameters that can change on an event-by-event basis
c$$$      
c$$$      call mad_setparam_uwer
c$$$      call compreal_uwer(kr_mad,rflav,amp2uwer)
c$$$
c$$$      tmp=abs(amp2uwer/amp2 -1d0)
c$$$
c$$$      if  (tmp.gt.tiny) then 
c$$$         write(*,*) rflav,red,'\n REAL: MUST BE EQUAL =====> ',
c$$$     $        amp2uwer,amp2,'\n RATIO: ', amp2uwer/amp2,reset
c$$$          
c$$$         write(*,*) "P:",((kr_mad(mu,ileg),mu=0,3),"\n",ileg=1,nlegreal)
c$$$         stop
c$$$      endif
c$$$
c$$$      tmp=(abs(amp2uwer-amp2)/amp2uwer)
c$$$
c$$$      if  (tmp.gt.tiny) then 
c$$$         write(*,*) rflav,red
c$$$     $        ,'\n REAL: % DIFFERENCE MUST BE 0 =====> ',abs(amp2uwer
c$$$     $        -amp2)/amp2uwer,reset
c$$$         stop
c$$$      endif
c$$$
c$$$c if checks are passed, use new results for real      
c$$$
c$$$      amp2=amp2uwer
c$$$      
c$$$      endif
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2*pi))     
      
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC    CHECKS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_DUW_paper_real
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'     
      integer rflav(nlegreal)
      double precision kr_mad(0:3,nlegreal),amp2mad
      double precision amp2uwer,uwer_paper_res,symm_fac
c$$$c     MADGRAPH routines in Uwer's library
c$$$      double precision SGG_TTBGG, SGG_TTBUUB,SUBDB_TTBUBDB,
c$$$     $     SUBUB_TTBUBUB, SUDB_TTBUDB,SUD_TTBUD,SUG_TTBUG,SUUB_TTBDDB
c$$$     $     ,SUUB_TTBGG,SUUB_TTBUUB,SUU_TTBUU,SGUB_TTBUBG
c$$$      external SGG_TTBGG, SGG_TTBUUB,SUBDB_TTBUBDB,SUBUB_TTBUBUB,
c$$$     $     SUDB_TTBUDB,SUD_TTBUD,SUG_TTBUG,SUUB_TTBDDB ,SUUB_TTBGG
c$$$     $     ,SUUB_TTBUUB,SUU_TTBUU,SGUB_TTBUBG     
      integer ileg,mu
      double precision tiny
      parameter (tiny=1d-5)
c colored output
      character*1 red(5)
      character*1 reset(4) 
      data red /' ','[','3','1','m'/ 
      data reset /' ','[','0','m'/ 
      red(1) = char(27)     ! Escape character (ASCII 27).
      reset(1) = char(27)
      print *,red,"###########################################",reset
      print *,red,"###########################################",reset
      print *,red,"    CHECK REALS AGAINST DUW PAPER          ",reset
      print *,red,"###########################################",reset
      print *,red,"###########################################",reset
      write(*,*) "P:"


c NUMBERS to compare with P.Uwer paper
      kr_mad(0,1)=2100d0              
      kr_mad(1,1)=-0d0                
      kr_mad(2,1)=-0d0              
      kr_mad(3,1)=2100d0              
c                                   
      kr_mad(0,2)=2800d0              
      kr_mad(1,2)=-0d0                
      kr_mad(2,2)=-0d0                
      kr_mad(3,2)=-2800d0             
c                                   
      kr_mad(0,3)=1581.118367308447d0 
      kr_mad(1,3)=1254.462316247655d0 
      kr_mad(2,3)=-766.9360998604944d0
      kr_mad(3,3)=-554.7905976902205d0
c                                   
      kr_mad(0,4)=1460.449317799282d0 
      kr_mad(1,4)=-975.9731477430979d0
      kr_mad(2,4)=-466.5314749495881d0
      kr_mad(3,4)=965.6402060944737d0 
c                                   
      kr_mad(0,5)=545.4084744819d0    
      kr_mad(1,5)=218.7220720302516d0 
      kr_mad(2,5)=472.0439121434804d0 
      kr_mad(3,5)=-163.7241712507502d0
c                  
      kr_mad(0,6)=1313.023840410371d0 
      kr_mad(1,6)=-497.2112405348086d0
      kr_mad(2,6)=761.423662666602d0  
      kr_mad(3,6)=-947.1254371535031d0
            
      write(*,*) ((kr_mad(mu,ileg),mu=0,3),'\n',ileg=1,nlegreal)
                  
      call checkmomzeronew(nlegreal,kr_mad)
      do ileg=1,nlegreal
         call checkmassnew(ileg,kr_mad)
      enddo

c$$$      call set_madgraph_parameters_uwer
c$$$      call mad_setparam_uwer

      rflav(1)=0
      rflav(2)=0
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=0
      rflav(6)=0

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=sgg_ttbgg(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3),kr_mad(0
c$$$     $     ,4),kr_mad(0,5),kr_mad(0,6)) 

      uwer_paper_res=3.12815868347843d-09
      symm_fac=4

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset
      print *,"MISSING FACTOR ",symm_fac

      
      rflav(1)=1
      rflav(2)=-1
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=0
      rflav(6)=0

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=suub_ttbgg(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3),kr_mad(0
c$$$     $     ,4),kr_mad(0,5),kr_mad(0,6))

      uwer_paper_res=4.48308845446477d-10
      symm_fac=4

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset
      print *,"MISSING FACTOR ",symm_fac
      
      
      rflav(1)=1
      rflav(2)=0
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=0
      rflav(6)=1

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=sug_ttbug(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3),kr_mad(0
c$$$     $     ,4),kr_mad(0,6),kr_mad(0,5))

      uwer_paper_res=1.10256509258713d-10
      symm_fac=4

      
      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset
      print *,"MISSING FACTOR ",symm_fac
      
      rflav(1)=-1
      rflav(2)=0
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=-1
      rflav(6)=0

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=sgub_ttbubg(kr_mad(0,2),kr_mad(0,1),kr_mad(0,3),kr_mad(0
c$$$     $     ,4),kr_mad(0,5),kr_mad(0,6))

      uwer_paper_res=1.384600673183816d-10
      symm_fac=4

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset
      print *,"MISSING FACTOR ",symm_fac
      
      rflav(1)=0
      rflav(2)=0
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=-1
      rflav(6)=1

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=sgg_ttbuub(kr_mad(0,2),kr_mad(0,1),kr_mad(0,3),kr_mad(0
c$$$     $     ,4),kr_mad(0,6),kr_mad(0,5))

      uwer_paper_res=2.42841040229558d-10
      symm_fac=4*5

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red,' MUST BE 1 =====> '
     $     ,symm_fac*amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $      ,reset
      print *,"MISSING FACTOR ",symm_fac
      
      
      rflav(1)=1
      rflav(2)=2
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=2
      rflav(6)=1

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=sud_ttbud(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3),kr_mad(0
c$$$     $     ,4),kr_mad(0,6),kr_mad(0,5))

      uwer_paper_res=4.44137855516180d-12
      symm_fac=4

      
      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red,' MUST BE 1 =====> '
     $     ,symm_fac*amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res 
     $     ,reset
      print *,"MISSING FACTOR ",symm_fac
      
  
      rflav(1)=-1
      rflav(2)=-2
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=-1
      rflav(6)=-2

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=subdb_ttbubdb(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3)
c$$$     $     ,kr_mad(0,4),kr_mad(0,5),kr_mad(0,6))

      uwer_paper_res=1.733763330485899d-11
      symm_fac=4

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset
      print *,"MISSING FACTOR ",symm_fac
      

      rflav(1)=1
      rflav(2)=-2
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=-2
      rflav(6)=1

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=sudb_ttbudb(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3)
c$$$     $     ,kr_mad(0,4),kr_mad(0,6),kr_mad(0,5))

      uwer_paper_res=4.796260245409952d-12
      symm_fac=4

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset 
      print *,"MISSING FACTOR ",symm_fac
      

      rflav(1)=1
      rflav(2)=-1
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=-2
      rflav(6)=2

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=suub_ttbddb(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3),kr_mad(0
c$$$     $ ,4),kr_mad(0,6),kr_mad(0,5))

      uwer_paper_res=6.13924303047741d-11
      symm_fac=16

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset
      print *,"MISSING FACTOR ",symm_fac

    

      rflav(1)=1
      rflav(2)=1
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=1
      rflav(6)=1

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=suu_ttbuu(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3),kr_mad(0
c$$$     $     ,4),kr_mad(0,5),kr_mad(0,6))

      uwer_paper_res=1.371477814148721d-11
      symm_fac=4

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset 
      print *,"MISSING FACTOR ",symm_fac
      
      

      rflav(1)=-1
      rflav(2)=-1
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=-1
      rflav(6)=-1

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=subub_ttbubub(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3)
c$$$     $     ,kr_mad(0,4),kr_mad(0,5),kr_mad(0,6))

      uwer_paper_res=1.411042000289490d-11
      symm_fac=4

      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset 
      print *,"MISSING FACTOR ",symm_fac
     

      rflav(1)=1
      rflav(2)=-1
      rflav(3)=6
      rflav(4)=-6
      rflav(5)=-1
      rflav(6)=1

      call compreal(kr_mad,rflav,amp2mad)

c$$$      amp2uwer=suub_ttbuub(kr_mad(0,1),kr_mad(0,2),kr_mad(0,3)
c$$$     $     ,kr_mad(0,4),kr_mad(0,6),kr_mad(0,5))

  
      uwer_paper_res=2.054843839960259d-11
      symm_fac=4

      
      print *,'as =',st_alpha
      print *,rflav
      print *,amp2mad*symm_fac,uwer_paper_res!,amp2uwer*symm_fac
      print *,red
     $     ,' MUST BE 1 =====> ',symm_fac
     $     *amp2mad/uwer_paper_res!,amp2uwer*symm_fac/uwer_paper_res
     $     ,reset 
      print *,"MISSING FACTOR ",symm_fac
      
      


      stop "COMPARISON WITH P.UWER NUMBERS FINISHED. PROGRAM STOPPED"
  
      
      end

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cccccc   ROUTINES TO INTERFACE TO UWER'S CODE       cccccccccc
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$
c$$$      subroutine set_madgraph_parameters_uwer
c$$$      implicit none
c$$$      include 'PhysPars.h'
c$$$      DOUBLE PRECISION   FMASS(12), FWIDTH(12)
c$$$      COMMON /FERMIONS/ FMASS,     FWIDTH
c$$$      integer i
c$$$      do i=1,12
c$$$         fmass(i)=0
c$$$         fwidth(i)=0
c$$$      enddo
c$$$      fmass(11)=ph_topmass
c$$$      end
c$$$
c$$$      subroutine mad_setparam_uwer
c$$$      implicit none
c$$$      include 'pwhg_st.h'
c$$$      include 'pwhg_math.h'
c$$$      DOUBLE PRECISION  GG(2),G
c$$$      COMMON /COUPQCD/ GG,    G
c$$$      G = DSQRT(4d0*PI*st_alpha) 
c$$$      GG(1) = -G
c$$$      GG(2) = -G     
c$$$      end
c$$$
c$$$
c$$$
c$$$      subroutine compreal_uwer(p,flav,amp2)
c$$$      implicit none
c$$$      real * 8 p(0:3,1: 6)
c$$$      integer flav( 6)
c$$$      real * 8 amp2
c$$$      real * 8 madp(0:3,1: 6)
c$$$      integer madflav( 6),perm( 6)
c$$$      logical flavequiv_perm
c$$$      external flavequiv_perm
c$$$      integer i, mu
c$$$c     MADGRAPH routines in Uwer's library
c$$$      double precision SGG_TTBGG, SGG_TTBUUB,SUBDB_TTBUBDB,
c$$$     $     SUBUB_TTBUBUB, SUDB_TTBUDB,SUD_TTBUD,SUG_TTBUG,SUUB_TTBDDB
c$$$     $     ,SUUB_TTBGG,SUUB_TTBUUB,SUU_TTBUU,SGUB_TTBUBG
c$$$      external SGG_TTBGG, SGG_TTBUUB,SUBDB_TTBUBDB,SUBUB_TTBUBUB,
c$$$     $     SUDB_TTBUDB,SUD_TTBUD,SUG_TTBUG,SUUB_TTBDDB ,SUUB_TTBGG
c$$$     $     ,SUUB_TTBUUB,SUU_TTBUU,SGUB_TTBUBG      
c$$$      madflav(1)=5
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUU_TTBUU(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=5
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-4
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBUB_TTBUBUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-5
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-4
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGG_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGG_TTBGG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGG_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGG_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGG_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGG_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=0
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUD_TTBUD(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUG_TTBUG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUU_TTBUU(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUU_TTBUU(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUU_TTBUU(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUU_TTBUU(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$         
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,1),madp(0,2),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=2
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=1
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=3
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=4
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-4
c$$$      madflav(6)=5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$         
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=-5
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-4
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$         
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUDB_TTBUDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,6),madp(0,5))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBDB_TTBUBDB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=0
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SGUB_TTBUBG(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=5
c$$$      madflav(6)=-5
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBDDB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=0
c$$$      madflav(6)=0
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBGG(madp(0,2),madp(0,1),madp(0,3),madp(0,4),madp(0
c$$$     $        ,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUUB_TTBUUB(madp(0,2),madp(0,1),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-2
c$$$      madflav(2)=-2
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-2
c$$$      madflav(6)=-2
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBUB_TTBUBUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-1
c$$$      madflav(2)=-1
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-1
c$$$      madflav(6)=-1
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBUB_TTBUBUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-3
c$$$      madflav(2)=-3
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-3
c$$$      madflav(6)=-3
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBUB_TTBUBUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      madflav(1)=-4
c$$$      madflav(2)=-4
c$$$      madflav(3)=6
c$$$      madflav(4)=-6
c$$$      madflav(5)=-4
c$$$      madflav(6)=-4
c$$$      if (flavequiv_perm( 6,flav,madflav,perm)) then
c$$$         do i=1, 6 
c$$$            do mu=0,3
c$$$               madp(mu,perm(i))=p(mu,i)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         amp2= SUBUB_TTBUBUB(madp(0,1),madp(0,2),madp(0,3),madp(0,4)
c$$$     $        ,madp(0,5),madp(0,6))
c$$$         return
c$$$      endif
c$$$
c$$$      write(*,*) 'ERROR: the flavour list', flav
c$$$      write(*,*) 'is not in the list of Uwer  MADGRAPH routines'
c$$$      call exit(1)
c$$$      end








