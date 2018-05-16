c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_br.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      double precision p(0:3,nlegs)
      integer vflav(nlegs)
      double precision virtual,powheginput
      external powheginput
      logical ini
      data ini/.true./ 
      save ini
      real *8 dotp
      external dotp
c     
      logical use_OLP_Interface
      common /colp/use_OLP_Interface
c
      integer OLP_code,jb,j
c     from born flavour structure to jborn (for simplicity, here
c     use only MASSLESS colored particles, sorted as:
c     id_plus,id_minus,id_final)
c     Particle 3 and 4 are always heavy colored particles
      integer bornflst2pwhgcode(-5:5,-5:5,-5:5)
      common/cbornflst2pwhgcode/bornflst2pwhgcode
c     from powheg born code to the OLP code
      integer pwhgcode2OLPcode(maxprocborn)
      common/cpwhgcode2OLPcode/pwhgcode2OLPcode
c     
      real * 8 c(-6:6),gamma(-6:6),gammap(-6:6)
      save c,gamma,gammap
      double precision bcut,largecorrfact
      logical checklargecorr
      common /ccheckvirtuals/bcut,largecorrfact,checklargecorr
      save /ccheckvirtuals/
      logical fftestflag
      common /cfftestflag/fftestflag
      save /cfftestflag/
      if (ini) then
         bcut=powheginput("#bcut") 
         largecorrfact=powheginput("#largecorrfact") 
         fftestflag=.false.
         if(powheginput("#ffltest").eq.1d0) fftestflag=.true.
c from 2.100 of FNO2007
         do j=-6,6
            if(j.eq.0) then
               c(j)=ca
               gamma(j)=(11*ca-4*tf*st_nlight)/6
               gammap(j)=(67d0/9-2*pi**2/3)*ca-23d0/9*tf*st_nlight
            else
               c(j)=cf
               gamma(j)=3d0/2*cf
               gammap(j)=(13d0/2-2*pi**2/3)*cf
            endif
         enddo
         if(use_OLP_interface) then
            call virtual_initialize_OLP
         else
c     initialize SM parameters and Uwer's routines
            call virtual_initialize
         endif
         ini=.false.
      endif

      if(.not.use_OLP_interface) then
         call virtual_evaluate(p,vflav,virtual)
      else
         OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $        ,vflav(5)))
         call virtual_OLP(p,OLP_code,virtual)
      endif


c     Add couplings
!     This is the finite part of the virtual contribution ===> Vfin
!     a factor as/(2pi) is missing as required by sigsoftvirt
!     The division by 2 is the leftover from having factorized out
!     as/(2pi) instead of as/(4pi) 
      virtual=(4*pi*st_alpha)**3 * virtual  /2d0

c     missing terms due to different definitions of Vfin in Uwer's code
c     first identify the flst label
      jb=bornflst2pwhgcode(vflav(1),vflav(2),vflav(5))
c     then subtract the missing term avoiding recalculating born
      virtual= virtual - 4d0*dotp(p(0,1),p(0,2)) * br_born(jb)
     $     *(c(flst_born(1,jb))+c(flst_born(2,jb))+c(flst_born(5,jb)))
     $     *pi*pi/6d0

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine virtual_evaluate(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_br.h'
      include 'pwhg_kn.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      double precision p(0:3,nlegs),pp(0:3,nlegs)
      integer vflav(nlegs),ileg,mu      
      double precision virtual,tmp
c     External functions
      double precision gggtt2virtfin
      external gggtt2virtfin
      double precision qqttg2virtfin
      external qqttg2virtfin
      double precision qgttq2virtfin
      external qgttq2virtfin
      double precision gqbttqb2virtfin
      external gqbttqb2virtfin
c     Variables needed for debugging
      logical debug
      parameter(debug=.false.)
      logical compare_with_DUW_paper
      parameter(compare_with_DUW_paper=.false.)
      double precision bcut,largecorrfact
      logical checklargecorr
      common /ccheckvirtuals/bcut,largecorrfact,checklargecorr
      save /ccheckvirtuals/
c     

      if((compare_with_DUW_paper).or.(debug)) then
         checklargecorr=.false.
      else 
         checklargecorr=.true. !!!! VERY IMPORTANT TO ALWAYS SET THIS
      endif       
c
c     set parameters that may change event-by-event
         call set_virtual_scales
c     Since the old definition of spinors in Uwer's library
c     does not allow for massless spinors parallel to z-axis
c     we apply a set of rotations x->y , y->z, z->x
         do ileg=1, nlegs
            pp(0,ileg)=p(0,ileg)
            pp(1,ileg)=p(3,ileg)
            pp(2,ileg)=p(1,ileg)
            pp(3,ileg)=p(2,ileg)
         enddo
c     to avoid bugs, restore exact masslessness of incoming partons
         pp(0,1)=dabs(pp(1,1))
         pp(0,2)=dabs(pp(1,2))

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CALL  DEBUGGING ROUTINES
         if(compare_with_DUW_paper) then
            call check_DUW_paper_virtual
         endif
         if(debug) then
c     Perform an extensive series of checks
c     really slow !!! Turn off for production
c     write(*,*) "P:",((p(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
c     write(*,*) "PP:",((pp(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
            call virtual_checks(pp,vflav)
         endif
c     ccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Start calculating virtuals
         if(vflav(1).ne.0.and.vflav(1)+vflav(2).eq.0) then
            if(vflav(1).gt.0) then
c     q qbar
               virtual=qqttg2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $              ,pp(0,4))/36d0
            else
c     qbar q
               virtual=qqttg2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $              ,pp(0,4))/36d0
            endif
         endif
         if(vflav(1)*vflav(2).eq.0) then
            if(vflav(1).gt.0) then
c     a quark and a gluon
               virtual=qgttq2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96d0
            elseif(vflav(1).lt.0) then
c     an anti-quark and a gluon
               virtual=gqbttqb2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96d0
            elseif(vflav(2).lt.0) then
c     a gluon and an anti-quark
               virtual=gqbttqb2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0 ,3)
     $              ,pp(0,4))/96d0
            elseif(vflav(2).gt.0) then
c     a gluon and a quark
               virtual=qgttq2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96d0  
            endif
         endif
c     two gluons
         if((vflav(1).eq.0).and.(vflav(2).eq.0))then      
            virtual=gggtt2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $           ,pp(0,4))/256d0

         endif
         
         end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC       OLP INTERFACE ROUTINES           CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine virtual_initialize_OLP
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_flst.h'
      logical debug
      parameter (debug=.false.)
c     from powheg born code to the OLP code
      integer pwhgcode2OLPcode(maxprocborn)
      common/cpwhgcode2OLPcode/pwhgcode2OLPcode


      character*11 filename
      integer iproc
      integer iun
      integer ios,j,k,l
      character *5 status
      integer code,foundbproc

      character * 100 line,line0

      write(*,*) 
      write(*,*) ' Checking contract file for external OLP '
      write(*,*) 
      filename="contract.lh"
      call OLP_Start(ios,filename//CHAR(0))
      if (ios.ne.1) then
         write(*,*) ' Error: OLP cannot handle contract file'
         call exit(1)
      endif
      
c     open the negotiation file, read it and
c     fill the array pwhgbcode2OLPcode properly
      call newunit(iun)
      open(unit=iun,file='contract.lh',status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*) 'cannot open contract.lh'
         call exit(1)
      endif

      foundbproc=0
      do l=1,maxprocborn
 111     continue
         line0=' '
         read(unit=iun,fmt='(a)',iostat=ios) line0
         if(debug) write(*,*) line0
         if(ios.ne.0.and.line0.eq.' ') goto 10
c     this means end of file...
         line=line0
         if(line(11:12).ne.'->') then
c     this means that current line is not a line with a 2->3 subprocess
            if(debug) write(*,*) 'Found a line without subprocess'
            goto 111
         endif
         do k=1,100
            if(line(k:k).eq.'#'.or.line(k:k).eq.'!') then
               print*, 'commented line'
               line(k:)=' '
               goto 123
            endif
            if(line(k:k).eq.'|') then
               line=line(k+1:)
               read(unit=line,fmt=*,iostat=ios) status,code
               if(debug) then
                  write(*,*) '\t line with info: '//
     $                 'status,OLPcode: ',status,code
               endif
               if(status.eq."OK") then
                  foundbproc=foundbproc+1
                  pwhgcode2OLPcode(foundbproc)=code
                  goto 123                  
               endif
            endif
         enddo
 123     continue
      enddo

 10   continue

      if(foundbproc.ne.flst_nborn) then
         write(*,*)'**********************************'
         write(*,*)
     $ 'ERROR: OLP and POWHEG-BOX have different number of'//
     $  ' born subprocesses: ',foundbproc,flst_nborn
         write(*,*)'**********************************'
         call exit(1)
      endif

      close(iun)

c     After the mapping between OLP and POWHEG code has been 
c     performed call OLP_Start again to initialize the virtual 
c     routines inside the OLP
      call OLP_Start(ios,filename//CHAR(0))
      if (ios.ne.-1) then
         write(*,*) ' Error: OLP cannot handle contract file'
         call exit(1)
      endif
      end


      subroutine virtual_OLP(p,code,virtual)
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_st.h'
      include '../include/pwhg_kn.h'
      real *8 p(0:3,nlegborn)
      integer code
      real *8 virtual
      real *8 p_olp(0:4,nlegborn)
      real *8 scales(4) ! these are: mu_ren,mu_fac,Q_ES,mu_reg
      real *8 virt_wgts(4) ! these are: pole2,pole1,pole0,born
      integer mu,ileg
      real *8 s
      real *8 dotp
      external dotp

c     Since the definition of spinors in Uwer's library
c     does not allow for massless spinors parallel to z-axis
c     we apply a set of rotations x->y , y->z, z->x
      do ileg=1, nlegborn
         p_olp(0,ileg)=p(0,ileg)
         p_olp(1,ileg)=p(3,ileg)
         p_olp(2,ileg)=p(1,ileg)
         p_olp(3,ileg)=p(2,ileg)
         p_olp(4,ileg)=0d0
      enddo
c     to avoid bugs, restore exact masslessness of incoming partons 
      p_olp(0,1)=dabs(p_olp(1,1))
      p_olp(0,2)=dabs(p_olp(1,2))
c     sets the masses of heavy particles
      p_olp(4,3)=kn_masses(3)
      p_olp(4,4)=kn_masses(4)

c     assign scales
      scales(1)=sqrt(st_muren2)
      scales(2)=sqrt(st_mufact2)
      scales(3)=sqrt(st_muren2)
      scales(4)=sqrt(st_muren2)

      call OLP_EvalSubProcess(code,p_olp,scales,1d0,virt_wgts)

      
      virtual=virt_wgts(3)

      end





CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC    ROUTINES USEFUL FOR  CHECKS         CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_DUW_paper_virtual
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_br.h'
      include 'pwhg_kn.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      double precision huge
      parameter (huge=1d100)
      double precision kr_mad(0:3,nlegs)
      double precision pp(0:3,nlegs)
      integer vflav(nlegs),ileg,mu      
      double precision virtual,tmp,born,pole2,pole1,pole0,born_uwer
      double precision uwer_res_cm2,uwer_res_cm1,uwer_res_cm0
      double precision iop_cm2,iop_cm1,iop_cm0
c     External virtual functions
      double precision gggtt2virtfin
      external gggtt2virtfin
      double precision qqttg2virtfin
      external qqttg2virtfin
      double precision qgttq2virtfin
      external qgttq2virtfin
      double precision gqbttqb2virtfin
      external gqbttqb2virtfin
c     External Born functions
      double precision gggtt2
      external gggtt2
      double precision qqttg2
      external qqttg2
      double precision qgttq2
      external qgttq2
      double precision gqbttqb2
      external gqbttqb2
c     colored output
      character*1 red(5)
      character*1 reset(4) 
      data red /' ','[','3','1','m'/ 
      data reset /' ','[','0','m'/ 
      red(1) = char(27)         ! Escape character (ASCII 27).
      reset(1) = char(27)
      
      print *,red,"###########################################",reset
      print *,red,"###########################################",reset
      print *,red," CHECK VIRTUAL ROUTINES AGAINST DUW PAPER ",reset
      print *,red,"###########################################",reset
      print *,red,"###########################################",reset

     
      write(*,*) "MOMENTA P:"
      kr_mad(0,1)=500d0
      kr_mad(1,1)=0d0
      kr_mad(2,1)=0d0
      kr_mad(3,1)=500d0
c     
      kr_mad(0,2)=500d0
      kr_mad(1,2)=0d0
      kr_mad(2,2)=0d0
      kr_mad(3,2)=-500d0
c     
      kr_mad(0,3)=458.5331753852783d0
      kr_mad(1,3)=207.0255169909440d0
      kr_mad(2,3)=0d0
      kr_mad(3,3)=370.2932732896167d0
c     
      kr_mad(0,4)=206.6000026080000d0
      kr_mad(1,4)=-10.65693677252589d0
      kr_mad(2,4)=42.52372780926147d0
      kr_mad(3,4)=-102.3998210421085d0
c     
      kr_mad(0,5)=334.8668220067217d0
      kr_mad(1,5)=-196.3685802184181d0
      kr_mad(2,5)=-42.52372780926147d0
      kr_mad(3,5)=-267.8934522475083d0

      write(*,'(4f22.13,a)') ((kr_mad(mu,ileg),mu=0,3),'\n',ileg=1,nlegs)

      write(*,*) "VIRTUAL: PERFORMING ROTATION OF MOMENTA  "
      write(*,*) "         TO COMPLY WITH SPINOR DEFINITION"
      write(*,*) "NEW MOMENTA P':"

c     Since the old definition of spinors in Uwer's library
c     does not allow for massless spinors parallel to z-axis
c     we apply a set of rotations x->y , y->z, z->x
      do ileg=1, nlegs
         pp(0,ileg)=kr_mad(0,ileg)
         pp(1,ileg)=kr_mad(3,ileg)
         pp(2,ileg)=kr_mad(1,ileg)
         pp(3,ileg)=kr_mad(2,ileg)
      enddo
      write(*,'(4f22.13,a)') ((pp(mu,ileg),mu=0,3),'\n',ileg=1,nlegs)

      call checkmomzeronew(nlegs,pp)
      do ileg=1,nlegs
         call checkmassnew(ileg,pp)
      enddo

      vflav(1)=0
      vflav(2)=0
      vflav(3)=6
      vflav(4)=-6
      vflav(5)=0
      

      uwer_res_cm2=-0.1540118420981379 
      uwer_res_cm1= 0.0731096895036588       	
      uwer_res_cm0=0.5295183452346090	
      
      
      iop_cm2=0.1540118421074569
      iop_cm1=-0.0731096894943435
      iop_cm0=-0.5280576886301999

c    output of Uwer's routine is the Born without couplings and before
c    the average on color and helicity. Corresponds to 4 Nc a_0 of formula A.2
      born=gggtt2(pp(0,1),pp(0,2),pp(0,5),pp(0,3),pp(0,4))


c     call MNR_ttbarj(kr_mad,vflav,born)
c     born=born*256d0/(4d0*pi*st_alpha)**3 ! remove couplings and averaging from MNR results to compare

      pole0=gggtt2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR2(huge)
      pole2=gggtt2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR2(0d0)
      call setDELTAIR1(huge)
      pole1=gggtt2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR1(0d0)

      pole2=pole2/huge/born*(st_alpha/(4d0*pi))
      pole1=pole1/huge/born*(st_alpha/(4d0*pi))
      pole0=pole0/born*(st_alpha/(4d0*pi))
     
      print *,'as =',st_alpha
      print *,vflav
      print *,' CM2, MYPOLE2'
      print *, uwer_res_cm2,pole2
      print *,red,' RATIO MUST BE 1 =====> ',pole2 /uwer_res_cm2 ,reset
      print *,' CM1, MYPOLE1'
      print *, uwer_res_cm1,pole1
      print *,red,' RATIO MUST BE 1 =====> ',pole1 /uwer_res_cm1 ,reset
      print *,' CM0, MYPOLE0'
      print *, uwer_res_cm0,pole0
      print *,red,' RATIO MUST BE 1 =====> ',pole0 /uwer_res_cm0 ,reset


      vflav(1)=1
      vflav(2)=-1
      vflav(3)=6
      vflav(4)=-6
      vflav(5)=0


      uwer_res_cm2=-0.0969704191047176
      uwer_res_cm1=-0.0126983208241891       	
      uwer_res_cm0=0.2435672439083931

      iop_cm2=0.0969704191046952
      iop_cm1=0.0126983208241661
      iop_cm0=-0.3992776407671517


c    output of Uwer's routine is the Born without couplings and before
c    the average on color and helicity. Corresponds to 4 Nc a_0 of formula A.2
      born=qqttg2(pp(0,1),pp(0,2),pp(0,5),pp(0,3),pp(0,4))

c     call MNR_ttbarj(kr_mad,vflav,born)
c     born=born*36d0/(4d0*pi*st_alpha)**3 ! remove couplings and averaging from MNR results to compare
  
      pole0=qqttg2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR2(huge)
      pole2=qqttg2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR2(0d0)
      call setDELTAIR1(huge)
      pole1=qqttg2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR1(0d0)

      pole2=(pole2/huge)/born*(st_alpha/(4d0*pi))
      pole1=(pole1/huge)/born*(st_alpha /(4d0*pi))
      pole0=pole0/born*(st_alpha/(4d0*pi))
     
      print *,'as =',st_alpha
      print *,vflav
      print *,' CM2, MYPOLE2'
      print *, uwer_res_cm2,pole2
      print *,red,' RATIO MUST BE 1 =====> ',pole2 /uwer_res_cm2 ,reset
      print *,' CM1, MYPOLE1'
      print *, uwer_res_cm1,pole1
      print *,red,' RATIO MUST BE 1 =====> ',pole1 /uwer_res_cm1 ,reset
      print *,' CM0, MYPOLE0'
      print *, uwer_res_cm0,pole0
      print *,red,' RATIO MUST BE 1 =====> ',pole0 /uwer_res_cm0 ,reset

      
      vflav(1)=1
      vflav(2)=0
      vflav(3)=6
      vflav(4)=-6
      vflav(5)=1

      uwer_res_cm2=-0.0969704191047088
      uwer_res_cm1=-0.0056430956994203
      uwer_res_cm0=0.4003849386477017
      
      iop_cm2=0.0969704191046950
      iop_cm1=0.0056430956994063
      iop_cm0=-0.4069645466913195


c    output of Uwer's routine is the Born without couplings and before
c    the average on color and helicity. Corresponds to 4 Nc a_0 of formula A.2
      born=qgttq2(pp(0,1),pp(0,2),pp(0,5),pp(0,3),pp(0,4))

c     call MNR_ttbarj(kr_mad,vflav,born)
c     born=born*96d0/(4d0*pi*st_alpha)**3 ! remove couplings and averaging from MNR results to compare

  
      pole0=qgttq2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR2(huge)
      pole2=qgttq2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR2(0d0)
      call setDELTAIR1(huge)
      pole1=qgttq2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
      call setDELTAIR1(0d0)

      pole2=(pole2/huge)/born*(st_alpha/(4d0*pi))
      pole1=(pole1/huge)/born*(st_alpha /(4d0*pi))
      pole0=pole0/born*(st_alpha/(4d0*pi))
     
      print *,'as =',st_alpha
      print *,vflav
      print *,' CM2, MYPOLE2'
      print *, uwer_res_cm2,pole2
      print *,red,' RATIO MUST BE 1 =====> ',pole2 /uwer_res_cm2 ,reset
      print *,' CM1, MYPOLE1'
      print *, uwer_res_cm1,pole1
      print *,red,' RATIO MUST BE 1 =====> ',pole1 /uwer_res_cm1 ,reset
      print *,' CM0, MYPOLE0'
      print *, uwer_res_cm0,pole0
      print *,red,' RATIO MUST BE 1 =====> ',pole0 /uwer_res_cm0 ,reset

      vflav(1)=0
      vflav(2)=-1
      vflav(3)=6
      vflav(4)=-6
      vflav(5)=-1



      uwer_res_cm2=-0.0969704191046802
      uwer_res_cm1=0.0833362739128030
      uwer_res_cm0=0.5384721403213878
      
      
      iop_cm2=0.0969704191046950
      iop_cm1=-0.0833362739127882
      iop_cm0=-0.3392937280293060

      
c    output of Uwer's routine is the Born without couplings and before
c    the average on color and helicity. Corresponds to 4 Nc a_0 of formula A.2
      born=gqbttqb2(pp(0,1),pp(0,2),pp(0,5),pp(0,3),pp(0,4))
  
c     call MNR_ttbarj(kr_mad,vflav,born)
c     born=born*96d0/(4d0*pi*st_alpha)**3 ! remove couplings and averaging from MNR results to compare

      pole0=gqbttqb2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) 
     $        ,pp(0,4))
      call setDELTAIR2(huge)
      pole2=gqbttqb2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) 
     $        ,pp(0,4))
      call setDELTAIR2(0d0)
      call setDELTAIR1(huge)
      pole1=gqbttqb2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) 
     $        ,pp(0,4))
      call setDELTAIR1(0d0)

      pole2=(pole2/huge)/born*(st_alpha/(4d0*pi))
      pole1=(pole1/huge)/born*(st_alpha /(4d0*pi))
      pole0=pole0/born*(st_alpha/(4d0*pi))
     
      print *,'as =',st_alpha
      print *,vflav
      print *,' CM2, MYPOLE2'
      print *, uwer_res_cm2,pole2
      print *,red,' RATIO MUST BE 1 =====> ',pole2 /uwer_res_cm2 ,reset
      print *,' CM1, MYPOLE1'
      print *, uwer_res_cm1,pole1
      print *,red,' RATIO MUST BE 1 =====> ',pole1 /uwer_res_cm1 ,reset
      print *,' CM0, MYPOLE0'
      print *, uwer_res_cm0,pole0
      print *,red,' RATIO MUST BE 1 =====> ',pole0 /uwer_res_cm0 ,reset


      print *
      print *,red, "COMPARISON WITH ARXIV 0810.0452 FINISHED",reset
      print *,red, "PROGRAM STOPPED",reset
      call exit(1)

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine virtual_checks(pp,vflav)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_br.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      integer nlegs,j
      parameter (nlegs=nlegborn)
      integer vflav(nlegs),ileg,mu,jleg,nu      
      double precision virtual,born,bornjk(nlegs,nlegs),bmunu(0:3,0:3
     $     ,nlegs)
      double precision Ioperator_gg,gggtt2virt
      external Ioperator_gg,gggtt2virt
      double precision Ioperator_qq,qqttg2virt
      external Ioperator_qq,qqttg2virt
      double precision Ioperator_qg,qgttq2virt
      external Ioperator_qg,qgttq2virt
      double precision Ioperator_gqb,gqbttqb2virt
      external Ioperator_gqb,gqbttqb2virt
      double precision VCS,ICS,tmp,temp
      double precision c(-6:6),gamma(-6:6),gammap(-6:6)
      double precision pole2,pole1,pole0
      double precision pp(0:3,nlegs),amp2MNR
      double precision tiny
      parameter (tiny=1d-2)
      double precision gggtt2,gggtt2virtfin
      external gggtt2,gggtt2virtfin
      double precision qqttg2,qqttg2virtfin
      external qqttg2,qqttg2virtfin
      double precision qgttq2,qgttq2virtfin
      external qgttq2,qgttq2virtfin
      double precision gqbttqb2,gqbttqb2virtfin
      external gqbttqb2,gqbttqb2virtfin
      logical pwhg_isfinite
      external pwhg_isfinite
      logical ini
      integer count
      data ini/.true./
      save ini,count
      integer nchecks
      parameter(nchecks=1000)
      double precision huge
      parameter (huge=1d100)
c     colored output
      character*1 red(5)
      character*1 reset(4) 
      data red /' ','[','3','1','m'/ 
      data reset /' ','[','0','m'/ 
      double precision dotp
      external dotp
      double precision betahk
      save c,gamma,gammap
      logical verbose
      parameter (verbose=.false.)
      double precision intDUWscale2
c     useful function
      betahk(ileg,jleg)=sqrt(1- (kn_masses(ileg)*kn_masses(jleg)
     $     /dotp(pp(0,ileg),pp(0,jleg)))**2)

      red(1) = char(27)         ! Escape character (ASCII 27).
      reset(1) = char(27)

      if (ini) then
         print *,red,"###########################################",reset
         print *,red,"###########################################",reset
         print *,red,"   PERFORMING CHECKS ON VIRTUAL ROUTINES   ",reset
         print *,red,"           tiny = ",tiny,"                 ",reset
         print *,red,"###########################################",reset
         print *,red,"###########################################",reset
          do j=-6,6
            if(j.eq.0) then
               c(j)=ca
               gamma(j)=(11*ca-4*tf*st_nlight)/6
               gammap(j)=(67d0/9-2*pi**2/3)*ca-23d0/9*tf*st_nlight
            else
               c(j)=cf
               gamma(j)=3d0/2*cf
               gammap(j)=(13d0/2-2*pi**2/3)*cf
            endif
         enddo
         ini=.false.
         count=0
      endif
      
      count=count+1

c     Checks kinematics

      call checkmomzeronew(nlegs,pp)
      call checkmassnew(1,pp)
      call checkmassnew(2,pp)
      call checkmassnew(3,pp)
      call checkmassnew(4,pp)
      call checkmassnew(5,pp)
         
      born=0d0
      virtual=0d0
      vcs=0d0
      ics=0d0

c     Calculate  everything
      if(vflav(1).ne.0.and.vflav(1)+vflav(2).eq.0) then
         if(vflav(1).gt.0) then
c     q qbar
            born=qqttg2(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0,4))
     $              /36.0
            virtual=qqttg2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0
     $           ,4))/36.0
            vcs=qqttg2virt(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0 ,4))
     $              /36.0
            ics=Ioperator_qq(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0
     $           ,4))/36.0
         else
c     qbar q
            born=qqttg2(pp(0,2),pp(0,1),pp(0,5),pp(0,3) ,pp(0,4))
     $           /36.0
            virtual=qqttg2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/36.0
            vcs=qqttg2virt(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/36.0
            ics=Ioperator_qq(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/36.0
            endif
      endif
      if(vflav(1)*vflav(2).eq.0) then
         if(vflav(1).gt.0) then
c     a quark and a gluon
            born=qgttq2(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
            virtual=qgttq2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
            vcs=qgttq2virt(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96.0
            ics=Ioperator_qg(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96.0
         elseif(vflav(1).lt.0) then
c     an anti-quark and a gluon
            born=gqbttqb2(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
            virtual=gqbttqb2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
            vcs=gqbttqb2virt(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96.0
            ics=Ioperator_gqb(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96.0
         elseif(vflav(2).lt.0) then
c     a gluon and an anti-quark
            born=gqbttqb2(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
            virtual=gqbttqb2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
            vcs=gqbttqb2virt(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96.0
            ics=Ioperator_gqb(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $              ,pp(0,4))/96.0
         elseif(vflav(2).gt.0) then
c     a gluon and a quark
            born=qgttq2(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0  
            virtual=qgttq2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0  
            vcs=qgttq2virt(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
            ics=Ioperator_qg(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
         endif
      endif
c     two gluons
      if((vflav(1).eq.0).and.(vflav(2).eq.0))then
         born=gggtt2(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0
     $        ,4))/256.0
         virtual=gggtt2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0
     $        ,4))/256.0
         vcs=gggtt2virt(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0
     $           ,4))/256.0
         ics=Ioperator_gg(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0
     $           ,4))/256.0
      endif


c     Check Born against MNR code

      call MNR_ttbarj(pp,vflav,amp2MNR)

      born=born*((4*pi*st_alpha)**3) 
      
      tmp=abs(amp2MNR/born -1d0)
 
      if(.not.pwhg_isfinite(tmp)) then
         STOP "Error in checking born in virtual=NaN"
      endif

        
      if  (tmp.gt.tiny) then 
         write(*,*) red,vflav
     $        ,'\n VIRTUAL_CHECKS: BORN MUST BE EQUAL =====> ',amp2MNR
     $        ,born,'\n RATIO: ', amp2MNR/born,reset
         
         write(*,*) "P:"
         write(*,'(4f22.15,a)')((pp(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
         stop
      endif
         
      tmp=(abs(amp2MNR-born)/amp2MNR)

      if(.not.pwhg_isfinite(tmp)) then
         STOP "Error in checking born % diff in virtual=NaN"
      endif

         
      if  (tmp.gt.tiny) then 
         write(*,*) red,vflav
     $  ,'\n VIRTUAL_CHECKS: BORN % DIFFERENCE MUST BE 0 =====> '
     $   ,abs(amp2MNR-born)/amp2MNR,reset
         write(*,*) "P:"
         write(*,'(4f22.15,a)')((pp(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
         stop
      endif

c     Uncomment the line below to skip the
c     check that (V-Vfin)/(2*Ioperator) = 1 
      goto 888
      tmp=abs((vcs-virtual)/ics/2d0 -1d0)

      if(.not.pwhg_isfinite(tmp)) then
         STOP "Error in checking I operator subtraction in virtual=NaN"
      endif

      if ((tmp.gt.tiny).and.(virtual.ne.0d0).and.(vcs.ne.0d0)) then 
         write(*,*) red,vflav
     $        ,'\n VIRTUAL_CHECKS: (V-Vfin)/2I MUST BE 1 ===> ',(vcs
     $        -virtual)/ics/2d0, "() ",vcs,virtual,2*ics,reset
         write(*,*) "P:"
         write(*,'(4f22.15,a)')((pp(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
         do j=1,10
         print *
         print *,(qqttg2virt(pp(0,2),pp(0,1),pp(0,5),pp(0,3) ,pp(0 ,4))
     $        -qqttg2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3) ,pp(0 ,4)))
     $        /Ioperator_qq(pp(0,2),pp(0,1),pp(0,5),pp(0 ,3) ,pp(0 ,4))
     $        /2.0,qqttg2virt(pp(0,2),pp(0,1),pp(0,5),pp(0,3) ,pp(0 ,4))
     $        ,qqttg2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3) ,pp(0 ,4))
     $        ,Ioperator_qq(pp(0,2),pp(0,1),pp(0,5),pp(0 ,3) ,pp(0 ,4))
         enddo
         stop
      endif
 888  continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c CHECK POLE STRUCTURE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      amp2MNR=born
c     evaluate color and spin correlated born
      call compborn(pp,vflav,born,bornjk,bmunu)

c re-check born against madgraph routines
      tmp=abs(amp2MNR/born -1d0)
 
      if(.not.pwhg_isfinite(tmp)) then
         STOP "Error in checking MAD born in virtual=NaN"
      endif

        
      if  (tmp.gt.tiny) then 
         write(*,*) red,vflav
     $        ,'\n VIRTUAL_CHECKS MAD: BORN MUST BE EQUAL =====> '
     $        ,amp2MNR ,born,'\n RATIO: ', amp2MNR/born,reset
         
         write(*,*) "P:"
         write(*,'(4f22.15,a)')((pp(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
         stop
      endif
         
         
      tmp=(abs(amp2MNR-born)/amp2MNR)

      if(.not.pwhg_isfinite(tmp)) then
         STOP "Error in checking MAD born % diff in virtual=NaN"
      endif

         
      if  (tmp.gt.tiny) then 
         write(*,*) red,vflav
     $  ,'\n VIRTUAL_CHECKS MAD: BORN % DIFFERENCE MUST BE 0 =====> '
     $   ,abs(amp2MNR-born)/amp2MNR,reset
         write(*,*) "P:"
         write(*,'(4f22.15,a)')((pp(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
         stop
      endif



      call setDELTAIR2(huge)

      if(vflav(1).ne.0.and.vflav(1)+vflav(2).eq.0) then
         if(vflav(1).gt.0) then
c     q qbar
            virtual=qqttg2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3),pp(0
     $           ,4))/36.0
         else
c     qbar q
            virtual=qqttg2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/36.0
         endif
      endif
      if(vflav(1)*vflav(2).eq.0) then
         if(vflav(1).gt.0) then
c     a quark and a gluon
            virtual=qgttq2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
          elseif(vflav(1).lt.0) then
c     an anti-quark and a gluon
             virtual=gqbttqb2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $            ,pp(0,4))/96.0
          elseif(vflav(2).lt.0) then
c     a gluon and an anti-quark
             virtual=gqbttqb2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $            ,pp(0,4))/96.0
          elseif(vflav(2).gt.0) then
c     a gluon and a quark
             virtual=qgttq2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $            ,pp(0,4))/96.0  
         endif
      endif
c     two gluons
      if((vflav(1).eq.0).and.(vflav(2).eq.0))then
         virtual=gggtt2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3),pp(0
     $        ,4))/256.0
      endif
      

      pole2= virtual/huge
      
      temp=0
      do ileg=1,nlegborn
         if (kn_masses(ileg).eq.0) then
            temp=temp-c(vflav(ileg))
         endif
      enddo

      pole2=pole2*(4d0*pi*st_alpha)**3/2d0
      tmp= pole2/temp/born
      
      if(.not.pwhg_isfinite(tmp)) then
         STOP "Error in checking double poles virtual=NaN"
      endif

      if (abs(tmp -1d0) .gt.tiny) then 
         write(*,*) red,vflav
     $        ,'\n VIRTUAL  EPS2 POLE CHECK: MUST BE 1 ===> '
     $        ,tmp," ()", pole0, virtual,reset
ccccccccccccccccccccccccccccccccccccccccccccc
         write(*,*) "P:"
         write(*,'(4f22.15,a)')((pp(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
         do j=1,10
         print *
         print *,gggtt2virt(pp(0,1),pp(0,2),pp(0,5),pp(0,3) ,pp(0 ,4))
         enddo
ccccccccccccccccccccccccccccccccccccccccccccc
         stop
      else 
         if(verbose) then
            write(*,*) "DOUBLE POLE CHECK: MUST BE 1 ==> ",tmp
         endif
      endif

c     Check single pole 

      call setDELTAIR2(0d0)
      call setDELTAIR1(huge)

      if(vflav(1).ne.0.and.vflav(1)+vflav(2).eq.0) then
         if(vflav(1).gt.0) then
c     q qbar
            virtual=qqttg2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3),pp(0
     $           ,4))/36.0
         else
c     qbar q
            virtual=qqttg2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/36.0
         endif
      endif
      if(vflav(1)*vflav(2).eq.0) then
         if(vflav(1).gt.0) then
c     a quark and a gluon
            virtual=qgttq2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
         elseif(vflav(1).lt.0) then
c     an anti-quark and a gluon
            virtual=gqbttqb2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
         elseif(vflav(2).lt.0) then
c     a gluon and an anti-quark
            virtual=gqbttqb2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0
         elseif(vflav(2).gt.0) then
c     a gluon and a quark
            virtual=qgttq2virtfin(pp(0,2),pp(0,1),pp(0,5),pp(0,3)
     $           ,pp(0,4))/96.0  
         endif
      endif
c     two gluons
      if((vflav(1).eq.0).and.(vflav(2).eq.0))then
         virtual=gggtt2virtfin(pp(0,1),pp(0,2),pp(0,5),pp(0,3),pp(0
     $        ,4))/256.0
      endif
      
      
c evaluate single pole residue      
c put here the same choiche for internal scale 
c as done in  SetVirtualScales inside wrap_virtual.cpp
      intDUWscale2= st_muren2
c
      temp=0
      do ileg=1,nlegborn
         do jleg=ileg+1,nlegborn
c     both particles are colored
            if (abs(vflav(ileg)).le.6.and.
     #              abs(vflav(jleg)).le.6) then
c     massless-massless case
               if (kn_masses(ileg).eq.0.and.kn_masses(jleg).eq.0) then
                  temp=temp+2*(log(2d0*dotp(pp(0,ileg),pp(0,jleg))
     $                 /intDUWscale2)* bornjk(ileg,jleg))
               endif
c     massless-massive case
               if (kn_masses(ileg).eq.0.and.kn_masses(jleg).gt.0) then
                  temp=temp+2*((log(2d0*dotp(pp(0,ileg),pp(0,jleg))
     $                 /intDUWscale2)-0.5*log((kn_masses(jleg)**2)
     $                 /intDUWscale2)) * bornjk(ileg,jleg))
               endif
c     massive-massless case
               if (kn_masses(ileg).gt.0.and.kn_masses(jleg).eq.0) then
                  temp=temp+2*((log(2d0*dotp(pp(0,ileg),pp(0,jleg))
     $                 /intDUWscale2)-0.5*log((kn_masses(ileg)**2)
     $                 /intDUWscale2)) * bornjk(jleg,ileg))
               endif
c     massive-massive case
               if (kn_masses(ileg).gt.0.and.kn_masses(jleg).gt.0) then
                  temp=temp+log((1d0+betahk(ileg,jleg))/((1-betahk(ileg
     $                 ,jleg))))/betahk(ileg,jleg)* bornjk(ileg,jleg)
               endif
            endif
         enddo
         if (kn_masses(ileg).ne.0) then
            temp=temp-c(vflav(ileg))*born
         endif
         if (kn_masses(ileg).eq.0) then
            temp=temp-gamma(vflav(ileg))*born
        endif

      enddo



      pole1=(virtual/huge)
      tmp=(((4*pi*st_alpha)**3 * pole1/2d0)/ temp)

      if(.not.pwhg_isfinite(tmp)) then
         STOP "Error in checking single poles virtual=NaN"
      endif


      if (abs(tmp -1d0) .gt.tiny) then 
         write(*,*) red,vflav
     $        ,'\n VIRTUAL  EPS1 POLE CHECK: MUST BE 1 ===> '
     $        ,tmp,reset
         write(*,*) "P:"
         write(*,'(4f22.15,a)')((pp(mu,ileg),mu=0,3),"\n",ileg=1 ,nlegs)
         stop
      else 
         if(verbose) then
            write(*,*) "SINGLE POLE CHECK: MUST BE 1 ==> ",tmp
         endif
      endif
           
      call setDELTAIR1(0d0)



      if(count.ge.nchecks) then
      print *,red,"###########################################",reset
      print *,red, "CHECKS PERFORMED FOR ",nchecks
     $     ," RANDOM PHASE SPACE POINTS",reset
      print *,red, " CHECKS PASSED!   PROGRAM STOPPED",reset
      print *,red,"###########################################",reset
             
      call exit(1)
      endif


      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  

     
