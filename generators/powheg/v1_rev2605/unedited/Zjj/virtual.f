c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !: when linking to an external OLP,
c     make sure that code lines just after the comments 
c     starting with !: should be uncommented.
cccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      double precision p(0:3,nlegborn)
      integer vflav(nlegborn)
      double precision virtual
      logical ini
      save ini
      data ini/.true./
      double precision powheginput
      external powheginput
      double precision born,dummy(0:3,0:3,6),bornjk(6,6)
      integer fakevirt
      save fakevirt 

cccccccccccccccccccccccccccccccccccccccccccccccc
      logical use_OLP_Interface
      common /colp/use_OLP_Interface
      data use_OLP_Interface/.false./
ccccccccccccccccccccccccccccccccccccccccccccccc

      if (ini) then
         if(powheginput("#use-OLP-interface").eq.1) then
            use_OLP_Interface=.true.
         endif

         fakevirt=powheginput("#fakevirt")
         if(fakevirt.eq.1) then
            write(*,*) 'WARNING: Using fakevirt !'
         endif

         if(use_OLP_interface) then
c     initialize OLP (a-la Binoth-Les-Houches)
            call virtual_initialize_OLP
         else
c     initialize virtuals (generic)
            call virtual_initialize
         endif
         ini=.false.
      endif

      if(fakevirt.eq.1) then
         call compborn(p,vflav,born,dummy,bornjk)
         virtual=0.2*born
      else
         if(use_OLP_interface) then
            call virtual_OLP(p,vflav,virtual)
         else
            call virtual_evaluate(p,vflav,virtual)
         endif
      endif
      end


cccccccccccccccccccccccccccccccccc
c     BLACKHAT
cccccccccccccccccccccccccccccccccc
      subroutine virtual_initialize_OLP
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      logical debug
      parameter (debug=.true.)
c     from powheg born code to the OLP code
      integer pwhgcode2OLPcode(maxprocborn)
      common/cpwhgcode2OLPcode/pwhgcode2OLPcode


      character*11 filename
      integer iun
      integer ios,k,l
      character *5 status
      integer code,foundbproc

      character * 100 line,line0

      write(*,*) 
      write(*,*) ' Checking contract file for external OLP '
      write(*,*) 
      filename="contract.lh"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
      call OLP_start(filename//CHAR(0),ios)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) ' OLP_init called '
      if (ios.ne.1) then
         write(*,*) ' Error: OLP cannot handle contract file'
         call exit(1)
      endif
      
c     open the negotiation file, read it and
c     fill the array pwhgbcode2OLPcode properly
      call newunit(iun)
      open(unit=iun,file='contract.lh',status='old',iostat=ios)
      if(ios.ne.0) then
c         write(*,*) 'cannot open contract.lh'
         write(*,*) ' WARNING: contract.lh file not found.'
         write(*,*) ' Run the OLP appropriate command to create it'
         write(*,*) ' (for BH, the command is LH_reader)'
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
c         if(debug) write(*,*) line(11:12)
         if(line(11:12).ne.'->') then !:
c     this means that current line is not a line with a subprocess
            if(debug) write(*,*) 'Found a line without subprocess'
            goto 111
         endif
         if(debug) write(*,*) 'Found a line with a subprocess string'
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
               if(status.eq."1") then !:
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

      end


      subroutine virtual_OLP(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      double precision p(0:3,nlegborn)
      integer vflav(nlegborn)
      double precision virtual

cccccccccccccccccccccccccccccc
      integer OLP_code

c     from born flavour structure to jborn (for simplicity, here
c     use only MASSLESS colored particles, sorted as:
c     id_plus,id_minus,id_final1,id_final2)
c     Particle 3 and 4 are the 2 leptons
      integer bornflst2pwhgcode(-5:5,-5:5,-5:5,-5:5)
      common/cbornflst2pwhgcode/bornflst2pwhgcode

c     from powheg born code to the OLP code
      integer pwhgcode2OLPcode(maxprocborn)
      common/cpwhgcode2OLPcode/pwhgcode2OLPcode
ccccccccccccccccccccccccccccccc

      double precision p_olp(0:4,nlegborn),couplings_olp(2)
      double precision ren_scale2,ren_scale
      double precision virt_wgts(4),virtual_0,virtual_1,virtual_2

      double precision born,bornjk(6,6),bornmunu(0:3,0:3,6)


      integer mu,ileg,jleg
      double precision dotp
      external dotp

cccccccccccccccccccccccccccccccccccccccccc
      double precision c(-6:6),gamma(-6:6),gammap(-6:6)
      double precision double_pole,single_pole,logtmp,OLP_tree
      logical polescheck
      double precision tiny
      parameter (polescheck=.false.,tiny=1d-4)
cccccccccccccccccccccccccccccccccccccccccc
      logical error
      double precision p_ref(0:3,nlegborn),born_ref,virtual_ref,bratio,
     $     vratio
      double precision pwhg_alphas,powheginput
      external pwhg_alphas,powheginput
      double precision tiny_ref
      parameter (tiny_ref=1d-7)
cccccccccccccccccccccccccccccccccccccccccc

      OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $     ,vflav(5),vflav(6)))
      do ileg=1,nlegborn
         do mu=0,3
            p_olp(mu,ileg)=p(mu,ileg)
         enddo
         p_olp(4,ileg)=0d0
      enddo
      ren_scale2=st_muren2
      ren_scale=sqrt(ren_scale2)
      couplings_olp(1)=0d0
      couplings_olp(2)=0d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
c      call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale2,0d0,0d0,virt_wgts)
      call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale,
     $     couplings_olp,virt_wgts)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      virtual_2=virt_wgts(1)
      virtual_1=virt_wgts(2)
      virtual_0=virt_wgts(3)

      call setborn(p,vflav,born,bornjk,bornmunu)

      virtual=virtual_0*born


ccccccccccccccccccccccccccccccccccccccc     
      if(polescheck) then
c     CHECKS

c     from 2.100 of FNO2007
         do jleg=-6,6
            if(jleg.eq.0) then
               c(jleg)=ca
               gamma(jleg)=(11*ca-4*tf*st_nlight)/6
               gammap(jleg)=(67d0/9-2*pi**2/3)*ca-23d0/9*tf*st_nlight
            else
               c(jleg)=cf
               gamma(jleg)=3d0/2*cf
               gammap(jleg)=(13d0/2-2*pi**2/3)*cf
            endif
         enddo

c     double pole
         double_pole=0d0
         do ileg=1,nlegborn
            if(abs(vflav(ileg)).le.6) then
               double_pole=double_pole - c(vflav(ileg))
            endif
         enddo
         if(abs(double_pole/virtual_2-1d0).gt.tiny) then
            write(*,*) 'vflav is ',vflav
            write(*,*) "double pole: th/OLP | th | OLP= ",double_pole/virtual_2,
     $           " | ",double_pole," | ",virtual_2
         endif

c     single_pole
         single_pole=0d0
         do ileg=1,nlegborn
            do jleg=ileg+1,nlegborn
               if((abs(vflav(ileg)).le.6).and.(abs(vflav(jleg)).le.6)) then
                  logtmp=log(2d0*dotp(p(0,ileg),p(0,jleg))/st_muren2)
                  single_pole=single_pole + 2d0*bornjk(ileg,jleg)*logtmp
               endif
            enddo
            if(abs(vflav(ileg)).le.6) single_pole=single_pole - gamma(vflav(ileg))*born
         enddo        
         single_pole=single_pole/born
         if(abs(single_pole/virtual_1-1d0).gt.tiny) then
            write(*,*) '-------'
            write(*,*) 'vflav is ',vflav
            write(*,*) 'momenta are ',p
            write(*,*) "single pole: th/OLP | th | OLP= ",single_pole/virtual_1,
     $           " | ",single_pole," | ",virtual_1
         endif

c     tree-level (not always returned by the OLP, BH doesn't)
c$$$         OLP_tree=virt_wgts(4)*(ph_unit_e**4 * (4*pi*st_alpha)**2)
c$$$         write(*,*) "tree level: POW/OLP | th | OLP= ",born/OLP_tree,
c$$$     $        " | ",born," | ",OLP_tree

      endif
ccccccccccccccccccccccccccccccccccccccc
      if(powheginput("#check_ref_amp").eq.1) then

         write(*,*) '======================================'
         write(*,*) ' checking born and virtual amplitudes '
         write(*,*) '======================================'

         p_ref(0,1)= 138.4784456174d0
         p_ref(1,1)= 0d0 
         p_ref(2,1)= 0d0 
         p_ref(3,1)= 138.4784456174d0

         p_ref(0,2)= 138.4784456174d0
         p_ref(1,2)= 0d0
         p_ref(2,2)= 0d0
         p_ref(3,2)=-138.4784456174d0

         p_ref(0,3)= 39.8695145736d0
         p_ref(1,3)=-30.2118421708d0
         p_ref(2,3)= 7.3821554172d0
         p_ref(3,3)=-24.9464740270d0

         p_ref(0,4)= 61.7949044527d0
         p_ref(1,4)= 60.7995170788d0
         p_ref(2,4)=-7.3821554172d0
         p_ref(3,4)=-8.2178294395d0

         p_ref(0,5)= 104.9831601822d0
         p_ref(1,5)=-38.2721470356d0
         p_ref(2,5)= 44.5354030781d0
         p_ref(3,5)= 87.0247353101d0

         p_ref(0,6)= 70.3093120261d0
         p_ref(1,6)= 7.6844721276d0  
         p_ref(2,6)=-44.5354030781d0
         p_ref(3,6)=-53.8604318436d0

         do ileg=1,nlegborn
            do mu=0,3
               p_olp(mu,ileg)=p_ref(mu,ileg)
            enddo
            p_olp(4,ileg)=0d0
         enddo

         ren_scale=91.1876d0
         st_muren2 = ren_scale**2
         st_alpha  = pwhg_alphas(st_muren2,st_lambda5MSB,st_nlight)
         couplings_olp(1)=0d0
         couplings_olp(2)=0d0

         error=.false.
         vflav(3)=11
         vflav(4)=-11
cccccccccccccccccccccc
c     u ubar g g
         born_ref=    1.0489864146 E-005 
         virtual_ref= 1.7673954550 E-004
         vflav(1)=2
         vflav(2)=-2
         vflav(5)=0
         vflav(6)=0
         OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $        ,vflav(5),vflav(6)))
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
         call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale,
     $        couplings_olp,virt_wgts)
         virtual_0=virt_wgts(3)
         call setborn(p_ref,vflav,born,bornjk,bornmunu)
         virtual=virtual_0*born
         bratio=born/born_ref
         vratio=virtual/virtual_ref
         write(*,*) vflav(1),vflav(2),vflav(5),vflav(6),bratio,vratio
         if((dabs(bratio-1).gt.tiny_ref).or.
     $        (dabs(vratio-1).gt.tiny_ref)) error=.true.
cccccccccccccccccccccc
c     d dbar g g
         born_ref=    1.3620189441 E-005
         virtual_ref= 2.4684280244 E-004
         vflav(1)=1
         vflav(2)=-1
         vflav(5)=0
         vflav(6)=0
         OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $        ,vflav(5),vflav(6)))
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
         call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale,
     $        couplings_olp,virt_wgts)
         virtual_0=virt_wgts(3)
         call setborn(p_ref,vflav,born,bornjk,bornmunu)
         virtual=virtual_0*born
         bratio=born/born_ref
         vratio=virtual/virtual_ref
         write(*,*) vflav(1),vflav(2),vflav(5),vflav(6),bratio,vratio
         if((dabs(bratio-1).gt.tiny_ref).or.
     $        (dabs(vratio-1).gt.tiny_ref)) error=.true.
cccccccccccccccccccccc
c     u ubar u ubar
         born_ref=    6.9320177545 E-006 
         virtual_ref= 2.1474734072 E-004
         vflav(1)=2
         vflav(2)=-2
         vflav(5)=2
         vflav(6)=-2
         OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $        ,vflav(5),vflav(6)))
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
         call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale,
     $        couplings_olp,virt_wgts)
         virtual_0=virt_wgts(3)
         call setborn(p_ref,vflav,born,bornjk,bornmunu)
         virtual=virtual_0*born
         bratio=born/born_ref
         vratio=virtual/virtual_ref
         write(*,*) vflav(1),vflav(2),vflav(5),vflav(6),bratio,vratio
         if((dabs(bratio-1).gt.tiny_ref).or.
     $        (dabs(vratio-1).gt.tiny_ref)) error=.true.
cccccccccccccccccccccc
c     u ubar c cbar
         born_ref=    1.8468108910 E-007
         virtual_ref= 1.8007297663 E-006
         vflav(1)=2
         vflav(2)=-2
         vflav(5)=4
         vflav(6)=-4
         OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $        ,vflav(5),vflav(6)))
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
         call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale,
     $        couplings_olp,virt_wgts)
         virtual_0=virt_wgts(3)
         call setborn(p_ref,vflav,born,bornjk,bornmunu)
         virtual=virtual_0*born
         bratio=born/born_ref
         vratio=virtual/virtual_ref
         write(*,*) vflav(1),vflav(2),vflav(5),vflav(6),bratio,vratio
         if((dabs(bratio-1).gt.tiny_ref).or.
     $        (dabs(vratio-1).gt.tiny_ref)) error=.true.
cccccccccccccccccccccc
c     d dbar d dbar
         born_ref=    1.0023268662 E-005
         virtual_ref= 3.1647659850 E-004
         vflav(1)=1
         vflav(2)=-1
         vflav(5)=1
         vflav(6)=-1
         OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $        ,vflav(5),vflav(6)))
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
         call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale,
     $        couplings_olp,virt_wgts)
         virtual_0=virt_wgts(3)
         call setborn(p_ref,vflav,born,bornjk,bornmunu)
         virtual=virtual_0*born
         bratio=born/born_ref
         vratio=virtual/virtual_ref
         write(*,*) vflav(1),vflav(2),vflav(5),vflav(6),bratio,vratio
         if((dabs(bratio-1).gt.tiny_ref).or.
     $        (dabs(vratio-1).gt.tiny_ref)) error=.true.
cccccccccccccccccccccc
c     d dbar s sbar
         born_ref=    2.3673784898 E-007
         virtual_ref= 2.4179197657 E-006
         vflav(1)=1
         vflav(2)=-1
         vflav(5)=3
         vflav(6)=-3
         OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $        ,vflav(5),vflav(6)))
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
         call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale,
     $        couplings_olp,virt_wgts)
         virtual_0=virt_wgts(3)
         call setborn(p_ref,vflav,born,bornjk,bornmunu)
         virtual=virtual_0*born
         bratio=born/born_ref
         vratio=virtual/virtual_ref
         write(*,*) vflav(1),vflav(2),vflav(5),vflav(6),bratio,vratio
         if((dabs(bratio-1).gt.tiny_ref).or.
     $        (dabs(vratio-1).gt.tiny_ref)) error=.true.
cccccccccccccccccccccc
c     u ubar d dbar
         born_ref=    6.6505600059 E-007
         virtual_ref= 1.3323314116 E-006
         vflav(1)=2
         vflav(2)=-2
         vflav(5)=1
         vflav(6)=-1
         OLP_code=pwhgcode2OLPcode(bornflst2pwhgcode(vflav(1),vflav(2)
     $        ,vflav(5),vflav(6)))
c     !: Uncomment this line when linking an external code to evaluate virtual
c     !: corrections using the Binoth-LH interface
         call OLP_EvalSubprocess(OLP_code,p_olp,ren_scale,
     $        couplings_olp,virt_wgts)
         virtual_0=virt_wgts(3)
         call setborn(p_ref,vflav,born,bornjk,bornmunu)
         virtual=virtual_0*born
         bratio=born/born_ref
         vratio=virtual/virtual_ref
         write(*,*) vflav(1),vflav(2),vflav(5),vflav(6),bratio,vratio
         if((dabs(bratio-1).gt.tiny_ref).or.
     $        (dabs(vratio-1).gt.tiny_ref)) error=.true.
cccccccccccccccccccccccc

         write(*,*) '======================================'
         write(*,*) ' born and virtual amplitudes checked '
         write(*,*) '======================================'

         if(error) then
            write(*,*) '   **** ERROR **** '
            write(*,*) '   Cannot reproduce born and virtual amplitudes'
            write(*,*) '   Likely something wrong in the B-LH interface'
         else
            write(*,*) '   **** TEST OK **** '
            write(*,*) '   Remove check_ref_amp flag from input card'
            write(*,*) '   and run with your own parameters'
         endif
         call exit(-1)

      endif




      end

      





ccccccccccccccccccccccccccccccccccccccccc
c     Generic virtual interface
ccccccccccccccccccccccccccccccccccccccccc
      subroutine virtual_initialize
      implicit none
      write(*,*) 'Error in virtual_initialize (virtual.f):'
      write(*,*) 'No virtual corrections present'
      write(*,*) 'An external OLP should be linked'
      call exit(-1)
      return
      end

      subroutine virtual_evaluate(p,vflav,virtual)
      implicit none
      double precision p(0:3,6)
      integer vflav(6)
      double precision virtual
      double precision born,dummy(0:3,0:3,6),bornjk(6,6)
      call compborn(p,vflav,born,dummy,bornjk)
      virtual=0.2*born
      write(*,*) 'Error in virtual_evaluate:'
      write(*,*) 'No virtual corrections present'
      call exit(-1)
      return
      end

