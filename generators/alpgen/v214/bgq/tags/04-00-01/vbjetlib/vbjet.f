c-------------------------------------------------------------------
      subroutine alsprc
c     assigns the hard process code
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      ihrd=5
      end

c-------------------------------------------------------------------
      subroutine alhset
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'vbjet.inc'
      integer i,il,jprocmax,jproc,ngrid,jgrid,maxsubproc
      real*8 hwidthlim
      parameter(hwidthlim= 10.d0)
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c     ngrid is the total number of grids allowd for P.S. variables.
c     jgrid(jproc) labels the grid associated to a given jproc.
c     jgrid= 1   -->   q qbar and qbar q (no quarks in final state)
c     jgrid= 2   -->   g q    and g qbar
c     jgrid= 3   -->   q g    and qbar g
c     jgrid= 4   -->   g g
c     jgrid= 5   -->   qq->qq
      data jgrid/100*0/
      character*2 numvb(8)
      data numvb/'  ',' 2',' 3',' 4',' 5',' 6',' 7',' 8'/
c     
      character*50 tmpstr
      integer j1,ntmp,res
      character*1 zd(maxpar-2),zr(5)
      common/zdec/zd
      data zr/'n','e','q','b','i'/
c     parameters for the gauge invariance prescription:
      idecay='y'
      winsize  = 2.d0/pi
      resonance= 'n'
      wmode    = 'yn'  
*     
      if (imode.le.1.and.nz+nw+nh+nph.eq.0) then
        print*,'nz+nw+nh+nph=  0 is not allowed!'
        stop
      endif
c     
c     process input parameters
      nprtns= njets
*     
c     initialise alpha
      if (apar(6).gt.10.d0) then
        print*,'Gamma(H) >',hwidthlim,' is not allowed!'
        stop
      endif
c     ngrid is the total number of grids allowd for P.S. variables.
c     jgrid(jproc) labels the grid associated to a given jproc.
*     
      if (mod(nw,2).eq.0) then
*     
c     jgrid= 1   -->   q qbar and qbar q (no quarks in final state)
c     jgrid= 2   -->   g q    and g qbar
c     jgrid= 3   -->   q g    and qbar g
c     jgrid= 4   -->   g g
c     jgrid= 5   -->   q q and qbar qbar
c     jgrid= 6   -->   q qbar and qbar q (+  quarks in final state)
*     
c     poi
        if(njets.eq.0) then     
          jprocmax= 1           
          ngrid   = 1
        elseif(njets.eq.1) then 
          jprocmax= 3 
          ngrid   = 3
        elseif(njets.eq.2) then 
          jprocmax = 16
          ngrid= 5
        elseif(njets.ge.3) then 
          jprocmax = 26          
          ngrid    =  5
        endif
        jgrid(1) =  1                        
        jgrid(2) =  2                        
        jgrid(3) =  3                        
        jgrid(4) =  4                        
        jgrid(5) =  5                        
        jgrid(6) =  5                        
        jgrid(7) =  5                        
        jgrid(8) =  5                        
        jgrid(9) =  5                        
        jgrid(10)=  5                        
        jgrid(11)=  5                        
        jgrid(12)=  5                        
        jgrid(13)=  5                        
        jgrid(14)=  5                        
        jgrid(15)=  5                        
        jgrid(16)=  5                        
        jgrid(17)=  2   
        jgrid(18)=  3   
        jgrid(19)=  2   
        jgrid(20)=  3   
        jgrid(21)=  2   
        jgrid(22)=  3   
        jgrid(23)=  2   
        jgrid(24)=  3   
        jgrid(25)=  2   
        jgrid(26)=  3
      else
        if(njets.eq.0) then     
          jprocmax= 1           
          ngrid   = 1
        elseif(njets.eq.1) then 
          jprocmax= 3           
          ngrid   = 3
        elseif(njets.eq.2) then 
          jprocmax=15          
          ngrid   = 4
        elseif(njets.ge.3) then 
          jprocmax=21          
          ngrid   = 5
        endif                   
        jgrid(1) = 1                         
        jgrid(2) = 2                         
        jgrid(3) = 3                         
        jgrid(4) = 4                         
        jgrid(5) = 4                         
        jgrid(6) = 4                         
        jgrid(7) = 4                         
        jgrid(8) = 4                         
        jgrid(9) = 4                         
        jgrid(10)= 4                         
        jgrid(11)= 4                         
        jgrid(12)= 4                         
        jgrid(13)= 4                         
        jgrid(14)= 4                         
        jgrid(15)= 4                         
        jgrid(16)= 5                         
        jgrid(17)= 5    
        jgrid(18)= 5    
        jgrid(19)= 5    
        jgrid(20)= 5    
        jgrid(21)= 5    
      endif
*     
c     protection:
*     
      if (jprocmax.gt.100) then
        write (*,*) 'increase maxsubproc!'
        stop
      endif
*     
      nvb  =nw+nz
      npart=2+nvb+nh+njets+nph
*     
c     masses to feed the generator
c     The ordering of the final state particles is
c     
c     jets, photons, Z's, W+ W- ... W+-      
c     
*     
      do i= 1,njets+2+nph
        p(5,i)=0
      enddo
      do i= 1,nz
        p(5,njets+2+i+nph)=mz
      enddo
      do i= 1,nw
        p(5,njets+2+nz+i+nph)=mw
      enddo
      do i= 1,nh
        p(5,njets+2+nz+nw+i+nph)=mh
      enddo
      do i=1,npart-2
        xm(i)=p(5,i+2)
      enddo
c     
c     run parameters:
      ntmp=zfstate
      j1= 0
      do while(ntmp.ne.0) 
        j1= j1+1
        res= mod(ntmp,10)
        zd(j1)= zr(res)
        ntmp= ntmp/10
      enddo
      if(j1.ne.nz) then
        write(*,*) 'wrong Z decay string, stop'
        stop
      endif
      i=1
      if(nw.gt.0) then
        tmpstr(i:i+4)=numvb(nw)//' W'
        i=i+4
      endif 
      if(nz.gt.0) then
        if(i.gt.1) then
          tmpstr(i:i+2)=' +'
          i=i+2
        endif 
        tmpstr(i:i+4)=numvb(nz)//' Z'
        i=i+4
      endif 
      if(nh.gt.0) then
        if(i.gt.1) then
          tmpstr(i:i+2)=' +'
          i=i+2
        endif 
        tmpstr(i:i+4)=numvb(nh)//' H'
        i=i+4
      endif 
c pippo
      if(nph.gt.0) then
        if(i.gt.1) then
          tmpstr(i:i+2)=' +'
          i=i+2
        endif 
        tmpstr(i:i+4)=numvb(nph)//'ph'
        i=i+4
      endif 
c pippo
      if(njets.gt.0) then
        if(i.gt.1) then
          tmpstr(i:i+2)=' +'
          i=i+2
        endif 
        if(njets.eq.1) then
          tmpstr(i:i+7)=numvb(njets)//' jet'
        else
          tmpstr(i:i+7)=numvb(njets)//' jets'
        endif
        i=i+7
      endif 
      do il=i,50
        tmpstr(il:il)=' '
      enddo
      write(niosta,*) tmpstr
      do il=1,i
        tmpstr(il:il)='='
      enddo
c     write(niosta,*) 'W -> ell nu,  Z -> ell+ ell-'
      write(niosta,*) tmpstr
      write(niosta,*)
     $     'Generation cuts for the partonic event sample:' 
      write(niosta,*) '     Light jets:'
      if (irapgap.eq.1) then
        write(niosta,*) 'irapgap=',irapgap
        write(niosta,*) 'etagap=',etagap
        if (njets.eq.3) write(niosta,*) 'ptcen=',ptcen
      endif
      write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $     ,' dR(j-j)>',drjmin 
c pippo
      if (nph.ge.1) then
        write(niosta,*) '     Photons:'
        write(niosta,*) 'ptmin(ph)=',ptphmin,' |etamax(ph)|=',etaphmax
        write(niosta,*) 'dR(ph-j)>',drphjmin
      endif
c pippo
      end

c-------------------------------------------------------------------
      subroutine selflav(jproc,xlum,afl)
c     
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4
c     
c     The quantity W stands for nvb= nw(odd) + nz Vector Bosons + nh
C     Higgs
c     The quantity Z stands for nvb= nw(even)+ nz Vector Bosons + nh
C     Higgs
c     
c     The ordering of the final state particles is 
c     
c     jets, Z's, W+ W- ... W+-, H's      
c     
c     jproc: (nw= odd) (Cabibbo structure taken into account)
c     
c     (jproc= 1     --> njets >= 0)
c     (jproc= 2:3   --> njets >= 1)
c     (jproc= 4:15  --> njets >= 2)
c     (jproc= 16:21 --> njets >= 3)
c     
ccc   To be completed: (jproc= 22:23 --> njets >= 4)
c     
c     1  q qbar' -> W  ,   qbar q' -> W                           
c     2  g q -> q' W    and    g qbar -> q'bar W                 
c     3  q g -> q' W    and  qbar g -> q'bar W                   
c     4  gg -> q qbar' W                                         
c     5  q qbar' -> W q'' qbar'' ,   qbar q' -> W q'' qbar''     
c     6  q(bar) q(bar)'' -> W q'(bar) q(bar)''                   
c     7  q''(bar) q(bar) -> W q'(bar) q(bar)''                   
c     8  q'' q''bar -> W q q'bar                                 
c     9  q qbar' -> W q qbar ,   qbar q' -> W qbar q             
c     10 qbar' q -> W q qbar ,   q' qbar -> W qbar q             
c     11 q qbar -> W q qbar'  + c.c.                             
c     12 sq qbar -> W q' qbar                                     
c     13 q q -> W q q'                                           
c     14 q q' -> W q q                                           
c     15 q q' -> W q' q'                                         
c     16 q g -> W q' q'' qbar'' ,   qbar g -> W qbar' q'' qbar'' 
c     17 g q -> W q' q'' qbar'' ,   g qbar -> W qbar' q'' qbar'' 
c     18 q g -> W q q qbar'                                      
c     19 q g -> W q' q qbar                                      
c     20 g q -> W q q qbar'                                      
c     21 g q -> W q' q qbar                                      
c     
ccc   To be completed:  22 gg -> q qbar' q'' qbar'' W
ccc   23 gg -> q qbar q qbar' W  + c.c.
c     
c     jproc: (nw= even)  (Cos_theta_cab= 1 approximation)
c     
c     (jproc= 1     --> njets >= 0)
c     (jproc= 2:3   --> njets >= 1)
c     (jproc= 4:16  --> njets >= 2)
c     (jproc= 17:26 --> njets >= 3)
c     (jproc= 17:26 --> njets >= 4)
c     
c     1  q qbar  -> Z                and complex conjugate             
c     2  g q     -> q Z              and c.c.          
c     3  q g     -> q Z              and c.c.        
c     4  gg      -> q qbar Z               
c     5  q q     -> q q Z            and c.c.
c     6  q q     -> q' q' Z          and c.c.
c     7  q q'    -> q  q' Z          and c.c. (q, q' same isodoublet)
c     8  q q''   -> q q'' Z          and c.c. (q, q''diff iso.)
c     9  q q'    -> q'' q''' Z       and c.c.    
c     10  q qbar  -> q qbar Z         and c.c.    
c     11  q qbar  -> q' qbar' Z       and c.c. (q, q' same iso.)
c     12  q qbar  -> q'' qbar'' Z     and c.c. (q, q''diff iso.)
c     13  q qbar' -> q qbar' or qbar q' Z and c.c. (q, q' same iso.)
c     14  q qbar''-> q qbar''         and c.c. (q, q''diff iso.)  
c     15  q qbar' -> q'' qbar''' Z    and c.c. (q, q' same iso.)  
c     16  q qbar''-> q' qbar''' Z     and c.c. (q, q''diff iso.) 
c     17  g q     -> q q qbar Z       and c.c.                    
c     18  q g     -> q q qbar Z       and c.c.                    
c     19  g q     -> q q' qbar' Z     and c.c. (q, q' same iso.)  
c     20  q g     -> q q' qbar' Z     and c.c. (q, q' same iso.)       
c     21  g q     -> q q'' qbar'' Z   and c.c. (q, q''diff iso.)  
c     22  q g     -> q q'' qbar'' Z   and c.c. (q, q''diff iso.)       
c     23  g q     -> q' q' qbar Z     and c.c. (q, q' same iso.)       
c     24  q g     -> q' q' qbar Z     and c.c. (q, q' same iso.)
C     
c     25  g q     -> q' q'' qbar''' Z and c.c.          
c     26  q g     -> q' q'' qbar''' Z and c.c. 
c     
c------------------------------------------------------------------------
C     -
      implicit none
      include 'alpgen.inc'
      include 'vbjet.inc'
      integer icconj
      common/hwconv/icconj
      integer mask(11),jproc1
      integer il1,il2,afl(maxpar)
      real tmp(100),slum,cwgt,swgt,rn,tmptot
      real fcab(2)
c     pippoq
      integer jj(10)
      integer nchan(1:6,0:4),jad(1:6,1:7),lrun
      real ipro
      integer icab(-4:4,2) 
c     cabibbo partner: (*,1)=cabibbo allowed, (*,2)-cabibbo suppressed
c     where 1=d 2=u 3=s 4=c 0=g and negatives are antiparticles
      data icab/-3,-4,-1,-2,0,2,1,4,3, 
     +     -1,-2,-3,-4,0,4,3,2,1/
      integer iqpp(2,2,4)
c     iqpp(k,j,i):  if q=i, q'=cabibbo partner of i (allowed for j=1,
c     suppressed for j=2) then iqpp are the two (k=1,2) flavours .ne. (q
C     .or.q')
      data iqpp/
c     du   dc   ud   us   sc   su   cs   cd   
     +     3,4, 2,3, 3,4, 1,4, 1,2, 1,4, 1,2, 2,3/
c     parton charges
      real dch,uch,dbch,ubch,wchrg,chrg(-6:6)
      parameter (dch=-0.333333333e0,uch=0.666666666e0)
      parameter (dbch=0.333333333e0,ubch=-0.666666666e0)
      data chrg/ubch,dbch,ubch,dbch,ubch,dbch,0e0,  
     +     dch, uch, dch, uch, dch, uch/
      integer nf
      real *8 xlum,xrn,rcho                                 
      integer i,j,k,l,ik,is,itmp,icount,init,jproc,ng,nlqp
      real cfac(0:5),ifact(0:8),invbosew(0:8)
      data cfac/8e0,5*3e0/,init/0/
      data invbosew/3*1e0,2e0,4e0,12e0,36e0,144e0,576e0/
      integer imap(-4:4),iswap
      data imap/-2,-1,-2,-1,0,1,2,1,2/
      data mask/5,7,8,9,10,11,12,13,14,15,16/
c     
      integer jp
      common/jpr/jp
c     overall efficiency for extraction of colour states:
c     #(non-zero color states) / 3**nq*8**ng
c     ccoef(i,j) for njets=j and i=#(light quark pairs)
      real*8 ccoef(3,0:6)       !ccoeff(number of q-pb pairs, number of gluons)
      data ccoef/  0.333333333,0.185185185,0.127572016,
     $     0.166666667,0.12037037,0.0936213992,
     $     0.114583333,0.0902777778,0.0743312757,
     $     0.0872395833,0.072337963,-1.d0,
     $     0.0704752604,0.0603841146,-1.d0,
     $     0.0591227214,-1d0,-1d0,
     $     0.0509236654,-1d0,-1d0/
c     2 quark pairs, ng=0,1,2,3:
c     15/81,0d0, 78/648,0d0, 468/5184, 93/729, 
c     3 quark pairs, ng=0,1,2,3:
c     3000/41472, 546/5832, 20034/331776, 3468/46656
c     1 quark pair, ng=0,1,2,3,4,5,6:
c     3/9, 12/72, 66/576, 402/4608, 2598/36864, 17436/294912, 120144
C     /2359296
c     
c     pippoa
      common/char/wchrg
C     
      real*8 qw(4),qf(4),qfb(4)
      integer j1,jwp,jwm,jz,jh,j2
      integer lblwm(maxpar),lblwp(maxpar),lblz(maxpar),lblh(maxpar)
      complex*16 eps(4,maxpar),aux(4)
      common/vdec/eps
      real*8 gv,ga,gvn,gve,gvu,gvd,gvv(4),gaa(4)
      real*8 zbr,zbree,zbruu,zbrnn,zbrdd,zbrcum(4),zqcum,zbrt
      character*1 zd(maxpar-2)
      common/zdec/zd
      integer exch(4)
      data exch/2,3,4,1/
C     
      save init,chrg,cfac,ifact,invbosew,iqpp,icab,fcab,zbrt
c     
      if(init.eq.0) then
        init=1
*     
c     protection:
*     
        fcab(1)=real(ccab2)
        fcab(2)=real(scab2)
        if (nz+nw+nh+njets+nph+2.gt.10) then
          print*,'nz+nw+nh+njets > 8 is not allowed!'
          stop
        endif
        if (njets.gt.3) then
          print*,'njets > 3 is not allowed!'
          stop
        endif
        ifact(0)=1e0
        do i=1,8
          ifact(i)=ifact(i-1)/real(i)
        enddo
c     
ccc modified
c        if(nz0.eq.0) then
        if(nz.eq.0) then
ccc modified
           zbrt= 1.d0
        else
           gvn=0.5d0
           gve=0.5d0*(-1.d0+4.d0*stw*stw)
           gvu=0.5d0*(1.d0-8.d0/3.d0*stw*stw)
           gvd=0.5d0*(-1.d0+4.d0/3.d0*stw*stw)
           gvv(1)=gvn
           gvv(2)=gve
           gvv(3)=gvu
           gvv(4)=gvd
           gaa(1)=-0.5d0
           gaa(2)=0.5d0
           gaa(3)=-0.5d0
           gaa(4)=0.5d0
           zbrnn=0.25d0+gvn**2
           zbree=0.25d0+gve**2
           zbruu=0.25d0+gvu**2
           zbrdd=0.25d0+gvd**2
           zbr=3.d0*(zbrnn+zbree+3.d0*zbrdd)+2.d0*3.d0*zbruu
           zbrnn=3.d0*zbrnn/zbr
           zbree=3.d0*zbree/zbr
           zbruu=2.d0*3.d0*zbruu/zbr
           zbrdd=3.d0*3.d0*zbrdd/zbr
           zbrcum(1)=zbrnn
           zbrcum(2)=zbree+zbrcum(1)
           zbrcum(3)=zbruu+zbrcum(2)
           zbrcum(4)=zbrdd+zbrcum(3)
           zqcum=zbruu/(zbruu+zbrdd)
           call setbr(nz,zbrcum,zd,zbrt)
        endif
c
      endif
*     
 1    call randa(xrn)
      rn=real(xrn)
      if(1e0-rn.lt.1e-6) goto 1
*     
      icount=0
      slum=0e0
      tmptot=0e0
      iswap=0
      do i=1,npart
        afl(i)=0
        ifl(i)=0
      enddo
*     
c     for nw= even goto 50
*     
      if (mod(nw,2).eq.0) goto 50 
*     
c     nw= odd
*     
c     start two light-quark processes 
*     
c     jproc=1   njets>=0
c     q qbar' -> W  ,   qbar q' -> W
*     
      if(jproc.eq.1) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(-icab(i,j))*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo  
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=-icab(i,j)
                nf= 2
                afl(1)= i
                afl(2)=-icab(i,1)
                goto 100
              endif 
            enddo
          endif
        enddo 
*     
c     jproc=2   
c     g q -> q' W    and    g qbar -> q'bar W
*     
      elseif(jproc.eq.2) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(0)*f2(i)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=0
                ifl(2)=i
                ifl(3)=icab(i,j)
                nf= 3
                afl(1)=0
                afl(2)=i
                afl(3)=icab(i,1)
                goto 100   
              endif 
            enddo
          endif
        enddo 
*     
c     jproc=3   
c     q g -> q' W    and  qbar g -> q'bar W
*     
      elseif(jproc.eq.3) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(0)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=0
                ifl(3)=icab(i,j)
                nf= 3
                afl(1)=i
                afl(2)=0
                afl(3)=icab(i,1)
                goto 100 
              endif 
            enddo
          endif
        enddo 
*     
c     jproc=4   
c     gg -> q qbar' W              
*     
      elseif(jproc.eq.4) then
        do i=1,4
          do j=1,2
            icount=icount+1
            tmp(icount)=fcab(j)
            slum=slum+tmp(icount)
          enddo
        enddo 
        icount=0
        rn=rn*slum
        slum=slum*f1(0)*f2(0)
        do i=1,4
          do j=1,2
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)=0
              ifl(2)=0
              ifl(3)=i
              ifl(4)=-icab(i,j)
              nf= 4
              afl(1)=0
              afl(2)=0
              afl(3)=i
              afl(4)=-icab(i,1)
              goto 100
            endif
          enddo
        enddo 
c     
c     start 4 light quark processes  (in order of increasing number of
c     minimum njets, different flavour first, like-flavours after
c     
*     
c     jproc=5  
c     q qbar' -> W q'' qbar'' ,   qbar q' -> W q'' qbar''
*     
      elseif(jproc.eq.5) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(-icab(i,j))*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
c     include factor of two for the two possible q'' flavours .ne. q or
c     q'
        slum=2e0*slum
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              do k=1,2
c     consider the two possible q'' flavours
                tmptot=tmptot+tmp(icount)
                if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=-icab(i,j)
                  ifl(3)=iqpp(k,j,abs(i))
                  ifl(4)=-ifl(3)
                  nf= 4
                  afl(1)=i
                  afl(2)=-icab(i,1)
                  afl(3)=iqpp(k,1,abs(i))
                  afl(4)=-afl(3)
                  goto 100
                endif 
              enddo
            enddo 
          endif
        enddo
*     
c     jproc=6 
c     q(bar) q(bar)'' -> W q'(bar) q(bar)''  
*     
      elseif(jproc.eq.6) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              is=1
              do k=1,4
                ik=k
                if(k.gt.2) then
                  ik=k-2
                  is=-1
                endif
                icount=icount+1
                tmp(icount)=f1(i)*f2(is*iqpp(ik,j,abs(i)))*fcab(j)
                slum=slum+tmp(icount)
              enddo
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              is=1
              do k=1,4
                icount=icount+1
                tmptot=tmptot+tmp(icount)
                if(tmptot.ge.rn) then 
                  ik=k
                  if(k.gt.2) then
                    ik=k-2
                    is=-1
                  endif
                  ifl(1)=i
                  ifl(2)=is*iqpp(ik,j,abs(i))
                  ifl(3)=icab(i,j)
                  ifl(4)=ifl(2)
                  nf= 4
                  afl(1)=i
                  afl(2)=is*iqpp(ik,1,abs(i))
                  afl(3)=icab(i,1)
                  afl(4)=afl(2)
                  goto 100
                endif 
              enddo
            enddo
          endif
        enddo 
*     
c     jproc=7  
c     q''(bar) q(bar) -> W q'(bar) q''(bar)
*     
      elseif(jproc.eq.7) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              is=1
              do k=1,4
                ik=k
                if(k.gt.2) then
                  ik=k-2
                  is=-1
                endif
                icount=icount+1
                tmp(icount)=f1(is*iqpp(ik,j,abs(i)))*f2(i)*fcab(j)
                slum=slum+tmp(icount)
              enddo
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              is=1
              do k=1,4
                icount=icount+1
                tmptot=tmptot+tmp(icount)
                if(tmptot.ge.rn) then 
                  ik=k
                  if(k.gt.2) then
                    ik=k-2
                    is=-1
                  endif
                  ifl(1)=is*iqpp(ik,j,abs(i))
                  ifl(2)=i
                  ifl(3)=icab(i,j)
                  ifl(4)=ifl(1)
                  nf= 4
                  afl(1)=is*iqpp(ik,1,abs(i))
                  afl(2)=i
                  afl(3)=icab(i,1)
                  afl(4)=afl(1)
                  goto 100
                endif 
              enddo
            enddo
          endif
        enddo 
*     
c     jproc=8 
c     q'' q''bar -> W q q'bar
*     
      elseif(jproc.eq.8) then
        do i=-4,4
          if(i.ne.0) then
            is=1
            do k=1,4
              if(k.ne.abs(i)) then
                do j=1,2
                  if(icab(k,j).ne.abs(i)) then
                    icount=icount+1
                    tmp(icount)=f1(i)*f2(-i)*fcab(j)
                    slum=slum+tmp(icount)
                  endif
                enddo
              endif
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,4
              if(k.ne.abs(i)) then
                do j=1,2
                  if(icab(k,j).ne.abs(i)) then
                    icount=icount+1
                    tmptot=tmptot+tmp(icount)
                    if(tmptot.ge.rn) then 
                      ifl(1)=i
                      ifl(2)=-i
                      ifl(3)=k
                      ifl(4)=-icab(k,j)
                      nf= 4
                      afl(1)=i
                      afl(2)=-i
                      afl(3)=k
                      afl(4)=-icab(k,1)
                      goto 100
                    endif
                  endif
                enddo
              endif
            enddo 
          endif
        enddo 
c     
c     then processes with like flavours
c     
*     
c     jproc=9  
c     q qbar' -> W q qbar ,   qbar q' -> W qbar q 
*     
      elseif(jproc.eq.9) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(-icab(i,j))*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=-icab(i,j)
                ifl(3)=i
                ifl(4)=-ifl(3)
                nf= 4
                afl(1)=i
                afl(2)=-icab(i,1)
                afl(3)=i
                afl(4)=-afl(3)
                goto 100
              endif 
            enddo
          endif
        enddo
*     
c     jproc=10 
c     qbar' q -> W q qbar ,   q' qbar -> W qbar q
*     
      elseif(jproc.eq.10) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(-icab(i,j))*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(2)=i
                ifl(1)=-icab(i,j)
                ifl(3)=i
                ifl(4)=-ifl(3)
                nf= 4
                afl(2)=i
                afl(1)=-icab(i,1)
                afl(3)=i
                afl(4)=-afl(3)
                goto 100
              endif 
            enddo 
          endif
        enddo
*     
c     jproc=11 
c     q qbar -> W q qbar'  + c.c.
*     
      elseif(jproc.eq.11) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(-i)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=-i
                ifl(3)=i
                ifl(4)=-icab(i,j)
                nf= 4
                afl(1)=i
                afl(2)=-i
                afl(3)=i
                afl(4)=-icab(i,1)
                goto 100
              endif
            enddo 
          endif
        enddo
*     
c     jproc=12 nejts>=4
c     q qbar -> W q' qbar
*     
      elseif(jproc.eq.12) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(-i)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=-i
                ifl(3)=icab(i,j)
                ifl(4)=-i
                nf= 4
                afl(1)=i
                afl(2)=-i
                afl(3)=icab(i,1)
                afl(4)=-i
                goto 100
              endif
            enddo 
          endif
        enddo
*     
c     jproc=13  
c     q q -> W q q'
*     
      elseif(jproc.eq.13) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(i)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=i
                ifl(3)=i
                ifl(4)=icab(i,j)
                nf= 4
                afl(1)=i
                afl(2)=i
                afl(3)=i
                afl(4)=icab(i,1)
                goto 100
              endif
            enddo 
          endif
        enddo
c     
c     jproc=14 
c     q q' -> W q q
c     
      elseif(jproc.eq.14) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(icab(i,j))*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=icab(i,j)
                ifl(3)=i
                ifl(4)=i
                nf= 4
                afl(1)=i
                afl(2)=icab(i,1)
                afl(3)=i
                afl(4)=i
                goto 100
              endif
            enddo 
          endif
        enddo
*     
c     jproc=15 
c     q q' -> W q' q'
*     
      elseif(jproc.eq.15) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(icab(i,j))*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=icab(i,j)
                ifl(3)=ifl(2)
                ifl(4)=ifl(2)
                nf= 4
                afl(1)=i
                afl(2)=icab(i,1)
                afl(3)=afl(2)
                afl(4)=afl(2)
                goto 100
              endif
            enddo 
          endif
        enddo
*     
c     jproc=16 
c     q g -> W q' q'' qbar'' ,   qbar g -> W qbar' q'' qbar''
*     
      elseif(jproc.eq.16) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(0)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
c     include factor of two for the two possible q'' flavours .ne. q or
c     q'
        slum=2e0*slum
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              do k=1,2
c     consider the two possible q'' flavours
                tmptot=tmptot+tmp(icount)
                if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=0
                  ifl(3)=icab(i,j)
                  ifl(4)=iqpp(k,j,abs(i))
                  ifl(5)=-ifl(4)
                  nf= 5
                  afl(1)=i
                  afl(2)=0
                  afl(3)=icab(i,1)
                  afl(4)=iqpp(k,1,abs(i))
                  afl(5)=-afl(4)
                  goto 100
                endif 
              enddo
            enddo
          endif
        enddo
*     
c     jproc=17  
c     g q -> W q' q'' qbar'' ,   g qbar -> W qbar' q'' qbar''
*     
      elseif(jproc.eq.17) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(0)*f2(i)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
c     include factor of two for the two possible q'' flavours .ne. q or
c     q'
        slum=2e0*slum
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              do k=1,2
c     consider the two possible q'' flavours
                tmptot=tmptot+tmp(icount)
                if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=i
                  ifl(3)=icab(i,j)
                  ifl(4)=iqpp(k,j,abs(i))
                  ifl(5)=-ifl(4)
                  nf= 5
                  afl(1)=0
                  afl(2)=i
                  afl(3)=icab(i,1)
                  afl(4)=iqpp(k,1,abs(i))
                  afl(5)=-afl(4)
                  goto 100
                endif 
              enddo
            enddo 
          endif
        enddo
*     
c     jproc=18  
c     q g -> W q q qbar'
*     
      elseif(jproc.eq.18) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(0)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=0
                ifl(3)=i
                ifl(4)=i
                ifl(5)=-icab(i,j)
                nf= 5
                afl(1)=i
                afl(2)=0
                afl(3)=i
                afl(4)=i
                afl(5)=-icab(i,1)
                goto 100
              endif 
            enddo
          endif
        enddo
*     
c     jproc=19  
c     q g -> W q' q qbar
*     
      elseif(jproc.eq.19) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(i)*f2(0)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=i
                ifl(2)=0
                ifl(3)=icab(i,j)
                ifl(4)=i
                ifl(5)=-i
                nf= 5
                afl(1)=i
                afl(2)=0
                afl(3)=icab(i,1)
                afl(4)=i
                afl(5)=-i
                goto 100
              endif 
            enddo
          endif
        enddo
*     
c     jproc=20  
c     g q -> W q q qbar' 
*     
      elseif(jproc.eq.20) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(0)*f2(i)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=0
                ifl(2)=i
                ifl(3)=i
                ifl(4)=i
                ifl(5)=-icab(i,j)
                nf= 5
                afl(1)=0
                afl(2)=i
                afl(3)=i
                afl(4)=i
                afl(5)=-icab(i,1)
                goto 100
              endif 
            enddo
          endif
        enddo
*     
c     jproc=21  
c     g q -> W q' q qbar
*     
      elseif(jproc.eq.21) then
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmp(icount)=f1(0)*f2(i)*fcab(j)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do j=1,2
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)=0
                ifl(2)=i
                ifl(3)=icab(i,j)
                ifl(4)=i
                ifl(5)=-i
                nf= 5
                afl(1)=0
                afl(2)=i
                afl(3)=icab(i,1)
                afl(4)=i
                afl(5)=-i
                goto 100
              endif 
            enddo
          endif
        enddo

      else
        write(*,*) 'jproc not defined, slum=0'
        slum=0
        stop
      endif
*     
      xlum=-1d0
      return
 50   continue 
*     
c     nw= even 
*     
c     jproc=1  
c     q qbar -> Z  + c.c.
*     
      if(jproc.eq.1) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)= f1(i)*f2(-i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)=  i
              ifl(2)= -i 
              nf= 2
              goto 100
            endif 
          endif
        enddo 
*     
c     jproc=2   
c     g q -> q  Z  + c.c
*     
      elseif(jproc.eq.2) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)=f1(0)*f2(i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= 0
              ifl(2)= i
              ifl(3)= ifl(2)
              nf= 3
              goto 100   
            endif 
          endif
        enddo 
*     
c     jproc=3   
c     q g -> q Z + c.c.
*     
      elseif(jproc.eq.3) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)=f1(i)*f2(0)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= i
              ifl(2)= 0
              ifl(3)= ifl(1)
              nf= 3
              goto 100 
            endif 
          endif
        enddo 
*     
c     jproc=4   
c     gg -> q qbar Z              
*     
      elseif(jproc.eq.4) then
        do i=1,4
          icount=icount+1
          tmp(icount)=f1(0)*f2(0)
          slum=slum+tmp(icount)
        enddo 
        icount=0
        rn=rn*slum
        do i=1,4
          icount=icount+1
          tmptot=tmptot+tmp(icount)
          if(tmptot.ge.rn) then 
            ifl(1)= 0
            ifl(2)= 0
            ifl(3)= i
            ifl(4)=-ifl(3)
            nf= 4
            goto 100
          endif
        enddo 
*     
c     jproc=5  
c     q q -> q q Z + c.c.
*     
      elseif(jproc.eq.5) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)=f1(i)*f2(i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)=i
              ifl(2)=i
              ifl(3)=ifl(1)
              ifl(4)=ifl(2)
              nf= 4
              goto 100
            endif
          endif 
        enddo
*     
c     jproc=6  
c     q q -> q' q' Z + c.c.
*     
      elseif(jproc.eq.6) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)=f1(i)*f2(i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= i
              ifl(2)= i
              ifl(3)= icab(ifl(1),1)
              ifl(4)= icab(ifl(2),1)
              nf= 4
              goto 100
            endif
          endif 
        enddo
*     
c     jproc=7  
c     q q' -> q q' Z + c.c. (q q' in the same isodoublet)
*     
      elseif(jproc.eq.7) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)=f1(i)*f2(icab(i,1))
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= i
              ifl(2)= icab(i,1)
              ifl(3)= ifl(1)
              ifl(4)= ifl(2)
              nf= 4
              goto 100
            endif
          endif 
        enddo
*     
c     jproc=8  
c     q q'' -> q q'' Z + c.c. (q q'' in different isodoublets)
*     
      elseif(jproc.eq.8) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(sign(iqpp(k,1,abs(i)),i))
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)= sign(iqpp(k,1,abs(i)),i)
                ifl(3)= ifl(1)
                ifl(4)= ifl(2)
                nf= 4
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=9  
c     q q' -> q'' q''' Z + c.c. 
*     
      elseif(jproc.eq.9) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(sign(iqpp(k,1,abs(i)),i))
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)= sign(iqpp(k,1,abs(i)),i)
                ifl(3)= icab(ifl(1),1)
                ifl(4)= icab(ifl(2),1)
                nf= 4
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=10  
c     q qbar -> q qbar Z + c.c.
*     
      elseif(jproc.eq.10) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)=f1(i)*f2(-i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= i
              ifl(2)=-i
              ifl(3)= ifl(1)
              ifl(4)= ifl(2)
              nf= 4
              goto 100
            endif
          endif 
        enddo
*     
c     jproc=11  
c     q qbar -> q' qbar' Z + c.c. (q q' in the same isodoublet)
*     
      elseif(jproc.eq.11) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)=f1(i)*f2(-i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= i
              ifl(2)=-i
              ifl(3)= icab(ifl(1),1)
              ifl(4)= icab(ifl(2),1)
              nf= 4
              goto 100
            endif
          endif 
        enddo
*     
c     jproc=12  
c     q qbar -> q'' qbar'' Z + c.c. (q q'' in different isodoublets)
*     
      elseif(jproc.eq.12) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(-i)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)=-i
                ifl(3)= iqpp(k,1,abs(i)) 
                ifl(4)=-ifl(3)
                nf= 4
                goto 100
              endif

            enddo
          endif 
        enddo
*     
c     jproc=13  
c     q qbar' -> q qbar' or qbar q' Z + c.c. (q q' in the same
C     isodoublet)
*     
      elseif(jproc.eq.13) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(-icab(i,1))
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)=-icab(i,1)
                if (k.eq.1) then
                  ifl(3)= ifl(1)
                  ifl(4)= ifl(2)
                else
                  ifl(3)=-ifl(2)
                  ifl(4)=-ifl(1) 
                endif
                nf= 4
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=14  
c     q qbar'' -> q qbar'' Z + c.c. (q q'' in different isodoublets)
*     
      elseif(jproc.eq.14) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(-sign(iqpp(k,1,abs(i)),i))
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)=-sign(iqpp(k,1,abs(i)),i)
                ifl(3)= ifl(1)
                ifl(4)= ifl(2)
                nf= 4
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=15  
c     q qbar' -> q'' qbar''' Z + c.c. (q q' in the same isodoublet)
*     
      elseif(jproc.eq.15) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(-icab(i,1))
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)=-icab(i,1)
                ifl(3)= sign(iqpp(k,1,abs(ifl(1))),ifl(1))
                ifl(4)=-icab(ifl(3),1)
                nf= 4
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=16  
c     q qbar'' -> q' qbar''' Z + c.c. (q q'' in different isodoublets)
*     
      elseif(jproc.eq.16) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(-sign(iqpp(k,1,abs(i)),i))
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)=-sign(iqpp(k,1,abs(i)),i)
                ifl(3)= icab(ifl(1),1)
                ifl(4)= icab(ifl(2),1)
                nf= 4
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=17  
c     g q -> q q qbar Z  + c.c.
*     
      elseif(jproc.eq.17) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)= f1(0)*f2(i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= 0
              ifl(2)= i 
              ifl(3)= ifl(2)
              ifl(4)= ifl(2)
              ifl(5)=-ifl(2)
              nf= 5
              goto 100
            endif 
          endif
        enddo 
*     
c     jproc=18  
c     q g -> q q qbar Z  + c.c.
*     
      elseif(jproc.eq.18) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)= f1(i)*f2(0)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= i
              ifl(2)= 0 
              ifl(3)= ifl(1)
              ifl(4)= ifl(1)
              ifl(5)=-ifl(1)
              nf= 5
              goto 100
            endif 
          endif
        enddo 
*     
c     jproc=19 
c     g q -> q q' qbar' Z  + c.c. (q q' in the same isodoublet)
*     
      elseif(jproc.eq.19) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)= f1(0)*f2(i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= 0
              ifl(2)= i 
              ifl(3)= ifl(2)
              ifl(4)= icab(ifl(2),1)
              ifl(5)=-ifl(4)
              nf= 5
              goto 100
            endif 
          endif
        enddo 
*     
c     jproc=20  
c     q g -> q q' qbar' Z  + c.c. (q q' in the same isodoublet)
*     
      elseif(jproc.eq.20) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)= f1(i)*f2(0)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= i
              ifl(2)= 0 
              ifl(3)= ifl(1)
              ifl(4)= icab(ifl(1),1)
              ifl(5)=-ifl(4)
              nf= 5
              goto 100
            endif 
          endif
        enddo 
*     
c     jproc=21 
c     g q -> q q'' qbar'' Z  + c.c. (q q'' in different isodoublets)
*     
      elseif(jproc.eq.21) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(0)*f2(i)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= 0
                ifl(2)= i
                ifl(3)= ifl(2)
                ifl(4)= sign(iqpp(k,1,abs(ifl(2))),ifl(2))
                ifl(5)=-ifl(4)
                nf= 5
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=22  
c     q g -> q q'' qbar'' Z  + c.c. (q q'' in different isodoublets)
*     
      elseif(jproc.eq.22) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(0)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)= 0
                ifl(3)= ifl(1)
                ifl(4)= sign(iqpp(k,1,abs(ifl(1))),ifl(1))
                ifl(5)=-ifl(4)
                nf= 5
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=23  
c     g q -> q' q' qbar Z  + c.c. (q q' in the same isodoublet)
*     
      elseif(jproc.eq.23) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)= f1(0)*f2(i)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= 0
              ifl(2)= i 
              ifl(3)= icab(ifl(2),1)
              ifl(4)= ifl(3)
              ifl(5)=-ifl(2)
              nf= 5
              goto 100
            endif 
          endif
        enddo 
*     
c     jproc=24  
c     q g -> q' q' qbar Z  + c.c. (q q' in the same isodoublet)
*     
      elseif(jproc.eq.24) then
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmp(icount)= f1(i)*f2(0)
            slum=slum+tmp(icount)
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
              ifl(1)= i
              ifl(2)= 0 
              ifl(3)= icab(ifl(1),1)
              ifl(4)= ifl(3)
              ifl(5)=-ifl(1)
              nf= 5
              goto 100
            endif 
          endif
        enddo 
*     
c     jproc=25 
c     g q -> q' q'' qbar''' Z  + c.c. 
*     
      elseif(jproc.eq.25) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(0)*f2(i)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= 0
                ifl(2)= i
                ifl(3)= icab(ifl(2),1)
                ifl(4)= sign(iqpp(k,1,abs(ifl(2))),ifl(2))
                ifl(5)=-icab(ifl(4),1)
                nf= 5
                goto 100
              endif
            enddo
          endif 
        enddo
*     
c     jproc=26 
c     q g -> q' q'' qbar''' Z  + c.c. 
*     
      elseif(jproc.eq.26) then
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmp(icount)=f1(i)*f2(0)
              slum=slum+tmp(icount)
            enddo
          endif
        enddo 
        icount=0
        rn=rn*slum
        do i=-4,4
          if(i.ne.0) then
            do k=1,2  
              icount=icount+1
              tmptot=tmptot+tmp(icount)
              if(tmptot.ge.rn) then 
                ifl(1)= i
                ifl(2)= 0
                ifl(3)= icab(ifl(1),1)
                ifl(4)= sign(iqpp(k,1,abs(ifl(1))),ifl(1))
                ifl(5)=-icab(ifl(4),1)
                nf=5
                goto 100
              endif
            enddo
          endif 
        enddo
      else
        print*,'ERROR in jproc when nw= even'
        stop
      endif
*     
      xlum=-1d0
      return
 100  continue
*     
      if (mod(nw,2).eq.0) then
        do l= 1,nf
          afl(l)= ifl(l)
        enddo
      endif
*     
c     Additional gluons
*     
      do l= nf+1,njets+2
        ifl(l)=0
      enddo 
*
c     gamma's
      do k= 1,nph               
        ifl(k+njets+2)= 22                 
      enddo                         
*     
c     Z's
*     
      do k= 1,nz               
        ifl(k+njets+nph+2)= 23                 
      enddo                         
*     
      if (nw.eq.0) then
        is= 0
        goto 110
      endif
*     
c     Assign charges
*     
      wchrg= chrg(ifl(1))+chrg(ifl(2))
      do i= 3,nf
        wchrg=wchrg-chrg(ifl(i))
      enddo
*     
c     W's 
*     
      if (mod(nw,2).eq.0) then
*     
c     when wchrg= 2, the order is  W+(1) W-(2) ... W+(nw-1) W+(nw) 
c     when wchrg= 0, the order is  W+(1) W-(2) ... W+(nw-1) W-(nw) 
c     when wchrg=-2, the order is  W+(1) W-(2) ... W-(nw-1) W-(nw) 
*     
        do k= 1,nw-2
          ifl(k+njets+2+nz+nph)= -24*(-1)**k
        enddo 
*     
        if(abs(wchrg-2.e0).lt.1.e-3) then
          is= 2
          ifl(nw-1+njets+2+nz+nph)= 24
          ifl(nw  +njets+2+nz+nph)= 24
        elseif(abs(wchrg).lt.1.e-3) then
          is= 0
          ifl(nw-1+njets+2+nz+nph)=  24
          ifl(nw  +njets+2+nz+nph)= -24
        elseif(abs(wchrg+2.e0).lt.1.e-3) then
          is=-2
          ifl(nw-1+njets+2+nz+nph)= -24
          ifl(nw  +njets+2+nz+nph)= -24
        else
          write(*,*) 'Something wrong with W charges, '
          write(*,*) 'stop and check charge conservation'
          stop
        endif
      else
*     
c     when wchrg= 1, the order is  W+(1) W-(2) ... W+(nw) 
c     when wchrg=-1, the order is  W+(1) W-(2) ... W-(nw) 
*     
        do k= 1,nw-1
          ifl(k+njets+2+nz+nph)= -24*(-1)**k
        enddo 
*     
        if (abs(wchrg-1.e0).lt.1.e-3) then
          is= 1
          ifl(nw+njets+2+nz+nph)= 24
        elseif (abs(wchrg+1.e0).lt.1.e-3) then
          is=-1
          ifl(nw+njets+2+nz+nph)= -24
        else
          write(*,*) 'Something wrong with W charges, '
          write(*,*) 'stop and check charge conservation'
          stop
        endif
      endif
*     
c     and then H's
*     

 110  do k= 1,nh               
        ifl(k+njets+2+nw+nz+nph)= 25                 
      enddo                         
c
c      do k=1,nph
c        ifl(k+njets+2+nw+nz+nh)= 22                 
c      enddo
*     
c     evaluate colour weight factors
*     
      do i=1,2
        if(ifl(i).eq.0) ifl(i)=21
      enddo
*     
      ng=0
      cwgt=1e0
      do i=3,njets+2
        cwgt=cwgt*cfac(abs(ifl(i)))
        if(ifl(i).eq.0) then
          ifl(i)=21
          ng=ng+1
        endif
      enddo
*     
c     evaluate spin weight factors
*     
      swgt=2e0
      swgt=swgt**(njets+nph)
      swgt=swgt*3.d0**(nz+nw)
*     
      do k= 1,nph               
        afl(k+njets+2)= ifl(k+njets+2)
      enddo                         
      do k= 1,nz               
        afl(k+njets+2+nph)= ifl(k+njets+2+nph)                  
      enddo                         
*     
      do k= 1,nw
        afl(k+njets+2+nz+nph)= ifl(k+njets+2+nz+nph)
      enddo 
*     
c     poi
      do k= 1,nh               
        afl(k+njets+2+nz+nph+nw)= ifl(k+njets+2+nz+nph+nw)                  
      enddo                         
*     
*     
      if(resc.ne.1) then
        write(*,*) 'resc.ne.1 is not suitable for this process, stop'
        stop
      endif 
      xlum= dble(slum*cwgt*swgt*ifact(ng)*ifact(nz)*ifact(nh)
     +     *ifact((nw-is)/2)*ifact((nw+is)/2)*ifact(nph))/resc**njets
c     rescaling for the strong coupling
c     xlum= dble(xlum*(4*pi*as)**njets)
c     
*     
      if (mod(nw,2).eq.0) then
!     identical quarks in final state
        if(jproc.eq.5.or.jproc.eq.6.or.jproc.eq.17.or.jproc.eq.18.or
     $       .jproc.eq.23.or.jproc.eq.24)  xlum=0.5d0*xlum 
        if(jproc.le.4) then
          nlqp=1
        elseif(jproc.le.26) then
          nlqp=2
        else
          write(*,*) 'jproc not valid'
          stop
        endif 
      else
        if(jproc.eq.14.or.jproc.eq.23.or.jproc.eq.20.or.jproc.eq.18.or
     $       .jproc.eq.15)  xlum=0.5d0*xlum ! identical quarks in final state
        if(jproc.le.4) then
          nlqp=1
        elseif(jproc.le.23) then
          nlqp=2
        else
          write(*,*) 'jproc not valid'
          stop
        endif 
      endif
c     
      xlum=xlum*ccoef(nlqp,2+njets-2*nlqp)
      jp=jproc
c     
      jwp=0
      do j1=1,npart
        lblwp(j1)=0
        if(ifl(j1).eq.24) then
          jwp=jwp+1
          lblwp(jwp)=j1
        endif
      enddo
      jwm=0
      do j1=1,npart
        lblwm(j1)=0
        if(ifl(j1).eq.-24) then
          jwm=jwm+1
          lblwm(jwm)=j1
        endif
      enddo
      jz=0
      do j1=1,npart
        lblz(j1)=0
        if(abs(ifl(j1)).eq.23) then
          jz=jz+1
          lblz(jz)=j1
        endif
      enddo
C     
      if(idecay.eq.'y') then
        do j2=1,jz
          qw(1)=p(4,lblz(j2))
          do j1=2,4
            qw(j1)=p(j1-1,lblz(j2))
          enddo
          if (zd(j2).eq.'n') then
            call randa(xrn)
            if(xrn.lt.0.3333333d0) then
              zfl(j2)= 12
            elseif(xrn.lt.0.66666666d0) then
              zfl(j2)= 14
            else
              zfl(j2)= 16
            endif
            gv=0.5d0
            ga=-0.5d0
          elseif (zd(j2).eq.'e') then
            call randa(xrn)
            if(xrn.lt.0.3333333d0) then
              zfl(j2)= 11
            elseif(xrn.lt.0.66666666d0) then
              zfl(j2)= 13
            else
              zfl(j2)= 15
            endif
            gv=gve
            ga=0.5d0
          elseif (zd(j2).eq.'q') then
            call randa(xrn)
            if(xrn.lt.zqcum) then
              call randa(xrn)
              if(xrn.lt.0.5d0) then
                zfl(j2)= 2
              else
                zfl(j2)= 4
              endif
              gv=gvu
              ga=-0.5d0            
            else
              call randa(xrn)
              if(xrn.lt.0.3333333d0) then
                zfl(j2)= 1
              elseif(xrn.lt.0.66666666d0) then
                zfl(j2)= 3
              else
                zfl(j2)= 5
              endif
              gv=gvd
              ga=0.5d0            
            endif
          elseif (zd(j2).eq.'b') then
            zfl(j2)= 5
            gv=gvd
            ga=0.5d0
          elseif (zd(j2).eq.'i') then
            call randa(xrn)
            j1= 1
            do i=1,3
              if(xrn.gt.zbrcum(i)) j1=i+1
            enddo
            gv=gvv(j1)
            ga=gaa(j1)
            call randa(xrn)
            if(j1.eq.3) then
              if(xrn.lt.0.5d0) then
                zfl(j2)= 2
              else
                zfl(j2)= 4
              endif
            else 
              if(j1.eq.1) j1= 12
              if(j1.eq.2) j1= 11
              if(j1.eq.4) j1= 1
              if(xrn.lt.0.3333333d0) then
                zfl(j2)= j1
              elseif(xrn.lt.0.66666666d0) then
                zfl(j2)= j1+2
              else
                zfl(j2)= j1+4
              endif
            endif
          endif
          call vpol(mz,qw,ga,gv,qf,qfb,aux)
          do i=1,4
            idec(i,1,lblz(j2)-2-njets)=qf(exch(i))
            idec(i,2,lblz(j2)-2-njets)=qfb(exch(i))
          enddo
          do j1=1,4
            eps(j1,j2)=aux(j1)
          enddo
        enddo

        do j2=1+jz,jwm+jz
          qw(1)=p(4,lblwm(j2-jz))
          do j1=2,4
            qw(j1)=p(j1-1,lblwm(j2-jz))
          enddo
c     subroutine Vpol(mw,pw,ga,gv,pf,pfb,eps)
          call vpol(mw,qw,-1.d0,1.d0,qf,qfb,aux)
          do i=1,4
            idec(i,1,lblwm(j2-jz)-2-njets)=qf(exch(i))
            idec(i,2,lblwm(j2-jz)-2-njets)=qfb(exch(i))
          enddo
          do j1=1,4
            eps(j1,j2)=aux(j1)
          enddo
        enddo
        do j2=1+jwm+jz,jwm+jwp+jz
          qw(1)=p(4,lblwp(j2-jz-jwm))
          do j1=2,4
            qw(j1)=p(j1-1,lblwp(j2-jz-jwm))
          enddo
c     subroutine Vpol(mw,pw,ga,gv,pf,pfb,eps)
          call vpol(mw,qw,-1.d0,1.d0,qf,qfb,aux)
          do i=1,4
            idec(i,1,lblwp(j2-jz-jwm)-2-njets)=qf(exch(i))
            idec(i,2,lblwp(j2-jz-jwm)-2-njets)=qfb(exch(i))
          enddo
          do j1=1,4
            eps(j1,j2)=aux(j1)
          enddo
        enddo
        xlum=xlum*zbrt 
      endif
c
      if(abs(wchrg).ne.2) xlum= xlum*1.d-10

C     
      end                                                
C     
C     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setbr(nz0,zbrcum,zd,zbrt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Input:  NZ0  number of Z bosons
c     ZBRCUM(I) partial sums of branching ratios
c     ZD(I) decay mode of the I-th Z, ZD(I)= n, e, q, b, i == 
c     neutrinos, charged leptons, hadrons, b-quarks, inclusive
c     Output: ZBRT  branching ratio for the chosen decay mode
C     
c     
      implicit none
c     
      include "alpgen.inc"
      integer nz0
      character*1 zd(maxpar-2)
      real*8 zbrcum(4),zbrt
c     
      integer zdprm(maxpar-2,40320) !40320= 8!
      integer zdn(maxpar-2),j1,j2,j3,jtmp,nmx,zdtmp(maxpar-2)
      integer ncases,ndecmod
      parameter(ncases=5,ndecmod=5)
      integer perm(maxpar-2),nperm,map(ncases,ndecmod)
      real*8 br(ncases),zbrtmp
      logical acc
c     
c     n     e     q     b     i
c     n        1     0     0     0     1
c     e        0     1     0     0     1
c     u,c      0     0     1     0     1
c     d,s      0     0     1     0     1
c     b        0     0     1     1     1
c     
      data map/1,0,0,0,0,
     + 0,1,0,0,0,
     + 0,0,1,1,1,
     + 0,0,0,0,1,
     + 1,1,1,1,1/
c     
      do j1=1,nz0
        if(zd(j1).eq.'n') zdn(j1)= 1
        if(zd(j1).eq.'e') zdn(j1)= 2
        if(zd(j1).eq.'q') zdn(j1)= 3
        if(zd(j1).eq.'b') zdn(j1)= 4
        if(zd(j1).eq.'i') zdn(j1)= 5
      enddo
      br(1)= zbrcum(1)
      do j1= 2,4
        br(j1)= zbrcum(j1)-zbrcum(j1-1)
      enddo
      br(4)= 2.d0/3.d0*br(4)
      br(5)= 0.5d0*br(4)
c     
      if (nz0.ne.0) then
        nmx= nz0**nz0-1
      else
        nmx= 0
      endif
      nperm= 0
      do j1=0,nmx
        jtmp= j1
        do j2=1,nz0
          perm(j2)= mod(jtmp,nz0)+1
          do j3=1,j2-1
            if(perm(j2).eq.perm(j3)) goto 10
          enddo
          jtmp= jtmp/nz0
        enddo
        do j2=1,nz0
          zdtmp(j2)= zdn(perm(j2))
        enddo
        do j2=1,nperm
          jtmp= 0
          do j3=1,nz0
            jtmp= zdtmp(j3)-zdprm(j3,j2)
            if(jtmp.ne.0) goto 11
          enddo
        enddo
        if (nperm.ne.0) goto 10
 11     nperm= nperm+1
        do j2=1,nz0
          zdprm(j2,nperm)= zdtmp(j2)
        enddo
 10   continue
      enddo
c     
      zbrt= 0.d0
      nmx= ncases**nz0-1
      do j1=0,nmx
        jtmp= j1
        do j2=1,nz0
          perm(j2)= mod(jtmp,ncases)+1
          jtmp= jtmp/ncases
        enddo
        acc= .false.
        do j2= 1,nperm
          do j3= 1,nz0
            if(map(perm(j3),zdprm(j3,j2)).eq.0) goto 20
          enddo
          acc= .true.
          goto 30
 20       continue
        enddo
 30     if(acc) then
          zbrtmp= 1.d0
          do j2=1,nz0
            zbrtmp= zbrtmp*br(perm(j2))
          enddo
          zbrt= zbrt+zbrtmp
        endif
      enddo
c     
      return
      end
c     
      subroutine setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
      implicit none
      include 'alpgen.inc'
      include 'vbjet.inc'
c     locals
      integer maxdec
      parameter (maxdec=40)
      integer ip, iw, ic, il, irn,id1,id2,i,j1,ires,itmp,ich
     $     ,wdmode,idecmode(8),idec0,rnint(8)
      integer iwfl(2,2)
      integer init
      data init/0/
      integer icab(-4:4,2) 
      double precision xrn,rangen2
      double precision pdec(5),m1,m2
c     cabibbo partner: (*,1)=cabibbo allowed, (*,2)-cabibbo suppressed
c     where 1=d 2=u 3=s 4=c 0=g and negatives are antiparticles
      data icab/-3,-4,-1,-2,0,2,1,4,3, 
     +     -1,-2,-3,-4,0,4,3,2,1/
      integer ilepc(3),inu(3),iup(2),idn(2)
      data ilepc/11,13,15/,inu/12,14,16/,iup/2,4/,idn/1,3/
c     BR's
c     double precision decbr,br(7)
c     data br/3*0.148148148,0.444444444,0.111111111,0.444444444,1d0/
      double precision wbr(6)
      data wbr/3*0.111111d0,0.3333333d0,0.66666666d0,1d0/
      double precision dmass(16)
      data dmass/16*0d0/
c     
      character *3 nth(8)
      data nth/'1st','2nd','3rd','4th','5th','6th','7th','8th'/
c     arguments
      integer nwrt,iflwrt(maxdec),icuwrt(2,maxdec)
      double precision pwrt(5,maxdec),decbr
      integer lblw(maxpar)
C
      integer decfl(maxpar)
c     lblw points the W positions inside ifl (order W+W-...)
c     
c     debug
      integer idbg
      double precision fcount
      common/fldbg/fcount(16),idbg
      data idbg/1/
c     
      save dmass
c     
      do ip=1,npart
        iflwrt(ip)=ifl(ip)
        do ic=1,2
          icuwrt(ic,ip)=icu(ic,ip)
        enddo
        do il=1,5
          pwrt(il,ip)=p(il,ip)
        enddo
      enddo
      nwrt=npart
      decbr=1d0
      if(idecay.ne.'y') return
c     
      if(init.eq.0) then
        if(nw.gt.0) then
c 1        write(*,*)
c     $ 'Input integer string with the decay modes of the individual Ws:'
c          write(*,*) '1: e nu'
c          write(*,*) '2: mu nu'
c          write(*,*) '3: tau nu'
c          write(*,*) '4: e/mu/tau nu '
c          write(*,*) '5: q qbar'' '
c          write(*,*) '6: fully inclusive'
c          write(*,*) 'E.g.: input 24 for WW -> mu l nu_mu nu_l'
c          write(*,*) 'E.g.: input 45 for WW -> l nu q qbar'' '
c          write(*,*)
c     $      'E.g.: input 126 for WWW -> e nue mu nu(mu) ff''    bar'
c          read(*,*) itmp
          itmp=iWdecmode
          j1= 0
          do while(itmp.ne.0) 
            j1= j1+1
            ires= mod(itmp,10)
            idecmode(j1)= ires
            itmp= itmp/10
          enddo
          if(j1.ne.nw) then
            write(*,*) 'W decay mode string incorrect, stop'
            stop
          endif
        endif
        dmass(4)=1.5
        dmass(5)=mb
        dmass(11)=0.5d-3
        dmass(13)=0.10566d0
        dmass(15)=1.777d0
        init=1
      endif

c     Z final states
      do ip=1,nz
        nwrt=nwrt+2
        do il=1,4
          pwrt(il,nwrt-1)=idec(il,1,ip+nph)
          pwrt(il,nwrt  )=idec(il,2,ip+nph)
          pdec(il)=pwrt(il,nwrt-1)+pwrt(il,nwrt)
        enddo
        iflwrt(nwrt-1)=zfl(ip)
        iflwrt(nwrt)=-zfl(ip)
        m1=dmass(abs(iflwrt(nwrt-1)))
        m2=dmass(abs(iflwrt(nwrt)))
        pwrt(5,nwrt-1)=0
        pwrt(5,nwrt)=0
        pdec(5)=mz
c     rescale momenta to incorporate mass effects in the decays
        call rescms(pdec,pwrt(1,nwrt-1),pwrt(1,nwrt),m1,m2)
        if(abs(iflwrt(nwrt-1)).le.6) then
          icuwrt(1,nwrt-1)=nwrt-1
          icuwrt(2,nwrt)=nwrt-1
        else
          icuwrt(1,nwrt-1)=0
          icuwrt(2,nwrt)=0
        endif
        icuwrt(2,nwrt-1)=0
        icuwrt(1,nwrt)=0
      enddo
c     
c     W final states
c      do i= nz+1,nw+nz
c         lblw(i-nz)= i+2
c      enddo

c     randomize W ordering
      if(nw.gt.0) then
!        call rndint(nw,rnint)
        call wdec(nw,idecmode,decfl,decbr)
        do ip=nph+nz+1,nph+nz+nw
          iw=ip-nz-nph
!          decbr=decbr*wbr(idecmode(iw))
          nwrt=nwrt+2
!          idec0=idecmode(iw)
          if(ifl(2+ip+njets).eq.24) then          
            if (decfl(iw).le.3) then
              iflwrt(nwrt-1)=10+2*decfl(iw)
              iflwrt(nwrt)=-iflwrt(nwrt-1)+1
            else
              iflwrt(nwrt-1)=2*(decfl(iw)-3)
              iflwrt(nwrt)=-iflwrt(nwrt-1)+1
            endif
!            iflwrt(nwrt-1)=12
!            iflwrt(nwrt)=-11
          elseif(ifl(2+ip+njets).eq.-24) then
            if (decfl(iw).le.3) then
              iflwrt(nwrt-1)=9+2*decfl(iw)
              iflwrt(nwrt)=-iflwrt(nwrt-1)-1
            else
              iflwrt(nwrt-1)=2*(decfl(iw)-3)-1
              iflwrt(nwrt)=-iflwrt(nwrt-1)-1
            endif
!             iflwrt(nwrt-1)=11
!            iflwrt(nwrt)=-12
          else
            write(*,*) 'something wrong with W code:',ifl(2+ip+njets)
          endif
          do j1=nwrt-1,nwrt
            if(abs(iflwrt(j1)).gt.10) then
              icuwrt(1,j1)=0
              icuwrt(2,j1)=0
            else
              if(iflwrt(j1).gt.0) then
                icuwrt(1,j1)=nwrt-1
                icuwrt(2,j1)=0
              else
                icuwrt(1,j1)=0
                icuwrt(2,j1)=nwrt-1
              endif
            endif
          enddo
!!!!!
          do il=1,4
            pwrt(il,nwrt-1)=idec(il,1,ip)
            pwrt(il,nwrt  )=idec(il,2,ip)
            pdec(il)=pwrt(il,nwrt-1)+pwrt(il,nwrt)
          enddo
          m1=dmass(abs(iflwrt(nwrt-1)))
          m2=dmass(abs(iflwrt(nwrt)))
          pwrt(5,nwrt-1)=0
          pwrt(5,nwrt)=0
          pdec(5)=mw
c     rescale momenta to incorporate mass effects in the decays
          call rescms(pdec,pwrt(1,nwrt-1),pwrt(1,nwrt),m1,m2)
        enddo
      endif
c$$$      do i=nwrt-1,nwrt
c$$$        do il=1,16
c$$$          if(abs(iflwrt(i)).eq.il) fcount(il)=fcount(il)+1d0
c$$$        enddo
c$$$      enddo
      end
*

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phspace(lnot,pswgt,djpd,djg)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'alpgen.inc'
      include 'vbjet.inc'
c      real *8 djbin,djbintot,dummy,djg
c      integer nvar,nbin,jprocmax,jproc,nv,ngrid,jgrid,maxsubproc
c      common/psopt/nvar,nv
c      parameter (maxsubproc= 100)
c      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c      integer nct,nx1,nx2,maxn
c      parameter (maxn= 100)
c      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      real *8 dummy,djg
      real *8 pswgt
      real *8 djpd,factor
      real *8 cutkin(12)
      real *8 wgt
      common/loccut/cutkin
      real *8 pl(maxpar),y(maxpar)
      real *8 p1(0:3,maxpar)
*     
c-    debugging variables
*     
      real *8 totpt
      integer nx,ninit,i,l,m,lnot
      parameter (nx= 20)
      data ninit/0/
      save ninit
      if(ninit.eq.0) then
*     
c-setup local generation cuts
*     
        cutkin(1)=ptjmin
        cutkin(3)=etajmax
        cutkin(5)=drjmin
        cutkin(6)=etagap
        cutkin(7)=ptcen
        cutkin(8)=ptphmin
        cutkin(9)=drphjmin
        cutkin(11)=drphmin
        cutkin(10)=etaphmax
        ninit=1
      endif
*     
c-    The generation starts
*     
      pswgt=0.d0     
*
      call momgen(njets,nph,nz,nw,nh,irapgap,roots,xm,x1,x2,
     +     p1,wgt,lnot)
      djg= 1.d0                 ! dummy variable
      if (lnot.eq.1) then
        pswgt= 0.d0
        goto 100
      endif
*     
c-    will write factor=factor0/(x1*x2), with factor0 function of
c-    njets etc.
*     
      factor= 1d0/(2.d0*pi)**(3*(nvb+nh+njets+nph)-4)/2.d0/s/x1/x2
*     
c-    Momenta in the LAB frame:
*     
      do l=1,4
        do m= npart,3,-1
          if (l.le.3) p(l,m)= p1(l,m-2)
          if (l.eq.4) p(l,m)= p1(0,m-2)
        enddo
      enddo
*     
      p(4,1)= roots/2.d0*x1   
      p(1,1)= 0.d0   
      p(2,1)= 0.d0   
      p(3,1)= roots/2.d0*x1   
*     
      p(4,2)= roots/2.d0*x2   
      p(1,2)= 0.d0   
      p(2,2)= 0.d0   
      p(3,2)=-roots/2.d0*x2   
comment
c      do l= 1,npart
c        print*,'osh',l,'=',p(4,l)**2-p(1,l)**2-p(2,l)**2-p(3,l)**2
c      enddo     
c
c      do l= 1,4
c        dummy= -p(l,1)-p(l,2)
c        do m= 3,npart
c          dummy= dummy+p(l,m)
c        enddo
c        print*,'summ',l,'=',dummy
c      enddo
c      print*,'   '
comment
*     
      if (njets+nph.eq.0) goto 50
*     
c-    Rapidities and pseudo-rapidities (in the LAB system), pt's,
c     deltar's:
*     
      do l= 3,njets+2+nph
        pt(l) = sqrt(p(1,l)**2+p(2,l)**2)
        pl(l) = p(3,l)
        eta(l)= -log(tan(0.5d0*atan2(pt(l),pl(l))))
        y(l)  = 0.5d0*log((p(4,l)+p(3,l))/(p(4,l)-p(3,l)))
      enddo
*     
c-    Calculates jet-jet distances and jet-photon:
*     
      do l= 3,njets+1+nph
        do m= l+1,njets+2+nph
          if(min(pt(l),pt(m)).gt.0d0) then
            dphi(l,m)= 
     +           (p(1,l)*p(1,m)+p(2,l)*p(2,m))
     +           /pt(l)/pt(m)
            if(dabs(dphi(l,m)).gt.1.d0) then
c$$$  write(*,*) 'partons',l,m,', cos(Dphi)=', dphi(l,m)
c$$$  .                 ,', set to +-1'
c$$$  write(*,*) 'pt(',l,')=',pt(l),'p(',l,')=',  (p(i,l)
c$$$  .                 ,i=1,4)
c$$$  write(*,*) 'pt(',m,')=',pt(m),'p(',m,')=',  (p(i,m)
c$$$  .                 ,i=1,4)
              if (dphi(l,m).gt.0.d0) dphi(l,m)= 1.d0
              if (dphi(l,m).lt.0.d0) dphi(l,m)=-1.d0
            endif
            dphi(l,m)=acos(dphi(l,m))
            dphi(m,l)=dphi(l,m)
          else
c***  fix***
c     dphi(l,m)=pi
c     dphi(m,l)=pi
            goto 100
          endif
        enddo
      enddo
*     
      do l= 3,njets+1+nph
        do m= l+1,njets+2+nph
          dr(l,m)= sqrt(dphi(l,m)**2+(eta(l)-eta(m))**2)
          dr(m,l)=dr(l,m)
        enddo
      enddo
*     
      call chkcut(lnot,pt,p,eta,dr,nvb,njets,nph)
      if (lnot.eq.1) then
        pswgt= 0.d0
        goto 100
      endif
*     
 50   continue
*     
      pswgt = factor*wgt
*     
c-    Evaluate q2: will include several possible options:
*     
      totpt= 0.d0
*     
c-    Sum of masses
*     
      do i= 1,nvb+nh
        totpt=totpt + xm(i+njets+nph)
      enddo
      if(iqopt.eq.0) then 
         qsq=1d0
      elseif(iqopt.eq.1) then
        qsq=totpt**2
         do l= 1,njets
           qsq=qsq+pt(l+2)**2
         enddo
      elseif(iqopt.eq.2) then
         qsq= (totpt)**2
      elseif(iqopt.eq.3) then
         qsq= x1*x2*roots**2
      endif
      qsq=qfac**2*qsq
 100  continue
      end
*
      subroutine chkcut(lnot,pt,p,eta,dr,nvb,njets,nph)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     c
c     Applies kinematical cuts to the final state during         c 
c     the phase-space generation                                 c
c     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ptjmin,drbmin,drjmin,etabmax,etajmax,ptbmin,ptphmin
      real*8 ptlmin,drlmin,etalmax,ptnmin,drbmax,etaphmax,drphjmin
     $     ,drphmin
      integer maxpar,ninit,njets,j,lnot,i,nvb,nph
      real*8 cutkin(12)
      common/loccut/cutkin
      data ninit/0/
      parameter (maxpar=10)
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),p(5,maxpar)
      real*8 pbjet(3),ptbjet,ptb,ptbb,etab,etabb,drbb,etabjet
      real *8 ib,ibb,sjj,sjjmin,etajmin
      real *8 mb2,dummy,ran(1:2),pa(0:3,5:6),pb(1:4,3:6),p56(0:3)
      integer l,m,ord(maxpar),jord(maxpar),nn,iopt,l3,l4,l5,l6
      real*8 pc(5,maxpar),wtdec,summ(4),dphic(maxpar,maxpar)
      real*8 ptc(maxpar),plc(maxpar),etac(maxpar),etc
      real*8 drc(maxpar,maxpar),mjj,mw,mz,s34,s35,s45,ran0
      common/ptanal/mjj,ptc,etac
      integer itmp
      data itmp/0/
      save ninit,ptjmin,ptbmin,etajmax,etabmax,drjmin,drbmin,ptlmin
     $     ,etalmax,drlmin,ptnmin 
*     
      if(ninit.eq.0) then
        ninit=1
        ptjmin = cutkin(1)
        etajmax= cutkin(3)
        drjmin = cutkin(5)
        ptphmin = cutkin(8)
        etaphmax= cutkin(10)
        drphjmin = cutkin(9)
        drphmin = cutkin(11)
      endif
*     
      lnot= 0 
*     
c     impose minimum pt and require eta within allowed range, for jets
*     
      do i= 3,njets+2
        if (pt(i).lt.ptjmin)          goto 10
        if (abs(eta(i)).gt.etajmax)   goto 10
      enddo 
      do i= njets+3,njets+2+nph
        if (pt(i).lt.ptphmin)          goto 10
        if (abs(eta(i)).gt.etaphmax)   goto 10
      enddo 
*
c     dR(phot-phot)<drphmin
*
      do i= njets+3,njets+nph+1
        do j= i+1,njets+nph+2
          if(dr(i,j).lt.drphmin)         goto 10
        enddo
      enddo
*
c     dR(phot-jet)<drphjmin
*
      do i= njets+3,njets+nph+2
        do j= 3,njets+2
          if(dr(i,j).lt.drphjmin)         goto 10
        enddo
      enddo
*     
c     dR(jet-jet)<drjmin
*     
      do i= 3,njets+1
        do j= i+1,njets+2
          if(dr(i,j).lt.drjmin)       goto 10
        enddo
      enddo
*
c     dR(jet-photon)<drphjmin
*
      do i= 3,njets+nph+1
        do j= i+1,njets+nph+2
          if (i.le.(njets+2).and.dr(i,j).lt.drphjmin) goto 10
        enddo
      enddo

      return 
*     
 10   lnot= 1
      return
      end
*     
      subroutine momgen(njets,nph,nz,nw,nh,irapgap,roots,xm1,x1,x2,p,
     +     wgt,lw)
*     
c-    Generator of np= njets+nz+nw+nh particles in the LAB frame.
*     
c     Here the ordering of the final state particles is 
c     
c     jet(1)     ... jet(njets) 
c     phot(1)... phot(nph)
c     Z(njets+1) ... Z(njets+nz) 
c     W+(njets+nz+1)   ... W+-(njets+nz+nw)      
c     H(njets+nz+nw+1)  ... H(njets+nz+nw+nh)      
c     
*     
      implicit none
      real*8 roots,x1,x2,wgt,wgt1,wgt2,cutkin,wjr,wtau
      real*8 etajmax,eta0min,eps,etacut,ptjmin,ptphmin
      real*8 ag,tau0,tau,etagap,ptcen,etaphmax
      real*8 ranram,sq,y0,yr,pt0lmax,pt0lsum
      real*8 xmsum,p,ranch0
      real*8 ptlim,s,en,pz
      integer lw,npar,mpar,np,njets,nz,nw,nh,nvb,j,lw1,nph
      parameter (npar= 100)
      real*8 pt0(npar),pt1(npar),eta0(npar),eta1(npar),
     +     xmr(npar)
      real*8 pm(0:4,npar)
      real*8 pt0_l(npar),pt1_l(npar),xm(npar),xm1(npar)
      dimension ranram(111)
      common/loccut/cutkin(12)
      dimension p(0:3,npar) 
      real *8 djbin,djbintot,djb0,djb1,dummy,djg,dj,apw(1:2)
      integer nvar,nbin,jprocmax,jproc,nv,ngrid,jgrid,maxsubproc
      common/psopt/nvar,nv
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      integer nct,nx1,nx2,maxn
      parameter (maxn= 100)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      integer mask,mmask,md,lmin,lmax,m,l,ndummy  
      real*8 peropt,wsymm,ran0 
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      data mpar/0/
*     
c-    irapgap=0: symmetrization NOT PERFORMED on jets 
c-    irapgap=1: symmetrization   PERFORMED   on jets
c-
c-    Note that the compensating statystical factor are computed
c-    for 
c-    NON OVERLAPPING REGIONS in eta_jets(i)
c-     
*     
      integer irapgap
      integer kp(10,6),jj
      data kp/1,2,3,4,5,6,7,8,9,10,
     +     2,3,1,4,5,6,7,8,9,10,
     +     3,1,2,4,5,6,7,8,9,10,
     +     1,3,2,4,5,6,7,8,9,10,
     +     3,2,1,4,5,6,7,8,9,10,
     +     2,1,3,4,5,6,7,8,9,10/
      save
*     
      if (mpar.eq.0) then
        mpar   = 1
        nvb    = nz+nw
        np     = nvb+nh+njets+nph
        eps    = 1.d-5
        etacut = 40.d0
        s      = roots*roots
*     
c-      here generation cuts can be inserted:
*     
        ptjmin = cutkin(1)
        etajmax= cutkin(3)
        etagap = cutkin(6)
        ptcen  = cutkin(7)
        ptphmin=cutkin(8)
        etaphmax=cutkin(10)
*     
c-      eta1 < eta < eta0
*     
        do j= 1,njets
          pt0_l(j) = ptjmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etajmax
          eta1(j)  =-etajmax
          xm(j)    = xm1(j)
        enddo
        do j= 1+njets,njets+nph
          pt0_l(j) = ptphmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etaphmax
          eta1(j)  =-etaphmax
          xm(j)    = xm1(j)
        enddo
        if (irapgap.eq.1) then
          if    (njets.eq.2) then
            eta0(1)  = etajmax
            eta1(1)  = etagap
            eta0(2)  =-etagap
            eta1(2)  =-etajmax
            wsymm    = 2.d0
          elseif(njets.eq.3) then
            eta0(1)  = etajmax
            eta1(1)  = etagap
            eta0(2)  =-etagap
            eta1(2)  =-etajmax
            eta0(3)  = etagap
            eta1(3)  =-etagap
            pt0_l(3) = ptcen
            wsymm    = 6.d0
          else
            print*,'OPTION DO NOT ALLLOWED'
            stop
          endif
        else
          wsymm= 1.d0
        endif 
*     
        do j= 1,nvb+nh
          pt0_l(njets+nph+j) = eps
          pt1_l(njets+nph+j) = roots/2.d0
          eta0(njets+nph+j)  = etacut
          eta1(njets+j+nph)  =-etacut
          xm   (njets+j+nph) = xm1(njets+j+nph)
        enddo
*     
c-      Parameters for the x1,x2 integration:
*     
        ag  = 0.98d0
*     
c-      tau0 is the lower cut on tau= x1*x2
*     
        xmsum= 0.d0
        pt0lmax= 0.d0
        pt0lsum= 0.d0
        eta0min= 100.d0
        do j= 1,np
          xmsum  = xmsum+xm(j)
          pt0lmax= max(pt0lmax,pt0_l(j))
          pt0lsum= pt0lsum+pt0_l(j)
          eta0min= min(eta0min,eta0(j))
        enddo
        tau0   = 1.d0/s*(xmsum)**2
        tau0   = max(tau0,4.d0*pt0lmax**2/s)  
        tau0   = max(tau0,pt0lsum**2/s)
*
c-      sets the ratios of calls mom/momr
*
        apw(1)= 0.9d0
        apw(2)= 1.d0-apw(1)
      endif
*
      lw= 0
      lmin= (jgrid(jproc)-1)*nv+2
      lmax= jgrid(jproc)*nv+1
*
      dj= 0.d0
      call rans(ran0)
      if (ran0.le.apw(1)) then
         md= 1
      else
         md= 2 
      endif
*
      djb1= 1.d0
      djb0= 1.d0
      m= 0
      do l= lmin,lmax
        m= m+1
        if (md.eq.1) then
          mmask(l)= 1
        else
          mmask(l)= mask(l)
        endif
        if (mmask(l).eq.1) then
          call onedimbin(1,nbin,djbin,l,ndummy,dummy)
          call grans(nbin,nx1(l),ranram(m))
          if (mask(l).eq.1) then
            djb1= djb1*djbin
          elseif (mask(l).eq.0) then
            djb0= djb0*djbin
          endif
        else
          ranram(m)= 0.d0
        endif
      enddo
      djbintot= djb1
      if (md.eq.1) djbintot= djbintot*djb0 
*     
c-    np.eq.1:
*     
      if (np.eq.1) then
        tau= xm(1)**2/roots/roots
        y0 =-0.5d0*log(tau) 
        call ttriangle(0,-y0,y0,-y0,y0,yr,wjr,ranram(1),lw1)
        sq = sqrt(tau)
        x1 = sq*exp(yr) 
        x2 = sq*exp(-yr) 
*     
        p(0,1)= roots/2.d0*(x1+x2)
        p(1,1)= 0.d0
        p(2,1)= 0.d0
        p(3,1)= roots/2.d0*(x1-x2)
        wgt= wjr/roots/roots/djbintot        
        return
      endif
*     
c-    Generates x1 and x2: 
*     
      call ppeaka(0,tau0,ag,tau0,1.d0,tau,wtau,ranram(1),lw1)
      sq = sqrt(tau)
      y0 = -0.5d0*log(tau) 
      call ttriangle(0,-y0,y0,-y0,y0,yr,wjr,ranram(2),lw1)

      x1 = sq*exp(yr) 
      x2 = sq*exp(-yr) 
*     
c-    Protection:
*     
      if (roots*sq.lt.xmsum) goto 100
*     
c-    Rescalings and transformations to feed mom:
*     
      do j= 1,np
        ptlim = s*s
     +       +(xm(j)**2-(xmsum-xm(j))**2)**2
     +       -2.d0*s*(xm(j)**2+(xmsum-xm(j))**2)   
        if (ptlim.le.0.d0) goto 100
        ptlim =  0.5d0/roots*sqrt(ptlim)
        pt0(j)= pt0_l(j)/roots
        pt1(j)= min(pt1_l(j),ptlim)/roots
*     
c-      Protection:
*     
        if (pt0(j).gt.pt1(j)) goto 100
        xmr(j)= xm(j)/roots
      enddo
*     
c-    momenta:
*     
      en    = 0.5d0*(x1+x2)     ! Rescaled initial total energy
      pz    = 0.5d0*(x1-x2)     ! Rescaled initial longitudinal momentum
*     
      if (md.eq.1) then
        call mom(0,np,pt0,pt1,eta0,eta1,xmr,en,pz,pm,wgt1,ranram,lw)
        if (lw.eq.0) then
          call momr(1,irapgap,njets,etagap,
     +              np,pt0,pt1,eta0,xmr,en,pz,pm,wgt2,lw)
          if (lw.eq.0) then
            if (wgt1.eq.0.d0.and.wgt2.eq.0.d0) goto 100
            wgt= wgt1*wgt2/
     +          (wgt2*apw(1)*djbintot+wgt1*apw(2)*djb1)
          else
            wgt= wgt1/(apw(1)*djbintot)
          endif
        else
          goto 100
        endif
      else
        call momr(0,irapgap,njets,etagap,
     +            np,pt0,pt1,eta0,xmr,en,pz,pm,wgt2,lw)
        if (lw.eq.0) then
          call mom(1,np,pt0,pt1,eta0,eta1,xmr,en,pz,pm,wgt1,ranram,lw)
          if (lw.eq.0) then
            m= 0
            do l= lmin,lmax
              m= m+1
              if (mask(l).eq.0) then
                nbin= min(nx1(l),int(1.d0+dfloat(nx1(l))*ranram(m)))
                call onedimbin(3,nbin,djbin,l,ndummy,dummy)
                djb0= djb0*djbin
              endif
            enddo
            if (wgt1.eq.0.d0.and.wgt2.eq.0.d0) goto 100
            wgt= wgt1*wgt2/
     +          (wgt1*apw(2)*djbintot+wgt2*apw(1)*djb0*djb1)
          else
            wgt= wgt2/(apw(2)*djbintot)
          endif
        else
          goto 100
        endif
      endif
*     
c-    from the x1,x2 integration
*     
      wgt= wgt*wtau*wjr     
*     
c-    Rescaling of momenta and weight and
c-    symmetrization performed if required
*     
      if (irapgap.eq.1) then
        call rans(ranch0)
        if (njets.eq.2) then
          if    (ranch0.lt.1.d0/2.d0) then
            jj= 1
          elseif(ranch0.lt.2.d0/2.d0) then
            jj= 6
          else
            print*,'ERROR IN MOMGEN'
            stop
          endif
        endif
        if (njets.eq.3) then
          if    (ranch0.lt.1.d0/6.d0) then
            jj= 1 
          elseif(ranch0.lt.2.d0/6.d0) then
            jj= 2 
          elseif(ranch0.lt.3.d0/6.d0) then
            jj= 3 
          elseif(ranch0.lt.4.d0/6.d0) then
            jj= 4 
          elseif(ranch0.lt.5.d0/6.d0) then
            jj= 5 
          elseif(ranch0.lt.6.d0/6.d0) then
            jj= 6 
          else
            print*,'ERROR IN MOMGEN'
            stop
          endif
        endif
      else
        jj= 1
      endif
      do j= 1,np
        p(0,j)= pm(0,kp(j,jj))*roots
        p(1,j)= pm(1,kp(j,jj))*roots
        p(2,j)= pm(2,kp(j,jj))*roots
        p(3,j)= pm(3,kp(j,jj))*roots
      enddo
*     
      wgt= wgt*s**(np-2)*wsymm
*     
      return
 100  lw= 1
      wgt= 0.d0
      return
      end
*     
      integer function ifact(n)
      ifact= 1
      do i= 1,n
        ifact= ifact*i
      enddo
      return
      end
*
c-------------------------------------------------------------------
      subroutine alsgrd
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'vbjet.inc'
      integer nout,nct1
      integer nvar,nch,n,nv
      integer init,j,k
      real*8 v1,ni
      real*8 al1,bet1
      integer nct,nx1,nx2
      integer maxn,ncmax
      parameter (maxn= 100,ncmax= 1000)  
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      common/ausil/init(ncmax),bet1(0:ncmax,maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/psopt/nvar,nv
      common/book/v1(ncmax),ni(maxn)
      data v1/ncmax*0d0/,ni/maxn*0/
      integer mask,mmask 
      real*8 peropt
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      data mask/maxn*0/,mmask/maxn*1/,peropt/maxn*0.9d0/
*     
c-    Protection:
*     
      if (nvb+nh+nph.eq.0) then
        print*,'nz+nw+nh+nph=  0 is not allowed!'
        stop
      endif
*     
c-    initialise size of grids:
*     
c-    for MC over jproc:

      nct(1)   = 0              ! first bin of variable 1
      nx1(1)   = jprocmax       ! number of jprocs
*     
c-    for the Phase-space reweighting:
*     
      if (nvb+njets+nh+nph.eq.1) then
        nv  = 1
      else 
        nv  = 2*(npart-2)-1
      endif
      nvar = nv*ngrid    
      do n= 2,nvar+1
        nch= 10
*     
        nct(n)   =  nct(n-1)+nx1(n-1)   
        nx1(n)   =  nch 
*
c-      set the channels shared by mom and momr: 
*
        if (mod(n-1,nv).eq.1)    mask(n)= 1
        if (mod(n-1,nv).eq.2)    mask(n)= 1
        if (mod(n-1,nv).eq.0.and.nv.eq.1) mask(n)= 1
*
c-      set the percentages of optimization: 
*
        if (mask(n).eq.1) peropt(n)= 0.9d0 
        if (mask(n).eq.0) peropt(n)= 1.d0-dfloat(nv-1)/200.d0
      enddo
      do n= 1,nvar+1
        nct1(n)= nx1(n)
        do j= 1,nct1(n)                                         

          init(nct(n)+j)= 0                                     

c-        the weights are initialized to be all equal 

          al1(j,n)= 1.d0/nct1(n)                                  
        enddo

        do j= 1,nct1(n)                                         
          bet1(j,n)= 0                                         
          do k= 1,j                                          
            bet1(j,n)= bet1(j,n)+al1(k,n)                          
          enddo                                              
        enddo                                                
      enddo
*     
c-    protection:
*     
      if (nct(nvar+1)+nct1(nvar+1).gt.ncmax) then
        print*,'INCREASE NCMAX',nct(nvar+1)+nct1(nvar+1)
        stop
      endif
*     
      if ((nvar+1).gt.maxn) then
        print*,'INCREASE MAXN'
        stop
      endif
      end
*
c-------------------------------------------------------------------
      subroutine dumpwgt
c-------------------------------------------------------------------
c     routine for the debugging of individual large-weight events
c     Not functional to the rnning of the code
      implicit none
      include 'alpgen.inc'
      include 'vbjet.inc'
      real *8 tmpmw,ptmp(4)
      integer i
c     write(*,*) 'x1,x2=',x1,x2
c     write(*,*) 'p(b)=',(p(i,3),i=1,3)
c     write(*,*) 'p(bbar)=',(p(i,4),i=1,3)
c     write(*,*) 'p(j1)=',(p(i,5),i=1,3)
c     write(*,*) 'p(j2)=',(p(i,6),i=1,3)
c     write(*,*) 'p(e)=',(p(i,7),i=1,3)
c     write(*,*) 'p(nu)=',(p(i,8),i=1,3)
      write(niosta,*) 'x1,x2=',x1,x2
      end

C***********************************************************************
      subroutine matrix(flvmlm,posi,impul,hel,rst,labcol,wgcol,colfl)
C***********************************************************************
C     
C     Input parameters: FLVMLM(NMAX) particles according
C     to michelangelo convention, POSI(NMAX) positions of
C     incoming particles IMPUL(4,8) particles four momenta,
C     HEL, array of particle helicities 
C     LABCOL to choose whether su(3) mode only (LABCOL=0) or
C     color flow is required (LABCOL=1) or dual mode only
C     (LABCOL=2), WGCOL random number for color flow unweighting
C     Output parameters:   RESULT squared matrix element (modulus),
C     COLFLOW 
C     string representing the selected coulor flow 
C     
C     
C     At present:
C     
C     the color flow is returned in the following way: each particle is
C     represented by a pair of integers (n1,n2); two particle (n1,n2),
C     (n3,n4)
C     are colour connected if either n2=n3 or n4=n1 (in this case the
C     coulor
C     flows from 2 to 1). COLFLOW will contain a string of pairs of
C     integers 
C     ordered as follows: d (or dbar) nu_ebar bbar ubar (or u) e^- b glu
C     glu
C     (with gluons ordered according to the ordering of momenta).
C     (COLFLOW=(n1,....,nn,....) only the first nn=2*particles_number
C     elements to be used)   
C     
C     
C     IMPUL(J,NPART) ===   J=1 Energy of the NPART-th particle
C     IMPUL(J,NPART) ===   J=2,3,4 x,y,z components  of the NPART-th
C     particle 
C     three momentum
C     the order of the momenta is: 
C     dbar nubar bbar u e- b ( glu glu glu)  (processes 1,2,3,4)
C     dbar nubar cbar bbar u e- c b (glu)  (processes > 5 )
C     where all particles are assumed outcoming and incoming gluons
C     must be before outcoming ones
C     
      implicit none
C     
      integer nmax              !maximum number of external particles, 
      parameter (nmax=10)
      integer nlb
      parameter (nlb=4)    
      real*8 impul(4,nmax),wgcol !the contribution to the integral
      real*8 rst
c     complex*16 rst
      integer labcol,colfl(2*nmax),flvmlm(nmax)
      integer posi(2)
      integer hel(nmax)         !particles helicities
C     
      integer flgdual           !dual (0) or su3 (1) amplitudes
      common/dual/flgdual
      integer color(2*nmax)     !coulor string
      integer colst(2*nmax)
      common/colore/color
      integer ncls,nprc,ant4(4,2),ant6a(6,6),ant2(2,1)
      integer nant,nantl,j3
      parameter (ncls=2)        !different class of processes
      integer nx,proc(nmax),antns(6),ndual,j1,j2
      parameter (ndual=41000)
      real*8 dualamp(ndual),colfac,damp,dampref,avgspin
      integer colaux(ndual,2*nmax),ndl,ndla(0:10,1:3),class
      integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
      common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
      integer flginit
      data flginit/0/
      integer fq(nmax),pq(nmax)
      data ndla/0,0,1,5,23,119,719,5039,40319,0,0,
     + 0,1,5,23,119,719,5039,0,0,0,0,
     + 0,2,11,59,359,6*0/

      save ndla,ant4,ant2,ant6a,flginit
      complex*16 result
      integer inpint(1000)
      common/initinter/inpint  
C     
      data ant2  /1,2/
      data ant4  /1,4,2,3,
     + 1,3,2,4/
      data ant6a  /1,5,2,6,3,4,
     + 1,6,2,4,3,5,
     + 1,6,2,5,3,4,
     + 1,4,2,6,3,5,
     + 1,4,2,5,3,6,
     + 1,5,2,4,3,6/         
c     
      data inpint/ 16,          ! N V-A 
     + 1,1,1,1,1,   1,1,2,1,2,   2,1,1,1,2,   3,1,2,1,1, ! zuu, zdd, w+ud, w-ud,
     + 10,1,1,1,1,  10,1,2,1,2,   1,2,1,2,1,   1,2,2,2,2, ! Auu, Add, zcc, zss
     + 2,2,1,2,2,   3,2,2,2,1,  10,2,1,2,1,  10,2,2,2,2, ! W+cs, W-sc, Acc, Ass
     + 11,1,1,1,1,  11,1,2,1,2,  11,2,1,2,1,  11,2,2,2,2, ! guu, gdd, gcc, gss
     + 0,                   ! N of yukawa
     + 15,                  ! N self-gauge
     + 1, 2, 3,  10, 2, 3,   4, 2, 3,   5, 1, 1, ! ZWW, AWW, auxiliary
     + 5, 1,10,   5,10,10,   5, 2, 3,   6, 1, 2, ! auxiliary
     + 6,10, 2,   7, 1, 3,   7,10, 3,   8, 2, 2, ! auxiliary
     + 9, 3, 3,  11,11,11,  12,11,11, ! auxiliary, ggg, Xgg
     + 4,                   ! N H-GAUGE
     + 1, 1, 1,   1, 2, 3,   3,1,1,   3,2,3, ! HZZ, HWW, HHZZ, HHWW
     + 3,                   ! N of self higgs
     + 1, 2, 3,             ! HHH, HHHH e aux per higgs quartici
     + 855*-100       /     ! EOF signal
c     
c     data inpint/ 6,
c     >               11,1,1,1,1,  11,1,2,1,2,  11,3,2,3,2,   2,1,1,1,2,  ! guu, gdd, gbb, w+ud
c     >                3,1,4,1,3,  11,2,1,2,1,                            ! w-en, gcc
c     >                0,                                                 ! N of yukawa
c     >                2,                                                 ! N self-gauge
c     >                11,11,11,  12,11,11,                               ! ggg Auxgg
c     >                961*-100/
c     
      if (labcol.eq.0) then
C     
        flgdual=1
        call matrix0(flvmlm,posi,impul,hel,result)
c     
      elseif (labcol.eq.1) then
c     
        j1=0
        do j3=1,nmax
          if (abs(flvmlm(j3)).gt.0.and.abs(flvmlm(j3)).le.6) then
            j1=j1+1
            fq(j1)=flvmlm(j3)
            pq(j1)=j3
          endif
        enddo
        if(j1.ne.nqrk) then
          write(*,*)'something wrong in matrix, NQRK not ok',nqrk,j3
        endif
c     
        flgdual=0
        do j1=1,2*nmax
          colst(j1)=color(j1)
        enddo 
c     
        do j1=1,nmax
          proc(j1)=abs(rep(j1))
        enddo
c     
        ndl=0
        if(nqrk.eq.2) then
          nantl=1
        elseif(nqrk.eq.4) then
          nantl=2
        elseif(nqrk.eq.6) then
          nantl=6
        else
          write(*,*)'Cannot deal with',nqrk,'quarks'
          stop
        endif
        do j3=1,nantl
          do j1=1,nqrk
            if(nqrk.eq.2) then
              antns(j1)=pq(ant2(j1,j3))
            elseif(nqrk.eq.4) then
              antns(j1)=pq(ant4(j1,j3))
            elseif(nqrk.eq.6) then
              antns(j1)=pq(ant6a(j1,j3))
            endif
          enddo
          j2=nqrk/2
          do j1=0,ndla(nglu,j2)
            nx=j1
            if(nqrk.eq.2) then
              call gendual2q(proc,antns,nglu,colst,nx,colfac) 
            elseif(nqrk.eq.4) then
              call gendual4q(proc,antns,nglu,colst,nx,colfac) 
            elseif(nqrk.eq.6) then
              call gendual6q(proc,antns,nglu,colst,nx,colfac) 
            endif 
            if(abs(colfac).gt.1.d-20) then
              ndl=ndl+1
              call matrix0(flvmlm,posi,impul,hel,result)
              dualamp(ndl)=abs(result*colfac)**2
              do j2=1,2*nmax
                colaux(ndl,j2)=color(j2)
              enddo
            endif
          enddo
        enddo
c     
        damp=0
        do j1=1,ndl
          damp=damp+dualamp(j1)
        enddo
        dampref=damp*wgcol
        j1=0
        damp=0.
        do while(damp.lt.dampref.and.j1.lt.ndl)
          j1=j1+1
          damp=damp+dualamp(j1)
        enddo
        ndl=j1
        if(ndl.eq.0) then
          do j1=1,2*nmax
            colfl(j1)=colst(j1)
          enddo
        else 
          do j1=1,2*nmax
            colfl(j1)=colaux(ndl,j1)
          enddo
        endif
C     
      else
        write(6,*)'wrong LABCOL in matrix',labcol
        stop
      endif
C     
      rst=abs(result)**2
c     rst=result
C     
      return
      end



      subroutine translation(jj,npar)
      implicit none
      integer i,j,jj,npar
      character*3 par,p
      dimension par(-4:4)
      dimension jj(10),p(10)
      data par/' cb',' sb',' ub',' db',' gl','  d','  u','  s','  c'/
*     
      do i= 1,npar
        p(i)= par(jj(i))
      enddo     
*     
      write(*,*) (p(j),j= 1,2),'  -> ',(p(j),j= 3,npar)
      return
      end
*     
      subroutine mom(lflag,np,pt0,pt1,eta0,eta1,xmr,en0,pz0,
     +               pm,wgt,ranram,lw)
      implicit none
      real*8 pi,wgt,rpr,phr,ran0,wj,en0,pz0,prx,pry
      real*8 etai1,etai2,etaj1,etaj2,xmri2,bxmri2,wjeta
      real*8 cn,wjaci,wjacj,phi,ranj,ref,alimp,alimm
      real*8 p0m,p0p 
      real*8 phj,alim,phip,phim,etap,etam
      integer npar,n,nri,nr,iter,i,j,k,lw,icont,lflag
      parameter (pi= 3.14159265358979323846264338327950d0)
      parameter (npar= 100)
      real*8 pt0(npar),pt1(npar),eta0(npar),eta1(npar),
     +       xmr(npar),bxmr(npar)
      real*8 pm(0:4,npar),pt(npar),eta(npar)
      integer jp(npar),init,np
      real*8 p0r(npar),p1r(npar),qcut(npar),pcut(npar),bpt0(npar)
      real*8 bet0(npar),ph(npar)
      real*8 ranram(111)
      real*8 en,pz,ga,de,p,q,a,b,cut,ausp,ausm,v,vmr,vpr,v2
      real*8 det,rx,x1m,x2m,pmod,qmod,den
      real*8 argsinh,aus,betp,alpp
      real*8 etalim,csi,eta0i,eps
      data eps/1d-8/
      data init/0/
      save
*
      if (init.eq.0) then
        init= 1
        do iter= 1,np
          jp(iter)= iter
        enddo
        n  = np-1
        nri= 3
        etalim= 100.d0
*
        do iter= 1,n
          bpt0(jp(iter))= 0.d0
          pcut(jp(iter))= 0.d0
          qcut(jp(iter))= 0.d0
          bxmr(jp(iter))= 0.d0
          bet0(jp(iter))= etalim
          do k= iter+1,np
              bxmr(jp(iter))= bxmr(jp(iter))+xmr(jp(k))
              bpt0(jp(iter))= bpt0(jp(iter))+pt0(jp(k))
              pcut(jp(iter))= max(pcut(jp(iter)),pt0(jp(k)))
              bet0(jp(iter))= min(bet0(jp(iter)),eta0(jp(k)))
          enddo
        enddo
        qcut(jp(n))= pt0(jp(np))
*
      endif
*
      lw  = 0
      nr  = nri
      rpr = 0.d0
      wgt = 1.d0
      en  = en0
      pz  = pz0
      if (lflag.eq.0) then
        phr = 0.d0
      elseif (lflag.eq.1) then
        prx = 0.d0
        pry = 0.d0
      else
        goto 101
      endif
*
      do iter= 1,n
        i= jp(iter)
        j= jp(iter+1)
*
c-      Generation of transverse momenta when lflag= 0
*
        ga = en+pz
        de = en-pz 
*
        if (ga.lt.0.d0) goto 100
        if (de.lt.0.d0) goto 100
        v  = sqrt(ga*de)
*
c-      Protection:
*
        vmr= v-rpr
        if (vmr.lt.0.d0) then
          vmr= dabs(vmr)
          v  = vmr+2.d0*rpr
          de = ga/v/v
        endif
        vpr= vmr+2.d0*rpr
        v2 = v*v
*
        eta0i= max(dabs(eta0(i)),dabs(eta1(i)))
        if (ga.ge.de) then
          alpp= 2.d0*ga*exp(-eta0i)
        else
          alpp= 2.d0*de*exp(-eta0i)
        endif
        betp= ((ga+de)*cosh(eta0i)-dabs(ga-de)*sinh(eta0i))
        csi = min(alpp,betp)
*
        if (lflag.eq.1) then
          pt(i)= sqrt(pm(1,i)**2+pm(2,i)**2)
          if (iter.ne.1) then
            phi= ((pm(1,i)*prx+pm(2,i)*pry)/(eps+pt(i))/rpr)
            if (phi.gt. 1.d0) phi= 1.d0
            if (phi.lt.-1.d0) phi=-1.d0
            phi= acos(phi)
          endif
          prx  = prx-pm(1,i)
          pry  = pry-pm(2,i)
          pt(j)= sqrt(prx*prx+pry*pry)
*
          eta(i)= -log(tan(0.5d0*atan2(pt(i),pm(3,i))))
        endif
*
        xmri2 = xmr(i)**2
        bxmri2= bxmr(i)**2
*
        aus   =  vpr*vmr 
        ausp  = (aus+xmri2-bxmri2)
        det   = (ausp**2-4.d0*xmri2*aus)
        if (det.lt.0.d0) then
          det= abs(det)
          if (det.gt.1.d-12) goto 100
        endif
        det   = sqrt(det)
        if (ausp.gt.0.d0) then
           p0p= (rpr*ausp+v*det)/2.d0/aus
           p0m=-(ausp**2-4.d0*v2*xmri2)/4.d0/aus/p0p
        else
           p0m= (rpr*ausp-v*det)/2.d0/aus
           p0p=-(ausp**2-4.d0*v2*xmri2)/4.d0/aus/p0m
        endif
        p0r(i)= max(pt0(i),p0m,qcut(i)-rpr)
        p1r(i)= min(p0p,v-bpt0(i))
        if (csi.gt.2.d0*v) then
           p1r(i)= min(p1r(i),(vpr*vmr)/(csi-2.d0*rpr))
        endif
        cn= 0.9d0
        if (p0r(i).eq.0.d0) cn= 0.9d0
        if (p0r(i).gt.p1r(i)) goto 100
        call ppeaka(lflag,max(pt0(i),xmr(i),1.d-2),
     +             cn,p0r(i),p1r(i),pt(i),wjaci,ranram(nr),lw)
        nr= nr+1
        if (lw.eq.1) goto 100
*
        if (lflag.eq.0) call rans(ran0)
*
        if (iter.eq.1) then
          if (lflag.eq.0) then
            phi= 2.d0*pi*ran0
            pt(j)= pt(i)
            phj= phi+pi
          endif
          wgt= wgt*2.d0*pi*pt(i)*wjaci
        else
          if (lflag.eq.0) then
            if (ran0.le.0.5d0) then
              ref = 1.d0
            else
              ref =-1.d0
            endif
            call rans(ranj)
          endif
*
          p0r(j)= max(dabs(pt(i)-rpr),qcut(i),pt(i)-v+2.d0*pcut(i))
          det= (v-sqrt(pt(i)**2+xmri2))**2-bxmri2
          if (det.lt.0.d0) then
             det= abs(det)
             if (det.gt.1.d-12) goto 100
          endif
          det   = sqrt(det)
          p1r(j)= min(pt(i)+rpr,det)
          if (csi.gt.2.d0*v) then
            p1r(j)= min(p1r(j),sqrt(pt(i)*(pt(i)-csi)+v2))
          endif
          if (p0r(j).gt.p1r(j)) goto 100
*
          alimp= (rpr**2+pt(i)**2-p1r(j)**2)/2.d0/rpr/(eps+pt(i))
          alimm= (rpr**2+pt(i)**2-p0r(j)**2)/2.d0/rpr/(eps+pt(i))
*
c-        Protections:
*
          if (alimp.gt. 1.d0) alimp= 1.d0
          if (alimp.lt.-1.d0) alimp=-1.d0 
          if (alimm.gt. 1.d0) alimm= 1.d0
          if (alimm.lt.-1.d0) alimm=-1.d0
*          
          phip= acos(alimp)
          phim= acos(alimm)
          if (phim.gt.phip) goto 100
*
          call fflat(lflag,phi,phim,phip,wjacj,ranj,lw)
          if (lw.eq.1) goto 100
          wgt= wgt*wjaci*wjacj*pt(i)*2.d0
*
          if (lflag.eq.0) then
            phi= phi*ref
            pt(j)= dabs(rpr**2+pt(i)**2-2.d0*rpr*pt(i)*cos(phi))
            pt(j)= sqrt(pt(j))
*
c-          Protections:
*
            if (rpr.eq.0.d0)   rpr  = 1.d-25
            if (pt(j).eq.0.d0) pt(j)= 1.d-25
            alim= (rpr**2+pt(j)**2-pt(i)**2)/2.d0/rpr/(eps+pt(j))
*
c-          Protections:
*
            if (alim.gt. 1.d0) alim= 1.d0
            if (alim.lt.-1.d0) alim=-1.d0
*
            phj=-acos(alim)*ref
*
          endif
        endif  
*
        if (lflag.eq.0) then
          phi= phi+phr
          phj= phj+phr
          pm(1,i)= pt(i)*sin(phi)
          pm(2,i)= pt(i)*cos(phi)
          pm(1,j)= pt(j)*sin(phj)
          pm(2,j)= pt(j)*cos(phj)
          ph(i)  = phi
          if (iter.eq.n) ph(j)  = phj
*
          phr= phj
        else
          ph(i)= acos(pm(2,i)/(eps+pt(i))) 
          if (pm(1,i).lt.0.d0) ph(i)= 2.d0*pi-ph(i)
          if (iter.eq.n) then
            ph(j)= acos(pm(2,j)/(eps+pt(j)) )
            if (pm(1,j).lt.0.d0) ph(j)= 2.d0*pi-ph(j)
          endif
        endif
*
        rpr= pt(j)
*
        p  = pt(i)
        q  = pt(j)
        a  = 1.d0+xmri2/p/p
        b  = max(1.d0+bxmri2/q/q,1.d0+4.d0/q/q*pcut(i)*
     +          (pcut(i)-q),(bpt0(i)/q)**2)
*
        cut= (sqrt(a)*p+sqrt(b)*q)**2-v2
        if (cut.gt.1.d-12) goto 100
        ausm= a*p*p-b*q*q
        det = (v2+ausm)**2-4.d0*v2*a*p*p
        if (det.lt.0.d0) then
           det= abs(det)
           if (det.gt.1.d-12) goto 100
        endif
        det = sqrt(det)
        rx  = v2+ausm
        x1m= (rx+det)/2.d0/p/de  
        x2m= ga*a/de/x1m                  ! x2m < x1m
        etai1= argsinh(0.5d0*(x1m-a/x1m))
        etai2= argsinh(0.5d0*(x2m-a/x2m))
        if (iter.le.n-1) then
*
c-        Generation of etas when lflag= 0
*
          etam= max(etai2, eta1(i))        
          etap= min(etai1, eta0(i))        
          if (etam.ge.etap) goto 100
          call ttriangle(lflag,etai2,etai1,etam,etap,eta(i),
     +                  wjeta,ranram(nr),lw)
          nr= nr+1
          if (lw.eq.1) goto 100
          wgt   = wgt*wjeta 
          if (lflag.eq.0) then
            pm(0,i)= sqrt(pt(i)**2*cosh(eta(i))**2+xmri2)   
            pm(3,i)= pt(i)*sinh(eta(i))                         
            pm(0,j)= en-pm(0,i)
            pm(3,j)= pz-pm(3,i)
          endif
          en     = en-pm(0,i)
          pz     = pz-pm(3,i)
          wgt    = wgt*sqrt(1.d0-(xmr(i)/pm(0,i))**2)
        else
*
c-        Cuts on the last 2 etas
*
          etaj1= argsinh((pz-pt(i)*sinh(etai1))/q)
          etaj2= argsinh((pz-pt(i)*sinh(etai2))/q)
*
c-        Look whether 0,1 or 2 solutions may contribute
*
          icont= 0                        
          if ((etai1.lt.eta0(i)).and.
     +        (etai1.gt.eta1(i)).and.
     +        (etaj1.lt.eta0(j)).and.
     +        (etaj1.gt.eta1(j))) icont= icont+1
          if ((etai2.lt.eta0(i)).and.
     +        (etai2.gt.eta1(i)).and.
     +        (etaj2.lt.eta0(j)).and.
     +        (etaj2.gt.eta1(j))) icont= icont+2
*
          if (icont.eq.0) then
             goto 100
          elseif (icont.eq.1) then
             if (lflag.eq.0) then
               eta(i)= etai1
               eta(j)= etaj1
             endif
             wj= 1.d0
          elseif (icont.eq.2) then
             if (lflag.eq.0) then
               eta(i)= etai2
               eta(j)= etaj2
             endif
             wj= 1.d0
          else
             if (lflag.eq.0) then
               call rans(ran0)
               if (ran0.lt.0.5d0) then
                 eta(i)= etai1
                 eta(j)= etaj1
               else
                 eta(i)= etai2
                 eta(j)= etaj2
               endif
             endif
             wj= 2.d0
          endif
          if (lflag.eq.0) then
             pm(0,i)= sqrt(pt(i)**2*cosh(eta(i))**2+xmri2)   
             pm(3,i)= pt(i)*sinh(eta(i))                         
             pm(0,j)= en-pm(0,i)
             pm(3,j)= pz-pm(3,i)                                  
          elseif (lflag.eq.1) then
             en     = en-pm(0,i)-pm(0,j)
             pz     = pz-pm(3,i)-pm(3,j) 
             eta(j) =-log(tan(0.5d0*atan2(pt(j),pm(3,j))))
          else
             goto 101
          endif
          pmod= pm(0,i)**2-xmri2
          qmod= pm(0,j)**2-xmr(j)*xmr(j)
          if (pmod.lt.0.d0) goto 100    
          if (qmod.lt.0.d0) goto 100
          pmod= sqrt(pmod)
          qmod= sqrt(qmod)
          den = dabs(pm(0,j)*pm(3,i)-pm(0,i)*pm(3,j))
          wgt = wgt*wj*
     +          pmod*qmod/p/q/cosh(eta(j))/cosh(eta(i))/den
        endif 
      enddo     
*
      wgt= wgt/2.d0**np      
      return
 100  wgt= 0.d0
      lw= 1
      return
 101  print*,'ERROR IN SUBROUTINE MOM'
      stop
      end
*
      subroutine momr(lflag,irapgap,njets,etagap,
     +                np,pt0,pt1,eta0,xmr,en0,pz0,pm,wgt,lw)
      implicit none
      integer lflag,np,npar,lw,j
      integer irapgap,njets   
      parameter (npar= 20)
      real*8 xmr(npar),en0,pz0,pm(0:4,npar),wgt,dj
      real*8 pt0(npar),pt1(npar),eta0(npar)
      real*8 pt(npar),pl(npar),eta(npar)
      real*8 etagap
      real*8 x1,x2,bvel,gvel,rootshr,pr(4,npar)
      integer kp(10,6),jj,lk
      data kp/1,2,3,4,5,6,7,8,9,10,
     +        2,3,1,4,5,6,7,8,9,10,
     +        3,1,2,4,5,6,7,8,9,10,
     +        1,3,2,4,5,6,7,8,9,10,
     +        3,2,1,4,5,6,7,8,9,10,
     +        2,1,3,4,5,6,7,8,9,10/
*
      lw= 0
      if (lflag.eq.1) then
        en0= 0.d0
        pz0= 0.d0
        do j= 1,np
          en0= en0+pm(0,j)          
          pz0= pz0+pm(3,j)
        enddo
      endif
      x1  = dabs(en0+pz0)
      x2  = dabs(en0-pz0)
      bvel= (x1-x2)/(x1+x2)
      gvel= 1.d0/sqrt(1.d0-bvel*bvel)
      rootshr= sqrt(x1*x2)
*
      if (lflag.eq.1) then
        do j= 1,np
          pr(4,j)= gvel*(pm(0,j)-bvel*pm(3,j))
          pr(1,j)= pm(1,j)
          pr(2,j)= pm(2,j)
          pr(3,j)= gvel*(pm(3,j)-bvel*pm(0,j))
        enddo 
      endif
      call rambo(lflag,np,rootshr,xmr,pr,dj)
      wgt= 1.d0/dj
      if (lflag.eq.1) goto 10
*
      do j= 1,np
        pm(0,j)= gvel*(pr(4,j)+bvel*pr(3,j))
        pm(1,j)= pr(1,j)
        pm(2,j)= pr(2,j)
        pm(3,j)= gvel*(pr(3,j)+bvel*pr(4,j))
      enddo
*
 10   if (irapgap.eq.1) then
        do j= 1,njets
          pt(j) = sqrt(pm(1,j)**2+pm(2,j)**2)
          pl(j) = pm(3,j)
          eta(j)= -log(tan(0.5d0*atan2(pt(j),pl(j))))
        enddo
        if    (njets.eq.2) then
           if (dabs(eta(1)).lt.etagap) goto 100
           if (dabs(eta(2)).lt.etagap) goto 100
           if (eta(1)*eta(2).gt.0.d0)  goto 100
        elseif(njets.eq.3) then
           lk= 0  
           do jj= 1,6
             if ((eta(kp(1,jj)).gt. etagap).and.
     +           (eta(kp(2,jj)).lt.-etagap).and.
     +      (dabs(eta(kp(3,jj))).lt.etagap).and.
     +           (pt(kp(3,jj)).gt.pt0(3))) lk= 1 
           enddo
           if (lk.eq.0) goto 100
        endif
      endif
      return
 100  lw= 1
      wgt= 0.d0
      end


C*********************************************************************
      subroutine wdec(nw,idecmode,decfl,decbr)
C*********************************************************************
C
C given the number of decaying W and the string IMODE labelling the
C chosen decay mode returns the flavour of final fermions in DECFL and
C the branching ratio of the chosen mode
C
      implicit none
C
      integer nw,idecmode(8)
      integer decfl(8)
      real*8 decbr
C
      integer nmx
      parameter (nmx=1000)
      integer tab(8,nmx),tab1(8,nmx),ncomb,wgt(nmx),wgt1(nmx)
      integer tabdec(3,5),sum(6),sumref(6),ncomb1,idecmode_tmp(8)
      integer decfl0(8),idx(8),j1,j2,j3,j4,tbtmp(8),fac(0:8),j5,nlp
      real*8 wgtot,decbr0,cum(nmx),rr
      logical init,found
      save init,tab1,wgt1,ncomb1,tabdec,decbr0,cum,fac
      data init/.true./
      data fac /1,1,2,6,24,120,720,5040,40320/
C imode           e  m  t  l  q  1
      data tabdec/1, 4, 6,             ! e
     $            2, 4, 6,             ! mu
     $            3, 4, 6,             ! tau 
     $            5, 6, 0,             ! up
     $            5, 6, 0 /            ! charm
C
      if (init) then
        init=.false.
        call finsw(nw,tab,ncomb,wgt)
        do j1=1,6
          sumref(j1)=0
        enddo
        do j1=1,nw
          sumref(idecmode(j1))=sumref(idecmode(j1))+1
        enddo
        nlp=sumref(1)+sumref(2)+sumref(3)+sumref(4)
        ncomb1=0
        wgtot=0.d0
        decbr0=0.d0
        do j1=1,ncomb
          wgtot=wgtot+wgt(j1)
          do j2=1,nw
            tbtmp(j2)=tab(j2,j1)
          enddo
          do j3=1,5
            sum(j3)=0
          enddo
          do j3=1,nw
            sum(tbtmp(j3))=sum(tbtmp(j3))+1
          enddo
          found=.true.
          do j3=1,3
            if(sum(j3).lt.sumref(j3)) found=.false.
          enddo
          if(sum(1)+sum(2)+sum(3).lt.nlp) found=.false.
          if(sum(4)+sum(5).lt.sumref(5)) found=.false. 
          if(found) then
            ncomb1=ncomb1+1
            do j4=1,nw
              tab1(j4,ncomb1)=tbtmp(j4)
            enddo              
            wgt1(ncomb1)=wgt(j1)
            decbr0=decbr0+wgt1(ncomb1)
            if(ncomb1.eq.1) then
              cum(ncomb1)=wgt1(ncomb1)
            else
              cum(ncomb1)=wgt1(ncomb1)+cum(ncomb1-1)
            endif
          endif
        enddo
        do j1=1,ncomb1
          cum(j1)=cum(j1)/decbr0
        enddo
        decbr0=decbr0/wgtot
      endif
C
      decbr=decbr0
      call randa(rr)
      do j1=1,ncomb1
        if(rr.lt.cum(j1)) goto 12
      enddo
 12   continue
C
      do j2=1,nw
        decfl0(j2)=tab1(j2,j1)
      enddo
      call randa(rr)
      j1=rr*fac(nw)+1
      j1=min(fac(nw),j1)
      call genprm(nw,j1,idx)      
      do j2=1,nw
        decfl(j2)=decfl0(idx(j2))
      enddo
C
      return
      end
C*********************************************************************
      subroutine finsw(nw,tab,ncomb,wgt)
C*********************************************************************
C
C Given the number of w bosons NW <= 8 returns in TAB(I,J) the J-th possible
C combination of the NW decay products (I=1,NW). 
C TAB(I,J) = 1(e),2(mu),3(tau),4(u),5(c)
C and TAB(I1,J)<=TAB(I2,J) if I1 < I2, WGT(J) returns the partial weight of 
C the jth configuration
C
      implicit none
C
      integer nw,nmx
      parameter (nmx=1000)
      integer tab(8,nmx),ncomb,wgt(nmx)
C
      integer j1,j2,nt,j3,tbtmp(8),j4,fac(0:8),j5,j6
      logical found
      save fac
      data fac /1,1,2,6,24,120,720,5040,40320/
C
      if(nw.gt.8.or.nw.lt.1) then 
        write(*,*)'nw out of range in finsw'
        stop
      endif
C
      nt=5**nw
      j1=0
      do j2=0,nt-1
        j3=j2
        do j4=1,nw
          tbtmp(j4)=mod(j3,5)+1
          j3=j3/5
        enddo
        found=.true.
        do j3=1,nw-1
          if(tbtmp(j3).gt.tbtmp(j3+1)) found =.false.
        enddo
        if(found) then
          j1=j1+1
          if(j1.gt.nmx) then
            write(*,*)'nmx too small in finsw'
            stop
          endif
          do j3=1,nw
            tab(j3,j1)=tbtmp(j3)
          enddo
          wgt(j1)=fac(nw)
          do j3=1,5
            j4=0
            do j5=1,nw
              if(tbtmp(j5).eq.j3) j4=j4+1
            enddo
            wgt(j1)=wgt(j1)/fac(j4)
            if(tbtmp(j3).gt.3)j6=j6+1
          enddo
          j6=0
          do j3=1,nw
            if(tbtmp(j3).gt.3)j6=j6+1
          enddo
          wgt(j1)=wgt(j1)*3**j6
        endif
      enddo
C
      ncomb=j1
C
      return
      end
