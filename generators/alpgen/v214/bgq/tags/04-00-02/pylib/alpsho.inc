c process code
      integer ihrd
      integer itopprc
c--event file data
      integer Nunit,NunitOut,NunitIni
      character*100 filename,tmpstr
c total number of input events
      common/AHio/Nunit,NunitOut,NunitIni,filename
c inputs for the matching algorithm
      integer iexc,npfst,nplst,nljets,njstart,njlast
     $     ,ickkw
      double precision etclus,rclus,etaclmax
      common/AHopts/etclus,rclus,etaclmax,iexc,npfst
     $     ,nplst,nljets,njstart,njlast,ickkw
c process and particles parameters
      double precision mc,mb,mt,mw,mz,mh
      double precision ebeam
      integer ndns,ih1,ih2       
      integer nw,nz,nh,nph
      integer ihvy,ihvy2
c pdf set type
      character pdftyp*25
c total number of partons
      integer npart
      common/AHppara/mc,mb,mt,mw,mz,mh,
     & ebeam,ih1,ih2,
     & ihrd,itopprc,
     & nw,nz,nh,nph,
     & ihvy,ihvy2,
     & npart,ndns,pdftyp
c weight information
      real *8 maxwgt,avgwgt,errwgt,totlum
      integer unwev
      common/AHwgts/maxwgt,avgwgt,errwgt,totlum,unwev
c general parameters
      integer nparam
      parameter (nparam=200)
      integer parlen,partyp
      character chpar*8,chpdes*70
      double precision parval
      common/AHpars/parval(nparam),chpar(nparam),chpdes(nparam)
     $     ,parlen(nparam),partyp(nparam)
c global event cuts
      double precision ptjmin,ptjmax,etajmax,drjmin,
     +     ptbmin,ptbmax,etabmax,drbmin,
     +     ptcmin,ptcmax,etacmax,drcmin,
     +     ptphmin,etaphmax,drphjmin,drphmin,drphlmin,
     +     ptlmin,etalmax,drlmin,metmin,mllmin,mllmax
      common/AHcuts/ptjmin,ptjmax,etajmax,drjmin,
     +     ptbmin,ptbmax,etabmax,drbmin,
     +        ptcmin,ptcmax,etacmax,drcmin,
     +        ptphmin,etaphmax,drphjmin,drphmin,drphlmin,
     +        ptlmin,etalmax,drlmin,metmin,mllmin,mllmax

