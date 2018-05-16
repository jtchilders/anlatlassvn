      subroutine runFH(mA,M2,M3,mstl,mstr,msbr,msff,Af,mu,tanb,
     $     mh_int,alpha)

      implicit none
#include "FHRecord.h"
#include "FHCouplings.h"
#include "SLHA.h"
      include 'PhysPars.h'
      double precision mA,M2,M3,mstl,mstr,msbr,msff,
     $     Af,mu,tanb,mh_int(2),alpha
c     FeynHiggs stuff
      integer error, nmfv,fast,key
      integer mssmpart, fieldren, tanbren, higgsmix, p2approx
      integer looplevel, runningMT, botResum, tlCplxApprox
      double precision MSf(2,4,3), MASf(6,4), MCha(2), MNeu(4)
      double complex USf(2,2,4,3),UASf(6,6,4)
      double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
      double complex DeltaMB
      double precision MGl
      double precision MHtree(4), SAtree
      double precision MHiggs(4)
      double complex SAeff, UHiggs(3,3), ZHiggs(3,3)
      double complex couplings(ncouplings), couplingsms(ncouplingsms)
      double precision gammas(ngammas), gammasms(ngammasms)
      character*50 filename,infilename,slhafilename
      double complex slhadata(nslhadata)
      double precision mst(2),msb(2)
      double complex Ust(2,2)
      double complex USb(2,2), Deltab
      integer FHslhaflag
      double precision powheginput
      external powheginput
      double precision invAlfa, AlfasMZ, GF
      double precision ME, upmass, MD, MM, MC, MS, ML, MB
      double precision MW, MZ
      double precision CKMlambda, CKMA, CKMrhobar, CKMetabar
      double precision Qtau,Qt,Qb,scalefactor,MHp
      double precision M1,M1imag,Afimag,Abimag,Atimag,muimag,M2imag,
     $M3imag
      RecordDecl(record)
      double complex M3cmplx,mucmplx,Abcmplx,Atcmplx
      double precision sqrts

      mssmpart = 4
      fieldren = 0
      tanbren = 0
      higgsmix = 2
      p2approx = 4
      looplevel = 2
      runningMT = 1
      botResum = 1
      tlCplxApprox = 0

      call FHSetFlags(error,
     $     mssmpart, fieldren, tanbren, higgsmix, p2approx,
     $     looplevel, runningMT, botResum, tlCplxApprox)

      call FHSetDebug(0)

      FHslhaflag = int(powheginput('FHslhaflag'))
      if(FHslhaflag.eq.0) then
        filename = 'powheg-fh.in'
      else
        filename = 'powheg-fh.slha'
      endif
      write(*,*) 'Reading from ',filename
      call FHReadRecord(error,record,slhadata,filename)
      if(error.eq.0) then
        write(*,*) 'Read succesfully input data from SLHA file'
      else if (error.eq.2) then
        write(*,*) 'Read succesfully input data from FH input file'
      else
        write(*,*) 'Error in reading ', filename, 'at line ',error
      endif
      call FHLoopRecord(error, record)
      call FHSetRecord(error,record)

!     Here we retrieve the parameter we need for the others subroutines

      call FHRetrieveSMPara(error,invAlfa, ph_asmz_nnlo, ph_GF,
     $     ME, upmass,MD,MM,MC,MS,ML,
     $     ph_mbmb, ph_Wmass, ph_Zmass,
     $     CKMlambda, CKMA, CKMrhobar, CKMetabar)

      call FHRetrievePara(error, scalefactor, ph_topmass, tanb, mA, MHp,
     $     msff, msff, mstl, mstr, msbr,
     $     msff, msff, msff, msff, msff,
     $     msff, msff, msff, msff, msff,
     $     mucmplx,
     $     dcmplx(Af,Afimag), Atcmplx, Abcmplx,
     $     dcmplx(Af,Afimag), dcmplx(Af,Afimag), dcmplx(Af,Afimag),
     $     dcmplx(Af,Afimag), dcmplx(Af,Afimag), dcmplx(Af,Afimag),
     $     dcmplx(M1,M1imag),dcmplx(M2,M2imag), m3cmplx,
     $     Qtau, Qt, Qb)

      At=DREAL(Atcmplx)
      Ab=DREAL(Abcmplx)
      mu=-DREAL(mucmplx)
      M3=DREAL(m3cmplx)

      call FHGetPara(error, nmfv, MSf, Usf, MASf, UASf,
     $     MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
     $     MHtree, SAtree)

      call FHGetTLPara(error, MSb, USb, Msbl2, Deltab)

      call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)

      mh_int(1) = MHiggs(1)
      mh_int(2) = MHiggs(2)
      alpha = asin(dreal(SAeff))

      mst(1) = MASf(3,3)
      mst(2) = MASf(6,3)
      Ust(1,1) = UASf(3,3,3)
      Ust(1,2) = UASf(3,6,3)
      Ust(2,1) = UASf(6,3,3)
      Ust(2,2) = UASf(6,6,3)

      write(*,*) 'Higgs and alpha in FH', mh_int(1), mh_int(2), alpha

      write(*,*) 'Stops in FeynHiggs:',MSf(1,3,3),MSf(2,3,3),
     $     dreal(USf(1,1,3,3)),dreal(USf(1,2,3,3))
      write(*,*) 'Sbots in FeynHiggs:',MSb(1),MSb(2),
     $     dreal(USb(1,1)),dreal(USb(1,2))

      ph_t1 = MSf(1,3,3)
      ph_t1_2 = ph_t1**2
      ph_t2 = MSf(2,3,3)
      ph_t2_2 = ph_t2**2

      ph_b1 = MSb(1)
      ph_b1_2 = ph_b1**2
      ph_b2 = MSb(2)
      ph_b2_2 = ph_b2**2

      ph_ctt = dreal(ust(1,1))
      ph_stt = dreal(ust(1,2))
      ph_ctb = dreal(usb(1,1))
      ph_stb = dreal(usb(1,2))

      ph_deltab = dreal(DeltaMB)
      write(*,*) 'deltab from FH', ph_deltab

c     We use the same indexing of FH for the Higgses.
      call FHCouplings(error, couplings, couplingsms, gammas,
     $ gammasms, fast)

      fh_Hwidth = GammaTot(ih)

      write(*,*) 'Writing FH data to powheg-fh.slha'
      slhafilename = 'powheg-fh-output.slha'
      call SLHAClear(slhadata)
      call FHOutputSLHA(error, slhadata, 255)
      call SLHAWrite(error, slhadata, slhafilename)


      return
      end
