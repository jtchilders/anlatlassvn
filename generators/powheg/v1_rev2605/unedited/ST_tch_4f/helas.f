

c
c ----------------------------------------------------------------------
c
      subroutine ficxxx(fi , fic)
c
c this subroutine charge conjugates a flowing-in fermion wavefunction.  
c                                                                       
c input:                                                                
c       complex fi(6)          : flowing-in fermion                 |fi>
c                                                                       
c output:                                                               
c       complex fic(6)         : charge conjugated fermion         <fic|
c
      implicit none
      double complex fi(6), fic(6)
c
      fic(1) = -fi(2)
      fic(2) =  fi(1)
      fic(3) =  fi(4)
      fic(4) = -fi(3)
      fic(5) = -fi(5)
      fic(6) = -fi(6)
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine focxxx(fo , foc)
c
c this subroutine charge conjugates a flowing-out fermion wavefunction. 
c                                                                       
c input:                                                                
c       complex fo(6)          : flowing-out fermion                <fo|
c                                                                       
c output:                                                               
c       complex foc(6)         : charge conjugated fermion         |foc>
c
      implicit none
      double complex fo(6), foc(6)
c
      foc(1) =  fo(2)
      foc(2) = -fo(1)
      foc(3) = -fo(4)
      foc(4) =  fo(3)
      foc(5) = -fo(5)
      foc(6) = -fo(6)
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine fsicxx(fic,sc,gc,fmass,fwidth , fsic)
c
c this subroutine computes an off-shell antifermion wavefunction from a 
c flowing-in external antifermion and a vector boson.                   
c                                                                       
c input:                                                                
c       complex fic(6)         : flow-in  antifermion           |fic>
c       complex sc(3)          : input    scalar                   s 
c       complex gc(2)          : coupling constants              gchf
c       real    fmass          : mass  of output antifermion fc'     
c       real    fwidth         : width of output antifermion fc'     
c                                                                       
c output:                                                               
c       complex fsic(6)        : off-shell fermion        |fc',s,fic>
c
      implicit none
      double complex fic(6),sc(3),fsic(6),gc(2),sl1,sl2,sr1,sr2,ds
      double precision pf(0:3),fmass,fwidth,pf2,p0p3,p0m3
c
      fsic(5) = fic(5)-sc(2)
      fsic(6) = fic(6)-sc(3)

      pf(0) = dble( fsic(5))
      pf(1) = dble( fsic(6))
      pf(2) = dimag(fsic(6))
      pf(3) = dimag(fsic(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      ds = -sc(1)/dcmplx( pf2-fmass**2, fmass*fwidth )
      p0p3 = pf(0)+pf(3)
      p0m3 = pf(0)-pf(3)
      sl1 = gc(1)*(p0p3*fic(1)+dconjg(fsic(6))*fic(2))
      sl2 = gc(1)*(p0m3*fic(2)       +fsic(6) *fic(1))
      sr1 = gc(2)*(p0m3*fic(3)-dconjg(fsic(6))*fic(4))
      sr2 = gc(2)*(p0p3*fic(4)       -fsic(6) *fic(3))

      fsic(1) = ( gc(1)*fmass*fic(1) + sr1 )*ds
      fsic(2) = ( gc(1)*fmass*fic(2) + sr2 )*ds
      fsic(3) = ( gc(2)*fmass*fic(3) + sl1 )*ds
      fsic(4) = ( gc(2)*fmass*fic(4) + sl2 )*ds
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine fsigld(fi,sc,gc,fmass,fwidth,smass,mNLSP,idecay , fsi)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-in external fermion and a scalar boson, for the NLSP-boson-
c Goldstino vertex. The h.c. of the NLSP decay is handled via the
c input parameter idecay.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex sc(3)          : input    scalar                      s
c       complex gc(2)          : coupling constants                  gsf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c       integer idecay         :  1 for NLSP decay to Goldstino
c                              : -1 for Goldstino to NLSP (h.c. of above)
c
c output:
c       complex fsi(6)         : off-shell fermion             |f',s,fi>
c
      implicit none
      double complex  fi(6), sc(3), gc(2), fsi(6), s1, s2, s3, s4, ds
      double complex  p14p, p14m, p23p, p23m
      double precision  pf(0:3), fmass, fwidth, mNLSP, smass, pf2
      integer idecay

      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fsi(5) = fi(5) - sc(2)
      fsi(6) = fi(6) - sc(3)

      pf(0) = dble( fsi(5))
      pf(1) = dble( fsi(6))
      pf(2) = dimag(fsi(6))
      pf(3) = dimag(fsi(5))
      pf2 = pf(0)**2 - pf(1)**2 - pf(2)**2 - pf(3)**2

      if ( idecay.ne.1 .or. idecay.ne.-1 ) then
         write(6,*) 'error in idecay of FSIGLD'
         stop
      end if

      p14p = dble(sc(2)) + dimag(sc(2))
      p14m = dble(sc(2)) - dimag(sc(2))
      p23p = dble(sc(3)) + dimag(sc(3))*ci
      p23m = dble(sc(3)) - dimag(sc(3))*ci

      ds = -sc(1)/dcmplx( pf2-fmass**2, fmass*fwidth )

      s1 = -idecay*gc(1)*fi(1)*smass**2
      s2 = -idecay*gc(1)*fi(2)*smass**2
      s3 = gc(1)*mNLSP*( fi(1)*p14p + fi(2)*p23m ) 
      s4 = gc(1)*mNLSP*( fi(1)*p23p + fi(2)*p14m )

      if ( gc(2).ne.cZero ) then
         s1 = s1 + gc(2)*mNLSP*( fi(3)*p14m - fi(4)*p23m )
         s2 = s2 + gc(2)*mNLSP*(-fi(3)*p23p + fi(4)*p14p )
         s3 = s3 - gc(2)*idecay*fi(3)*smass**2
         s4 = s4 - gc(2)*idecay*fi(4)*smass**2
      end if

      fsi(1) = ( (pf(0)-pf(3))*s3 - dconjg(fsi(6))*s4 + fmass*s1 )*ds
      fsi(2) = (       -fsi(6)*s3 +  (pf(0)+pf(3))*s4 + fmass*s2 )*ds
      fsi(3) = ( (pf(0)+pf(3))*s1 + dconjg(fsi(6))*s2 + fmass*s3 )*ds
      fsi(4) = (        fsi(6)*s1 +  (pf(0)-pf(3))*s2 + fmass*s4 )*ds
c
      return          
      end
c
c ----------------------------------------------------------------------
c
      subroutine fsocxx(foc,sc,gc,fmass,fwidth , fsoc)
c
c this subroutine computes an off-shell antifermion wavefunction from a 
c flowing-out external antifermion and a vector boson.                  
c                                                                       
c input:                                                                
c       complex foc(6)         : flow-out fermion               <foc|
c       complex sc(6)          : input    scalar                   s 
c       complex gc(2)          : coupling constants              gchf
c       real     fmass         : mass  of output antifermion fc'     
c       real     fwidth        : width of output antifermion fc'     
c                                                                       
c output:                                                               
c       complex fsoc(6)        : off-shell fermion         <fo,s,fc'|
c
      implicit none
      double complex foc(6),sc(6),fsoc(6),gc(2),sl1,sl2,sr1,sr2,ds
      double precision pf(0:3),fmass,fwidth,pf2,p0p3,p0m3
c
      fsoc(5) = foc(5)+sc(2)
      fsoc(6) = foc(6)+sc(3)

      pf(0) = dble( fsoc(5))
      pf(1) = dble( fsoc(6))
      pf(2) = dimag(fsoc(6))
      pf(3) = dimag(fsoc(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      ds = -sc(1)/dcmplx( pf2-fmass**2, fmass*fwidth )
      p0p3 = pf(0)+pf(3)
      p0m3 = pf(0)-pf(3)
      sl1 = gc(2)*(p0p3*foc(3)       +fsoc(6) *foc(4))
      sl2 = gc(2)*(p0m3*foc(4)+dconjg(fsoc(6))*foc(3))
      sr1 = gc(1)*(p0m3*foc(1)       -fsoc(6) *foc(2))
      sr2 = gc(1)*(p0p3*foc(2)-dconjg(fsoc(6))*foc(1))

      fsoc(1) = ( gc(1)*fmass*foc(1) + sl1 )*ds
      fsoc(2) = ( gc(1)*fmass*foc(2) + sl2 )*ds
      fsoc(3) = ( gc(2)*fmass*foc(3) + sr1 )*ds
      fsoc(4) = ( gc(2)*fmass*foc(4) + sr2 )*ds
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine fsogld(fo,sc,gc,fmass,fwidth,smass,mNLSP,idecay , fso)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-out external fermion and a scalar boson, for the NLSP-boson-
c Goldstino vertex. The h.c. of the NLSP decay is handled via the
c input parameter idecay.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex sc(3)          : input    scalar                      s
c       complex gc(2)          : coupling constants                  gsf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c       integer idecay         :  1 for NLSP decay to Goldstino
c                              : -1 for Goldstino to NLSP (h.c. of above)
c
c output:
c       complex fso(6)         : off-shell fermion             <fo,s,f'|
c
      implicit none
      double complex  fo(6), sc(3), gc(2), fso(6), s1, s2, s3, s4, ds
      double precision  pf(0:3), fmass, fwidth, mNLSP, smass, pf2
      double precision  p14p, p14m, p23p, p23m
      integer idecay

      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fso(5) = fo(5) + sc(2)
      fso(6) = fo(6) + sc(3)

      pf(0) = dble( fso(5))
      pf(1) = dble( fso(6))
      pf(2) = dimag(fso(6))
      pf(3) = dimag(fso(5))
      pf2 = pf(0)**2 - pf(1)**2 - pf(2)**2 - pf(3)**2

      if ( idecay.ne.1 .or. idecay.ne.-1 ) then
         write(6,*) 'error in idecay of FSOGLD'
         stop
      end if

      p14p = dble(sc(2)) + dimag(sc(2))
      p14m = dble(sc(2)) - dimag(sc(2))
      p23p = dble(sc(3)) + dimag(sc(3))*ci
      p23m = dble(sc(3)) - dimag(sc(3))*ci

      ds = -sc(1)/dcmplx( pf2-fmass**2, fmass*fwidth )

      s1 = -idecay*gc(1)*fo(1)*smass**2
      s2 = -idecay*gc(1)*fo(2)*smass**2
      s3 = gc(1)*mNLSP*( fo(1)*p14m - fo(2)*p23p )
      s4 = gc(1)*mNLSP*(-fo(1)*p23m + fo(2)*p14p )

      if ( gc(2).ne.cZero ) then
         s1 = s1 + gc(2)*mNLSP*( fo(3)*p14p + fo(4)*p23p )
         s2 = s2 + gc(2)*mNLSP*( fo(3)*p23m + fo(4)*p14m )
         s3 = s3 - gc(2)*idecay*fo(3)*smass**2
         s4 = s4 - gc(2)*idecay*fo(4)*smass**2
      end if

      fso(1) = (  (pf(0)+pf(3))*s3 +         fso(6)*s4 + fmass*s1 )*ds
      fso(2) = ( dconjg(fso(6))*s3 +  (pf(0)-pf(3))*s4 + fmass*s2 )*ds
      fso(3) = (  (pf(0)-pf(3))*s1 -         fso(6)*s2 + fmass*s3 )*ds
      fso(4) = (-dconjg(fso(6))*s1 +  (pf(0)+pf(3))*s2 + fmass*s4 )*ds
c
      return          
      end
c
c ----------------------------------------------------------------------
c
      subroutine fvicxx(fic,vc,gc,fmass,fwidth , fvic)
c
c this subroutine computes an off-shell antifermion wavefunction from a 
c flowing-in external antifermion and a vector boson.                   
c                                                                       
c input:                                                                
c       complex fic(6)         : flow-in  antifermion              |fic>
c       complex vc(6)          : input    vector                      v 
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output antifermion f'         
c       real    fwidth         : width of output antifermion f'         
c                                                                       
c output:                                                               
c       complex fvic(6)        : off-shell antifermion       |fc',v,fic>
c
      implicit none
      double complex fic(6),vc(6),gc(2),fvic(6),sl1,sl2,sr1,sr2,d
      double precision pf(0:3),fmass,fwidth,pf2

      double precision rOne
      parameter( rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fvic(5) = fic(5)-vc(5)
      fvic(6) = fic(6)-vc(6)

      pf(0) = dble( fvic(5))
      pf(1) = dble( fvic(6))
      pf(2) = dimag(fvic(6))
      pf(3) = dimag(fvic(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d = rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
      sr1 =   (vc(1)-      vc(4))*fic(3)
     &      - (vc(2)-cImag*vc(3))*fic(4)
      sr2 = - (vc(2)+cImag*vc(3))*fic(3)
     &      + (vc(1)+      vc(4))*fic(4)

      if ( gc(2).ne.cZero ) then
         sl1 =   (vc(1)+      vc(4))*fic(1)
     &         + (vc(2)-cImag*vc(3))*fic(2)
         sl2 =   (vc(2)+cImag*vc(3))*fic(1)
     &         + (vc(1)-      vc(4))*fic(2)

         fvic(1) = ( gc(2)*((pf(0)-pf(3))*sl1 - dconjg(fvic(6))*sl2)
     &              +gc(1)*fmass*sr1 )*d
         fvic(2) = ( gc(2)*(     -fvic(6)*sl1 +   (pf(0)+pf(3))*sl2)
     &              +gc(1)*fmass*sr2 )*d
         fvic(3) = ( gc(1)*((pf(0)+pf(3))*sr1 + dconjg(fvic(6))*sr2)
     &              +gc(2)*fmass*sl1 )*d
         fvic(4) = ( gc(1)*(      fvic(6)*sr1 +   (pf(0)-pf(3))*sr2)
     &              +gc(2)*fmass*sl2 )*d
      else
         d = d * gc(1)
         fvic(1) = fmass*sr1*d
         fvic(2) = fmass*sr2*d
         fvic(3) = ((pf(0)+pf(3))*sr1 + dconjg(fvic(6))*sr2)*d
         fvic(4) = (      fvic(6)*sr1 +   (pf(0)-pf(3))*sr2)*d
      end if
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine fvigld(fi,vc,gc,fmass,fwidth,idecay , fvi)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-in external fermion and a vector boson, for the NLSP-boson-
c Goldstino vertex. The h.c. of the NLSP decay is handled via the
c input parameter idecay (picks out correct Goldstino momentum).
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c       integer idecay         :  1 for NLSP decay to Goldstino
c                              : -1 for Goldstino to NLSP (h.c. of above)
c
c output:
c       complex fvi(6)         : off-shell fermion             |f',v,fi>
c
      implicit none
      double complex  fi(6), vc(6), gc(2), fvi(6), sl1, sl2, sr1, sr2, d
      double complex  p14p, p14m, p23p, p23m, A14p, A14m, A23p, A23m
      double complex  AdotpG
      double precision  fmass, fwidth
      double precision  pf(0:3), pv(4), pf2, pdotpG
      integer idecay

      double precision rOne
      parameter( rOne = 1.0d0 )
      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fvi(5) = fi(5) - vc(5)
      fvi(6) = fi(6) - vc(6)

      pv(1) = dble( vc(5))
      pv(2) = dble( vc(6))
      pv(3) = dimag(vc(6))
      pv(4) = dimag(vc(5))

      pf(0) = dble( fvi(5))
      pf(1) = dble( fvi(6))
      pf(2) = dimag(fvi(6))
      pf(3) = dimag(fvi(5))
      pf2 = pf(0)**2 - pf(1)**2 - pf(2)**2 - pf(3)**2

      if ( idecay.eq.1 ) then
         pdotpG = pv(1)*pf(0) - pv(2)*pf(1) - pv(3)*pf(2) - pv(4)*pf(3)
         AdotpG = vc(1)*pf(0) - vc(2)*pf(1) - vc(3)*pf(2) - vc(4)*pf(3)
      else if ( idecay.eq.-1 ) then
         pdotpG =  pv(1)*dble( fi(5)) - pv(2)*dble( fi(6))
     &           - pv(3)*dimag(fi(6)) - pv(4)*dimag(fi(5))
         AdotpG =  vc(1)*dble( fi(5)) - vc(2)*dble( fi(6))
     &           - vc(3)*dimag(fi(6)) - vc(4)*dimag(fi(5))
      else
         write(6,*) 'error in idecay of FVIGLD'
         stop
      end if

      p14p = dble(vc(5)) + dimag(vc(5))
      p14m = dble(vc(5)) - dimag(vc(5))
      p23p = vc(6)
      p23m = dconjg(vc(6))
      A14p = vc(1) + vc(4)
      A14m = vc(1) - vc(4)
      A23p = vc(2) + vc(3)*ci
      A23m = vc(2) - vc(3)*ci

      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
      d = d*idecay

      sl1 =  (p14p*AdotpG - A14p*pdotpG)*fi(1)
     &      +(p23m*AdotpG - A23m*pdotpG)*fi(2)
      sl2 =  (p23p*AdotpG - A23p*pdotpG)*fi(1)
     &      +(p14m*AdotpG - A14m*pdotpG)*fi(2)

      if ( gc(2).ne.cZero ) then
         sr1 =  (p14m*AdotpG - A14m*pdotpG)*fi(3)
     &         -(p23m*AdotpG - A23m*pdotpG)*fi(4)
         sr2 = -(p23p*AdotpG - A23p*pdotpG)*fi(3)
     &         +(p14p*AdotpG - A14p*pdotpG)*fi(4)

         fvi(1) = ( gc(1)*((pf(0)-pf(3))*sl1 - dconjg(fvi(6))*sl2 )
     &             +gc(2)*fmass*sr1 )*d
         fvi(2) = ( gc(1)*(      -fvi(6)*sl1 +  (pf(0)+pf(3))*sl2 )
     &             +gc(2)*fmass*sr2 )*d
         fvi(3) = ( gc(2)*((pf(0)+pf(3))*sr1 + dconjg(fvi(6))*sr2 )
     &             +gc(1)*fmass*sl1 )*d
         fvi(4) = ( gc(2)*(       fvi(6)*sr1 +  (pf(0)-pf(3))*sr2 )
     &             +gc(1)*fmass*sl2 )*d

      else
         d = d*gc(1)
         fvi(1) = d*((pf(0)-pf(3))*sl1 - dconjg(fvi(6))*sl2)
         fvi(2) = d*(      -fvi(6)*sl1 +  (pf(0)+pf(3))*sl2)
         fvi(3) = d*fmass*sl1
         fvi(4) = d*fmass*sl2
      end if
c
      return          
      end
      subroutine fvigox(fi,vc,gc,fmass,fwidth , fvi)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN external fermion and a vector boson in the case when 
c they are all color octets.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvi(6)         : off-shell fermion             |f',v,fi>
c
      implicit none
      double complex fi(6),vc(6),gc(2),fvi(6),sl1,sl2,sr1,sr2,d
      double precision pf(0:3),fmass,fwidth,pf2
      integer i
      
      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      fvi(5) = fi(5)-vc(5)
      fvi(6) = fi(6)-vc(6)

      pf(0) = dble( fvi(5))
      pf(1) = dble( fvi(6))
      pf(2) = dimag(fvi(6))
      pf(3) = dimag(fvi(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
      sl1 =   (vc(1)+      vc(4))*fi(1)
     &      + (vc(2)-cImag*vc(3))*fi(2)
      sl2 =   (vc(2)+cImag*vc(3))*fi(1)
     &      + (vc(1)-      vc(4))*fi(2)

      if ( gc(2).ne.cZero ) then
         sr1 =   (vc(1)-      vc(4))*fi(3)
     &         - (vc(2)-cImag*vc(3))*fi(4)
         sr2 = - (vc(2)+cImag*vc(3))*fi(3)
     &         + (vc(1)+      vc(4))*fi(4)

         fvi(1) = ( gc(1)*((pf(0)-pf(3))*sl1 - dconjg(fvi(6))*sl2)
     &             +gc(2)*fmass*sr1 )*d
         fvi(2) = ( gc(1)*(      -fvi(6)*sl1 +  (pf(0)+pf(3))*sl2)
     &             +gc(2)*fmass*sr2 )*d
         fvi(3) = ( gc(2)*((pf(0)+pf(3))*sr1 + dconjg(fvi(6))*sr2)
     &             +gc(1)*fmass*sl1 )*d
         fvi(4) = ( gc(2)*(       fvi(6)*sr1 +  (pf(0)-pf(3))*sr2)
     &             +gc(1)*fmass*sl2 )*d

      else
         d = d * gc(1)
         fvi(1) = ((pf(0)-pf(3))*sl1 - dconjg(fvi(6))*sl2)*d
         fvi(2) = (      -fvi(6)*sl1 +  (pf(0)+pf(3))*sl2)*d
         fvi(3) = fmass*sl1*d
         fvi(4) = fmass*sl2*d
      end if
c
c    Need sign flip for wavefunction to account for 
c    difference between Helas convention and color
      do i=1,4
        fvi(i)=-fvi(i)
      enddo
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine fvocxx(foc,vc,gc,fmass,fwidth , fvoc)
c
c this subroutine computes an off-shell antifermion wavefunction from a 
c flowing-out external antifermion and a vector boson.                  
c                                                                       
c input:                                                                
c       complex foc(6)         : flow-out antifermion              <foc|
c       complex vc(6)          : input    vector                      v 
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output antifermion f'         
c       real    fwidth         : width of output antifermion f'         
c                                                                       
c output:                                                               
c       complex fvoc(6)        : off-shell antifermion       <foc,v,fc'|
c
      implicit none
      double complex foc(6),vc(6),gc(2),fvoc(6),sl1,sl2,sr1,sr2,d
      double precision pf(0:3),fmass,fwidth,pf2

      double precision rOne
      parameter( rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fvoc(5) = foc(5)+vc(5)
      fvoc(6) = foc(6)+vc(6)

      pf(0) = dble( fvoc(5))
      pf(1) = dble( fvoc(6))
      pf(2) = dimag(fvoc(6))
      pf(3) = dimag(fvoc(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d = rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
      sr1 =   (vc(1)-      vc(4))*foc(1)
     &      - (vc(2)+cImag*vc(3))*foc(2)
      sr2 = - (vc(2)-cImag*vc(3))*foc(1)
     &      + (vc(1)+      vc(4))*foc(2)

      if ( gc(2).ne.cZero ) then
         sl1 =   (vc(1)+      vc(4))*foc(3)
     &         + (vc(2)+cImag*vc(3))*foc(4)
         sl2 =   (vc(2)-cImag*vc(3))*foc(3)
     &         + (vc(1)-      vc(4))*foc(4)

         fvoc(1) = ( gc(1)*( (pf(0)+pf(3))*sr1   +       fvoc(6)*sr2)
     &              +gc(2)*fmass*sl1 )*d
         fvoc(2) = ( gc(1)*( dconjg(fvoc(6))*sr1 + (pf(0)-pf(3))*sr2)
     &              +gc(2)*fmass*sl2 )*d
         fvoc(3) = ( gc(2)*( (pf(0)-pf(3))*sl1   -       fvoc(6)*sl2)
     &              +gc(1)*fmass*sr1 )*d
         fvoc(4) = ( gc(2)*(-dconjg(fvoc(6))*sl1 + (pf(0)+pf(3))*sl2)
     &              +gc(1)*fmass*sr2 )*d
      else
         d = d * gc(1)
         fvoc(1) = ((pf(0)+pf(3))*sr1          +fvoc(6)*sr2)*d
         fvoc(2) = (dconjg(fvoc(6))*sr1 + (pf(0)-pf(3))*sr2)*d
         fvoc(3) = fmass*sr1*d
         fvoc(4) = fmass*sr2*d
      end if
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine fvogld(fo,vc,gc,fmass,fwidth,idecay , fvo)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-out external fermion and a vector boson, for the NLSP-boson-
c Goldstino vertex. The h.c. of the NLSP decay is handled via the
c input parameter idecay (picks out correct Goldstino momentum).
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c       integer idecay         :  1 for NLSP decay to Goldstino
c                              : -1 for Goldstino to NLSP (h.c. of above)
c
c output:
c       complex fvo(6)         : off-shell fermion             <fo,v,f'|
c
      implicit none
      double complex  fo(6), vc(6), gc(2), fvo(6), sl1, sl2, sr1, sr2, d
      double complex  p14p, p14m, p23p, p23m, A14p, A14m, A23p, A23m
      double complex  AdotpG
      double precision  fmass, fwidth
      double precision  pf(0:3), pv(4), pf2, pdotpG
      integer idecay

      double precision rOne
      parameter( rOne = 1.0d0 )
      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fvo(5) = fo(5) + vc(5)
      fvo(6) = fo(6) + vc(6)

      pv(1) = dble( vc(5))
      pv(2) = dble( vc(6))
      pv(3) = dimag(vc(6))
      pv(4) = dimag(vc(5))

      pf(0) = dble( fvo(5))
      pf(1) = dble( fvo(6))
      pf(2) = dimag(fvo(6))
      pf(3) = dimag(fvo(5))
      pf2 = pf(0)**2 - pf(1)**2 - pf(2)**2 - pf(3)**2

      if ( idecay.eq.1 ) then
         pdotpG = pv(1)*pf(0) - pv(2)*pf(1) - pv(3)*pf(2) - pv(4)*pf(3)
         AdotpG = vc(1)*pf(0) - vc(2)*pf(1) - vc(3)*pf(2) - vc(4)*pf(3)
      else if ( idecay.eq.-1 ) then
         pdotpG =  pv(1)*dble( fo(5)) - pv(2)*dble( fo(6))
     &           - pv(3)*dimag(fo(6)) - pv(4)*dimag(fo(5))
         AdotpG =  vc(1)*dble( fo(5)) - vc(2)*dble( fo(6))
     &           - vc(3)*dimag(fo(6)) - vc(4)*dimag(fo(5))
      else
         write(6,*) 'error in idecay of FVOGLD'
         stop
      end if

      p14p = dble(vc(5)) + dimag(vc(5))
      p14m = dble(vc(5)) - dimag(vc(5))
      p23p = vc(6)
      p23m = dconjg(vc(6))
      A14p = vc(1) + vc(4)
      A14m = vc(1) - vc(4)
      A23p = vc(2) + vc(3)*ci
      A23m = vc(2) - vc(3)*ci

      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
      d = d*idecay

      sl1 =  (p14p*AdotpG - A14p*pdotpG)*fo(3)
     &      +(p23p*AdotpG - A23p*pdotpG)*fo(4)
      sl2 =  (p23m*AdotpG - A23m*pdotpG)*fo(3)
     &      +(p14m*AdotpG - A14m*pdotpG)*fo(4)

      if ( gc(2).ne.cZero ) then
         sr1 =  (p14m*AdotpG - A14m*pdotpG)*fo(1)
     &         -(p23p*AdotpG - A23p*pdotpG)*fo(2)
         sr2 = -(p23m*AdotpG - A23m*pdotpG)*fo(1)
     &         +(p14p*AdotpG - A14p*pdotpG)*fo(2)

         fvo(1) = ( gc(2)*(  (pf(0)+pf(3))*sr1         +fvo(6)*sr2 )
     &             +gc(1)*fmass*sl1 )*d
         fvo(2) = ( gc(2)*( dconjg(fvo(6))*sr1 + (pf(0)-pf(3))*sr2 )
     &             +gc(1)*fmass*sl2 )*d
         fvo(3) = ( gc(1)*(  (pf(0)-pf(3))*sl1         -fvo(6)*sl2 )
     &             +gc(2)*fmass*sr1 )*d
         fvo(4) = ( gc(1)*(-dconjg(fvo(6))*sl1 + (pf(0)+pf(3))*sl2 )
     &             +gc(2)*fmass*sr2 )*d

      else
         d = d*gc(1)
         fvo(1) = d*fmass*sl1
         fvo(2) = d*fmass*sl2
         fvo(3) = d*(  (pf(0)-pf(3))*sl1         -fvo(6)*sl2)
         fvo(4) = d*(-dconjg(fvo(6))*sl1 + (pf(0)+pf(3))*sl2)
      end if
c
      return          
      end
      subroutine fvogox(fo,vc,gc,fmass,fwidth , fvo)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external fermion and a vector boson in the case when 
c they are all color octets.
c
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvo(6)         : off-shell fermion             <fo,v,f'|
c
      implicit none
      double complex fo(6),vc(6),gc(2),fvo(6),sl1,sl2,sr1,sr2,d
      double precision pf(0:3),fmass,fwidth,pf2

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      fvo(5) = fo(5)+vc(5)
      fvo(6) = fo(6)+vc(6)

      pf(0) = dble( fvo(5))
      pf(1) = dble( fvo(6))
      pf(2) = dimag(fvo(6))
      pf(3) = dimag(fvo(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
      sl1 =   (vc(1)+      vc(4))*fo(3)
     &      + (vc(2)+cImag*vc(3))*fo(4)
      sl2 =   (vc(2)-cImag*vc(3))*fo(3)
     &      + (vc(1)-      vc(4))*fo(4)

      if ( gc(2).ne.cZero ) then
         sr1 =   (vc(1)-      vc(4))*fo(1)
     &         - (vc(2)+cImag*vc(3))*fo(2)
         sr2 = - (vc(2)-cImag*vc(3))*fo(1)
     &         + (vc(1)+      vc(4))*fo(2)

         fvo(1) = ( gc(2)*( (pf(0)+pf(3))*sr1  +        fvo(6)*sr2)
     &             +gc(1)*fmass*sl1 )*d
         fvo(2) = ( gc(2)*( dconjg(fvo(6))*sr1 + (pf(0)-pf(3))*sr2)
     &             +gc(1)*fmass*sl2 )*d
         fvo(3) = ( gc(1)*( (pf(0)-pf(3))*sl1  -        fvo(6)*sl2)
     &             +gc(2)*fmass*sr1 )*d
         fvo(4) = ( gc(1)*(-dconjg(fvo(6))*sl1 + (pf(0)+pf(3))*sl2)
     &             +gc(2)*fmass*sr2 )*d

      else
         d = d * gc(1)
         fvo(1) = fmass*sl1*d
         fvo(2) = fmass*sl2*d
         fvo(3) = (  (pf(0)-pf(3))*sl1 -        fvo(6)*sl2)*d
         fvo(4) = (-dconjg(fvo(6))*sl1 + (pf(0)+pf(3))*sl2)*d
      end if
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine hiocxx(fic,foc,gc,smass,swidth , hioc)                
c                                                                       
c this subroutine computes an off-shell scalar current from an external 
c antifermion pair.                                                     
c                                                                       
c input:                                                                
c       complex fic(6)         : flow-in  antifermion              |fic>
c       complex foc(6)         : flow-out antifermion              <foc|
c       complex gc(2)          : coupling constants                 gchf
c       real    smass          : mass  of output scalar s               
c       real    swidth         : width of output scalar s               
c                                                                       
c output:                                                               
c       complex hioc(3)        : scalar current           j(<fic|s|foc>)
c                                                                       
      implicit none
      double complex fic(6),foc(6),hioc(3),gc(2),dn
      double precision q(0:3),smass,swidth,q2
c
      hioc(2) = foc(5)-fic(5)
      hioc(3) = foc(6)-fic(6)

      q(0) = dble( hioc(2))
      q(1) = dble( hioc(3))
      q(2) = dimag(hioc(3))
      q(3) = dimag(hioc(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)

      dn = -dcmplx( q2-smass**2, smass*swidth )

      hioc(1) = (  gc(1)*(foc(1)*fic(1)+foc(2)*fic(2))
     &           + gc(2)*(foc(3)*fic(3)+foc(4)*fic(4)) )/dn
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine hiogld(fi,fo,gc,smass,swidth,mNLSP,idecay , hio)
c
c This subroutine computes an off-shell scalar current for the NLSP-
c Goldstino vertex from the external fermion pair. The h.c. of the NLSP
c decay is handled via the input parameter idecay.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gc(2)          : coupling constants                  gsf
c       real    smass          : mass  of output scalar s
c       real    swidth         : width of output scalar s
c       real    mNLSP          : mass of NLSP
c       integer idecay         :  1 for NLSP decay to Goldstino
c                              : -1 for Goldstino to NLSP (h.c. of above)
c
c output:
c       complex hio(3)         : scalar current          j^mu(<fo|s|fi>)
c
      implicit none
      double complex fi(6), fo(6), gc(2), hio(3)
      double complex dn, p14p, p14m, p23p, p23m
      double precision q(0:3), smass, swidth, mNLSP, q2
      double precision pG(1:4)
      integer idecay

      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      hio(2) = -fi(5) + fo(5)
      hio(3) = -fi(6) + fo(6)

      if ( idecay.ne.1 .or. idecay.ne.-1 ) then
         write(6,*) 'error in idecay of HIOGLD'
         stop
      end if

      q(0) = dble( hio(2))
      q(1) = dble( hio(3))
      q(2) = dimag(hio(3))
      q(3) = dimag(hio(2))
      q2 = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2

      p14p = q(0) + q(3)
      p14m = q(0) - q(3)
      p23p = q(1) + q(2)*ci
      p23m = q(1) - q(2)*ci

      dn = -dcmplx( q2-smass**2, smass*swidth )

      hio(1) = gc(1)*( ( ( fo(3)*p14p + fo(4)*p23p )*fi(1)
     &                  +( fo(3)*p23m + fo(4)*p14m )*fi(2) )*mNLSP
     &                -( fo(1)*fi(1) + fo(2)*fi(2) )*idecay*smass**2 )

      if ( gc(2).ne.cZero ) then
         hio(1) = hio(1) + gc(2) *
     &            ( ( ( fo(1)*p14m - fo(2)*p23p )*fi(3)
     &               -( fo(1)*p23m - fo(2)*p14p )*fi(4) )*mNLSP
     &             -( fo(3)*fi(3) + fo(4)*fi(4) )*idecay*smass**2 )
      end if

      hio(1) = hio(1)/dn
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine ioscxx(fic,foc,sc,gc , vertex)
c
c This subroutine computes an amplitude of the antifermion-antifermion- 
c scalar coupling.                                                      
c                                                                       
c input:                                                                
c       complex fic(6)         : flow-in  antifermion              |fic>
c       complex foc(6)         : flow-out antifermion              <foc|
c       complex sc(3)          : input    scalar                      s 
c       complex gc(2)          : coupling constants                 gchf
c                                                                       
c output:                                                               
c       complex vertex         : amplitude                   <foc|s|fic>
c
      implicit none
      double complex fic(6),foc(6),sc(3),gc(2),vertex
c
      vertex = sc(1)*( gc(1)*(fic(1)*foc(1)+fic(2)*foc(2))
     &                +gc(2)*(fic(3)*foc(3)+fic(4)*foc(4)) )
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine iosgld(fi,fo,sc,gc,smass,mNLSP,idecay , vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-scalar
c SUSY Goldstino coupling. In this routine, the NLSP is decaying to a
c boson and a Goldstino. The h.c. of the NLSP decay is handled via the
c input parameter idecay.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex sc(3)          : input    scalar                      s
c       complex gc(2)          : coupling constants                  gsf
c       real    mNLSP          : mass of the NLSP
c       real    smass          : mass of the scalar boson
c       integer idecay         :  1 for NLSP decay to Goldstino
c                              : -1 for Goldstino to NLSP (h.c. of above)
c
c output:
c       complex vertex         : amplitude                     <fo|s|fi>
c
      implicit none
      double complex  fi(6), fo(6), gc(2), sc(3), vertex
      double complex  p14p, p14m, p23p, p23m
      double precision  mNLSP, smass
      integer idecay

      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      if ( idecay.ne.1 .or. idecay.ne.-1 ) then
         write(6,*) 'error in idecay of IOSGLD'
         stop
      end if

      p14p = dble(sc(2)) + dimag(sc(2))
      p14m = dble(sc(2)) - dimag(sc(2))
      p23p = dble(sc(3)) + dimag(sc(3))*ci
      p23m = dble(sc(3)) - dimag(sc(3))*ci

      vertex = gc(1) *
     &         ( ( ( fo(3)*p14p + fo(4)*p23p )*fi(1)
     &            +( fo(3)*p23m + fo(4)*p14m )*fi(2) )*mNLSP
     &          -( fo(1)*fi(1) + fo(2)*fi(2) )*idecay*smass**2 )

      if ( gc(2).ne.cZero ) then
         vertex = vertex + gc(2) *
     &            ( ( ( fo(1)*p14m - fo(2)*p23p )*fi(3)
     &               -( fo(1)*p23m - fo(2)*p14p )*fi(4) )*mNLSP
     &             -( fo(3)*fi(3) + fo(4)*fi(4) )*idecay*smass**2 )
      end if

      vertex = vertex * sc(1)
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine iovcxx(fic,foc,vc,gc , vertex)
c
c this subroutine computes an amplitude of the antifermion-antifermion- 
c vector coupling.                                                      
c                                                                       
c input:                                                                
c       complex fic(6)         : flow-in  antifermion              |fic>
c       complex foc(6)         : flow-out antifermion              <foc|
c       complex vc(6)          : input    vector                      v 
c       complex gc(2)          : coupling constants                  gvf
c                                                                       
c output:                                                               
c       complex vertex         : amplitude                   <foc|v|fic>
c
      implicit none
      double complex fic(6),foc(6),vc(6),gc(2),vertex

      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      vertex = - gc(1)*( (foc(1)*fic(3)+foc(2)*fic(4))*vc(1)
     &                  -(foc(1)*fic(4)+foc(2)*fic(3))*vc(2)
     &                  +(foc(1)*fic(4)-foc(2)*fic(3))*vc(3)*cImag
     &                  -(foc(1)*fic(3)-foc(2)*fic(4))*vc(4)       )

      if ( gc(2).ne.cZero ) then
         vertex = vertex
     &            - gc(2)*( (foc(3)*fic(1)+foc(4)*fic(2))*vc(1)
     &                     +(foc(3)*fic(2)+foc(4)*fic(1))*vc(2)
     &                     -(foc(3)*fic(2)-foc(4)*fic(1))*vc(3)*cImag
     &                     +(foc(3)*fic(1)-foc(4)*fic(2))*vc(4)       )
      end if
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine iovgld(fi,fo,vc,gc,idecay , vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c SUSY Goldstino coupling. In this routine, the NLSP is decaying to a
c boson and a Goldstino. The h.c. of the NLSP decay is handled via the
c input parameter idecay (picks out correct Goldstino momentum).
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       integer idecay         :  1 for NLSP decay to Goldstino
c                              : -1 for Goldstino to NLSP (h.c. of above)
c
c output:
c       complex vertex         : amplitude                     <fo|v|fi>
c
      implicit none
      double complex  fi(6), fo(6), gc(2), vc(6), vertex
      double complex  AdotpG, A14p, A14m, A23p, A23m
      double complex  p14p, p14m, p23p, p23m
      double precision  pdotpG
      integer idecay

      double precision rOne
      parameter( rOne = 1.0d0 )
      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      if ( idecay.eq.1 ) then
         pdotpG =  dble( vc(5))*dble( fo(5))
     &           - dble( vc(6))*dble( fo(6))
     &           - dimag(vc(6))*dimag(fo(6))
     &           - dimag(vc(5))*dimag(fo(5))
         AdotpG =  vc(1)*dble( fo(5)) - vc(2)*dble( fo(6))
     &           - vc(3)*dimag(fo(6)) - vc(4)*dimag(fo(5))
      else if ( idecay.eq.-1 ) then
         pdotpG =  dble( vc(5))*dble( fi(5))
     &           - dble( vc(6))*dble( fi(6))
     &           - dimag(vc(6))*dimag(fi(6))
     &           - dimag(vc(5))*dimag(fi(5))
         AdotpG =  vc(1)*dble( fi(5)) - vc(2)*dble( fi(6))
     &           - vc(3)*dimag(fi(6)) - vc(4)*dimag(fi(5))
      else
         write(6,*) 'error in idecay of IOVGLD'
         stop
      end if

      p14p = dble(vc(5)) + dimag(vc(5))
      p14m = dble(vc(5)) - dimag(vc(5))
      p23p = vc(6)
      p23m = dconjg(vc(6))
      A14p = vc(1) + vc(4)
      A14m = vc(1) - vc(4)
      A23p = vc(2) + vc(3)*ci
      A23m = vc(2) - vc(3)*ci

      vertex = gc(1)*( ( ( fo(3)*p14p + fo(4)*p23p )*fi(1)
     &                  +( fo(3)*p23m + fo(4)*p14m )*fi(2) )*AdotpG
     &                -( ( fo(3)*A14p + fo(4)*A23p )*fi(1)
     &                  +( fo(3)*A23m + fo(4)*A14m )*fi(2) )*pdotpG )

      if ( gc(2).ne.cZero ) then
         vertex = vertex
     &          + gc(2)*( ( (fo(1)*p14m - fo(2)*p23p )*fi(3)
     &                     -(fo(1)*p23m - fo(2)*p14p )*fi(4) )*AdotpG
     &                   -( (fo(1)*A14m - fo(2)*A23p )*fi(3)
     &                     -(fo(1)*A23m - fo(2)*A14p )*fi(4) )*pdotpG )
      end if

      vertex = vertex * idecay
c
      return
      end
      subroutine iovgox(fi,fo,vc,gc , vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c coupling in the case when they are all color octets.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c
c output:
c       complex vertex         : amplitude                     <fo|v|fi>
c     
      implicit none
      double complex fi(6),fo(6),gc(2),vc(6),vertex

      double precision rZero, rOne
      parameter( rZero = 0.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      vertex =  gc(1)*( (fo(3)*fi(1)+fo(4)*fi(2))*vc(1)
     &                 +(fo(3)*fi(2)+fo(4)*fi(1))*vc(2)
     &                 -(fo(3)*fi(2)-fo(4)*fi(1))*vc(3)*cImag
     &                 +(fo(3)*fi(1)-fo(4)*fi(2))*vc(4)        )

      if ( gc(2).ne.cZero ) then
         vertex = vertex
     &          + gc(2)*( (fo(1)*fi(3)+fo(2)*fi(4))*vc(1)
     &                   -(fo(1)*fi(4)+fo(2)*fi(3))*vc(2)
     &                   +(fo(1)*fi(4)-fo(2)*fi(3))*vc(3)*cImag
     &                   -(fo(1)*fi(3)-fo(2)*fi(4))*vc(4)        )
      end if
c
      return
      end
c ----------------------------------------------------------------------
c
      subroutine jiocxx(fic,foc,gc,vmass,vwidth , jioc)
c
c This subroutine computes an off-shell vector current from an external 
c antifermion pair. The vector boson propagator is given in Feynman     
c gauge for a massless vector and in unitary gauge for a massive vector.
c                                                                       
c input:                                                                
c       complex fic(6)         : flow-in  antifermion              |fic>
c       complex foc(6)         : flow-out antifermion              <foc|
c       complex gc(2)          : coupling constants                  gvf
c       real    vmass          : mass  of output vector v               
c       real    vwidth         : width of output vector v               
c                                                                       
c output:                                                               
c       complex jioc(6)        : vector current        j^mu(<foc|v|fic>)
c
      implicit none
      double complex fic(6),foc(6),gc(2),jioc(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2,dd
      double complex cm2 ! mass**2- I Gamma mass (Benj&Claude)

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      jioc(5) = foc(5)-fic(5)
      jioc(6) = foc(6)-fic(6)

      q(0) = dble( jioc(5))
      q(1) = dble( jioc(6))
      q(2) = dimag(jioc(6))
      q(3) = dimag(jioc(5))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

      if ( vmass.ne.rZero ) then

         d = -rOne/dcmplx( q2-vm2, vmass*vwidth )
c  for the running width, use below instead of the above d.
c         d = -rOne/dcmplx( q2-vm2, max(vwidth*q2/vmass,rZero) )

         if ( gc(2).ne.cZero ) then
            c0=  gc(2)*( foc(3)*fic(1)+foc(4)*fic(2))
     &          +gc(1)*( foc(1)*fic(3)+foc(2)*fic(4))
            c1= -gc(2)*( foc(3)*fic(2)+foc(4)*fic(1))
     &          +gc(1)*( foc(1)*fic(4)+foc(2)*fic(3))
            c2=( gc(2)*( foc(3)*fic(2)-foc(4)*fic(1)) 
     &          +gc(1)*(-foc(1)*fic(4)+foc(2)*fic(3)))*cImag
            c3=  gc(2)*(-foc(3)*fic(1)+foc(4)*fic(2))
     &          +gc(1)*( foc(1)*fic(3)-foc(2)*fic(4))
         else
            d = d*gc(1)
            c0 =  foc(1)*fic(3)+foc(2)*fic(4)
            c1 =  foc(1)*fic(4)+foc(2)*fic(3)
            c2 =(-foc(1)*fic(4)+foc(2)*fic(3))*cImag
            c3 =  foc(1)*fic(3)-foc(2)*fic(4)
         end if
	 
c     Fabio's implementation of the fixed width (Benj&Claude)
         cm2=dcmplx( vm2, -vmass*vwidth )
c         cs = (q(0)*c0-q(1)*c1-q(2)*c2-q(3)*c3)/vm2
         cs = (q(0)*c0-q(1)*c1-q(2)*c2-q(3)*c3)/cm2
	 
         jioc(1) = (c0-cs*q(0))*d
         jioc(2) = (c1-cs*q(1))*d
         jioc(3) = (c2-cs*q(2))*d
         jioc(4) = (c3-cs*q(3))*d

      else

         d = dcmplx( -rOne/q2, rZero )
         if ( gc(2).ne.cZero ) then
            jioc(1) = ( gc(2)*( foc(3)*fic(1)+foc(4)*fic(2))
     &                 +gc(1)*( foc(1)*fic(3)+foc(2)*fic(4)) )*d
            jioc(2) = (-gc(2)*( foc(3)*fic(2)+foc(4)*fic(1))
     &                 +gc(1)*( foc(1)*fic(4)+foc(2)*fic(3)) )*d
            jioc(3) = ( gc(2)*( foc(3)*fic(2)-foc(4)*fic(1))
     &                 +gc(1)*(-foc(1)*fic(4)+foc(2)*fic(3)) )
     &                *d*cImag
            jioc(4) = ( gc(2)*(-foc(3)*fic(1)+foc(4)*fic(2))
     &                 +gc(1)*( foc(1)*fic(3)-foc(2)*fic(4)) )*d
         else
            d = d*gc(1)
            jioc(1) = ( foc(1)*fic(3)+foc(2)*fic(4))*d
            jioc(2) = ( foc(1)*fic(4)+foc(2)*fic(3))*d
            jioc(3) = (-foc(1)*fic(4)+foc(2)*fic(3))*d*cImag
            jioc(4) = ( foc(1)*fic(3)-foc(2)*fic(4))*d
         end if

      end if
c
      return
      end
c
c ----------------------------------------------------------------------
c
      subroutine jiogld(fi,fo,gc,vmass,vwidth,idecay , jio)
c
c This subroutine computes an off-shell vector current for the NLSP-
c Goldstino vertex from an external fermion pair. The vector boson 
c propagator is given in feynman gauge for a massless vector and in 
c unitary gauge for a massive vector. The h.c. of the NLSP decay is
c handled via the input parameter idecay (picks out correct
c Goldstino momentum).
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gc(2)          : coupling constants                  gvf
c       real    vmass          : mass  of output vector v
c       real    vwidth         : width of output vector v
c       integer idecay         :  1 for NLSP decay to Goldstino
c                              : -1 for Goldstino to NLSP (h.c. of above)
c
c output:
c       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
c
      implicit none
      double complex  fi(6), fo(6), gc(2), jio(6), c0, c1, c2, c3, cs
      double complex  d, dum, p14p, p14m, p23p, p23m
      double precision  q(0:3), vmass, vwidth, q2, vm2, dd
      double precision  pG(1:4), pdotpG
      integer idecay

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      if ( idecay.eq.1 ) then
         pG(1) =  dble(fo(5))
         pG(2) =  dble(fo(6))
         pG(3) = dimag(fo(6))
         pG(4) = dimag(fo(5))
      else if ( idecay.eq.-1 ) then
         pG(1) =  dble(fi(5))
         pG(2) =  dble(fi(6))
         pG(3) = dimag(fi(6))
         pG(4) = dimag(fi(5))
      else
         write(6,*) 'error in idecay of JIOGLD'
         stop
      end if

      jio(5) = fo(5) - fi(5)
      jio(6) = fo(6) - fi(6)

      q(0) = dble( jio(5))
      q(1) = dble( jio(6))
      q(2) = dimag(jio(6))
      q(3) = dimag(jio(5))
      q2  = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2
      vm2 = vmass**2

      pdotpG = q(0)*pG(1) - q(1)*pG(2) - q(2)*pG(3) - q(3)*pG(4)

      p14p = q(0) + q(3)
      p14m = q(0) - q(3)
      p23p = jio(6)
      p23m = dconjg(jio(6))

      if ( vmass.ne.rZero ) then

         d = rOne/dcmplx( q2-vm2, vmass*vwidth )
         d = d*idecay
c  for the running width, use below instead of the above d.
c         d = rOne/dcmplx( q2-vm2, max(vwidth*q2/vmass,rZero) )

         if ( gc(2).ne.cZero ) then
            dum =  ( (fo(3)*p14p + fo(4)*p23p)*fi(1)
     &              +(fo(3)*p23m + fo(4)*p14m)*fi(2) )*gc(1)
     &           + ( (fo(1)*p14m - fo(2)*p23p)*fi(3)
     &              -(fo(1)*p23m - fo(2)*p14p)*fi(4) )*gc(2)

            c0 =  dum*pG(1)
     &           -pdotpG*( gc(1)*( fo(3)*fi(1) + fo(4)*fi(2) )
     &                    +gc(2)*( fo(1)*fi(3) + fo(2)*fi(4) ) )
            c1 =  dum*pG(2)
     &           -pdotpG*(-gc(1)*( fo(4)*fi(1) + fo(3)*fi(2) )
     &                    +gc(2)*( fo(2)*fi(3) + fo(1)*fi(4) ) )
            c2 =  dum*pG(3)
     &           -pdotpG*( gc(1)*(-fo(4)*fi(1) + fo(3)*fi(2) )
     &                    +gc(2)*( fo(2)*fi(3) - fo(1)*fi(4) ) )*ci
            c3 =  dum*pG(4)
     &           -pdotpG*( gc(1)*(-fo(3)*fi(1) + fo(4)*fi(2) )
     &                    +gc(2)*( fo(1)*fi(3) - fo(2)*fi(4) ) )
         else
            d = d*gc(1)
            dum =  (fo(3)*p14p + fo(4)*p23p)*fi(1)
     &            +(fo(3)*p23m + fo(4)*p14m)*fi(2)

            c0 = dum*pG(1) - ( fo(3)*fi(1) + fo(4)*fi(2) )*pdotpG
            c1 = dum*pG(2) + ( fo(4)*fi(1) + fo(3)*fi(2) )*pdotpG
            c2 = dum*pG(3) + ( fo(4)*fi(1) - fo(3)*fi(2) )*pdotpG*ci
            c3 = dum*pG(4) + ( fo(3)*fi(1) - fo(4)*fi(2) )*pdotpG
         end if

         cs = (q(0)*c0 - q(1)*c1 - q(2)*c2 - q(3)*c3) / vm2

         jio(1) = (c0-cs*q(0))*d
         jio(2) = (c1-cs*q(1))*d
         jio(3) = (c2-cs*q(2))*d
         jio(4) = (c3-cs*q(3))*d

      else
         dd = idecay*rOne/q2

         if ( gc(2).ne.cZero ) then
            dum =  ( (fo(3)*p14p + fo(4)*p23p)*fi(1)
     &              +(fo(3)*p23m + fo(4)*p14m)*fi(2) )*gc(1)
     &           + ( (fo(1)*p14m - fo(2)*p23p)*fi(3)
     &              -(fo(1)*p23m - fo(2)*p14p)*fi(4) )*gc(2)

            jio(1) = ( dum*pG(1) - pdotpG*(
     &                 gc(1)*( fo(3)*fi(1) + fo(4)*fi(2) )
     &                +gc(2)*( fo(1)*fi(3) + fo(2)*fi(4) ) ) )*dd
            jio(2) = ( dum*pG(2) - pdotpG*(
     &                -gc(1)*( fo(4)*fi(1) + fo(3)*fi(2) )
     &                +gc(2)*( fo(2)*fi(3) + fo(1)*fi(4) ) ) )*dd
            jio(3) = ( dum*pG(3) - pdotpG*ci*(
     &                 gc(1)*(-fo(4)*fi(1) + fo(3)*fi(2) )
     &                +gc(2)*( fo(2)*fi(3) - fo(1)*fi(4) ) ) )*dd
            jio(4) = ( dum*pG(4) - pdotpG*(
     &                 gc(1)*(-fo(3)*fi(1) + fo(4)*fi(2) )
     &                +gc(2)*( fo(1)*fi(3) - fo(2)*fi(4) ) ) )*dd

         else
            dd = dd*gc(1)
            dum =  (fo(3)*p14p + fo(4)*p23p)*fi(1)
     &            +(fo(3)*p23m + fo(4)*p14m)*fi(2)

            jio(1)=dd*(dum*pG(1) - pdotpG*(fo(3)*fi(1) + fo(4)*fi(2)))
            jio(2)=dd*(dum*pG(2) + pdotpG*(fo(4)*fi(1) + fo(3)*fi(2)))
            jio(3)=dd*(dum*pG(3) + ci*pdotpG*(fo(4)*fi(1)-fo(3)*fi(2)))
            jio(4)=dd*(dum*pG(4) + pdotpG*(fo(3)*fi(1) - fo(4)*fi(2)))
         end if
      end if
c
      return
      end
      subroutine jiogox(fi,fo,gc,vmass,vwidth , jio)
c
c This subroutine computes an off-shell vector current from an external
c fermion pair.  The vector boson propagator is given in Feynman gauge
c for a massless vector and in unitary gauge for a massive vector in the 
c case when they are all color octets.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gc(2)          : coupling constants                  gvf
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
c     
      implicit none
      double complex fi(6),fo(6),gc(2),jio(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)


      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      jio(5) = fo(5)-fi(5)
      jio(6) = fo(6)-fi(6)

      q(0) = dble( jio(5))
      q(1) = dble( jio(6))
      q(2) = dimag(jio(6))
      q(3) = dimag(jio(5))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

      if ( vmass.ne.rZero ) then

         d = rOne/dcmplx( q2-vm2, vmass*vwidth )
c     For the running width, use below instead of the above d.
c     d = rOne/dcmplx( q2-vm2, max(vwidth*q2/vmass,rZero) )


         if ( gc(2).ne.cZero ) then
            c0 =  gc(1)*( fo(3)*fi(1)+fo(4)*fi(2))
     &           +gc(2)*( fo(1)*fi(3)+fo(2)*fi(4))
            c1 = -gc(1)*( fo(3)*fi(2)+fo(4)*fi(1))
     &           +gc(2)*( fo(1)*fi(4)+fo(2)*fi(3))
            c2 =( gc(1)*( fo(3)*fi(2)-fo(4)*fi(1))
     &           +gc(2)*(-fo(1)*fi(4)+fo(2)*fi(3)))*cImag
            c3 =  gc(1)*(-fo(3)*fi(1)+fo(4)*fi(2))
     &           +gc(2)*( fo(1)*fi(3)-fo(2)*fi(4))
         else
            d = d*gc(1)
            c0 =   fo(3)*fi(1)+fo(4)*fi(2)
            c1 =  -fo(3)*fi(2)-fo(4)*fi(1)
            c2 = ( fo(3)*fi(2)-fo(4)*fi(1))*cImag
            c3 =  -fo(3)*fi(1)+fo(4)*fi(2)
         end if

c     Fabio's implementation of the fixed width
         cm2=dcmplx( vm2, -vmass*vwidth )
c     cs = (q(0)*c0-q(1)*c1-q(2)*c2-q(3)*c3)/vm2
         cs = (q(0)*c0-q(1)*c1-q(2)*c2-q(3)*c3)/cm2
         jio(1) = (c0-cs*q(0))*d
         jio(2) = (c1-cs*q(1))*d
         jio(3) = (c2-cs*q(2))*d
         jio(4) = (c3-cs*q(3))*d

      else

         d = dcmplx( rOne/q2, rZero )
         if ( gc(2).ne.cZero ) then
            jio(1) = ( gc(1)*( fo(3)*fi(1)+fo(4)*fi(2))
     &                +gc(2)*( fo(1)*fi(3)+fo(2)*fi(4)) )*d
            jio(2) = (-gc(1)*( fo(3)*fi(2)+fo(4)*fi(1))
     &                +gc(2)*( fo(1)*fi(4)+fo(2)*fi(3)) )*d
            jio(3) = ( gc(1)*( fo(3)*fi(2)-fo(4)*fi(1))
     &                +gc(2)*(-fo(1)*fi(4)+fo(2)*fi(3)))
     &               *d*cImag
            jio(4) = ( gc(1)*(-fo(3)*fi(1)+fo(4)*fi(2))
     &                +gc(2)*( fo(1)*fi(3)-fo(2)*fi(4)) )*d
         else
            d = d*gc(1)
            jio(1) =  ( fo(3)*fi(1)+fo(4)*fi(2))*d
            jio(2) = -( fo(3)*fi(2)+fo(4)*fi(1))*d
            jio(3) =  ( fo(3)*fi(2)-fo(4)*fi(1))*d*cImag
            jio(4) =  (-fo(3)*fi(1)+fo(4)*fi(2))*d
         end if

      end if
c
      return
      end
      subroutine t2xxxx(p,tmass,nhel,nst ,ft)
c
c This subroutine computes k^mu*e^nu where e is delta(i,nhel).
c It is used to test for gauge invariance of the tensor routines.
c
c input:
c       real    p(0:3)             : four-momentum of tensor boson
c       real    tmass              : mass          of tensor boson
c       integer nhel = =-2..2        : construction of e^nu
c       integer nst  = -1 or 1     : +1 for final, -1 for initial
c
c output:
c       complex tc(6,4)            : tensor wavefunction       epsilon^mu^nu(t)
c     
      implicit none
      double complex ft(18),tc(6,4), ep(4), em(4)
      double precision p(0:3), tmass
      integer nhel, nst, i, j

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
c
      tc(5,1) = dcmplx(p(0),p(3))*nst
      tc(6,1) = dcmplx(p(1),p(2))*nst

      ep(1) = dcmplx( p(0), rZero )
      ep(2) = dcmplx( p(1), rZero )
      ep(3) = dcmplx( p(2), rZero )
      ep(4) = dcmplx( p(3), rZero )

      if ( nhel.eq.1 ) then
         em(1) = dcmplx( rOne , rZero )
         em(2) = dcmplx( rZero, rZero )
         em(3) = dcmplx( rZero, rZero )
         em(4) = dcmplx( rZero, rZero )
      else if ( nhel.eq.-1 ) then
         em(1) = dcmplx( rZero, rZero )
         em(2) = dcmplx( rOne , rZero )
         em(3) = dcmplx( rZero, rZero )
         em(4) = dcmplx( rZero, rZero )
      else if ( nhel.eq.2 ) then
         em(1) = dcmplx( rZero, rZero )
         em(2) = dcmplx( rZero, rZero )
         em(3) = dcmplx( rOne , rZero )
         em(4) = dcmplx( rZero, rZero )
      else if ( nhel.eq.-2.or.nhel.eq.0 ) then
         em(1) = dcmplx( rZero, rZero )
         em(2) = dcmplx( rZero, rZero )
         em(3) = dcmplx( rZero, rZero )
         em(4) = dcmplx( rOne , rZero )
      end if

      do j = 1,4
         do i = 1,4
            tc(i,j) = ep(i)*em(j)
         end do
      end do
c
      ft(1)=tc(1,1)
      ft(2)=tc(1,2)
      ft(3)=tc(1,3)
      ft(4)=tc(1,4)
      ft(5)=tc(2,1)
      ft(6)=tc(2,2)
      ft(7)=tc(2,3)
      ft(8)=tc(2,4)
      ft(9)=tc(3,1)
      ft(10)=tc(3,2)
      ft(11)=tc(3,3)
      ft(12)=tc(3,4)
      ft(13)=tc(4,1)
      ft(14)=tc(4,2)
      ft(15)=tc(4,3)
      ft(16)=tc(4,4)
      ft(17)=tc(5,1)
      ft(18)=tc(6,1)

      return
      end
      subroutine txxxx2(p,tmass,nhel,nst , tc)
c
c This subroutine computes k^mu*e^nu where e is delta(i,nhel).
c It is used to test for gauge invariance of the tensor routines.
c
c input:
c       real    p(0:3)             : four-momentum of tensor boson
c       real    tmass              : mass          of tensor boson
c       integer nhel = 1..4        : construction of e^nu
c       integer nst  = -1 or 1     : +1 for final, -1 for initial
c
c output:
c       complex tc(6,4)            : tensor wavefunction       epsilon^mu^nu(t)
c     
      implicit none
      double complex tc(6,4), ep(4), em(4)
      double precision p(0:3), tmass
      integer nhel, nst, i, j

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
c
      tc(5,1) = dcmplx(p(0),p(3))*nst
      tc(6,1) = dcmplx(p(1),p(2))*nst

      ep(1) = dcmplx( p(0), rZero )
      ep(2) = dcmplx( p(1), rZero )
      ep(3) = dcmplx( p(2), rZero )
      ep(4) = dcmplx( p(3), rZero )

      if ( nhel.eq.1 ) then
         em(1) = dcmplx( rOne , rZero )
         em(2) = dcmplx( rZero, rZero )
         em(3) = dcmplx( rZero, rZero )
         em(4) = dcmplx( rZero, rZero )
      else if ( nhel.eq.2 ) then
         em(1) = dcmplx( rZero, rZero )
         em(2) = dcmplx( rOne , rZero )
         em(3) = dcmplx( rZero, rZero )
         em(4) = dcmplx( rZero, rZero )
      else if ( nhel.eq.3 ) then
         em(1) = dcmplx( rZero, rZero )
         em(2) = dcmplx( rZero, rZero )
         em(3) = dcmplx( rOne , rZero )
         em(4) = dcmplx( rZero, rZero )
      else if ( nhel.eq.4 ) then
         em(1) = dcmplx( rZero, rZero )
         em(2) = dcmplx( rZero, rZero )
         em(3) = dcmplx( rZero, rZero )
         em(4) = dcmplx( rOne , rZero )
      end if

      do j = 1,4
         do i = 1,4
            tc(i,j) = ep(i)*em(j)
         end do
      end do
c
      return
      end
      subroutine txxxxx(p,tmass,nhel,nst , tc)
c
c This subroutine computes a TENSOR wavefunction.
c
c input:
c       real    p(0:3)         : four-momentum of tensor boson
c       real    tmass          : mass          of tensor boson
c       integer nhel           : helicity      of tensor boson
c                = -2,-1,0,1,2 : (0 is forbidden if tmass=0.0)
c       integer nst  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex tc(18)         : tensor wavefunction    epsilon^mu^nu(t)
c     
      implicit none
      double precision p(0:3), tmass
      integer nhel, nst
      double complex tc(18)

      double complex ft(6,4), ep(4), em(4), e0(4)
      double precision pt, pt2, pp, pzpt, emp, sqh, sqs
      integer i, j

      double precision rZero, rHalf, rOne, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0 )

      integer stdo
      parameter( stdo = 6 )

      sqh = sqrt(rHalf)
      sqs = sqrt(rHalf/3.d0)

      pt2 = p(1)**2 + p(2)**2
      pp = min(p(0),sqrt(pt2+p(3)**2))
      pt = min(pp,sqrt(pt2))

      ft(5,1) = dcmplx(p(0),p(3))*nst
      ft(6,1) = dcmplx(p(1),p(2))*nst

      if ( nhel.ge.0 ) then 
c construct eps+
         if ( pp.eq.rZero ) then
            ep(1) = dcmplx( rZero )
            ep(2) = dcmplx( -sqh )
            ep(3) = dcmplx( rZero , nst*sqh )
            ep(4) = dcmplx( rZero )
         else
            ep(1) = dcmplx( rZero )
            ep(4) = dcmplx( pt/pp*sqh )
            if ( pt.ne.rZero ) then
               pzpt = p(3)/(pp*pt)*sqh
               ep(2) = dcmplx( -p(1)*pzpt , -nst*p(2)/pt*sqh )
               ep(3) = dcmplx( -p(2)*pzpt ,  nst*p(1)/pt*sqh )
            else
               ep(2) = dcmplx( -sqh )
               ep(3) = dcmplx( rZero , nst*sign(sqh,p(3)) )
            endif
         endif
      end if

      if ( nhel.le.0 ) then 
c construct eps-
         if ( pp.eq.rZero ) then
            em(1) = dcmplx( rZero )
            em(2) = dcmplx( sqh )
            em(3) = dcmplx( rZero , nst*sqh )
            em(4) = dcmplx( rZero )
         else
            em(1) = dcmplx( rZero )
            em(4) = dcmplx( -pt/pp*sqh )
            if ( pt.ne.rZero ) then
               pzpt = -p(3)/(pp*pt)*sqh
               em(2) = dcmplx( -p(1)*pzpt , -nst*p(2)/pt*sqh )
               em(3) = dcmplx( -p(2)*pzpt ,  nst*p(1)/pt*sqh )
            else
               em(2) = dcmplx( sqh )
               em(3) = dcmplx( rZero , nst*sign(sqh,p(3)) )
            endif
         endif
      end if

      if ( abs(nhel).le.1 ) then  
c construct eps0
         if ( pp.eq.rZero ) then
            e0(1) = dcmplx( rZero )
            e0(2) = dcmplx( rZero )
            e0(3) = dcmplx( rZero )
            e0(4) = dcmplx( rOne )
         else
            emp = p(0)/(tmass*pp)
            e0(1) = dcmplx( pp/tmass )
            e0(4) = dcmplx( p(3)*emp )
            if ( pt.ne.rZero ) then
               e0(2) = dcmplx( p(1)*emp )
               e0(3) = dcmplx( p(2)*emp )
            else
               e0(2) = dcmplx( rZero )
               e0(3) = dcmplx( rZero )
            endif
         end if
      end if

      if ( nhel.eq.2 ) then
         do j = 1,4
            do i = 1,4
               ft(i,j) = ep(i)*ep(j)
            end do
         end do
      else if ( nhel.eq.-2 ) then
         do j = 1,4
            do i = 1,4
               ft(i,j) = em(i)*em(j)
            end do
         end do
      else if (tmass.eq.0) then
         do j = 1,4
            do i = 1,4
               ft(i,j) = 0
            end do
         end do
      else if (tmass.ne.0) then
        if  ( nhel.eq.1 ) then
           do j = 1,4
              do i = 1,4
                 ft(i,j) = sqh*( ep(i)*e0(j) + e0(i)*ep(j) )
              end do
           end do
        else if ( nhel.eq.0 ) then
           do j = 1,4
              do i = 1,4
                 ft(i,j) = sqs*( ep(i)*em(j) + em(i)*ep(j)
     &                                + rTwo*e0(i)*e0(j) )
              end do
           end do
        else if ( nhel.eq.-1 ) then
           do j = 1,4
              do i = 1,4
                 ft(i,j) = sqh*( em(i)*e0(j) + e0(i)*em(j) )
              end do
           end do
        else
           write(stdo,*) 'invalid helicity in TXXXXX'
           stop
        end if
      end if

      tc(1) = ft(1,1)
      tc(2) = ft(1,2)
      tc(3) = ft(1,3)
      tc(4) = ft(1,4)
      tc(5) = ft(2,1)
      tc(6) = ft(2,2)
      tc(7) = ft(2,3)
      tc(8) = ft(2,4)
      tc(9) = ft(3,1)
      tc(10) = ft(3,2)
      tc(11) = ft(3,3)
      tc(12) = ft(3,4)
      tc(13) = ft(4,1)
      tc(14) = ft(4,2)
      tc(15) = ft(4,3)
      tc(16) = ft(4,4)
      tc(17) = ft(5,1)
      tc(18) = ft(6,1)

      return
      end
      subroutine v2xxxx(p,vmass,nhel,nsv , vc)
c
c This subroutine computes a VECTOR wavefunction.
c
c input:
c       real    p(0:3)         : four-momentum of vector boson
c       real    vmass          : mass          of vector boson
c       integer nhel = -1, 0, 1: helicity      of vector boson
c                                (0 is forbidden if vmass=0.0)
c       integer nsv  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex vc(6)          : vector wavefunction       epsilon^mu(v)
c     
      implicit none
      double complex vc(6)
      double precision p(0:3),vmass,hel,hel0,pt,pt2,pp,pzpt,emp,sqh
      integer nhel,nsv,nsvahl

      double precision rZero, rHalf, rOne, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0 )
      
      
      vc(5) = dcmplx(p(0),p(3))*nsv
      vc(6) = dcmplx(p(1),p(2))*nsv

      if( vmass.eq.rZero ) then
            vc(1) = rOne
            vc(2) = p(1)/p(0)
            vc(3) = p(2)/p(0)
            vc(4) = p(3)/p(0)
         else
            vc(1) = p(0)/vmass
            vc(2) = p(1)/vmass
            vc(3) = p(2)/vmass
            vc(4) = p(3)/vmass
         endif
      return
      end
      subroutine boostx(p,q , pboost)
c
c This subroutine performs the Lorentz boost of a four-momentum.  The
c momentum p is assumed to be given in the rest frame of q.  pboost is
c the momentum p boosted to the frame in which q is given.  q must be a
c timelike momentum.
c
c input:
c       real    p(0:3)         : four-momentum p in the q rest  frame
c       real    q(0:3)         : four-momentum q in the boosted frame
c
c output:
c       real    pboost(0:3)    : four-momentum p in the boosted frame
c
      implicit none
      double precision p(0:3),q(0:3),pboost(0:3),pq,qq,m,lf

      double precision rZero
      parameter( rZero = 0.0d0 )






c
      qq = q(1)**2+q(2)**2+q(3)**2


      if ( qq.ne.rZero ) then
         pq = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
         m = sqrt(q(0)**2-qq)
         lf = ((q(0)-m)*pq/qq+p(0))/m
         pboost(0) = (p(0)*q(0)+pq)/m
         pboost(1) =  p(1)+q(1)*lf
         pboost(2) =  p(2)+q(2)*lf
         pboost(3) =  p(3)+q(3)*lf
      else
         pboost(0) = p(0)
         pboost(1) = p(1)
         pboost(2) = p(2)
         pboost(3) = p(3)
      endif
c
      return
      end
      subroutine eaixxx(eb,ea,shlf,chlf,phi,nhe,nha , eai)
c
c This subroutine computes an off-shell electron wavefunction after
c emitting a photon from the electron beam, with a special care for the
c small angle region.  The momenta are measured in the laboratory frame,
c where the e- beam is along the positive z axis.
c
c input:
c       real    eb             : energy (GeV)    of beam  e-
c       real    ea             : energy (GeV)    of final photon
c       real    shlf           : sin(theta/2)    of final photon
c       real    chlf           : cos(theta/2)    of final photon
c       real    phi            : azimuthal angle of final photon
c       integer nhe  = -1 or 1 : helicity        of beam  e-
c       integer nha  = -1 or 1 : helicity        of final photon
c
c output:
c       complex eai(6)         : off-shell electron             |e',A,e>
c
      implicit none
      double complex eai(6),phs
      double precision eb,ea,shlf,chlf,phi,alpha,gal,rnhe,x,c,s,d
      double precision coeff,xnnp,xnnm,snp,csp
      integer nhe,nha,nn

      double precision rHalf, rOne, rTwo, rFour, rOte
      double precision rPi, rIalph
      parameter( rHalf = 0.5d0, rOne = 1.0d0, rTwo = 2.0d0 )
      parameter( rFour = 4.0d0, rOte = 128.9d0 )
      parameter( rPi = 3.14159265358979323846d0 )
      parameter( rIalph = 137.0359895d0 )

      double precision me
      parameter( me = 0.510998902d-3 )

c

      alpha = rOne/rOte
      gal = sqrt(alpha*rFour*rPi)

      nn = nha*nhe
      rnhe = nhe
      x = ea/eb
      c = (chlf+shlf)*(chlf-shlf)
      s = rTwo*chlf*shlf
      d = -rOne/(ea*eb*(rFour*shlf**2+(me/eb)**2*c))
      coeff = -nn*gal*sqrt(eb)*d
      xnnp = x*(1+nn)
      xnnm = x*(1-nn)
      snp = sin(phi)
      csp = cos(phi)
      phs = dcmplx( csp, rnhe*snp )

      eai((5-3*nhe)/2) = -rnhe*coeff*me*s*(rOne+xnnp*rHalf)
      eai((5-nhe)/2)   =  xnnp*coeff*me*chlf**2*phs
      eai((5+nhe)/2)   =  rnhe*coeff*eb*s*(-rTwo+xnnm)
      eai((5+3*nhe)/2) =  xnnm*coeff*eb*shlf**2*phs*rTwo

      eai(5) =  eb*dcmplx( rOne-x, rOne-x*c )
      eai(6) = -eb*x*s*dcmplx( csp, snp )
c
      return
      end
      subroutine eaoxxx(eb,ea,shlf,chlf,phi,nhe,nha , eao)
c
c This subroutine computes an off-shell positron wavefunction after
c emitting a photon from the positron beam, with a special care for the
c small angle region.  The momenta are measured in the laboratory frame,
c where the e+ beam is along the negative z axis.
c
c input:
c       real    eb             : energy (GeV)    of beam  e+
c       real    ea             : energy (GeV)    of final photon
c       real    shlf           : sin(theta/2)    of final photon
c       real    chlf           : cos(theta/2)    of final photon
c       real    phi            : azimuthal angle of final photon
c       integer nhe  = -1 or 1 : helicity        of beam  e+
c       integer nha  = -1 or 1 : helicity        of final photon
c
c output:
c       complex eao(6)         : off-shell positron             <e,A,e'|
c
      implicit none
      double complex eao(6),phs
      double precision eb,ea,shlf,chlf,phi,alpha,gal,rnhe,x,c,s,d
      double precision coeff,xnnp,xnnm,snp,csp
      integer nhe,nha,nn

      double precision rHalf, rOne, rTwo, rFour, rOte
      double precision rPi, rIalph
      parameter( rHalf = 0.5d0, rOne = 1.0d0, rTwo = 2.0d0 )
      parameter( rFour = 4.0d0, rOte = 128.9d0 )
      parameter( rPi = 3.14159265358979323846d0 )
      parameter( rIalph = 137.0359895d0 )

      double precision me
      parameter( me = 0.510998902d-3 )

c

      alpha = rOne/rOte
      gal = sqrt(alpha*rFour*rPi)

      nn = nha*nhe
      rnhe = nhe
      x = ea/eb
      c = (chlf+shlf)*(chlf-shlf)
      s = rTwo*chlf*shlf
      d = -rOne/(ea*eb*(rFour*chlf**2-(me/eb)**2*c))
      coeff = nn*gal*sqrt(eb)*d
      xnnp = x*(1+nn)
      xnnm = x*(1-nn)
      snp = sin(phi)
      csp = cos(phi)
      phs = dcmplx( csp, -rnhe*snp )

      eao((5-3*nhe)/2) =               coeff*me*s*(rOne+xnnp*rHalf)
      eao((5-nhe)/2)   = rnhe*xnnp    *coeff*me*shlf**2*phs
      eao((5+nhe)/2)   =               coeff*eb*s*(-rTwo+xnnm)
      eao((5+3*nhe)/2) = real(nha-nhe)*coeff*eb*x*chlf**2*phs*rTwo

      eao(5) = eb*dcmplx( x-rOne, x*c+rOne )
      eao(6) = eb*x*s*dcmplx( csp, snp )
c
      return
      end
      subroutine fsixxx(fi,sc,gc,fmass,fwidth , fsi)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN external fermion and a vector boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex sc(3)          : input    scalar                      s
c       complex gc(2)          : coupling constants                 gchf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fsi(6)         : off-shell fermion             |f',s,fi>
c
      implicit none
      double complex fi(6),sc(3),fsi(6),gc(2),sl1,sl2,sr1,sr2,ds
      double precision pf(0:3),fmass,fwidth,pf2,p0p3,p0m3

c

      fsi(5) = fi(5)-sc(2)
      fsi(6) = fi(6)-sc(3)

      pf(0) = dble( fsi(5))
      pf(1) = dble( fsi(6))
      pf(2) = dimag(fsi(6))
      pf(3) = dimag(fsi(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)


      ds = -sc(1)/dcmplx( pf2-fmass**2, fmass*fwidth )
      p0p3 = pf(0)+pf(3)
      p0m3 = pf(0)-pf(3)
      sl1 = gc(1)*(p0p3*fi(1)+dconjg(fsi(6))*fi(2))
      sl2 = gc(1)*(p0m3*fi(2)       +fsi(6) *fi(1))
      sr1 = gc(2)*(p0m3*fi(3)-dconjg(fsi(6))*fi(4))
      sr2 = gc(2)*(p0p3*fi(4)       -fsi(6) *fi(3))

      fsi(1) = ( gc(1)*fmass*fi(1) + sr1 )*ds
      fsi(2) = ( gc(1)*fmass*fi(2) + sr2 )*ds
      fsi(3) = ( gc(2)*fmass*fi(3) + sl1 )*ds
      fsi(4) = ( gc(2)*fmass*fi(4) + sl2 )*ds
c
      return
      end
      subroutine fsoxxx(fo,sc,gc,fmass,fwidth , fso)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external fermion and a vector boson.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex sc(6)          : input    scalar                      s
c       complex gc(2)          : coupling constants                 gchf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fso(6)         : off-shell fermion             <fo,s,f'|
c
      implicit none
      double complex fo(6),sc(6),fso(6),gc(2),sl1,sl2,sr1,sr2,ds
      double precision pf(0:3),fmass,fwidth,pf2,p0p3,p0m3

c

      fso(5) = fo(5)+sc(2)
      fso(6) = fo(6)+sc(3)

      pf(0) = dble( fso(5))
      pf(1) = dble( fso(6))
      pf(2) = dimag(fso(6))
      pf(3) = dimag(fso(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)


      ds = -sc(1)/dcmplx( pf2-fmass**2, fmass*fwidth )
      p0p3 = pf(0)+pf(3)
      p0m3 = pf(0)-pf(3)
      sl1 = gc(2)*(p0p3*fo(3)       +fso(6) *fo(4))
      sl2 = gc(2)*(p0m3*fo(4)+dconjg(fso(6))*fo(3))
      sr1 = gc(1)*(p0m3*fo(1)       -fso(6) *fo(2))
      sr2 = gc(1)*(p0p3*fo(2)-dconjg(fso(6))*fo(1))

      fso(1) = ( gc(1)*fmass*fo(1) + sl1 )*ds
      fso(2) = ( gc(1)*fmass*fo(2) + sl2 )*ds
      fso(3) = ( gc(2)*fmass*fo(3) + sr1 )*ds
      fso(4) = ( gc(2)*fmass*fo(4) + sr2 )*ds
c
      return
      end
      subroutine ftixkk_1(fi,tc,g,fmass,fwidth , fti)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN external fermion and a KK tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex tc(6,4)        : input    tensor                      t
c       real    g              : coupling constant                   gtf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fti(6)         : off-shell fermion             |f',t,fi>
c
      implicit none
      double complex fi(6), tc(6,4), fti(6)
      double precision g, fmass, fwidth

      double complex k14p, k14m, k23, k23s, p14p, p14m
      double complex D1, D2, D3, D4, Tii, mTii, d
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double precision pf(4), k(4), pf2, m2
      integer i
      
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
c
      m2 = rTwo*fmass

      fti(5) = fi(5) - tc(5,1)
      fti(6) = fi(6) - tc(6,1)

      pf(1) = dreal(fti(5))
      pf(2) = dreal(fti(6))
      pf(3) = dimag(fti(6))
      pf(4) = dimag(fti(5))
      pf2 = pf(1)**2 - pf(2)**2 - pf(3)**2 - pf(4)**2

      k(1) = dreal(fi(5)) + pf(1)
      k(2) = dreal(fi(6)) + pf(2)
      k(3) = dimag(fi(6)) + pf(3)
      k(4) = dimag(fi(5)) + pf(4)

      k14p = dcmplx( k(1)+k(4), rZero )
      k14m = dcmplx( k(1)-k(4), rZero )
      k23  = dcmplx( k(2), k(3) )
      k23s = dconjg( k23 )
      p14p = dcmplx( pf(1)+pf(4), rZero )
      p14m = dcmplx( pf(1)-pf(4), rZero )

      T11 = rTwo*tc(1,1)
      T22 = rTwo*tc(2,2)
      T33 = rTwo*tc(3,3)
      T44 = rTwo*tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      D1 =   k(1)*(T11-T14) - k(2)*(T12-T24)
     &     - k(3)*(T13-T34) - k(4)*(T14-T44)

      D2 = - k(1)*(T12-ci*T13) + k(2)*(T22-ci*T23)
     &     + k(3)*(T23-ci*T33) + k(4)*(T24-ci*T34)

      D3 = - k(1)*(T12+ci*T13) + k(2)*(T22+ci*T23)
     &     + k(3)*(T23+ci*T33) + k(4)*(T24+ci*T34)

      D4 =   k(1)*(T11+T14) - k(2)*(T12+T24)
     &     - k(3)*(T13+T34) - k(4)*(T14+T44)

      Tii  = T11 - T22 - T33 - T44
      mTii = fmass*Tii

      if ( fmass.gt.rZero ) then
         d = - g/dcmplx( pf2-fmass**2, fmass*fwidth )
      else
         d = - g/dcmplx( pf2, rZero )
      end if

      fti(1) =   fi(1)*(p14m*D4 + dconjg(fti(6))*D3)
     &         - fi(2)*(p14m*D2 + dconjg(fti(6))*D1)
     &         + Tii*(  fi(1)*(dconjg(fti(6))*k23  - p14m*k14p)
     &                + fi(2)*(dconjg(fti(6))*k14m - p14m*k23s) )

      fti(2) = - fi(1)*(fti(6)*D4 + p14p*D3)
     &         + fi(2)*(fti(6)*D2 + p14p*D1)
     &         + Tii*(  fi(1)*(fti(6)*k14p - p14p*k23 )
     &                + fi(2)*(fti(6)*k23s - p14p*k14m) )

      fti(3) =   fi(3)*(p14p*D1 + dconjg(fti(6))*D3)
     &         + fi(4)*(p14p*D2 + dconjg(fti(6))*D4)
     &         + Tii*(  fi(3)*(dconjg(fti(6))*k23  - p14p*k14m)
     &                - fi(4)*(dconjg(fti(6))*k14p - p14p*k23s) )

      fti(4) =   fi(3)*(fti(6)*D1 + p14m*D3)
     &         + fi(4)*(fti(6)*D2 + p14m*D4)
     &         + Tii*(  fi(3)*(p14m*k23  - fti(6)*k14m)
     &                - fi(4)*(p14m*k14p - fti(6)*k23s) )

      if ( fmass.gt.rZero ) then
         fti(1) = fti(1) + fmass*( D1*fi(3) + D2*fi(4) )
         fti(2) = fti(2) + fmass*( D3*fi(3) + D4*fi(4) )
         fti(3) = fti(3) + fmass*( D4*fi(1) - D2*fi(2) )
         fti(4) = fti(4) + fmass*(-D3*fi(1) + D1*fi(2) )
         do i = 1,4
            fti(i) = fti(i) + mTii*m2*fi(i)
         end do
         fti(1) = fti(1) + mTii*(  fi(3)*(rTwo*p14m - k14m)
     &                           + fi(4)*(k23 - rTwo*dconjg(fti(6))) )
         fti(2) = fti(2) + mTii*(  fi(3)*(k23 - rTwo*fti(6))
     &                           + fi(4)*(rTwo*p14p - k14p) )
         fti(3) = fti(3) + mTii*(  fi(1)*(rTwo*p14p - k14p)
     &                           + fi(2)*(rTwo*dconjg(fti(6)) - k23s) )
         fti(4) = fti(4) + mTii*(  fi(1)*(rTwo*fti(6) - k23)
     &                           + fi(2)*(rTwo*p14m - k14m) )
      end if

      do i = 1,4
         fti(i) = fti(i)*d
      end do
c
      return
      end
      subroutine ftixxx(fi,tc,gt,fmass,fwidth , fti)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN external fermion and a tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex tc(18)         : input    tensor                       T
c       complex gt             : coupling constant       gtf=-1/Lambda/4
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fti(6)         : off-shell fermion             |f':T,fi>
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fi(6), tc(18), gt, fti(6)
      double precision fmass, fwidth

      double complex ft(6,4)
      double complex k14p, k14m, k23, k23s, p14p, p14m
      double complex D1, D2, D3, D4, Tii, mTii, d
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double precision pf(4), k(4), pf2, m2
      integer i
      
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )


      m2 = rTwo*fmass
      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      fti(5) = fi(5) - ft(5,1)
      fti(6) = fi(6) - ft(6,1)

      pf(1) = dreal(fti(5))
      pf(2) = dreal(fti(6))
      pf(3) = dimag(fti(6))
      pf(4) = dimag(fti(5))
      pf2 = pf(1)**2 - pf(2)**2 - pf(3)**2 - pf(4)**2

      k(1) = dreal(fi(5)) + pf(1)
      k(2) = dreal(fi(6)) + pf(2)
      k(3) = dimag(fi(6)) + pf(3)
      k(4) = dimag(fi(5)) + pf(4)

      k14p = dcmplx( k(1)+k(4), rZero )
      k14m = dcmplx( k(1)-k(4), rZero )
      k23  = dcmplx( k(2), k(3) )
      k23s = dconjg( k23 )
      p14p = dcmplx( pf(1)+pf(4), rZero )
      p14m = dcmplx( pf(1)-pf(4), rZero )

      T11 = rTwo*ft(1,1)
      T22 = rTwo*ft(2,2)
      T33 = rTwo*ft(3,3)
      T44 = rTwo*ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      D1 =   k(1)*(T11-T14) - k(2)*(T12-T24)
     &     - k(3)*(T13-T34) - k(4)*(T14-T44)

      D2 = - k(1)*(T12-ci*T13) + k(2)*(T22-ci*T23)
     &     + k(3)*(T23-ci*T33) + k(4)*(T24-ci*T34)

      D3 = - k(1)*(T12+ci*T13) + k(2)*(T22+ci*T23)
     &     + k(3)*(T23+ci*T33) + k(4)*(T24+ci*T34)

      D4 =   k(1)*(T11+T14) - k(2)*(T12+T24)
     &     - k(3)*(T13+T34) - k(4)*(T14+T44)

      Tii  = T11 - T22 - T33 - T44
      mTii = fmass*Tii

      if ( fmass.gt.rZero ) then
         d = - gt/dcmplx( pf2-fmass**2, fmass*fwidth )
      else
         d = - gt/dcmplx( pf2, rZero )
      end if

      fti(1) =   fi(1)*(p14m*D4 + dconjg(fti(6))*D3)
     &         - fi(2)*(p14m*D2 + dconjg(fti(6))*D1)
     &         + Tii*(  fi(1)*(dconjg(fti(6))*k23  - p14m*k14p)
     &                + fi(2)*(dconjg(fti(6))*k14m - p14m*k23s) )

      fti(2) = - fi(1)*(fti(6)*D4 + p14p*D3)
     &         + fi(2)*(fti(6)*D2 + p14p*D1)
     &         + Tii*(  fi(1)*(fti(6)*k14p - p14p*k23 )
     &                + fi(2)*(fti(6)*k23s - p14p*k14m) )

      fti(3) =   fi(3)*(p14p*D1 + dconjg(fti(6))*D3)
     &         + fi(4)*(p14p*D2 + dconjg(fti(6))*D4)
     &         + Tii*(  fi(3)*(dconjg(fti(6))*k23  - p14p*k14m)
     &                - fi(4)*(dconjg(fti(6))*k14p - p14p*k23s) )

      fti(4) =   fi(3)*(fti(6)*D1 + p14m*D3)
     &         + fi(4)*(fti(6)*D2 + p14m*D4)
     &         + Tii*(  fi(3)*(p14m*k23  - fti(6)*k14m)
     &                - fi(4)*(p14m*k14p - fti(6)*k23s) )

      if ( fmass.gt.rZero ) then
         fti(1) = fti(1) + fmass*( D1*fi(3) + D2*fi(4) )
         fti(2) = fti(2) + fmass*( D3*fi(3) + D4*fi(4) )
         fti(3) = fti(3) + fmass*( D4*fi(1) - D2*fi(2) )
         fti(4) = fti(4) + fmass*(-D3*fi(1) + D1*fi(2) )
         do i = 1,4
            fti(i) = fti(i) + mTii*m2*fi(i)
         end do
         fti(1) = fti(1) + mTii*(  fi(3)*(rTwo*p14m - k14m)
     &                           + fi(4)*(k23 - rTwo*dconjg(fti(6))) )
         fti(2) = fti(2) + mTii*(  fi(3)*(k23 - rTwo*fti(6))
     &                           + fi(4)*(rTwo*p14p - k14p) )
         fti(3) = fti(3) + mTii*(  fi(1)*(rTwo*p14p - k14p)
     &                           + fi(2)*(rTwo*dconjg(fti(6)) - k23s) )
         fti(4) = fti(4) + mTii*(  fi(1)*(rTwo*fti(6) - k23)
     &                           + fi(2)*(rTwo*p14m - k14m) )
      end if

      do i = 1,4
         fti(i) = fti(i)*d
      end do

      return
      end
      subroutine ftoxkk_1(fo,tc,g,fmass,fwidth , fto)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external fermion and a KK tensor boson.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex tc(6,4)        : input    tensor                      t
c       real    g              : coupling constant                   gtf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvo(6)         : off-shell fermion             <fo,t,f'|
c
      implicit none
      double complex fo(6), tc(6,4), fto(6)
      double precision g, fmass, fwidth

      double complex k14p, k14m, k23, k23s, p14p, p14m
      double complex D1, D2, D3, D4, Tii, mTii, d
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double precision pf(4), k(4), pf2, m2
      integer i

      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
c
      m2 = rTwo*fmass

      fto(5) = fo(5) + tc(5,1)
      fto(6) = fo(6) + tc(6,1)

      pf(1) = dreal(fto(5))
      pf(2) = dreal(fto(6))
      pf(3) = dimag(fto(6))
      pf(4) = dimag(fto(5))
      pf2 = pf(1)**2 - pf(2)**2 - pf(3)**2 - pf(4)**2

      k(1) = dreal(fo(5)) + pf(1)
      k(2) = dreal(fo(6)) + pf(2)
      k(3) = dimag(fo(6)) + pf(3)
      k(4) = dimag(fo(5)) + pf(4)

      k14p = dcmplx( k(1)+k(4), rZero )
      k14m = dcmplx( k(1)-k(4), rZero )
      k23  = dcmplx( k(2), k(3) )
      k23s = dconjg( k23 )
      p14p = dcmplx( pf(1)+pf(4), rZero )
      p14m = dcmplx( pf(1)-pf(4), rZero )

      T11 = rTwo*tc(1,1)
      T22 = rTwo*tc(2,2)
      T33 = rTwo*tc(3,3)
      T44 = rTwo*tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      D1 =   k(1)*(T11-T14) + k(2)*(T24-T12)
     &     + k(3)*(T34-T13) + k(4)*(T44-T14)

      D2 =   k(1)*(ci*T13-T12) + k(2)*(T22-ci*T23)
     &     + k(3)*(T23-ci*T33) + k(4)*(T24-ci*T34)

      D3 = - k(1)*(T12+ci*T13) + k(2)*(T22+ci*T23)
     &     + k(3)*(T23+ci*T33) + k(4)*(T24+ci*T34)

      D4 =   k(1)*(T11+T14) - k(2)*(T12+T24)
     &     - k(3)*(T13+T34) - k(4)*(T14+T44)

      Tii  = T11 - T22 - T33 - T44
      mTii = fmass*Tii

      if ( fmass.gt.rZero ) then
         d = - g/dcmplx( pf2-fmass**2, fmass*fwidth )
      else
         d = - g/dcmplx( pf2, rZero )
      end if

      fto(1) =   fo(1)*(p14p*D1 + fto(6)*D2)
     &         + fo(2)*(p14p*D3 + fto(6)*D4)
     &         + Tii*(- fo(1)*(p14p*k14m - fto(6)*k23s)
     &                + fo(2)*(p14p*k23  - fto(6)*k14p) )

      fto(2) =   fo(1)*(dconjg(fto(6))*D1 + p14m*D2)
     &         + fo(2)*(dconjg(fto(6))*D3 + p14m*D4)
     &         + Tii*(  fo(1)*(p14m*k23s - dconjg(fto(6))*k14m)
     &                - fo(2)*(p14m*k14p - dconjg(fto(6))*k23 ) )

      fto(3) =   fo(3)*(p14m*D4 + fto(6)*D2)
     &         + fo(4)*(p14m*D3 - fto(6)*D1)
     &         + Tii*(  fo(3)*(fto(6)*k23s - p14m*k14p)
     &                + fo(4)*(fto(6)*k14m - p14m*k23 ) )

      fto(4) = - fo(3)*(dconjg(fto(6))*D4 + p14p*D2)
     &         + fo(4)*(dconjg(fto(6))*D3 - p14p*D1)
     &         + Tii*(  fo(3)*(dconjg(fto(6))*k14p - p14p*k23s)
     &                + fo(4)*(dconjg(fto(6))*k23  - p14p*k14m) )

      if ( fmass.gt.rZero ) then
         fto(1) = fto(1) + fmass*( D4*fo(3) - D3*fo(4) )
         fto(2) = fto(2) + fmass*(-D2*fo(3) + D1*fo(4) )
         fto(3) = fto(3) + fmass*( D1*fo(1) + D3*fo(2) )
         fto(4) = fto(4) + fmass*( D2*fo(1) + D4*fo(2) )
         do i = 1,4
            fto(i) = fto(i) + mTii*m2*fo(i)
         end do
         fto(1) = fto(1) + mTii*(  fo(3)*(rTwo*p14p - k14p)
     &                           + fo(4)*(rTwo*fto(6) - k23) )
         fto(2) = fto(2) + mTii*(  fo(3)*(rTwo*dconjg(fto(6)) - k23s)
     &                           + fo(4)*(rTwo*p14m - k14m) )
         fto(3) = fto(3) + mTii*(  fo(1)*(rTwo*p14m - k14m)
     &                           + fo(2)*(k23 - rTwo*fto(6)) )
         fto(4) = fto(4) + mTii*(  fo(1)*(k23s - rTwo*dconjg(fto(6)))
     &                           + fo(2)*(rTwo*p14p - k14p) )
      end if

      do i = 1,4
         fto(i) = fto(i)*d
      end do
c
      return
      end
      subroutine ftoxxx(fo,tc,gt,fmass,fwidth , fto)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external fermion and a tensor boson.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex tc(18)         : input    tensor                       T
c       complex gt             : coupling constant       gtf=-1/Lambda/4
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fto(6)         : off-shell fermion             <fo,T:f'|
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fo(6), tc(18), gt, fto(6)
      double precision fmass, fwidth

      double complex ft(6,4)
      double complex Tii,d
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double precision p(4), k(4), p2
      integer i

      double precision rZero, r2,rtwo
      parameter( rZero = 0.0d0, r2 = 2.0d0,rtwo=2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )


      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      fto(5) = fo(5) + ft(5,1)
      fto(6) = fo(6) + ft(6,1)

      p(1) = dreal(fto(5))
      p(2) = dreal(fto(6))
      p(3) = dimag(fto(6))
      p(4) = dimag(fto(5))
      p2 = p(1)**2 - p(2)**2 - p(3)**2 - p(4)**2

      k(1) = dreal(fo(5)) + p(1)
      k(2) = dreal(fo(6)) + p(2)
      k(3) = dimag(fo(6)) + p(3)
      k(4) = dimag(fo(5)) + p(4)

      T11 = ft(1,1)
      T22 = ft(2,2)
      T33 = ft(3,3)
      T44 = ft(4,4)
      Tii = T11-T22-T33-T44
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)


      if ( fmass.gt.rZero ) then
         d = - gt/dcmplx( p2-fmass**2, fmass*fwidth )
      else
         d = - gt/dcmplx( p2, rZero )
      end if

      fto(1) =fmass*(fmass*fo(1)*r2**2*Tii 
     &	- fo(4)*r2*Tii*(k(2) + ci*k(3))
     & - fo(3)*r2*Tii*(k(1) + k(4))) + 
     &    fmass*(fo(4)*(T12*k(1) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) - k(3)) - T24*k(4) - ci*T34*k(4)) + 
     &       fo(3)*(-(T12*k(2)) - T24*k(2) - T13*k(3)
     & - T34*k(3) - T14*(-k(1) + k(4)))) + 
     &    (fmass*fo(4)*r2**2*Tii - fo(1)*r2*Tii*(-k(2) 
     &+ ci*k(3)) - fo(2)*r2*Tii*(k(1) + k(4)))*(p(2) + ci*p(3)) + 
     &    (fo(1)*(-(T12*k(1)) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) + k(3)) + T24*k(4) - ci*T34*k(4)) + 
     &       fo(2)*(-(T12*k(2)) - T24*k(2) - T13*k(3)
     & - T34*k(3) - T14*(-k(1) + k(4))))*(p(2) + ci*p(3)) + 
     &    (fmass*fo(3)*r2**2*Tii - fo(2)*r2*Tii*(-k(2) 
     &- ci*k(3)) - fo(1)*r2*Tii*(k(1) - k(4)))*(p(1) + p(4)) + 
     &    (fo(2)*(-(T12*k(1)) - ci*T13*k(1) 
     &+ T23*(ci*k(2) + k(3)) + T24*k(4) + ci*T34*k(4)) + 
     &       fo(1)*(-(T12*k(2)) + T24*k(2) 
     &- T13*k(3) + T34*k(3) - T14*(k(1) + k(4))))*(p(1) + p(4)) + 
     &    rtwo*(fmass*(fo(4)*(-(T22*k(2)) 
     &- ci*T33*k(3)) + fo(3)*(T11*k(1) - T44*k(4))) + 
     &       (fo(1)*(T22*k(2) - ci*T33*k(3)) 
     &+ fo(2)*(T11*k(1) - T44*k(4)))*(p(2) + ci*p(3)) + 
     &       (fo(2)*(T22*k(2) + ci*T33*k(3))
     & + fo(1)*(T11*k(1) + T44*k(4)))*(p(1) + p(4)))

      fto(2) = fmass*(fmass*fo(2)*r2**2*Tii 
     &	- fo(3)*r2*Tii*(k(2) - ci*k(3))
     & - fo(4)*r2*Tii*(k(1) - k(4))) + 
     &    fmass*(fo(3)*(T12*k(1) - ci*T13*k(1)
     & + T23*(ci*k(2) - k(3)) - T24*k(4) + ci*T34*k(4)) + 
     &       fo(4)*(-(T12*k(2)) + T24*k(2) 
     &- T13*k(3) + T34*k(3) - T14*(k(1) + k(4)))) + 
     &    (fmass*fo(3)*r2**2*Tii - fo(2)*r2*Tii*(-k(2) 
     &- ci*k(3)) - fo(1)*r2*Tii*(k(1) - k(4)))*(p(2) - ci*p(3)) + 
     &    (fo(2)*(-(T12*k(1)) - ci*T13*k(1) 
     &+ T23*(ci*k(2) + k(3)) + T24*k(4) + ci*T34*k(4)) + 
     &       fo(1)*(-(T12*k(2)) + T24*k(2) - T13*k(3)
     & + T34*k(3) - T14*(k(1) + k(4))))*(p(2) - ci*p(3)) + 
     &    rtwo*(fmass*(fo(3)*(-(T22*k(2)) + ci*T33*k(3)) 
     &+ fo(4)*(T11*k(1) + T44*k(4))) + 
     &       (fo(2)*(T22*k(2) + ci*T33*k(3)) 
     &+ fo(1)*(T11*k(1) + T44*k(4)))*(p(2) - ci*p(3)) + 
     &       (fo(1)*(T22*k(2) - ci*T33*k(3)) 
     &+ fo(2)*(T11*k(1) - T44*k(4)))*(p(1) - p(4))) + 
     &    (fmass*fo(4)*r2**2*Tii - fo(1)*r2*Tii*(-k(2) 
     &+ ci*k(3)) - fo(2)*r2*Tii*(k(1) + k(4)))*(p(1) - p(4)) + 
     &    (fo(1)*(-(T12*k(1)) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) + k(3)) + T24*k(4) - ci*T34*k(4)) + 
     &       fo(2)*(-(T12*k(2)) - T24*k(2) - T13*k(3)
     & - T34*k(3) - T14*(-k(1) + k(4))))*(p(1) - p(4))

      fto(3) = fmass*(fmass*fo(3)*r2**2*Tii 
     &	- fo(2)*r2*Tii*(-k(2) - ci*k(3)) 
     &- fo(1)*r2*Tii*(k(1) - k(4))) + 
     &    fmass*(fo(2)*(-(T12*k(1)) - ci*T13*k(1) 
     &+ T23*(ci*k(2) + k(3)) + T24*k(4) + ci*T34*k(4)) + 
     &       fo(1)*(-(T12*k(2)) + T24*k(2) - T13*k(3)
     & + T34*k(3) - T14*(k(1) + k(4)))) + 
     &    (fmass*fo(2)*r2**2*Tii - fo(3)*r2*Tii*(k(2) 
     &- ci*k(3)) - fo(4)*r2*Tii*(k(1) - k(4)))*(-p(2) - ci*p(3)) + 
     &    (fo(3)*(T12*k(1) - ci*T13*k(1) + T23*(ci*k(2)
     & - k(3)) - T24*k(4) + ci*T34*k(4)) + 
     &       fo(4)*(-(T12*k(2)) + T24*k(2) - T13*k(3) 
     &+ T34*k(3) - T14*(k(1) + k(4))))*(-p(2) - ci*p(3)) + 
     &    rtwo*(fmass*(fo(2)*(T22*k(2) + ci*T33*k(3)) 
     &+ fo(1)*(T11*k(1) + T44*k(4))) + 
     &       (fo(3)*(-(T22*k(2)) + ci*T33*k(3)) 
     &+ fo(4)*(T11*k(1) + T44*k(4)))*(-p(2) - ci*p(3)) + 
     &       (fo(4)*(-(T22*k(2)) - ci*T33*k(3))
     & + fo(3)*(T11*k(1) - T44*k(4)))*(p(1) - p(4))) + 
     &    (fmass*fo(1)*r2**2*Tii - fo(4)*r2*Tii*(k(2) 
     &+ ci*k(3)) - fo(3)*r2*Tii*(k(1) + k(4)))*(p(1) - p(4)) + 
     &    (fo(4)*(T12*k(1) + ci*T13*k(1) +
     & T23*(-(ci*k(2)) - k(3)) - T24*k(4) - ci*T34*k(4)) + 
     &       fo(3)*(-(T12*k(2)) - T24*k(2) - T13*k(3) 
     &- T34*k(3) - T14*(-k(1) + k(4))))*(p(1) - p(4))

      fto(4) =fmass*(fmass*fo(4)*r2**2*Tii
     &	 - fo(1)*r2*Tii*(-k(2) + ci*k(3))
     & - fo(2)*r2*Tii*(k(1) + k(4))) + 
     &    fmass*(fo(1)*(-(T12*k(1)) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) + k(3)) + T24*k(4) - ci*T34*k(4)) + 
     &       fo(2)*(-(T12*k(2)) - T24*k(2) - T13*k(3) 
     &- T34*k(3) - T14*(-k(1) + k(4)))) + 
     &    (fmass*fo(1)*r2**2*Tii - fo(4)*r2*Tii*(k(2) 
     &+ ci*k(3)) - fo(3)*r2*Tii*(k(1) + k(4)))*(-p(2) + ci*p(3)) + 
     &    (fo(4)*(T12*k(1) + ci*T13*k(1) 
     &+ T23*(-(ci*k(2)) - k(3)) - T24*k(4) - ci*T34*k(4)) + 
     &       fo(3)*(-(T12*k(2)) - T24*k(2) - T13*k(3)
     & - T34*k(3) - T14*(-k(1) + k(4))))*(-p(2) + ci*p(3)) + 
     &    (fmass*fo(2)*r2**2*Tii - fo(3)*r2*Tii*(k(2) 
     &- ci*k(3)) - fo(4)*r2*Tii*(k(1) - k(4)))*(p(1) + p(4)) + 
     &    (fo(3)*(T12*k(1) - ci*T13*k(1) + T23*(ci*k(2)
     & - k(3)) - T24*k(4) + ci*T34*k(4)) + 
     &       fo(4)*(-(T12*k(2)) + T24*k(2) - T13*k(3) 
     &+ T34*k(3) - T14*(k(1) + k(4))))*(p(1) + p(4)) + 
     &    rtwo*(fmass*(fo(1)*(T22*k(2) - ci*T33*k(3)) 
     &+ fo(2)*(T11*k(1) - T44*k(4))) + 
     &       (fo(4)*(-(T22*k(2)) - ci*T33*k(3))
     & + fo(3)*(T11*k(1) - T44*k(4)))*(-p(2) + ci*p(3)) + 
     &       (fo(3)*(-(T22*k(2)) + ci*T33*k(3))
     & + fo(4)*(T11*k(1) + T44*k(4)))*(p(1) + p(4)))
     
      do i = 1,4
         fto(i) = fto(i)*d
      end do

      return
      end
      subroutine fvidmx(fi,vc,gc,fmass,fwidth , fvi)
c
c This subroutine computes a dipole moment off-shell fermion
c wavefunction from a flowing-IN external fermion and a vector boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex vc(6)          : input    vector                      v
c       complex gc(2,2)        : coupling constants                  gvf
c                              : first index is L,R as normal
c                              : second index is EDM,-MDM
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvi(6)         : off-shell fermion             |f',v,fi>
c
      implicit none
      double complex fi(6), vc(6), fvi(6), sl1, sl2, sr1, sr2, d
      double complex gc(2,2), gL, gR
      double precision pf(0:3), fmass, fwidth, pf2

      double complex kvc21, kvc31, kvc41, kvc32, kvc42, kvc43
      double precision k(1:4)
      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      gL = -gc(1,1) + cImag*gc(1,2)
      gR =  gc(2,1) + cImag*gc(2,2)

c k in vertex formula defined as (pi - po)
      k(1) = dble( vc(5))
      k(2) = dble( vc(6))
      k(3) = dimag(vc(6))
      k(4) = dimag(vc(5))

      kvc21 = (k(2)*vc(1) - k(1)*vc(2))*cImag
      kvc31 =  k(3)*vc(1) - k(1)*vc(3)
      kvc41 = (k(4)*vc(1) - k(1)*vc(4))*cImag
      kvc32 =  k(3)*vc(2) - k(2)*vc(3)
      kvc42 = (k(4)*vc(2) - k(2)*vc(4))*cImag
      kvc43 =  k(4)*vc(3) - k(3)*vc(4)

      fvi(5) = fi(5) - vc(5)
      fvi(6) = fi(6) - vc(6)

      pf(0) = dble( fvi(5))
      pf(1) = dble( fvi(6))
      pf(2) = dimag(fvi(6))
      pf(3) = dimag(fvi(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)


      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )

      sl1 = gL*(  fi(1)*(kvc41 + kvc32)
     &          + fi(2)*(kvc42 + kvc21 + kvc43 + kvc31) )
      sl2 = gL*(- fi(1)*(kvc42 - kvc21 - kvc43 + kvc31)
     &          - fi(2)*(kvc41 + kvc32)                 )

      if ( gc(2,1).ne.cZero .or.
     &     gc(2,2).ne.cZero      ) then
         sr1 = gR*(- fi(3)*(kvc41 - kvc32)
     &             + fi(4)*(kvc42 - kvc21 + kvc43 - kvc31) )
         sr2 = gR*(- fi(3)*(kvc42 + kvc21 - kvc43 - kvc31)
     &             + fi(4)*(kvc41 - kvc32)                 )

         fvi(1) = ( (pf(0)-pf(3))*sr1 - dconjg(fvi(6))*sr2
     &             + fmass*sl1                             )*d
         fvi(2) = (      - fvi(6)*sr1 +  (pf(0)+pf(3))*sr2
     &             + fmass*sl2                             )*d
         fvi(3) = ( (pf(0)+pf(3))*sl1 + dconjg(fvi(6))*sl2
     &             + fmass*sr1                             )*d
         fvi(4) = (        fvi(6)*sl1 +  (pf(0)-pf(3))*sl2
     &             + fmass*sr2                             )*d

      else
         fvi(1) = fmass*sl1*d
         fvi(2) = fmass*sl2*d
         fvi(3) = ( (pf(0)+pf(3))*sl1 + dconjg(fvi(6))*sl2 )*d
         fvi(4) = (        fvi(6)*sl1 +  (pf(0)-pf(3))*sl2 )*d
      end if
c
      return
      end
      subroutine fvixxx(fi,vc,gc,fmass,fwidth , fvi)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN external fermion and a vector boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvi(6)         : off-shell fermion             |f',v,fi>
c
      implicit none
      double complex fi(6),vc(6),gc(2),fvi(6),sl1,sl2,sr1,sr2,d
      double precision pf(0:3),fmass,fwidth,pf2
      
      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      fvi(5) = fi(5)-vc(5)
      fvi(6) = fi(6)-vc(6)

      pf(0) = dble( fvi(5))
      pf(1) = dble( fvi(6))
      pf(2) = dimag(fvi(6))
      pf(3) = dimag(fvi(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)


      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
      sl1 =   (vc(1)+      vc(4))*fi(1)
     &      + (vc(2)-cImag*vc(3))*fi(2)
      sl2 =   (vc(2)+cImag*vc(3))*fi(1)
     &      + (vc(1)-      vc(4))*fi(2)

      if ( gc(2).ne.cZero ) then
         sr1 =   (vc(1)-      vc(4))*fi(3)
     &         - (vc(2)-cImag*vc(3))*fi(4)
         sr2 = - (vc(2)+cImag*vc(3))*fi(3)
     &         + (vc(1)+      vc(4))*fi(4)

         fvi(1) = ( gc(1)*((pf(0)-pf(3))*sl1 - dconjg(fvi(6))*sl2)
     &             +gc(2)*fmass*sr1 )*d
         fvi(2) = ( gc(1)*(      -fvi(6)*sl1 +  (pf(0)+pf(3))*sl2)
     &             +gc(2)*fmass*sr2 )*d
         fvi(3) = ( gc(2)*((pf(0)+pf(3))*sr1 + dconjg(fvi(6))*sr2)
     &             +gc(1)*fmass*sl1 )*d
         fvi(4) = ( gc(2)*(       fvi(6)*sr1 +  (pf(0)-pf(3))*sr2)
     &             +gc(1)*fmass*sl2 )*d

      else
         d = d * gc(1)
         fvi(1) = ((pf(0)-pf(3))*sl1 - dconjg(fvi(6))*sl2)*d
         fvi(2) = (      -fvi(6)*sl1 +  (pf(0)+pf(3))*sl2)*d
         fvi(3) = fmass*sl1*d
         fvi(4) = fmass*sl2*d
      end if
c
      return
      end
      subroutine fvodmx(fo,vc,gc,fmass,fwidth , fvo)
c
c This subroutine computes a dipole moment off-shell fermion
c wavefunction from a flowing-OUT external fermion and a vector boson.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2,2)        : coupling constants                  gvf
c                              : first index is L,R as normal
c                              : second index is EDM,-MDM
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvo(6)         : off-shell fermion             <fo,v,f'|
c
      implicit none
      double complex fo(6), vc(6), fvo(6), sl1, sl2, sr1, sr2, d
      double complex gc(2,2), gL, gR
      double precision  pf(0:3), fmass, fwidth, pf2

      double complex kvc21, kvc31, kvc41, kvc32, kvc42, kvc43
      double precision k(1:4)
      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      gL = -gc(1,1) + cImag*gc(1,2)
      gR =  gc(2,1) + cImag*gc(2,2)

c k in vertex formula defined as (pi - po)
      k(1) = dble( vc(5))
      k(2) = dble( vc(6))
      k(3) = dimag(vc(6))
      k(4) = dimag(vc(5))

      kvc21 = (k(2)*vc(1) - k(1)*vc(2))*cImag
      kvc31 =  k(3)*vc(1) - k(1)*vc(3)
      kvc41 = (k(4)*vc(1) - k(1)*vc(4))*cImag
      kvc32 =  k(3)*vc(2) - k(2)*vc(3)
      kvc42 = (k(4)*vc(2) - k(2)*vc(4))*cImag
      kvc43 =  k(4)*vc(3) - k(3)*vc(4)

      fvo(5) = fo(5) + vc(5)
      fvo(6) = fo(6) + vc(6)

      pf(0) = dble( fvo(5))
      pf(1) = dble( fvo(6))
      pf(2) = dimag(fvo(6))
      pf(3) = dimag(fvo(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)


      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )

      sl1 = gL*(  fo(1)*(kvc41 + kvc32)
     &          - fo(2)*(kvc42 - kvc21 - kvc43 + kvc31) )
      sl2 = gL*(  fo(1)*(kvc42 + kvc21 + kvc43 + kvc31)
     &          - fo(2)*(kvc41 + kvc32)                 )

      if ( gc(2,1).ne.cZero .or.
     &     gc(2,2).ne.cZero      ) then
         sr1 = gR*(- fo(3)*(kvc41 - kvc32)
     &             - fo(4)*(kvc42 + kvc21 - kvc43 - kvc31) )
         sr2 = gR*(  fo(3)*(kvc42 - kvc21 + kvc43 - kvc31)
     &             + fo(4)*(kvc41 - kvc32)                 )

         fvo(1) = (  (pf(0)+pf(3))*sr1       + fvo(6)*sr2
     &             + fmass*sl1                            )*d
         fvo(2) = ( dconjg(fvo(6))*sr1 +(pf(0)-pf(3))*sr2
     &             + fmass*sl2                            )*d
         fvo(3) = (  (pf(0)-pf(3))*sl1       - fvo(6)*sl2
     &             + fmass*sr1                            )*d
         fvo(4) = (-dconjg(fvo(6))*sl1 +(pf(0)+pf(3))*sl2
     &             + fmass*sr2                            )*d

      else
         fvo(1) = fmass*sl1*d
         fvo(2) = fmass*sl2*d
         fvo(3) = (  (pf(0)-pf(3))*sl1        - fvo(6)*sl2)*d
         fvo(4) = (-dconjg(fvo(6))*sl1 + (pf(0)+pf(3))*sl2)*d
      end if
c
      return
      end
      subroutine fvoxxx(fo,vc,gc,fmass,fwidth , fvo)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external fermion and a vector boson.
c
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvo(6)         : off-shell fermion             <fo,v,f'|
c
      implicit none
      double complex fo(6),vc(6),gc(2),fvo(6),sl1,sl2,sr1,sr2,d
      double precision pf(0:3),fmass,fwidth,pf2

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      fvo(5) = fo(5)+vc(5)
      fvo(6) = fo(6)+vc(6)

      pf(0) = dble( fvo(5))
      pf(1) = dble( fvo(6))
      pf(2) = dimag(fvo(6))
      pf(3) = dimag(fvo(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)


      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
      sl1 =   (vc(1)+      vc(4))*fo(3)
     &      + (vc(2)+cImag*vc(3))*fo(4)
      sl2 =   (vc(2)-cImag*vc(3))*fo(3)
     &      + (vc(1)-      vc(4))*fo(4)

      if ( gc(2).ne.cZero ) then
         sr1 =   (vc(1)-      vc(4))*fo(1)
     &         - (vc(2)+cImag*vc(3))*fo(2)
         sr2 = - (vc(2)-cImag*vc(3))*fo(1)
     &         + (vc(1)+      vc(4))*fo(2)

         fvo(1) = ( gc(2)*( (pf(0)+pf(3))*sr1  +        fvo(6)*sr2)
     &             +gc(1)*fmass*sl1 )*d
         fvo(2) = ( gc(2)*( dconjg(fvo(6))*sr1 + (pf(0)-pf(3))*sr2)
     &             +gc(1)*fmass*sl2 )*d
         fvo(3) = ( gc(1)*( (pf(0)-pf(3))*sl1  -        fvo(6)*sl2)
     &             +gc(2)*fmass*sr1 )*d
         fvo(4) = ( gc(1)*(-dconjg(fvo(6))*sl1 + (pf(0)+pf(3))*sl2)
     &             +gc(2)*fmass*sr2 )*d

      else
         d = d * gc(1)
         fvo(1) = fmass*sl1*d
         fvo(2) = fmass*sl2*d
         fvo(3) = (  (pf(0)-pf(3))*sl1 -        fvo(6)*sl2)*d
         fvo(4) = (-dconjg(fvo(6))*sl1 + (pf(0)+pf(3))*sl2)*d
      end if
c
      return
      end
      subroutine ftixkk_2(fi,tc,g,fmass,fwidth , fti)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN external fermion and a KK tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex tc(6,4)        : input    tensor                      t
c       real    g              : coupling constant                   gtf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fti(6)         : off-shell fermion             |f',t,fi>
c
      implicit none
      double complex fi(6), tc(6,4), fti(6), d
      double complex k14p, k14m, k23, k23s, p14p, p14m
      double complex D1, D2, D3, D4, Tii, mTii
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double precision g, fmass, fwidth, pf(4), k(4), pf2, m2
      integer i
      
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
c
      m2 = rTwo*fmass

      fti(5) = fi(5) - tc(5,1)
      fti(6) = fi(6) - tc(6,1)

      pf(1) = dreal(fti(5))
      pf(2) = dreal(fti(6))
      pf(3) = dimag(fti(6))
      pf(4) = dimag(fti(5))
      pf2 = pf(1)**2 - pf(2)**2 - pf(3)**2 - pf(4)**2

      k(1) = dreal(fi(5)) + pf(1)
      k(2) = dreal(fi(6)) + pf(2)
      k(3) = dimag(fi(6)) + pf(3)
      k(4) = dimag(fi(5)) + pf(4)

      k14p = dcmplx( k(1)+k(4), rZero )
      k14m = dcmplx( k(1)-k(4), rZero )
      k23  = dcmplx( k(2), k(3) )
      k23s = dconjg( k23 )
      p14p = dcmplx( pf(1)+pf(4), rZero )
      p14m = dcmplx( pf(1)-pf(4), rZero )

      T11 = rTwo*tc(1,1)
      T22 = rTwo*tc(2,2)
      T33 = rTwo*tc(3,3)
      T44 = rTwo*tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      D1 =   k(1)*(T11-T14) - k(2)*(T12-T24)
     &     - k(3)*(T13-T34) - k(4)*(T14-T44)

      D2 = - k(1)*(T12-ci*T13) + k(2)*(T22-ci*T23)
     &     + k(3)*(T23-ci*T33) + k(4)*(T24-ci*T34)

      D3 = - k(1)*(T12+ci*T13) + k(2)*(T22+ci*T23)
     &     + k(3)*(T23+ci*T33) + k(4)*(T24+ci*T34)

      D4 =   k(1)*(T11+T14) - k(2)*(T12+T24)
     &     - k(3)*(T13+T34) - k(4)*(T14+T44)

      Tii  = 16.d0*(T11 - T22 - T33 - T44)
      mTii = fmass*Tii

      if ( fmass.gt.rZero ) then
         d = - g/dcmplx( pf2-fmass**2, fmass*fwidth )
      else
         d = - g/dcmplx( pf2, rZero )
      end if

      fti(1) =   fi(1)*(p14m*D4 + dconjg(fti(6))*D3)
     &         - fi(2)*(p14m*D2 + dconjg(fti(6))*D1)
     &         + Tii*(  fi(1)*(dconjg(fti(6))*k23  - p14m*k14p)
     &                + fi(2)*(dconjg(fti(6))*k14m - p14m*k23s) )

      fti(2) = - fi(1)*(fti(6)*D4 + p14p*D3)
     &         + fi(2)*(fti(6)*D2 + p14p*D1)
     &         + Tii*(  fi(1)*(fti(6)*k14p - p14p*k23 )
     &                + fi(2)*(fti(6)*k23s - p14p*k14m) )

      fti(3) =   fi(3)*(p14p*D1 + dconjg(fti(6))*D3)
     &         + fi(4)*(p14p*D2 + dconjg(fti(6))*D4)
     &         + Tii*(  fi(3)*(dconjg(fti(6))*k23  - p14p*k14m)
     &                - fi(4)*(dconjg(fti(6))*k14p - p14p*k23s) )

      fti(4) =   fi(3)*(fti(6)*D1 + p14m*D3)
     &         + fi(4)*(fti(6)*D2 + p14m*D4)
     &         + Tii*(  fi(3)*(p14m*k23  - fti(6)*k14m)
     &                - fi(4)*(p14m*k14p - fti(6)*k23s) )

      if ( fmass.gt.rZero ) then
         fti(1) = fti(1) + fmass*( D1*fi(3) + D2*fi(4) )
         fti(2) = fti(2) + fmass*( D3*fi(3) + D4*fi(4) )
         fti(3) = fti(3) + fmass*( D4*fi(1) - D2*fi(2) )
         fti(4) = fti(4) + fmass*(-D3*fi(1) + D1*fi(2) )
         do i = 1,4
            fti(i) = fti(i) + mTii*m2*fi(i)
         end do
         fti(1) = fti(1) + mTii*(  fi(3)*(rTwo*p14m - k14m)
     &                           + fi(4)*(k23 - rTwo*dconjg(fti(6))) )
         fti(2) = fti(2) + mTii*(  fi(3)*(k23 - rTwo*fti(6))
     &                           + fi(4)*(rTwo*p14p - k14p) )
         fti(3) = fti(3) + mTii*(  fi(1)*(rTwo*p14p - k14p)
     &                           + fi(2)*(rTwo*dconjg(fti(6)) - k23s) )
         fti(4) = fti(4) + mTii*(  fi(1)*(rTwo*fti(6) - k23)
     &                           + fi(2)*(rTwo*p14m - k14m) )
      end if

      do i = 1,4
         fti(i) = fti(i)*d
      end do
c
      return
      end
      subroutine fvtixx(fi,vc,tc,gc,gt,fmass,fwidth , fvti)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN external fermion, a gauge boson and a tensor boson.
c
c input:
c       complex fi(6)          : flow-in fermion                    |fi>
c       complex vc(6)          : input   vector                        v
c       complex tc(18)         : input   tensor                        T
c       complex gc(2)          : coupling constants                  gvf
c       complex gt             : coupling constant      gtfv=-1/Lambda/2
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvti(6)        : off-shell fermion           |f':v,T,fi>
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fi(6), vc(6), tc(18), gc(2), gt, fvti(6)
      double precision  fmass, fwidth

      double complex ft(6,4)
      double complex d, T00, T12, T13, T14, T23, T24, T34
      double precision po(4), po2
      integer i
      
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex cone
      parameter( cone = ( 0.0d0, 1.0d0 ))

      
      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      fvti(5) = fi(5) - ft(5,1) -vc(5)
      fvti(6) = fi(6) - ft(6,1) -vc(6)

      po(1) = dreal(fvti(5))
      po(2) = dreal(fvti(6))
      po(3) = dimag(fvti(6))
      po(4) = dimag(fvti(5))
 
      po2=po(1)**2-po(2)**2-po(3)**2-po(4)**2
     
      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)
   
      if ( fmass.gt.rZero ) then
         d =  -1.0d0/dcmplx( po2-fmass**2, fmass*fwidth )
      else
         d =  -1.0d0/dcmplx( po2, rZero )
      end if

      fvti(1) =  fmass*fi(4)*gc(2)*(T12*vc(1) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) + vc(3)) + 2.d0*T00*(-vc(2) + cone*vc(3))
     & - T24*vc(4) + cone*T34*vc(4)) + fmass*fi(3)*gc(2)*
     &     (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3)
     &+ 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))) + 
     &    fi(1)*gc(1)*((-po(2) + cone*po(3))*(-(T12*vc(1))
     & - cone*T13*vc(1) - T23*(-(cone*vc(2)) - vc(3)) + 
     &          2.d0*T00*(vc(2) + cone*vc(3)) + T24*vc(4) 
     &+ cone*T34*vc(4)) + 
     &       (po(1) - po(4))*(T12*vc(2) + T24*vc(2) 
     &+ T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4)) 
     &+ 2.d0*T00*(vc(1) + vc(4)))) + 
     &    fi(2)*gc(1)*((po(1) - po(4))*(-(T12*vc(1))
     & + cone*T13*vc(1) - T23*(cone*vc(2) - vc(3)) 
     &+ 2.d0*T00*(vc(2) - cone*vc(3)) + 
     &          T24*vc(4) - cone*T34*vc(4)) 
     &+ (-po(2) + cone*po(3))*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) 
     &- T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))

      fvti(2) = fmass*fi(3)*gc(2)*(T12*vc(1) 
     &+ cone*T13*vc(1) - T23*(cone*vc(2) + vc(3))
     & + 2.d0*T00*(-vc(2) - cone*vc(3)) - 
     &       T24*vc(4) - cone*T34*vc(4)) + fmass*fi(4)*gc(2)*
     &     (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3)
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4))) + 
     &    fi(1)*gc(1)*((po(1) + po(4))*(-(T12*vc(1)) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) - vc(3)) + 
     &          2.d0*T00*(vc(2) + cone*vc(3))
     & + T24*vc(4) + cone*T34*vc(4)) + 
     &       (-po(2) - cone*po(3))*(T12*vc(2) + T24*vc(2)
     & + T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4)) + 
     &          2.d0*T00*(vc(1) + vc(4)))) + fi(2)*gc(1)*
     &     ((-po(2) - cone*po(3))*(-(T12*vc(1)) + cone*T13*vc(1)
     & - T23*(cone*vc(2) - vc(3)) + 2.d0*T00*(vc(2) - cone*vc(3)) + 
     &          T24*vc(4) - cone*T34*vc(4)) + (po(1) + po(4))*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3) 
     &+ 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))

      fvti(3) = fmass*fi(2)*gc(1)*(-(T12*vc(1)) + cone*T13*vc(1)
     & - T23*(cone*vc(2) - vc(3)) + 2.d0*T00*(vc(2) - cone*vc(3)) + 
     &       T24*vc(4) - cone*T34*vc(4)) + fmass*fi(1)*gc(1)*
     &     (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3)
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4))) + 
     &    fi(4)*gc(2)*((po(1) + po(4))*(T12*vc(1) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) + vc(3)) + 2.d0*T00*(-vc(2) + cone*vc(3))  
     &       -T24*vc(4) + cone*T34*vc(4)) + (po(2) - cone*po(3))*
     &        (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3)
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4)))) + 
     &    fi(3)*gc(2)*((po(2) - cone*po(3))*(T12*vc(1)
     & + cone*T13*vc(1) - T23*(cone*vc(2) + vc(3)) + 
     &          2.d0*T00*(-vc(2) - cone*vc(3))
     & - T24*vc(4) - cone*T34*vc(4)) + 
     &       (po(1) + po(4))*(T12*vc(2) - T24*vc(2)
     & + T13*vc(3) - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) 
     & + T14*(vc(1) + vc(4))))

      fvti(4) = fmass*fi(1)*gc(1)*(-(T12*vc(1)) - cone*T13*vc(1)
     &	 - T23*(-(cone*vc(2)) - vc(3)) 
     & + 2.d0*T00*(vc(2) + cone*vc(3)) + 
     &       T24*vc(4) + cone*T34*vc(4)) + fmass*fi(2)*gc(1)*
     &     (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3)
     & + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))) + 
     &    fi(4)*gc(2)*((po(2) + cone*po(3))*(T12*vc(1) 
     & - cone*T13*vc(1) - T23*(-(cone*vc(2)) + vc(3)) + 
     &          2.d0*T00*(-vc(2) + cone*vc(3)) - T24*vc(4)
     & + cone*T34*vc(4)) + 
     &       (po(1) - po(4))*(T12*vc(2) + T24*vc(2)
     & + T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4))
     & + 2.d0*T00*(vc(1) + vc(4)))) + 
     &    fi(3)*gc(2)*((po(1) - po(4))*(T12*vc(1) 
     &+ cone*T13*vc(1) - T23*(cone*vc(2) + vc(3))
     & + 2.d0*T00*(-vc(2) - cone*vc(3)) - 
     &          T24*vc(4) - cone*T34*vc(4))
     & + (po(2) + cone*po(3))*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) 
     &- T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))  

      fvti(1) = fvti(1)
     &	-2*fmass*fi(4)*gc(2)*(ft(2,2)*vc(2) - cone*ft(3,3)*vc(3))
     & - 2*fmass*fi(3)*gc(2)*(ft(1,1)*vc(1) + ft(4,4)*vc(4)) + 
     &    fi(1)*gc(1)*(-2*(-po(2) + cone*po(3))*(-(ft(2,2)*vc(2))
     & - cone*ft(3,3)*vc(3)) - 
     &       2*(po(1) - po(4))*(ft(1,1)*vc(1) - ft(4,4)*vc(4))) + 
     &    fi(2)*gc(1)*(-2*(po(1) - po(4))*(-(ft(2,2)*vc(2)) 
     &+ cone*ft(3,3)*vc(3)) - 
     &       2*(-po(2) + cone*po(3))*(ft(1,1)*vc(1) + ft(4,4)*vc(4)))

      fvti(2) = fvti(2)
     &-2*fmass*fi(3)*gc(2)*(ft(2,2)*vc(2) + cone*ft(3,3)*vc(3))
     & - 2*fmass*fi(4)*gc(2)*(ft(1,1)*vc(1) - ft(4,4)*vc(4)) + 
     &   fi(1)*gc(1)*(-2*(po(1) + po(4))*(-(ft(2,2)*vc(2))
     & - cone*ft(3,3)*vc(3)) - 
     &       2*(-po(2) - cone*po(3))*(ft(1,1)*vc(1) - ft(4,4)*vc(4))) + 
     &    fi(2)*gc(1)*(-2*(-po(2) - cone*po(3))*(-(ft(2,2)*vc(2))
     & + cone*ft(3,3)*vc(3)) - 
     &       2*(po(1) + po(4))*(ft(1,1)*vc(1) + ft(4,4)*vc(4)))
     	
      fvti(3) = fvti(3)
     &-2*fmass*fi(2)*gc(1)*(-(ft(2,2)*vc(2)) + cone*ft(3,3)*vc(3))
     & - 2*fmass*fi(1)*gc(1)*(ft(1,1)*vc(1) - ft(4,4)*vc(4)) + 
     &    fi(4)*gc(2)*(-2*(po(1) + po(4))*(ft(2,2)*vc(2) 
     &- cone*ft(3,3)*vc(3)) - 
     &       2*(po(2) - cone*po(3))*(ft(1,1)*vc(1) - ft(4,4)*vc(4))) + 
     &    fi(3)*gc(2)*(-2*(po(2) - cone*po(3))*(ft(2,2)*vc(2) 
     &+ cone*ft(3,3)*vc(3)) - 
     &       2*(po(1) + po(4))*(ft(1,1)*vc(1) + ft(4,4)*vc(4)))	
	
      fvti(4) = fvti(4)
     &-2*fmass*fi(1)*gc(1)*(-(ft(2,2)*vc(2)) - cone*ft(3,3)*vc(3))
     & - 2*fmass*fi(2)*gc(1)*(ft(1,1)*vc(1) + ft(4,4)*vc(4)) + 
     &    fi(4)*gc(2)*(-2*(po(2) + cone*po(3))*(ft(2,2)*vc(2)
     & - cone*ft(3,3)*vc(3)) - 
     &       2*(po(1) - po(4))*(ft(1,1)*vc(1) - ft(4,4)*vc(4))) + 
     &    fi(3)*gc(2)*(-2*(po(1) - po(4))*(ft(2,2)*vc(2)
     & + cone*ft(3,3)*vc(3)) - 
     &       2*(po(2) + cone*po(3))*(ft(1,1)*vc(1) + ft(4,4)*vc(4)))	

      do i = 1,4
         fvti(i) = -fvti(i)*d*gt
      end do

      return
      end
      subroutine ftoxkk_2(fo,tc,g,fmass,fwidth , fto)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external fermion and a KK tensor boson.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex tc(6,4)        : input    tensor                      t
c       real    g              : coupling constant                   gtf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvo(6)         : off-shell fermion             <fo,t,f'|
c
      implicit none
      double complex fo(6), tc(6,4), fto(6), d
      double complex k14p, k14m, k23, k23s, p14p, p14m
      double complex D1, D2, D3, D4, Tii, mTii
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double precision g, fmass, fwidth, pf(4), k(4), pf2, m2
      integer i

      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
c
      m2 = rTwo*fmass

      fto(5) = fo(5) + tc(5,1)
      fto(6) = fo(6) + tc(6,1)

      pf(1) = dreal(fto(5))
      pf(2) = dreal(fto(6))
      pf(3) = dimag(fto(6))
      pf(4) = dimag(fto(5))
      pf2 = pf(1)**2 - pf(2)**2 - pf(3)**2 - pf(4)**2

      k(1) = dreal(fo(5)) + pf(1)
      k(2) = dreal(fo(6)) + pf(2)
      k(3) = dimag(fo(6)) + pf(3)
      k(4) = dimag(fo(5)) + pf(4)

      k14p = dcmplx( k(1)+k(4), rZero )
      k14m = dcmplx( k(1)-k(4), rZero )
      k23  = dcmplx( k(2), k(3) )
      k23s = dconjg( k23 )
      p14p = dcmplx( pf(1)+pf(4), rZero )
      p14m = dcmplx( pf(1)-pf(4), rZero )

      T11 = rTwo*tc(1,1)
      T22 = rTwo*tc(2,2)
      T33 = rTwo*tc(3,3)
      T44 = rTwo*tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      D1 =   k(1)*(T11-T14) + k(2)*(T24-T12)
     &     + k(3)*(T34-T13) + k(4)*(T44-T14)

      D2 =   k(1)*(ci*T13-T12) + k(2)*(T22-ci*T23)
     &     + k(3)*(T23-ci*T33) + k(4)*(T24-ci*T34)

      D3 = - k(1)*(T12+ci*T13) + k(2)*(T22+ci*T23)
     &     + k(3)*(T23+ci*T33) + k(4)*(T24+ci*T34)

      D4 =   k(1)*(T11+T14) - k(2)*(T12+T24)
     &     - k(3)*(T13+T34) - k(4)*(T14+T44)

      Tii  = 16.d0*(T11 - T22 - T33 - T44)
      mTii = fmass*Tii

      if ( fmass.gt.rZero ) then
         d = - g/dcmplx( pf2-fmass**2, fmass*fwidth )
      else
         d = - g/dcmplx( pf2, rZero )
      end if

      fto(1) =   fo(1)*(p14p*D1 + fto(6)*D2)
     &         + fo(2)*(p14p*D3 + fto(6)*D4)
     &         + Tii*(- fo(1)*(p14p*k14m - fto(6)*k23s)
     &                + fo(2)*(p14p*k23  - fto(6)*k14p) )

      fto(2) =   fo(1)*(dconjg(fto(6))*D1 + p14m*D2)
     &         + fo(2)*(dconjg(fto(6))*D3 + p14m*D4)
     &         + Tii*(  fo(1)*(p14m*k23s - dconjg(fto(6))*k14m)
     &                - fo(2)*(p14m*k14p - dconjg(fto(6))*k23 ) )

      fto(3) =   fo(3)*(p14m*D4 + fto(6)*D2)
     &         + fo(4)*(p14m*D3 - fto(6)*D1)
     &         + Tii*(  fo(3)*(fto(6)*k23s - p14m*k14p)
     &                + fo(4)*(fto(6)*k14m - p14m*k23 ) )

      fto(4) = - fo(3)*(dconjg(fto(6))*D4 + p14p*D2)
     &         + fo(4)*(dconjg(fto(6))*D3 - p14p*D1)
     &         + Tii*(  fo(3)*(dconjg(fto(6))*k14p - p14p*k23s)
     &                + fo(4)*(dconjg(fto(6))*k23  - p14p*k14m) )

      if ( fmass.gt.rZero ) then
         fto(1) = fto(1) + fmass*( D4*fo(3) - D3*fo(4) )
         fto(2) = fto(2) + fmass*(-D2*fo(3) + D1*fo(4) )
         fto(3) = fto(3) + fmass*( D1*fo(1) + D3*fo(2) )
         fto(4) = fto(4) + fmass*( D2*fo(1) + D4*fo(2) )
         do i = 1,4
            fto(i) = fto(i) + mTii*m2*fo(i)
         end do
         fto(1) = fto(1) + mTii*(  fo(3)*(rTwo*p14p - k14p)
     &                           + fo(4)*(rTwo*fto(6) - k23) )
         fto(2) = fto(2) + mTii*(  fo(3)*(rTwo*dconjg(fto(6)) - k23s)
     &                           + fo(4)*(rTwo*p14m - k14m) )
         fto(3) = fto(3) + mTii*(  fo(1)*(rTwo*p14m - k14m)
     &                           + fo(2)*(k23 - rTwo*fto(6)) )
         fto(4) = fto(4) + mTii*(  fo(1)*(k23s - rTwo*dconjg(fto(6)))
     &                           + fo(2)*(rTwo*p14p - k14p) )
      end if

      do i = 1,4
         fto(i) = fto(i)*d
      end do
c
      return
      end
      subroutine fvtoxx(fo,vc,tc,gc,gt,fmass,fwidth , fvto)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-out external fermion, a gauge boson and a tensor boson.
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                       v
c       complex tc(18)         : input    tensor                       T
c       complex gc(2)          : coupling constants                  gvf
c       complex gt             : coupling constant      gtfv=-1/Lambda/2
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvto(6)        : off-shell fermion           <fo,v,T:f'|
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fo(6), vc(6), tc(18), gc(2),gt, fvto(6)
      double precision  fmass, fwidth

      double complex ft(6,4)
      double complex d, T00, T12, T13, T14, T23, T24, T34
      double precision pi(4), pi2
      integer i
      
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex cone
      parameter( cone = ( 0.0d0, 1.0d0 ))


      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      fvto(5) = fo(5) + ft(5,1) +vc(5)
      fvto(6) = fo(6) + ft(6,1) +vc(6)

      pi(1) = dreal(fvto(5))
      pi(2) = dreal(fvto(6))
      pi(3) = dimag(fvto(6))
      pi(4) = dimag(fvto(5))
 
      pi2 = pi(1)**2-pi(2)**2-pi(3)**2-pi(4)**2
     

      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      if ( fmass.gt.rZero ) then
         d =  -1.0d0/dcmplx( pi2-fmass**2, fmass*fwidth )
      else
         d =  -1.0d0/dcmplx( pi2, rZero )
      end if

      fvto(1) = gc(2)*(pi(2) + cone*pi(3))*(fo(1)*(T12*vc(1)
     &	 - cone*T13*vc(1) - T23*(-(cone*vc(2)) + vc(3)) + 
     &          2.d0*T00*(-vc(2) + cone*vc(3)) - T24*vc(4)
     & + cone*T34*vc(4)) + 
     &       fo(2)*(T12*vc(2) + T24*vc(2) + T13*vc(3) 
     &+ T34*vc(3) + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4))))  
     &   + fmass*gc(1)*(fo(4)*(-(T12*vc(1)) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) - vc(3)) + 2.d0*T00*(vc(2) + cone*vc(3)) + 
     &          T24*vc(4) + cone*T34*vc(4)) + fo(3)*
     &        (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3)
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4)))) + 
     &    gc(2)*(pi(1) + pi(4))*(fo(2)*(T12*vc(1) + cone*T13*vc(1)
     &- T23*(cone*vc(2) + vc(3)) + 2.d0*T00*(-vc(2) - cone*vc(3)) - 
     &          T24*vc(4) - cone*T34*vc(4)) + fo(1)*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3) 
     &+ 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))


      fvto(2) = gc(2)*(pi(1) - pi(4))*(fo(1)*(T12*vc(1) - cone*T13*vc(1)
     & - T23*(-(cone*vc(2)) + vc(3)) + 2.d0*T00*(-vc(2) + cone*vc(3)) 
     &        -   T24*vc(4) + cone*T34*vc(4)) + fo(2)*
     &        (T12*vc(2) + T24*vc(2) + T13*vc(3) + T34*vc(3) 
     & + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4)))) + 
     &    gc(2)*(pi(2) - cone*pi(3))*(fo(2)*(T12*vc(1) 
     & + cone*T13*vc(1) - T23*(cone*vc(2) + vc(3)) + 
     &          2.d0*T00*(-vc(2) - cone*vc(3)) - T24*vc(4)
     & - cone*T34*vc(4)) + 
     &       fo(1)*(T12*vc(2) - T24*vc(2) + T13*vc(3) 
     & - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))  
     &    +fmass*gc(1)*(fo(3)*(-(T12*vc(1)) + cone*T13*vc(1)
     & - T23*(cone*vc(2) - vc(3)) + 2.d0*T00*(vc(2) - cone*vc(3)) + 
     &          T24*vc(4) - cone*T34*vc(4)) + fo(4)*
     &        (T12*vc(2) - T24*vc(2) + T13*vc(3) - T34*vc(3)
     & + 2.d0*T00*(vc(1) - vc(4)) + T14*(vc(1) + vc(4))))



      fvto(3) = gc(1)*(pi(1) - pi(4))*(fo(4)*(-(T12*vc(1))
     &	 - cone*T13*vc(1) - T23*(-(cone*vc(2)) - vc(3)) + 
     &    2.d0*T00*(vc(2) + cone*vc(3)) + T24*vc(4) + cone*T34*vc(4))  
     &       +fo(3)*(T12*vc(2) + T24*vc(2) + T13*vc(3) 
     & + T34*vc(3) + T14*(-vc(1) + vc(4)) + 2.d0*T00*(vc(1) + vc(4))))  
     &   + fmass*gc(2)*(fo(2)*(T12*vc(1) + cone*T13*vc(1) 
     &- T23*(cone*vc(2) + vc(3)) + 2.d0*T00*(-vc(2) - cone*vc(3))
     & - T24*vc(4) - 
     &          cone*T34*vc(4)) + fo(1)*(T12*vc(2) 
     & - T24*vc(2) + T13*vc(3) - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4))  
     &         + T14*(vc(1) + vc(4)))) + gc(1)*(-pi(2) - cone*pi(3))*
     &     (fo(3)*(-(T12*vc(1)) + cone*T13*vc(1) 
     &- T23*(cone*vc(2) - vc(3))
     & + 2.d0*T00*(vc(2) - cone*vc(3)) + T24*vc(4) - 
     &          cone*T34*vc(4)) + fo(4)*(T12*vc(2) - T24*vc(2)
     & + T13*vc(3) - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + 
     &          T14*(vc(1) + vc(4))))

      fvto(4) = fmass*gc(2)*(fo(1)*
     &        (T12*vc(1) - cone*T13*vc(1) - T23*(-(cone*vc(2)) + vc(3))
     & + 2.d0*T00*(-vc(2) + cone*vc(3)) - T24*vc(4) + 
     &          cone*T34*vc(4)) + fo(2)*(T12*vc(2) + T24*vc(2) 
     & + T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4)) + 
     &   2.d0*T00*(vc(1) + vc(4)))) + gc(1)*(-pi(2) + cone*pi(3))*
     &     (fo(4)*(-(T12*vc(1)) - cone*T13*vc(1) - T23*(-(cone*vc(2))
     & - vc(3)) + 2.d0*T00*(vc(2) + cone*vc(3)) + T24*vc(4) + 
     &          cone*T34*vc(4)) + fo(3)*(T12*vc(2) + T24*vc(2)
     & + T13*vc(3) + T34*vc(3) + T14*(-vc(1) + vc(4)) + 
     &          2.d0*T00*(vc(1) + vc(4)))) + gc(1)*(pi(1) + pi(4))*
     &     (fo(3)*(-(T12*vc(1)) + cone*T13*vc(1) - T23*(cone*vc(2) 
     & - vc(3)) + 2.d0*T00*(vc(2) - cone*vc(3)) + T24*vc(4) - 
     &          cone*T34*vc(4)) + fo(4)*(T12*vc(2) - T24*vc(2)
     & + T13*vc(3) - T34*vc(3) + 2.d0*T00*(vc(1) - vc(4)) + 
     &          T14*(vc(1) + vc(4))))

      fvto(1) = fvto(1)+
     &gc(2)*(pi(2) + cone*pi(3))*(-2*fo(1)*(ft(2,2)*vc(2)
     & - cone*ft(3,3)*vc(3)) - 
     &       2*fo(2)*(ft(1,1)*vc(1) - ft(4,4)*vc(4))) + 
     &    fmass*gc(1)*(-2*fo(4)*(-(ft(2,2)*vc(2)) - cone*ft(3,3)*vc(3))
     & - 2*fo(3)*(ft(1,1)*vc(1) - ft(4,4)*vc(4))) + 
     &    gc(2)*(pi(1) + pi(4))*(-2*fo(2)*(ft(2,2)*vc(2) 
     &+ cone*ft(3,3)*vc(3)) - 2*fo(1)*(ft(1,1)*vc(1) + ft(4,4)*vc(4)))

      fvto(2) = fvto(2)+
     &gc(2)*(pi(1) - pi(4))*(-2*fo(1)*(ft(2,2)*vc(2)-cone*ft(3,3)*vc(3))
     & - 2*fo(2)*(ft(1,1)*vc(1) - ft(4,4)*vc(4))) + 
     &    gc(2)*(pi(2)-cone*pi(3))*(-2*fo(2)*(ft(2,2)*vc(2)
     & + cone*ft(3,3)*vc(3))-2*fo(1)*(ft(1,1)*vc(1) + ft(4,4)*vc(4))) + 
     &    fmass*gc(1)*(-2*fo(3)*(-(ft(2,2)*vc(2)) + cone*ft(3,3)*vc(3))
     & - 2*fo(4)*(ft(1,1)*vc(1) + ft(4,4)*vc(4)))
     	
      fvto(3) = fvto(3)+
     &gc(1)*(pi(1) - pi(4))*(-2*fo(4)*(-(ft(2,2)*vc(2))
     & - cone*ft(3,3)*vc(3)) - 2*fo(3)*(ft(1,1)*vc(1)-ft(4,4)*vc(4)))+ 
     &    fmass*gc(2)*(-2*fo(2)*(ft(2,2)*vc(2) + cone*ft(3,3)*vc(3))
     & - 2*fo(1)*(ft(1,1)*vc(1) + ft(4,4)*vc(4))) + 
     &    gc(1)*(-pi(2) - cone*pi(3))*(-2*fo(3)*(-(ft(2,2)*vc(2))
     & + cone*ft(3,3)*vc(3)) - 2*fo(4)*(ft(1,1)*vc(1) + ft(4,4)*vc(4)))
	
      fvto(4) = fvto(4)+
     &fmass*gc(2)*(-2*fo(1)*(ft(2,2)*vc(2) - cone*ft(3,3)*vc(3))
     & - 2*fo(2)*(ft(1,1)*vc(1) - ft(4,4)*vc(4))) + 
     &    gc(1)*(-pi(2) + cone*pi(3))*(-2*fo(4)*(-(ft(2,2)*vc(2))
     & - cone*ft(3,3)*vc(3)) - 2*fo(3)*(ft(1,1)*vc(1)-ft(4,4)*vc(4)))+ 
     &    gc(1)*(pi(1) + pi(4))*(-2*fo(3)*(-(ft(2,2)*vc(2)) 
     &+ cone*ft(3,3)*vc(3)) - 2*fo(4)*(ft(1,1)*vc(1) + ft(4,4)*vc(4)))

      do i = 1,4
         fvto(i) = -fvto(i)*d*gt
      end do

      return
      end
      subroutine ggggtx(va,vb,vc,vd,tc,gc,gt , vertex)
c      
c This subroutine computes the portion of the amplitude of the five-point 
c coupling of a tensor boson with 4 massless color octet gauge bosons
c (gluons) corresponding to the color structure f^{a,b,e} f{c,d,e}. 
c
c To obtain the complete amplitude, this coupling must be called three
c times (once for each color structure) with the following permutations:
c     call ggggtx(va,vb,vc,vd,tc,gc,gt , vertex1)
c     call ggggtx(va,vc,vd,vb,tc,gc,gt , vertex2)
c     call ggggtx(va,vd,vb,vc,tc,gc,gt , vertex3)
c corresponding to
c	f^{a,b,e} f^{c,d,e}
c	f^{a,c,e} f^{d,b,e}
c	f^{a,d,e} f^{b,c,e}
c                                                                       
c input:                                                                
c       complex va(6)          : boson with adjoint color index a     va
c       complex vb(6)          : boson with adjoint color index b     vb
c       complex vc(6)          : boson with adjoint color index c     vc
c       complex vd(6)          : boson with adjoint color index d     vd
c       complex tc(18)         : input tensor                          T
c       real    gc             : coupling constant                    gs
c       complex gt             : coupling constant         gtv=-1/Lambda
c
c output:
c       complex vertex         : amplitude          gamma(va,vb,vc,vd,T)
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex va(6), vb(6), vc(6), vd(6), tc(18), gt, vertex
      double precision gc

      double complex vab,vac,vad,vbc,vbd,vcd,ft(6,4),dvertx
      double complex T00, T12, T13, T14, T23, T24, T34	
      double complex TV24,TV23,TV14,TV13
      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      TV14 = rtwo*(ft(1,1)*va(1)*vd(1)+ft(2,2)*va(2)*vd(2)
     &+ft(3,3)*va(3)*vd(3)+ft(4,4)*va(4)*vd(4))

      TV13 = rtwo*(ft(1,1)*va(1)*vc(1)+ft(2,2)*va(2)*vc(2)
     &+ft(3,3)*va(3)*vc(3)+ft(4,4)*va(4)*vc(4))

      TV24 = rtwo*(ft(1,1)*vb(1)*vd(1)+ft(2,2)*vb(2)*vd(2)
     &+ft(3,3)*vb(3)*vd(3)+ft(4,4)*vb(4)*vd(4))

      TV23 = rtwo*(ft(1,1)*vb(1)*vc(1)+ft(2,2)*vb(2)*vc(2)
     &+ft(3,3)*vb(3)*vc(3)+ft(4,4)*vb(4)*vc(4))

	
      TV14 = TV14- T12*(va(1)*vd(2) + va(2)*vd(1))
     &          - T13*(va(1)*vd(3) + va(3)*vd(1))
     &          - T14*(va(1)*vd(4) + va(4)*vd(1))
     &          + T23*(va(2)*vd(3) + va(3)*vd(2))
     &          + T24*(va(2)*vd(4) + va(4)*vd(2))
     &          + T34*(va(3)*vd(4) + va(4)*vd(3))
      
      TV13 = TV13 - T12*(va(1)*vc(2) + va(2)*vc(1))
     &          - T13*(va(1)*vc(3) + va(3)*vc(1))
     &          - T14*(va(1)*vc(4) + va(4)*vc(1))
     &          + T23*(va(2)*vc(3) + va(3)*vc(2))
     &          + T24*(va(2)*vc(4) + va(4)*vc(2))
     &          + T34*(va(3)*vc(4) + va(4)*vc(3))

      TV24 = TV24 - T12*(vb(1)*vd(2) + vb(2)*vd(1))
     &          - T13*(vb(1)*vd(3) + vb(3)*vd(1))
     &          - T14*(vb(1)*vd(4) + vb(4)*vd(1))
     &          + T23*(vb(2)*vd(3) + vb(3)*vd(2))
     &          + T24*(vb(2)*vd(4) + vb(4)*vd(2))
     &          + T34*(vb(3)*vd(4) + vb(4)*vd(3))

      TV23 = TV23 - T12*(vb(1)*vc(2) + vb(2)*vc(1))
     &          - T13*(vb(1)*vc(3) + vb(3)*vc(1))
     &          - T14*(vb(1)*vc(4) + vb(4)*vc(1))
     &          + T23*(vb(2)*vc(3) + vb(3)*vc(2))
     &          + T24*(vb(2)*vc(4) + vb(4)*vc(2))
     &          + T34*(vb(3)*vc(4) + vb(4)*vc(3))
     	

      vab = va(1)*vb(1)-va(2)*vb(2)-va(3)*vb(3)-va(4)*vb(4)
      vac = va(1)*vc(1)-va(2)*vc(2)-va(3)*vc(3)-va(4)*vc(4)
      vad = va(1)*vd(1)-va(2)*vd(2)-va(3)*vd(3)-va(4)*vd(4)
      vbc = vb(1)*vc(1)-vb(2)*vc(2)-vb(3)*vc(3)-vb(4)*vc(4)
      vbd = vb(1)*vd(1)-vb(2)*vd(2)-vb(3)*vd(3)-vb(4)*vd(4)
      vcd = vc(1)*vd(1)-vc(2)*vd(2)-vc(3)*vd(3)-vc(4)*vd(4)

      dvertx = -TV13*vbd-TV24*vac+TV23*vad+TV14*vbc 
     &+vbd*vac*T00-vad*vbc*T00

      vertex = -dvertx * gc*gc*gt

      return
      end
      subroutine ggggxx(ga,gb,gc,gd,g, vertex)
c
c This subroutine computes the portion of the amplitude of the four-point 
c coupling of 4 massless color octet gauge bosons (gluons) corresponding 
c to the color structure f^{a,b,e} f{c,d,e}. 
c To optain the complete amplitude, this coupling must be called three
c times (once for each color structure) with the following permutations:
c	call ggggxx(ga,gb,gc,gd,g,v1)
c       call ggggxx(ga,gc,gd,gb,g,v2)
c       call ggggxx(ga,gd,gb,gc,g,v3)
c
c	f^{a,b,e} f{c,d,e}*v1+
c	f^{a,c,e} f{d,b,e}*v2+
c	f^{a,d,e} f{b,c,e}*v3
c (See 2.9.1 of the manual for more information).
c                                                                       
c input:                                                                
c       complex ga(0:3)        : Boson with adjoint color index a 
c       complex gb(0:3)        : Boson with adjoint color index b
c       complex gc(0:3)        : Boson with adjoint color index c 
c       complex gd(0:3)        : Boson with adjoint color index d
c       real    g              : coupling of w31 with w-/w+             
c
      implicit none
      double complex ga(6),gb(6),gc(6),gd(6),vertex
      double complex dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),
     &     dvertx,v12,v13,v14,v23,v24,v34
      double precision pga(0:3),pgb(0:3),pgc(0:3),pgd(0:3),g

      save dv1,dv2,dv3, dv4
c      save dv1,dv2,dv3,dv4,dvertx,v12,v13,v14,v23,v24,v34

c

      dv1(0) = dcmplx(ga(1))
      dv1(1) = dcmplx(ga(2))
      dv1(2) = dcmplx(ga(3))
      dv1(3) = dcmplx(ga(4))
      dv2(0) = dcmplx(gb(1))
      dv2(1) = dcmplx(gb(2))
      dv2(2) = dcmplx(gb(3))
      dv2(3) = dcmplx(gb(4))
      dv3(0) = dcmplx(gc(1))
      dv3(1) = dcmplx(gc(2))
      dv3(2) = dcmplx(gc(3))
      dv3(3) = dcmplx(gc(4))
      dv4(0) = dcmplx(gd(1))
      dv4(1) = dcmplx(gd(2))
      dv4(2) = dcmplx(gd(3))
      dv4(3) = dcmplx(gd(4))

      v12 = dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
      v13 = dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
      v14 = dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
      v23 = dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
      v24 = dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
      v34 = dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)

      dvertx = v14*v23 -v13*v24

      vertex = dcmplx( dvertx ) * (g*g)

c      if (abs(dvertx) .gt. 1d40) then
c         write(*,*) 'Careful',abs(dvertx)
c         write(*,*) v12,v13,v14
c         write(*,*) v23,v24,v34
c      endif
c
      return
      end
      subroutine gggxxx(wm,wp,w3,g , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c the gauge bosons.
c
c input:
c       complex wm(6)          : vector               flow-out W-
c       complex wp(6)          : vector               flow-out W+
c       complex w3(6)          : vector               j3 or A    or Z
c       real    g              : coupling constant    gw or gwwa or gwwz
c
c output:
c       complex vertex         : amplitude               gamma(wm,wp,w3)
c
      implicit none
      double complex wm(6),wp(6),w3(6),vertex
      double complex xv1,xv2,xv3,v12,v23,v31
      double complex p12,p13,p21,p23,p31,p32
      double precision pwm(0:3),pwp(0:3),pw3(0:3),g

      double precision rZero, rTenth
      parameter( rZero = 0.0d0, rTenth = 0.1d0 )

c
      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pw3(0) = dble( w3(5))
      pw3(1) = dble( w3(6))
      pw3(2) = dimag(w3(6))
      pw3(3) = dimag(w3(5))


      v12 = wm(1)*wp(1)-wm(2)*wp(2)-wm(3)*wp(3)-wm(4)*wp(4)
      v23 = wp(1)*w3(1)-wp(2)*w3(2)-wp(3)*w3(3)-wp(4)*w3(4)
      v31 = w3(1)*wm(1)-w3(2)*wm(2)-w3(3)*wm(3)-w3(4)*wm(4)
      xv1 = rZero
      xv2 = rZero
      xv3 = rZero
      if ( abs(wm(1)).ne.rZero ) then
         if ( abs(wm(1)).ge.max(abs(wm(2)),abs(wm(3)),abs(wm(4)))
     &        *rTenth )
     &      xv1 = pwm(0)/wm(1)
      endif
      if ( abs(wp(1)).ne.rZero ) then
         if ( abs(wp(1)).ge.max(abs(wp(2)),abs(wp(3)),abs(wp(4)))
     &        *rTenth )
     &      xv2 = pwp(0)/wp(1)
      endif
      if ( abs(w3(1)).ne.rZero ) then
         if ( abs(w3(1)).ge.max(abs(w3(2)),abs(w3(3)),abs(w3(4)))
     &        *rTenth )
     &      xv3 = pw3(0)/w3(1)
      endif

      p12 = (pwm(0)-xv1*wm(1))*wp(1)-(pwm(1)-xv1*wm(2))*wp(2)
     &     -(pwm(2)-xv1*wm(3))*wp(3)-(pwm(3)-xv1*wm(4))*wp(4)
      p13 = (pwm(0)-xv1*wm(1))*w3(1)-(pwm(1)-xv1*wm(2))*w3(2)
     &     -(pwm(2)-xv1*wm(3))*w3(3)-(pwm(3)-xv1*wm(4))*w3(4)
      p21 = (pwp(0)-xv2*wp(1))*wm(1)-(pwp(1)-xv2*wp(2))*wm(2)
     &     -(pwp(2)-xv2*wp(3))*wm(3)-(pwp(3)-xv2*wp(4))*wm(4)
      p23 = (pwp(0)-xv2*wp(1))*w3(1)-(pwp(1)-xv2*wp(2))*w3(2)
     &     -(pwp(2)-xv2*wp(3))*w3(3)-(pwp(3)-xv2*wp(4))*w3(4)
      p31 = (pw3(0)-xv3*w3(1))*wm(1)-(pw3(1)-xv3*w3(2))*wm(2)
     &     -(pw3(2)-xv3*w3(3))*wm(3)-(pw3(3)-xv3*w3(4))*wm(4)
      p32 = (pw3(0)-xv3*w3(1))*wp(1)-(pw3(1)-xv3*w3(2))*wp(2)
     &     -(pw3(2)-xv3*w3(3))*wp(3)-(pw3(3)-xv3*w3(4))*wp(4)

      vertex = -(v12*(p13-p23)+v23*(p21-p31)+v31*(p32-p12))*g
c
      return
      end
      subroutine hioxxx(fi,fo,gc,smass,swidth , hio)
c
c This subroutine computes an off-shell scalar current from an external
c fermion pair.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gc(2)          : coupling constants                 gchf
c       real    smass          : mass  of OUTPUT scalar s
c       real    swidth         : width of OUTPUT scalar s
c
c output:
c       complex hio(3)         : scalar current             j(<fi|s|fo>)
c
      implicit none
      double complex fi(6),fo(6),hio(3),gc(2),dn
      double precision q(0:3),smass,swidth,q2

c

      hio(2) = fo(5)-fi(5)
      hio(3) = fo(6)-fi(6)

      q(0) = dble( hio(2))
      q(1) = dble( hio(3))
      q(2) = dimag(hio(3))
      q(3) = dimag(hio(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)


      dn = -dcmplx( q2-smass**2, smass*swidth )

      hio(1) = ( gc(1)*(fo(1)*fi(1)+fo(2)*fi(2))
     &          +gc(2)*(fo(3)*fi(3)+fo(4)*fi(4)) )/dn
c
      return
      end
      subroutine hsssxx(s1,s2,s3,gc,smass,swidth , hsss)
c
c This subroutine computes an off-shell scalar current from the four-
c scalar coupling.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex s3(3)          : third  scalar                        s3
c       complex gc             : coupling constant                 ghhhh
c       real    smass          : mass  of OUTPUT scalar s'
c       real    swidth         : width of OUTPUT scalar s'
c
c output:
c       complex hsss(3)        : scalar current           j(s':s1,s2,s3)
c     
      implicit none
      double complex s1(3),s2(3),s3(3),gc,hsss(3),dg
      double precision q(0:3),smass,swidth,q2

c

      hsss(2) = s1(2)+s2(2)+s3(2)
      hsss(3) = s1(3)+s2(3)+s3(3)

      q(0) = dble( hsss(2))
      q(1) = dble( hsss(3))
      q(2) = dimag(hsss(3))
      q(3) = dimag(hsss(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)


      dg = -gc/dcmplx( q2-smass**2, smass*swidth )

      hsss(1) = dg * s1(1)*s2(1)*s3(1)
c
      return
      end
      subroutine hsstxx(tc1,sc,gt,xm,xw,jts)
c
c- by RF - Feb. 2006 
c
c This subroutine computes an off-shell tensor current from the tts coupling.
c
c     input:
c          complex tc1               : Incoming tensor particle
c          complex sc                : Incoming scalar particle (Higgs)
c          real    gt                : coupling constant for the tts vertex
c
c     output:
c          complex jts               : off-shell tensor current
c
c     not used:
c          xm, xw
c

      implicit none
c--   dimension of the current set to arbitrary length
      INTEGER DIM
      PARAMETER(DIM=18)
c      include 'dim.inc'
      double complex tc1(DIM),jts(DIM),sc(DIM)
      double precision gt, xm, xw

c The outgoing tensor current is the same as the incoming multiplied by the
c coupling constant and the scalar particle.
c Note that the diagonal tensor terms are always zero because
c the tensor particle is anti-symmetric.

      jts( 1) = 0 !gt * sc(1) * tc1( 1)
      jts( 2) =  gt * sc(1) * tc1( 2)
      jts( 3) =  gt * sc(1) * tc1( 3)
      jts( 4) =  gt * sc(1) * tc1( 4)

      jts( 5) =  gt * sc(1) * tc1( 5)
      jts( 6) = 0 !gt * sc(1) * tc1( 6)
      jts( 7) =  gt * sc(1) * tc1( 7)
      jts( 8) =  gt * sc(1) * tc1( 8)

      jts( 9) =  gt * sc(1) * tc1( 9)
      jts(10) =  gt * sc(1) * tc1(10)
      jts(11) = 0 !gt * sc(1) * tc1(11)
      jts(12) =  gt * sc(1) * tc1(12)

      jts(13) =  gt * sc(1) * tc1(13)
      jts(14) =  gt * sc(1) * tc1(14)
      jts(15) =  gt * sc(1) * tc1(15)
      jts(16) = 0 !gt * sc(1) * tc1(16)

      jts(17) = sc(2) + tc1(17)
      jts(18) = sc(3) + tc1(18)

      return
      end
      subroutine hssxxx(s1,s2,gc,smass,swidth , hss)
c
c This subroutine computes an off-shell scalar current from the three-
c scalar coupling.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex gc             : coupling constant                  ghhh
c       real    smass          : mass  of OUTPUT scalar s'
c       real    swidth         : width of OUTPUT scalar s'
c
c output:
c       complex hss(3)         : scalar current              j(s':s1,s2)
c     
      implicit none
      double complex s1(3),s2(3),gc,hss(3),dg
      double precision q(0:3),smass,swidth,q2

c

      hss(2) = s1(2)+s2(2)
      hss(3) = s1(3)+s2(3)

      q(0) = dble( hss(2))
      q(1) = dble( hss(3))
      q(2) = dimag(hss(3))
      q(3) = dimag(hss(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)


      dg = -gc/dcmplx( q2-smass**2, smass*swidth )

      hss(1) = dg*s1(1)*s2(1)
c
      return
      end
      subroutine hstlxx(s1,t2,gc,smass,swidth , hst)
c- by RF - Mar. 2006
c
c This subroutine computes an off-shell scalar current from the three-
c scalar coupling.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex t2(3)          : internal particle (scalar)           s2
c       complex gc             : coupling constant                  ghhh
c       real    smass          : mass  of OUTPUT scalar s'
c       real    swidth         : width of OUTPUT scalar s'
c
c output:
c       complex hst(3)         : scalar current              j(s':s1,s2)
c     
      implicit none
      include "dimension.inc"

      double complex s1(DIM),t2(DIM),hst(DIM),dg
      double precision q(0:3),smass,swidth,q2,gc

      hst(2) = s1(2)+t2(2)
      hst(3) = s1(3)+t2(3)

      q(0) = -dble( hst(2))
      q(1) = -dble( hst(3))
      q(2) = -dimag(hst(3))
      q(3) = -dimag(hst(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)


      dg = gc /dcmplx( q2-smass**2, smass*swidth )
      

      hst(1) = dg*s1(1)*t2(1)

      return
      end
      subroutine hstxxx(tc,sc,gt,smass,swidth , hst)
c
c This subroutine computes an off-shell scalar current from
c the scalar-scalar-tensor boson coupling.
c
c input:
c       complex tc(18)         : input tensor                          T
c       complex sc(3)          : input scalar                          s
c       complex gt             : coupling constant         gts=-1/Lambda
c       real    smass          : mass  of output scalar s'
c       real    swidth         : width of output scalar s'
c
c output:
c       complex hst(3)         : scalar current                j(s':T,s)     
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex tc(18), sc(3), hst(3)
      double precision smass, swidth
      double complex gt

      double complex ft(6,4)
      double complex T12, T13, T14, T23, T24, T34
      double complex TKK
      double precision ps1(4), ps2(4)
      integer i
      double complex cZero, d
      double precision rZero, rTwo, pf2
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      hst(2) = sc(2)+ft(5,1)
      hst(3) = sc(3)+ft(6,1)

      ps1(1) = dreal(sc(2))
      ps1(2) = dreal(sc(3))
      ps1(3) = dimag(sc(3))
      ps1(4) = dimag(sc(2))

      ps2(1) = dreal(hst(2))
      ps2(2) = dreal(hst(3))
      ps2(3) = dimag(hst(3))
      ps2(4) = dimag(hst(2))

      pf2 = ps2(1)**2 - ps2(2)**2 - ps2(3)**2 - ps2(4)**2
      
      if ( smass.gt.rZero ) then
         d = - gt/dcmplx( pf2-smass**2, smass*swidth )
      else
         d = - gt/dcmplx( pf2, rZero )
      end if

      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      TKK   = cZero
    
      do i = 1,4
         TKK=TKK+ft(i,i)*ps1(i)*ps2(i)
      end do

      TKK   = rTwo*TKK
    
      TKK = TKK - T12*(ps1(1)*ps2(2) + ps1(2)*ps2(1))
     &          - T13*(ps1(1)*ps2(3) + ps1(3)*ps2(1))
     &          - T14*(ps1(1)*ps2(4) + ps1(4)*ps2(1))
     &          + T23*(ps1(2)*ps2(3) + ps1(3)*ps2(2))
     &          + T24*(ps1(2)*ps2(4) + ps1(4)*ps2(2))
     &          + T34*(ps1(3)*ps2(4) + ps1(4)*ps2(3))

      hst(1) = TKK+(ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4))
     &	*(smass**2-ps1(1)*ps2(1)+ps1(2)*ps2(2)
     &      +ps1(3)*ps2(3)+ps1(4)*ps2(4))

      hst(1) = hst(1) * d*sc(1)

      return
      end
      subroutine httaxx(tc1,tc2,gc,mass,width,jsc)
c
c- by RF - Mar. 2006 
c
c This subroutine computes an off-shell scalar from the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex tc2(18)           : Incoming tensor particle
c          complex gc(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c          real    mass              : mass of the outgoing scalar
c          real    width             : width of the outgoing scalar
c
c     output:
c          complex sc(3)             : Incoming scalar particle
c
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),tc2(DIM),jsc(DIM)
      double complex vertex,dj,gc(2)
      double precision mass,width,q2,q(4)

c Take the inner product between the tensor particles. The 
c Note that the tensor particle is antisymmetric, thus all diagonal terms
c are zero.

      jsc(2)=tc1(17)+tc2(17)
      jsc(3)=tc1(18)+tc2(18)

      if (gc(1).NE.(0D0,0D0)) then

      q(1) = -dble( jsc(2))
      q(2) = -dble( jsc(3))
      q(3) = -dimag(jsc(3))
      q(4) = -dimag(jsc(2))

      q2 = q(1)**2 - q(2)**2 - q(3)**2 - q(4)**2

      dj = gc(1) /dcmplx( q2-mass**2, mass*width )

      jsc(1) = dj* (
c     &                       + tc1( 1) * tc2( 1)
     &                       - tc1( 2) * tc2( 2)
     &                       - tc1( 3) * tc2( 3)
     &                       - tc1( 4) * tc2( 4)

     &                       - tc1( 5) * tc2( 5)
c     &                       + tc1( 6) * tc2( 6)
     &                       + tc1( 7) * tc2( 7)
     &                       + tc1( 8) * tc2( 8)

     &                       - tc1( 9) * tc2( 9)
     &                       + tc1(10) * tc2(10)
c     &                       + tc1(11) * tc2(11)
     &                       + tc1(12) * tc2(12)

     &                       - tc1(13) * tc2(13)
     &                       + tc1(14) * tc2(14)
     &                       + tc1(15) * tc2(15)
c     &                       + tc1(16) * tc2(16)
     &                                           )

      else
         jsc( 1)=(0D0,0D0)
      endif

      return
      end
      subroutine httcxx(tc1,tc2,gc,mass,width,jsc)
c
c- by RF - Mar. 2006 
c  CP3  Modified Nov. 2009 

c This subroutine computes an off-shell scalar from the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex tc2(18)           : Incoming tensor particle
c          complex gc(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c          real    mass              : mass of the outgoing scalar
c          real    width             : width of the outgoing scalar
c
c     output:
c          complex sc(3)             : Incoming scalar particle
c
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),tc2(DIM),jsc(DIM),T1(6,4),T2(6,4)
      double complex vertex,dj,gc(2)
      double precision mass,width,q2,q(4)

c Take the inner product between the tensor particles. The 
c Note that the tensor particle is antisymmetric, thus all diagonal terms
c are zero.

      T1(1,1) = tc1(1)
      T1(1,2) = tc1(2)
      T1(1,3) = tc1(3)
      T1(1,4) = tc1(4)
      T1(2,1) = tc1(5)
      T1(2,2) = tc1(6)
      T1(2,3) = tc1(7)
      T1(2,4) = tc1(8)
      T1(3,1) = tc1(9)
      T1(3,2) = tc1(10)
      T1(3,3) = tc1(11)
      T1(3,4) = tc1(12)
      T1(4,1) = tc1(13)
      T1(4,2) = tc1(14)
      T1(4,3) = tc1(15)
      T1(4,4) = tc1(16)
      T1(5,1) = tc1(17)
      T1(6,1) = tc1(18)

      T2(1,1) = tc2(1)
      T2(1,2) = tc2(2)
      T2(1,3) = tc2(3)
      T2(1,4) = tc2(4)
      T2(2,1) = tc2(5)
      T2(2,2) = tc2(6)
      T2(2,3) = tc2(7)
      T2(2,4) = tc2(8)
      T2(3,1) = tc2(9)
      T2(3,2) = tc2(10)
      T2(3,3) = tc2(11)
      T2(3,4) = tc2(12)
      T2(4,1) = tc2(13)
      T2(4,2) = tc2(14)
      T2(4,3) = tc2(15)
      T2(4,4) = tc2(16)
      T2(5,1) = tc2(17)
      T2(6,1) = tc2(18)

      jsc(2)=tc1(17)+tc2(17)
      jsc(3)=tc1(18)+tc2(18)

      if (gc(1).NE.(0D0,0D0)) then

      q(1) = -dble( jsc(2))
      q(2) = -dble( jsc(3))
      q(3) = -dimag(jsc(3))
      q(4) = -dimag(jsc(2))

      q2 = q(1)**2 - q(2)**2 - q(3)**2 - q(4)**2

      dj = -gc(1) /dcmplx( q2-mass**2, mass*width )

      jsc(1) = dj* ( T1(1,2)*T2(1,2) - T1(2,1)*T2(1,2) 
     -  + T1(1,3)*T2(1,3) - 
     -  T1(3,1)*T2(1,3) + T1(1,4)*T2(1,4) - T1(4,1)*T2(1,4) - 
     -  T1(1,2)*T2(2,1) + T1(2,1)*T2(2,1) - T1(2,3)*T2(2,3) + 
     -  T1(3,2)*T2(2,3) - T1(2,4)*T2(2,4) + T1(4,2)*T2(2,4) - 
     -  T1(1,3)*T2(3,1) + T1(3,1)*T2(3,1) + T1(2,3)*T2(3,2) - 
     -  T1(3,2)*T2(3,2) - T1(3,4)*T2(3,4) + T1(4,3)*T2(3,4) - 
     -  T1(1,4)*T2(4,1) + T1(4,1)*T2(4,1) + T1(2,4)*T2(4,2) - 
     -  T1(4,2)*T2(4,2) + T1(3,4)*T2(4,3) - T1(4,3)*T2(4,3)
     &                                           )

      else
         jsc( 1)=(0D0,0D0)
      endif

      return
      end
      subroutine httsxx(tc1,tc2,sc,g1,g2,mass,width,htts)
c
c- by RF - Mar. 2006 
c
c This subroutine computes an off-shell tensor current from the ttss coupling.
c
c     input:
c          complex tc1(18)           : first incoming tensor particle
c          complex tc2(18)           : second incoming tensor particle
c          complex sc(3)             : Incoming scalar particle
c          complex g1(2)             : coupling constant (Higgs effc. theor)
c          real    g2                : coupling constant (include extra Higgs)
c          real    mass              : mass of the outgoing scalar
c          real    width             : width of the outgoing scalar
c
c     output:
c          complex htts              : off-shell tensor current
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),htts(DIM),sc(DIM),tc2(DIM)
      double complex dg,g1(2)

      double precision g2,mass,width,q2,q(0:3)

c The outgoing tensor current is the same as the incoming multiplied by the
c coupling constants and the scalar particles.
c Note that the diagonal tensor terms are always zero because
c the tensor particle is anti-symmetric.

      
      htts(2) = sc(2) + tc2(17) + tc1(17)
      htts(3) = sc(3) + tc2(18) + tc1(18)

 
      if (g1(1).NE.(0D0,0D0)) then

      q(0) = -dble( htts(2))
      q(1) = -dble( htts(3))
      q(2) = -dimag(htts(3))
      q(3) = -dimag(htts(2))

      q2 = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2

      dg = - g1(1)*g2/dcmplx( q2-mass**2, mass*width )

      htts(1)= dg*sc(1)*(
c     &                   + tc1( 1) * tc2( 1)
     &                   - tc1( 2) * tc2( 2)
     &                   - tc1( 3) * tc2( 3)
     &                   - tc1( 4) * tc2( 4)

     &                   - tc1( 5) * tc2( 5)
c     &                   + tc1( 6) * tc2( 6)
     &                   + tc1( 7) * tc2( 7)
     &                   + tc1( 8) * tc2( 8)

     &                   - tc1( 9) * tc2( 9)
     &                   + tc1(10) * tc2(10)
c     &                   + tc1(11) * tc2(11)
     &                   + tc1(12) * tc2(12)

     &                   - tc1(13) * tc2(13)
     &                   + tc1(14) * tc2(14)
     &                   + tc1(15) * tc2(15)
c     &                   + tc1(16) * tc2(16)
     &                                       )


      else
         htts( 1)=(0D0,0D0)
      endif


      return
      end
      subroutine httxxx(tc1,tc2,gc,mass,width,jsc)
c
c- by RF - Mar. 2006 
c
c This subroutine computes an off-shell scalar from the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex tc2(18)           : Incoming tensor particle
c          complex gc(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c          real    mass              : mass of the outgoing scalar
c          real    width             : width of the outgoing scalar
c
c     output:
c          complex sc(3)             : Incoming scalar particle
c
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),tc2(DIM),jsc(DIM)
      double complex vertex,dj,gc(2)
      double precision mass,width,q2,q(4)

c Take the inner product between the tensor particles. The 
c Note that the tensor particle is antisymmetric, thus all diagonal terms
c are zero.

      jsc(2)=tc1(17)+tc2(17)
      jsc(3)=tc1(18)+tc2(18)

      if (gc(1).NE.(0D0,0D0)) then

      q(1) = -dble( jsc(2))
      q(2) = -dble( jsc(3))
      q(3) = -dimag(jsc(3))
      q(4) = -dimag(jsc(2))

      q2 = q(1)**2 - q(2)**2 - q(3)**2 - q(4)**2

      dj = gc(1) /dcmplx( q2-mass**2, mass*width )

      jsc(1) = dj* (
c     &                       + tc1( 1) * tc2( 1)
     &                       - tc1( 2) * tc2( 2)
     &                       - tc1( 3) * tc2( 3)
     &                       - tc1( 4) * tc2( 4)

     &                       - tc1( 5) * tc2( 5)
c     &                       + tc1( 6) * tc2( 6)
     &                       + tc1( 7) * tc2( 7)
     &                       + tc1( 8) * tc2( 8)

     &                       - tc1( 9) * tc2( 9)
     &                       + tc1(10) * tc2(10)
c     &                       + tc1(11) * tc2(11)
     &                       + tc1(12) * tc2(12)

     &                       - tc1(13) * tc2(13)
     &                       + tc1(14) * tc2(14)
     &                       + tc1(15) * tc2(15)
c     &                       + tc1(16) * tc2(16)
     &                                           )

      else
         jsc( 1)=(0D0,0D0)
      endif

      return
      end
      subroutine hvsxxx(vc,sc,gc,smass,swidth , hvs)
c
c This subroutine computes an off-shell scalar current from the vector-
c scalar-scalar coupling.  The coupling is absent in the minimal SM in
c unitary gauge.
c
c input:
c       complex vc(6)          : input vector                          v
c       complex sc(3)          : input scalar                          s
c       complex gc             : coupling constant (s charge)
c       real    smass          : mass  of OUTPUT scalar s'
c       real    swidth         : width of OUTPUT scalar s'
c
c examples of the coupling constant gc for susy particles are as follows:
c   -----------------------------------------------------------
c   |    s1    | (q,i3) of s1  ||   v=A   |   v=Z   |   v=W   |
c   -----------------------------------------------------------
c   | nu~_l    | (  0  , +1/2) ||   ---   |  gzn(1) |  gwf(1) |
c   | e~_l     | ( -1  , -1/2) ||  gal(1) |  gzl(1) |  gwf(1) |
c   | u~_l     | (+2/3 , +1/2) ||  gau(1) |  gzu(1) |  gwf(1) |
c   | d~_l     | (-1/3 , -1/2) ||  gad(1) |  gzd(1) |  gwf(1) |
c   -----------------------------------------------------------
c   | e~_r-bar | ( +1  ,  0  ) || -gal(2) | -gzl(2) | -gwf(2) |
c   | u~_r-bar | (-2/3 ,  0  ) || -gau(2) | -gzu(2) | -gwf(2) |
c   | d~_r-bar | (+1/3 ,  0  ) || -gad(2) | -gzd(2) | -gwf(2) |
c   -----------------------------------------------------------
c where the sc charge is defined by the flowing-OUT quantum number.
c
c output:
c       complex hvs(3)         : scalar current                j(s':v,s)
c     
      implicit none
      double complex vc(6),sc(3),hvs(3),dg,qvv,qpv,gc
      double precision qv(0:3),qp(0:3),qa(0:3),smass,swidth,q2

      double precision rTwo
      parameter( rTwo = 2.0d0 )

c

      hvs(2) = vc(5)+sc(2)
      hvs(3) = vc(6)+sc(3)

      qv(0) = dble(  vc(5))
      qv(1) = dble(  vc(6))
      qv(2) = dimag( vc(6))
      qv(3) = dimag( vc(5))
      qp(0) = dble(  sc(2))
      qp(1) = dble(  sc(3))
      qp(2) = dimag( sc(3))
      qp(3) = dimag( sc(2))
      qa(0) = dble( hvs(2))
      qa(1) = dble( hvs(3))
      qa(2) = dimag(hvs(3))
      qa(3) = dimag(hvs(2))
      q2 = qa(0)**2-(qa(1)**2+qa(2)**2+qa(3)**2)


      dg = -gc/dcmplx( q2-smass**2, smass*swidth )
      qvv = qv(0)*vc(1)-qv(1)*vc(2)-qv(2)*vc(3)-qv(3)*vc(4)
      qpv = qp(0)*vc(1)-qp(1)*vc(2)-qp(2)*vc(3)-qp(3)*vc(4)

      hvs(1) = dg*(rTwo*qpv+qvv)*sc(1)
c
      return
      end
      subroutine hvvhxx(v1,v2,gc,smass,swidth , jsvv)
c
c- by RF - Mar. 2006
c
c
c This subroutine computes an off-shell (pseudo-)scalar current from
c two incoming vectors using the scalar effective coupling.
c
c input:
c       complex v1(6)          : first  vector
c       complex v2(6)          : second vector
c       complex gc(2)          : coupling constant: gc(1) scalar
c                                                   gc(2) pseudo-scalar
c       real smass             : mass of the outgoing (pseudo-)scalar
c       real swidth            : width of the outgoing (pseudo-)scalar
c
c output:
c       complex jsvv(3)        : (pseudo-)scalar current
c     
c
      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex v1(DIM),v2(DIM),jsvv(DIM),vertex1,vertex2,dj
      double complex v12,p2v1,p1v2	
      double complex v13,v14,v23,v24,v34
      double precision p12,p13,p14,p23,p24,p34
      double precision p1(0:3),p2(0:3),q(4),q2
      double precision smass,swidth
      double complex gc(2)

      p1(0) = dble( v1(5))
      p1(1) = dble( v1(6))
      p1(2) = dimag(v1(6))
      p1(3) = dimag(v1(5))

      p2(0) = dble( v2(5))
      p2(1) = dble( v2(6))
      p2(2) = dimag(v2(6))
      p2(3) = dimag(v2(5))

      vertex1 = (0D0,0D0)
      vertex2 = (0D0,0D0)
      
      jsvv(2) = v1(5) + v2(5)
      jsvv(3) = v1(6) + v2(6)

      q(1) = -dble( jsvv(2))
      q(2) = -dble( jsvv(3))
      q(3) = -dimag(jsvv(3))
      q(4) = -dimag(jsvv(2))

      q2 = q(1)**2 - q(2)**2 - q(3)**2 - q(4)**2

      dj = dcmplx( q2-smass**2, smass*swidth )


      if (gc(1).NE.(0D0,0D0)) then

         v12  = v1(1)*v2(1) - v1(2)*v2(2) - v1(3)*v2(3) - v1(4)*v2(4)
         p12  = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
         p2v1 = v1(1)*p2(0) - v1(2)*p2(1) - v1(3)*p2(2) - v1(4)*p2(3)
         p1v2 = p1(0)*v2(1) - p1(1)*v2(2) - p1(2)*v2(3) - p1(3)*v2(4)	

         vertex1 = gc(1)*(v12*p12 - p2v1*p1v2)
      endif

      if (gc(2).NE.(0D0,0D0)) then
          p12 = p1(0)*p2(1) - p1(1)*p2(0)
          p13 = p1(0)*p2(2) - p1(2)*p2(0)
          p14 = p1(0)*p2(3) - p1(3)*p2(0)
          p23 = p1(1)*p2(2) - p1(2)*p2(1)
          p24 = p1(1)*p2(3) - p1(3)*p2(1)
          p34 = p1(2)*p2(3) - p1(3)*p2(2)

          v12 = v1(1)*v2(2) - v1(2)*v2(1)
          v13 = v1(1)*v2(3) - v1(3)*v2(1)
          v14 = v1(1)*v2(4) - v1(4)*v2(1)
          v23 = v1(2)*v2(3) - v1(3)*v2(2)
          v24 = v1(2)*v2(4) - v1(4)*v2(2)
          v34 = v1(3)*v2(4) - v1(4)*v2(3)

          vertex2 = - gc(2)*( v12*p34 - v13*p24 + v14*p23
     &                       +v23*p14 - v24*p13 + v34*p12 )
      endif
       
      jsvv(1) = (vertex1 + vertex2) /dj


      return
      end
      subroutine hvvshx(v1,v2,sc,g1,mass,width, hvvsh)
c
c- by RF - Mar. 2006
c
c
c This subroutine computes an amplitude of the vector-vector-Higgs-Higgs 
c effective coupling.
c
c input:
c       complex v1(6)          : first  vector                        
c       complex v2(6)          : second vector                        
c       complex sc1(3)         : first  scalar                        
c       complex g1(2)          : first coupling constant                 
c       real    mass           : mass of the outgoing scalar
c       real    width          : width of the outgoing scalar
c
c output:
c       complex hvvsh(3)       : scalar current
c     
      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex v1(DIM),v2(DIM),sc(DIM),hvvsh(DIM)
      double complex hvvsh1,hvvsh2,dg,g1(2)
      double complex v12,p2v1,p1v2,v13,v14,v23,v24,v34
      double precision p12,p13,p14,p23,p24,p34
      
      double precision p1(0:3),p2(0:3),mass,width,q2,q(0:3)


      p1(0) = dble( v1(5))
      p1(1) = dble( v1(6))
      p1(2) = dimag(v1(6))
      p1(3) = dimag(v1(5))

      p2(0) = dble( v2(5))
      p2(1) = dble( v2(6))
      p2(2) = dimag(v2(6))
      p2(3) = dimag(v2(5))

      hvvsh(2) = v1(5)+v2(5)+sc(2)
      hvvsh(3) = v1(6)+v2(6)+sc(3)

      q(0) = -dble( hvvsh(2))
      q(1) = -dble( hvvsh(3))
      q(2) = -dimag(hvvsh(3))
      q(3) = -dimag(hvvsh(2))

      q2 = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2

      dg = dcmplx( q2-mass**2, mass*width )

      hvvsh1 = (0D0,0D0)
      hvvsh2 = (0D0,0D0)

      if (g1(1).NE.(0D0,0D0)) then

         v12  = v1(1)*v2(1) - v1(2)*v2(2) - v1(3)*v2(3) - v1(4)*v2(4)
         p12  = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
         p2v1 = v1(1)*p2(0) - v1(2)*p2(1) - v1(3)*p2(2) - v1(4)*p2(3)
         p1v2 = p1(0)*v2(1) - p1(1)*v2(2) - p1(2)*v2(3) - p1(3)*v2(4)	

         hvvsh1 = - g1(1) *(v12*p12 - p2v1*p1v2)
      endif

      if (g1(2).NE.(0D0,0D0)) then
          p12 = p1(0)*p2(1) - p1(1)*p2(0)
          p13 = p1(0)*p2(2) - p1(2)*p2(0)
          p14 = p1(0)*p2(3) - p1(3)*p2(0)
          p23 = p1(1)*p2(2) - p1(2)*p2(1)
          p24 = p1(1)*p2(3) - p1(3)*p2(1)
          p34 = p1(2)*p2(3) - p1(3)*p2(2)

          v12 = v1(1)*v2(2) - v1(2)*v2(1)
          v13 = v1(1)*v2(3) - v1(3)*v2(1)
          v14 = v1(1)*v2(4) - v1(4)*v2(1)
          v23 = v1(2)*v2(3) - v1(3)*v2(2)
          v24 = v1(2)*v2(4) - v1(4)*v2(2)
          v34 = v1(3)*v2(4) - v1(4)*v2(3)

          hvvsh2 = g1(2)*( v12*p34 - v13*p24 + v14*p23
     &                    +v23*p14 - v24*p13 + v34*p12 )
      endif
       
      hvvsh(1) = sc(1)*(hvvsh1 + hvvsh2) /dg

      return
      end
      subroutine hvvsxx(v1,v2,sc,gc,smass,swidth , hvvs)
c
c This subroutine computes an off-shell scalar current of the vector-
c vector-scalar-scalar coupling.
c
c input:
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       complex sc(3)          : input  scalar                        s
c       complex gc             : coupling constant                 gvvhh
c       real    smass          : mass  of OUTPUT scalar s'
c       real    swidth         : width of OUTPUT scalar s'
c
c output:
c       complex hvvs(3)        : scalar current            j(s':v1,v2,s)
c     
      implicit none
      double complex v1(6),v2(6),sc(3),gc,hvvs(3),dg
      double precision q(0:3),smass,swidth,q2

c

      hvvs(2) = v1(5)+v2(5)+sc(2)
      hvvs(3) = v1(6)+v2(6)+sc(3)

      q(0) = dble( hvvs(2))
      q(1) = dble( hvvs(3))
      q(2) = dimag(hvvs(3))
      q(3) = dimag(hvvs(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)


      dg = -gc/dcmplx( q2-smass**2, smass*swidth )

      hvvs(1) = dg*sc(1)
     &         *(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
c
      return
      end
      subroutine hvvtxx(w1,w2,g,xm,xw,jt)
c
c- by RF - Feb. 2006
c
c This subroutine computes the portion of the off-shell current
c for the color-octect tensor t3 in terms of w1 and w2 
c
c input:
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       real    g              : first  coupling constant
c       real    xm             : not used
c       real    xw             : not used
c
c output:
c       complex jt(18)        : tensor current  j^(mu,nu)(w':w1,w2,w3)
c
      implicit none

c dimension of the current set to arbitrary length
      integer DIM
      parameter (DIM=18)
c      include "dimension.inc"
      double complex w1(DIM),w2(DIM),jt(DIM)

      double precision xm,xw,g,s2g

      double precision sqrTwo
      parameter( sqrTwo = 1.41421356237309514547462185873882845044d0 )

      s2g = g * sqrTwo
      
      jt( 1) = 0 ! dcmplx( s2g * (w1(1)*w2(1)-w1(1)*w2(1)) )
      jt( 2) =  dcmplx( s2g * (w1(1)*w2(2)-w1(2)*w2(1)) )
      jt( 3) =  dcmplx( s2g * (w1(1)*w2(3)-w1(3)*w2(1)) )
      jt( 4) =  dcmplx( s2g * (w1(1)*w2(4)-w1(4)*w2(1)) )

      jt( 5) =  dcmplx( s2g * (w1(2)*w2(1)-w1(1)*w2(2)) )
      jt( 6) = 0 ! dcmplx( s2g * (w1(2)*w2(2)-w1(2)*w2(2)) )
      jt( 7) =  dcmplx( s2g * (w1(2)*w2(3)-w1(3)*w2(2)) )
      jt( 8) =  dcmplx( s2g * (w1(2)*w2(4)-w1(4)*w2(2)) )

      jt( 9) =  dcmplx( s2g * (w1(3)*w2(1)-w1(1)*w2(3)) )
      jt(10) =  dcmplx( s2g * (w1(3)*w2(2)-w1(2)*w2(3)) )
      jt(11) = 0 ! dcmplx( s2g * (w1(3)*w2(3)-w1(3)*w2(3)) )
      jt(12) =  dcmplx( s2g * (w1(3)*w2(4)-w1(4)*w2(3)) )

      jt(13) =  dcmplx( s2g * (w1(4)*w2(1)-w1(1)*w2(4)) )
      jt(14) =  dcmplx( s2g * (w1(4)*w2(2)-w1(2)*w2(4)) )
      jt(15) =  dcmplx( s2g * (w1(4)*w2(3)-w1(3)*w2(4)) )
      jt(16) = 0 ! dcmplx( s2g * (w1(4)*w2(4)-w1(4)*w2(4)) )

      jt(17) = w1(5) + w2(5)
      jt(18) = w1(6) + w2(6)

      return
      end
      subroutine hvvvxx(ga,gb,gc,g1,g2,mass,width,jhvvv)
c
c- by RF - Mar. 2006
c
c This subroutine computes an off-shell (pseudo-)scalar current
c from the coupling of three gauge bosons
c
c input:
c       complex ga(6)          : first  incoming vector   (gluon)
c       complex gb(6)          : second incoming vector   (gluon)
c       complex gc(6)          : third  incoming vector   (gluon)
c       real    g1             : coupling constant        (QCD)
c       complex g2             : coupling constant: g2(1) scalar
c                                                   g2(2) pseudo-scalar
c       real mass              : mass of the outgoing scalar
c       real width             : width of the outgoing scalar
c
c output:
c       complex jhvvv(3)       : output (pseudo-)scalar current
c

      implicit none

c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex ga(DIM),gb(DIM),gc(DIM),jhvvv(DIM)

      double complex dvertx, vertex, vertex1, vertex2,dj
      double complex vab, vbc, vca, v123, v124, v134, v234
      double complex pgagb, pgagc, pgbga, pgbgc, pgcga, pgcgb
      double precision pga(0:3),pgb(0:3),pgc(0:3),pabc(4)
      double precision g1,mass, width, q2, q(4)
      double complex g2(2)

      pga(0) = dble( ga(5))
      pga(1) = dble( ga(6))
      pga(2) = dimag(ga(6))
      pga(3) = dimag(ga(5))

      pgb(0) = dble( gb(5))
      pgb(1) = dble( gb(6))
      pgb(2) = dimag(gb(6))
      pgb(3) = dimag(gb(5))

      pgc(0) = dble( gc(5))
      pgc(1) = dble( gc(6))
      pgc(2) = dimag(gc(6))
      pgc(3) = dimag(gc(5))

      vertex1 = (0D0,0D0)
      vertex2 = (0D0,0D0)


      jhvvv(2) = ga(5) + gb(5) + gc(5)
      jhvvv(3) = ga(6) + gb(6) + gc(6)


      q(1) = -dble( jhvvv(2))
      q(2) = -dble( jhvvv(3))
      q(3) = -dimag(jhvvv(3))
      q(4) = -dimag(jhvvv(2))

      q2 = q(1)**2 - q(2)**2 - q(3)**2 - q(4)**2

      dj = g1 /dcmplx( q2-mass**2, mass*width )


      if (g2(1).NE.(0D0,0D0)) then
      vab = ga(1)*gb(1)-ga(2)*gb(2)-ga(3)*gb(3)-ga(4)*gb(4)
      vbc = gb(1)*gc(1)-gb(2)*gc(2)-gb(3)*gc(3)-gb(4)*gc(4)
      vca = gc(1)*ga(1)-gc(2)*ga(2)-gc(3)*ga(3)-gc(4)*ga(4)

      pgagb = pga(0)*gb(1) - pga(1)*gb(2) - pga(2)*gb(3) - pga(3)*gb(4)
      pgagc = pga(0)*gc(1) - pga(1)*gc(2) - pga(2)*gc(3) - pga(3)*gc(4)
      pgbga = pgb(0)*ga(1) - pgb(1)*ga(2) - pgb(2)*ga(3) - pgb(3)*ga(4)
      pgbgc = pgb(0)*gc(1) - pgb(1)*gc(2) - pgb(2)*gc(3) - pgb(3)*gc(4)
      pgcga = pgc(0)*ga(1) - pgc(1)*ga(2) - pgc(2)*ga(3) - pgc(3)*ga(4)
      pgcgb = pgc(0)*gb(1) - pgc(1)*gb(2) - pgc(2)*gb(3) - pgc(3)*gb(4)

      dvertx = vab*(pgagc-pgbgc) + vbc*(pgbga-pgcga) + vca*(pgcgb-pgagb)
      vertex1= - dvertx * g2(1)
      endif

      if (g2(2).NE.(0D0,0D0)) then
      pabc(1) = pga(0) + pgb(0) + pgc(0)
      pabc(2) = pga(1) + pgb(1) + pgc(1)
      pabc(3) = pga(2) + pgb(2) + pgc(2)
      pabc(4) = pga(3) + pgb(3) + pgc(3)

      v123 =   ga(1)*gb(2)*gc(3) - ga(1)*gb(3)*gc(2) - ga(2)*gb(1)*gc(3)
     &       + ga(2)*gb(3)*gc(1) + ga(3)*gb(1)*gc(2) - ga(3)*gb(2)*gc(1)
      v124 = - ga(1)*gb(2)*gc(4) + ga(1)*gb(4)*gc(2) + ga(2)*gb(1)*gc(4)
     &       - ga(2)*gb(4)*gc(1) - ga(4)*gb(1)*gc(2) + ga(4)*gb(2)*gc(1)
      v134 =   ga(1)*gb(3)*gc(4) - ga(1)*gb(4)*gc(3) - ga(3)*gb(1)*gc(4)
     &       + ga(3)*gb(4)*gc(1) + ga(4)*gb(1)*gc(3) - ga(4)*gb(3)*gc(1)
      v234 = - ga(2)*gb(3)*gc(4) + ga(2)*gb(4)*gc(3) + ga(3)*gb(2)*gc(4)
     &       - ga(3)*gb(4)*gc(2) - ga(4)*gb(2)*gc(3) + ga(4)*gb(3)*gc(2)


      vertex2= - g2(2) * (  v123*pabc(4) + v124*pabc(3)
     &                    + v134*pabc(2) + v234*pabc(1) )
      endif

      jhvvv(1) = dj * (vertex1 + vertex2)

      return
      end
      subroutine hvvxxx(v1,v2,gc,smass,swidth , hvv)
c
c This subroutine computes an off-shell scalar current from the vector-
c vector-scalar coupling.
c
c input:
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       complex gc             : coupling constant                  gvvh
c       real    smass          : mass  of OUTPUT scalar s
c       real    swidth         : width of OUTPUT scalar s
c
c output:
c       complex hvv(3)         : off-shell scalar current     j(s:v1,v2)
c     
      implicit none
      double complex v1(6),v2(6),gc,hvv(3),dg
      double precision q(0:3),smass,swidth,q2

      double precision rZero
      parameter( rZero = 0.0d0 )

c

      hvv(2) = v1(5)+v2(5)
      hvv(3) = v1(6)+v2(6)

      q(0) = dble( hvv(2))
      q(1) = dble( hvv(3))
      q(2) = dimag(hvv(3))
      q(3) = dimag(hvv(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)


      dg = -gc/dcmplx( q2-smass**2, smass*swidth )

      hvv(1) = dg*(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
c
      return
      end
      subroutine iosxxx(fi,fo,sc,gc , vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-scalar
c coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex sc(3)          : input    scalar                      s
c       complex gc(2)          : coupling constants                 gchf
c
c output:
c       complex vertex         : amplitude                     <fo|s|fi>
c     
      implicit none
      double complex fi(6),fo(6),sc(3),gc(2),vertex

c

      vertex = sc(1)*( gc(1)*(fi(1)*fo(1)+fi(2)*fo(2))
     &                +gc(2)*(fi(3)*fo(3)+fi(4)*fo(4)) )
c
      return
      end
      subroutine iotxkk(fi,fo,tc,g,fmass , vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex tc(6,4)        : input    tensor                      t
c       real    g              : coupling constant                 -kappa/8
c       real    fmass          : fermion mass                        m_f
c
c output:
c       complex vertex         : amplitude                        <fo|t|fi>
c     
      implicit none
      double complex fi(6), fo(6), tc(6,4), vertex
      double precision g, fmass

      double complex k23, k23s, D1, D2, D3, D4, Tii
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double complex f13, f14, f23, f24, f31, f32, f41, f42
      double precision k(4), k14p, k14m, m2

      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
c
      m2 = rTwo*fmass

      k(1) = dreal(fi(5)+fo(5))
      k(2) = dreal(fi(6)+fo(6))
      k(3) = dimag(fi(6)+fo(6))
      k(4) = dimag(fi(5)+fo(5))
      k23  = dcmplx( k(2),k(3) )
      k23s = dconjg( k23 )
      k14p = k(1) + k(4)
      k14m = k(1) - k(4)

      f13 = fo(1)*fi(3)
      f14 = fo(1)*fi(4)
      f23 = fo(2)*fi(3)
      f24 = fo(2)*fi(4)
      f31 = fo(3)*fi(1)
      f32 = fo(3)*fi(2)
      f41 = fo(4)*fi(1)
      f42 = fo(4)*fi(2)

      T11 = rTwo*tc(1,1)
      T22 = rTwo*tc(2,2)
      T33 = rTwo*tc(3,3)
      T44 = rTwo*tc(4,4)
      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      D1 =   k(1)*(T11-T14) - k(2)*(T12-T24)
     &     - k(3)*(T13-T34) - k(4)*(T14-T44)

      D2 = - k(1)*(T12-ci*T13) + k(2)*(T22-ci*T23)
     &     + k(3)*(T23-ci*T33) + k(4)*(T24-ci*T34)

      D3 = - k(1)*(T12+ci*T13) + k(2)*(T22+ci*T23)
     &     + k(3)*(T23+ci*T33) + k(4)*(T24+ci*T34)

      D4 =   k(1)*(T11+T14) - k(2)*(T12+T24)
     &     - k(3)*(T13+T34) - k(4)*(T14+T44)

      Tii = T11 - T22 - T33 - T44

      vertex = D1*(f13+f42) + D2*(f14-f32) + D3*(f23-f41) + D4*(f24+f31)

      vertex = vertex + Tii*( - k14p*(f24+f31) - k14m*(f13+f42)
     &                        +  k23*(f23-f41) + k23s*(f14-f32) )

      if ( fmass.ne.rZero ) then
         vertex = vertex + m2*Tii*(  fo(1)*fi(1) + fo(2)*fi(2)
     &                             + fo(3)*fi(3) + fo(4)*fi(4) )
      end if

      vertex = vertex * g
c
      return
      end
      subroutine iotxxx(fi,fo,tc,gt,fmass , vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-tensor
c coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex tc(18)         : input    tensor                       T
c       complex gt             : coupling constant       gtf=-1/Lambda/4
c       real    fmass          : fermion mass                        m_f
c
c output:
c       complex vertex         : amplitude                     <fo|T|fi>
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fi(6), fo(6), tc(18), gt, vertex
      double precision fmass

      double complex ft(6,4)
      double complex k23, k23s, D1, D2, D3, D4, Tii
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double complex f13, f14, f23, f24, f31, f32, f41, f42
      double precision k(4), k14p, k14m, m2

      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )


      m2 = rTwo*fmass

      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      k(1) = dreal(fi(5)+fo(5))
      k(2) = dreal(fi(6)+fo(6))
      k(3) = dimag(fi(6)+fo(6))
      k(4) = dimag(fi(5)+fo(5))
      k23  = dcmplx( k(2),k(3) )
      k23s = dconjg( k23 )
      k14p = k(1) + k(4)
      k14m = k(1) - k(4)

      f13 = fo(1)*fi(3)
      f14 = fo(1)*fi(4)
      f23 = fo(2)*fi(3)
      f24 = fo(2)*fi(4)
      f31 = fo(3)*fi(1)
      f32 = fo(3)*fi(2)
      f41 = fo(4)*fi(1)
      f42 = fo(4)*fi(2)

      T11 = rTwo*ft(1,1)
      T22 = rTwo*ft(2,2)
      T33 = rTwo*ft(3,3)
      T44 = rTwo*ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      D1 =   k(1)*(T11-T14) - k(2)*(T12-T24)
     &     - k(3)*(T13-T34) - k(4)*(T14-T44)

      D2 = - k(1)*(T12-ci*T13) + k(2)*(T22-ci*T23)
     &     + k(3)*(T23-ci*T33) + k(4)*(T24-ci*T34)

      D3 = - k(1)*(T12+ci*T13) + k(2)*(T22+ci*T23)
     &     + k(3)*(T23+ci*T33) + k(4)*(T24+ci*T34)

      D4 =   k(1)*(T11+T14) - k(2)*(T12+T24)
     &     - k(3)*(T13+T34) - k(4)*(T14+T44)

      Tii = T11 - T22 - T33 - T44

      vertex = D1*(f13+f42) + D2*(f14-f32) + D3*(f23-f41) + D4*(f24+f31)

      vertex = vertex + Tii*( - k14p*(f24+f31) - k14m*(f13+f42)
     &                        +  k23*(f23-f41) + k23s*(f14-f32) )

      if ( fmass.ne.rZero ) then
         vertex = vertex + m2*Tii*(  fo(1)*fi(1) + fo(2)*fi(2)
     &                             + fo(3)*fi(3) + fo(4)*fi(4) )
      end if

      vertex = vertex * gt

      return
      end
      subroutine iovdmx(fi,fo,vc,gc, vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c dipole moment (non-renormalizable) coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2,2)        : coupling constants                  gvf
c                              : first index is L,R as normal
c                              : second index is EDM,-MDM
c
c output:
c       complex vertex         : amplitude                     <fo|v|fi>
c
      implicit none
      double complex fi(6), fo(6), vc(6), vertex, gc(2,2)

      double complex q(5:6), dum1, dum2
      double complex f1122, f12, f21, f3344, f34, f43
      double complex f12p21, f12m21, f34p43, f34m43
      double complex kvc21, kvc31, kvc41, kvc32, kvc42, kvc43
      double precision  rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      q(5) = fi(5) - fo(5)
      q(6) = fi(6) - fo(6)

      f1122  = fo(1)*fi(1) - fo(2)*fi(2)
      f12    = fo(1)*fi(2)
      f21    = fo(2)*fi(1)
      f12p21 = f12 + f21
      f12m21 = f12 - f21

      kvc21 = ( dble(q(6))*vc(1) -  dble(q(5))*vc(2))*cImag
      kvc31 =  dimag(q(6))*vc(1) -  dble(q(5))*vc(3)
      kvc41 = (dimag(q(5))*vc(1) -  dble(q(5))*vc(4))*cImag
      kvc32 =  dimag(q(6))*vc(2) -  dble(q(6))*vc(3)
      kvc42 = (dimag(q(5))*vc(2) -  dble(q(6))*vc(4))*cImag
      kvc43 =  dimag(q(5))*vc(3) - dimag(q(6))*vc(4)

      dum1 =   ( kvc31 + kvc42 )*f12m21
     &       + ( kvc32 + kvc41 )*f1122
     &       + ( kvc43 + kvc21 )*f12p21

c     (-) from gamma^5 in EDM only
      vertex = ( -gc(1,1) + cImag*gc(1,2) )*dum1    

      if ( gc(2,1).ne.cZero .or.
     &     gc(2,2).ne.cZero      ) then
         f3344  = fo(3)*fi(3) - fo(4)*fi(4)
         f34    = fo(3)*fi(4)
         f43    = fo(4)*fi(3)
         f34p43 = f34 + f43
         f34m43 = f34 - f43
         dum2 =   (-kvc31 + kvc42 )*f34m43
     &          + ( kvc32 - kvc41 )*f3344
     &          + ( kvc43 - kvc21 )*f34p43
         vertex = vertex + ( gc(2,1) + cImag*gc(2,2) )*dum2
      end if
c
      return
      end
      subroutine iovkxx(fi,fo,tc,g, vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex tc(6,4)        : input    tensor                      t
c       complex g(1)           : coupling constant                 -kappa/8
c       real    g(2)           : fermion mass                        m_f
c
c output:
c       complex vertex         : amplitude                        <fo|t|fi>
c     
      implicit none
      double complex fi(18), fo(18), tc(18), vertex,g(2)
      double precision fmass

      double complex k23, k23s, D1, D2, D3, D4, Tii
      double complex T11, T22, T33, T44, T12, T13, T14, T23, T24, T34
      double complex f13, f14, f23, f24, f31, f32, f41, f42
      double precision k(4), k14p, k14m, m2

      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
c
      fmass = dreal(g(2))
      m2 = rTwo*fmass

      k(1) = dreal(fi(5)+fo(5))
      k(2) = dreal(fi(6)+fo(6))
      k(3) = dimag(fi(6)+fo(6))
      k(4) = dimag(fi(5)+fo(5))
      k23  = dcmplx( k(2),k(3) )
      k23s = dconjg( k23 )
      k14p = k(1) + k(4)
      k14m = k(1) - k(4)

c      write (*,*) k

      f13 = fo(1)*fi(3)
      f14 = fo(1)*fi(4)
      f23 = fo(2)*fi(3)
      f24 = fo(2)*fi(4)
      f31 = fo(3)*fi(1)
      f32 = fo(3)*fi(2)
      f41 = fo(4)*fi(1)
      f42 = fo(4)*fi(2)

      T11 = rTwo*tc(1)
      T22 = rTwo*tc(6)
      T33 = rTwo*tc(11)
      T44 = rTwo*tc(16)
      T12 = tc( 5) + tc( 2)
      T13 = tc( 9) + tc( 3)
      T14 = tc(13) + tc( 4)
      T23 = tc(10) + tc( 7)
      T24 = tc(14) + tc( 8)
      T34 = tc(15) + tc(12)

      D1 =   k(1)*(T11-T14) - k(2)*(T12-T24)
     &     - k(3)*(T13-T34) - k(4)*(T14-T44)

      D2 = - k(1)*(T12-ci*T13) + k(2)*(T22-ci*T23)
     &     + k(3)*(T23-ci*T33) + k(4)*(T24-ci*T34)

      D3 = - k(1)*(T12+ci*T13) + k(2)*(T22+ci*T23)
     &     + k(3)*(T23+ci*T33) + k(4)*(T24+ci*T34)

      D4 =   k(1)*(T11+T14) - k(2)*(T12+T24)
     &     - k(3)*(T13+T34) - k(4)*(T14+T44)

      Tii = T11 - T22 - T33 - T44

      vertex = D1*(f13+f42) + D2*(f14-f32) + D3*(f23-f41) + D4*(f24+f31)

      vertex = vertex + Tii*( - k14p*(f24+f31) - k14m*(f13+f42)
     &                        +  k23*(f23-f41) + k23s*(f14-f32) )

      if ( fmass.ne.rZero ) then
         vertex = vertex + m2*Tii*(  fo(1)*fi(1) + fo(2)*fi(2)
     &                             + fo(3)*fi(3) + fo(4)*fi(4) )
      end if

      vertex = vertex * g(1) 
c
      return
      end
      subroutine iovtkk(fi,fo,vc,tc,g , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two gauge bosons and a Kaluza-Klein tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion     SM |fi>
c       complex fo(6)          : flow-out fermion     SM <fo|
c       complex vc(6)          : vector               SM   v
c       complex tc(6,4)        : tensor               KK   t
c       real    g(2)           : coupling constant    -g(L,R)*kappa/4
c
c output:
c       complex vertex         : amplitude            gamma(fi,fo,vc,tc)
c     
      implicit none
      double complex fi(6), fo(6), vc(6), tc(6,4), vertex
      double precision g(2)

      double complex f13, f14, f23, f24, f31, f32, f41, f42
      double complex fs1L, fs1R, fs2L, fs2R, fs3L, fs3R, fs4L, fs4R
      double complex T12, T13, T14, T23, T24, T34

      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
c
      f31 = fo(3)*fi(1)*g(1)
      f32 = fo(3)*fi(2)*g(1)
      f41 = fo(4)*fi(1)*g(1)
      f42 = fo(4)*fi(2)*g(1)

      fs1L =  f31 + f42
      fs2L = -f32 - f41
      fs3L = (f32 - f41)*ci
      fs4L = -f31 + f42

      if ( g(2).ne.rZero ) then
         f14 = fo(1)*fi(4)*g(2)
         f13 = fo(1)*fi(3)*g(2)
         f23 = fo(2)*fi(3)*g(2)
         f24 = fo(2)*fi(4)*g(2)
         fs1R =  f13 + f24
         fs2R =  f23 + f14
         fs3R = (f23 - f14)*ci
         fs4R =  f13 - f24
      end if

      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      if ( g(2).ne.rZero ) then

         vertex =  (fs1L + fs1R)*(  vc(1)*rTwo*tc(1,1)
     &                            - vc(2)*T12 - vc(3)*T13 - vc(4)*T14 )

     &           + (fs2L + fs2R)*(  vc(2)*rTwo*tc(2,2)
     &                            - vc(1)*T12 + vc(3)*T23 + vc(4)*T24 )

     &           + (fs3L + fs3R)*(  vc(3)*rTwo*tc(3,3)
     &                            - vc(1)*T13 + vc(2)*T23 + vc(4)*T34 )

     &           + (fs4L + fs4R)*(  vc(4)*rTwo*tc(4,4)
     &                            - vc(1)*T14 + vc(2)*T24 + vc(3)*T34 )

         vertex = vertex - rTwo*( tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4) )
     &                         *(  (vc(1)+      vc(4))*(f31+f24)
     &                           + (vc(1)-      vc(4))*(f13+f42)
     &                           + (vc(2)+ci*vc(3))*(f41-f23)
     &                           + (vc(2)-ci*vc(3))*(f32-f14) )

      else

         vertex =  fs1L*(  vc(1)*rTwo*tc(1,1)
     &                   - vc(2)*T12 - vc(3)*T13 - vc(4)*T14 )

     &           + fs2L*(  vc(2)*rTwo*tc(2,2)
     &                   - vc(1)*T12 + vc(3)*T23 + vc(4)*T24 )

     &           + fs3L*(  vc(3)*rTwo*tc(3,3)
     &                   - vc(1)*T13 + vc(2)*T23 + vc(4)*T34 )

     &           + fs4L*(  vc(4)*rTwo*tc(4,4)
     &                   - vc(1)*T14 + vc(2)*T24 + vc(3)*T34 )

         vertex = vertex - rTwo*( tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4) )
     &                         *(  (vc(1)+      vc(4))*f31
     &                           + (vc(1)-      vc(4))*f42
     &                           + (vc(2)+ci*vc(3))*f41
     &                           + (vc(2)-ci*vc(3))*f32 )

      end if
c
      return
      end
      subroutine iovtxx(fi,fo,vc,tc,gc,gt , vertex)
c
c This subroutine computes an amplitude of the four-point coupling of
c a vector boson, two fermions and a tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                       v
c       complex tc(18)         : input    tensor                       T
c       complex gc(2)          : coupling constants                  gvf
c       complex gt             : coupling constant      gtfv=-1/Lambda/2
c
c output:
c       complex vertex         : amplitude                   <fo|v,T|fi>
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fi(6), fo(6), vc(6), tc(18), gc(2), gt, vertex

      double complex ft(6,4)
      double complex T00,T12, T13, T14, T23, T24, T34

      double precision rZero, r2
      parameter( rZero = 0.0d0, r2 = 2.0d0 )
      double complex cone
      parameter( cone = ( 0.0d0, 1.0d0 ) )

      
      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)


      vertex =-gt*(fi(4)*gc(2)*(fo(1)*(T12*vc(1) - cone*T13*vc(1)
     &	 - r2*ft(2,2)*vc(2) + cone*r2*ft(3,3)*vc(3) - 
     &          T23*(-(cone*vc(2)) + vc(3)) + r2*T00*(-vc(2) 
     &+ cone*vc(3)) - T24*vc(4) + cone*T34*vc(4)) + 
     &       fo(2)*(-(r2*ft(1,1)*vc(1)) + T12*vc(2) + T24*vc(2)
     & + T13*vc(3) + T34*vc(3) + r2*ft(4,4)*vc(4) + 
     &          T14*(-vc(1) + vc(4)) + r2*T00*(vc(1) + vc(4)))) + 
     &    fi(1)*gc(1)*(fo(4)*(-(T12*vc(1)) - cone*T13*vc(1) 
     &+ r2*ft(2,2)*vc(2) - T23*(-(cone*vc(2)) - vc(3)) + 
     &          cone*r2*ft(3,3)*vc(3) + r2*T00*(vc(2) 
     &+ cone*vc(3)) + T24*vc(4) + cone*T34*vc(4)) + 
     &       fo(3)*(-(r2*ft(1,1)*vc(1)) + T12*vc(2) 
     &+ T24*vc(2) + T13*vc(3) + T34*vc(3) + r2*ft(4,4)*vc(4) + 
     &          T14*(-vc(1) + vc(4)) + r2*T00*(vc(1) + vc(4)))) + 
     &    fi(3)*gc(2)*(fo(2)*(T12*vc(1) + cone*T13*vc(1)
     & - r2*ft(2,2)*vc(2) - cone*r2*ft(3,3)*vc(3)
     & - T23*(cone*vc(2) + vc(3)) + 
     &          r2*T00*(-vc(2) - cone*vc(3)) - T24*vc(4) 
     &- cone*T34*vc(4)) + 
     &       fo(1)*(-(r2*ft(1,1)*vc(1)) + T12*vc(2) 
     &- T24*vc(2) + T13*vc(3) - T34*vc(3) + r2*T00*(vc(1) - vc(4)) - 
     &          r2*ft(4,4)*vc(4) + T14*(vc(1) + vc(4)))) + 
     &    fi(2)*gc(1)*(fo(3)*(-(T12*vc(1)) + cone*T13*vc(1)
     & + r2*ft(2,2)*vc(2) - T23*(cone*vc(2) - vc(3)) - 
     &          cone*r2*ft(3,3)*vc(3) + r2*T00*(vc(2) 
     &- cone*vc(3)) + T24*vc(4) - cone*T34*vc(4)) + 
     &       fo(4)*(-(r2*ft(1,1)*vc(1)) + T12*vc(2) 
     &- T24*vc(2) + T13*vc(3) - T34*vc(3) + r2*T00*(vc(1) - vc(4)) - 
     &          r2*ft(4,4)*vc(4) + T14*(vc(1) + vc(4)))))

      return
      end
      subroutine iovxxx(fi,fo,vc,gc , vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c
c output:
c       complex vertex         : amplitude                     <fo|v|fi>
c     
      implicit none
      double complex fi(6),fo(6),gc(2),vc(6),vertex

      double precision rZero, rOne
      parameter( rZero = 0.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      vertex =  gc(1)*( (fo(3)*fi(1)+fo(4)*fi(2))*vc(1)
     &                 +(fo(3)*fi(2)+fo(4)*fi(1))*vc(2)
     &                 -(fo(3)*fi(2)-fo(4)*fi(1))*vc(3)*cImag
     &                 +(fo(3)*fi(1)-fo(4)*fi(2))*vc(4)        )

      if ( gc(2).ne.cZero ) then
         vertex = vertex
     &          + gc(2)*( (fo(1)*fi(3)+fo(2)*fi(4))*vc(1)
     &                   -(fo(1)*fi(4)+fo(2)*fi(3))*vc(2)
     &                   +(fo(1)*fi(4)-fo(2)*fi(3))*vc(3)*cImag
     &                   -(fo(1)*fi(3)-fo(2)*fi(4))*vc(4)        )
      end if
c
      return
      end
      subroutine ixxxxx(p,fmass,nhel,nsf , fi)
c
c This subroutine computes a fermion wavefunction with the flowing-IN
c fermion number.
c
c input:
c       real    p(0:3)         : four-momentum of fermion
c       real    fmass          : mass          of fermion
c       integer nhel = -1 or 1 : helicity      of fermion
c       integer nsf  = -1 or 1 : +1 for particle, -1 for anti-particle
c
c output:
c       complex fi(6)          : fermion wavefunction               |fi>
c     
      implicit none
      double complex fi(6),chi(2)
      double precision p(0:3),sf(2),sfomeg(2),omega(2),fmass,
     &     pp,pp3,sqp0p3,sqm(0:1)
      integer nhel,nsf,ip,im,nh

      double precision rZero, rHalf, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0, rTwo = 2.0d0 )
      
c

      fi(5) = dcmplx(p(0),p(3))*nsf
      fi(6) = dcmplx(p(1),p(2))*nsf

      nh = nhel*nsf

      if ( fmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))
         
         if ( pp.eq.rZero ) then
            
            sqm(0) = dsqrt(abs(fmass)) ! possibility of negative fermion masses
            sqm(1) = sign(sqm(0),fmass) ! possibility of negative fermion masses
            ip = (1+nh)/2
            im = (1-nh)/2
            
            fi(1) = ip     * sqm(ip)
            fi(2) = im*nsf * sqm(ip)
            fi(3) = ip*nsf * sqm(im)
            fi(4) = im     * sqm(im)

         else

            sf(1) = dble(1+nsf+(1-nsf)*nh)*rHalf
            sf(2) = dble(1+nsf-(1-nsf)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = fmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , p(2) )/dsqrt(rTwo*pp*pp3)
            endif
            
            fi(1) = sfomeg(1)*chi(im)
            fi(2) = sfomeg(1)*chi(ip)
            fi(3) = sfomeg(2)*chi(im)
            fi(4) = sfomeg(2)*chi(ip)
            
         endif
         
      else
         
         if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
            sqp0p3 = 0d0
         else
            sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsf
         end if
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fi(1) = dcmplx( rZero )
            fi(2) = dcmplx( rZero )
            fi(3) = chi(1)
            fi(4) = chi(2)
         else
            fi(1) = chi(2)
            fi(2) = chi(1)
            fi(3) = dcmplx( rZero )
            fi(4) = dcmplx( rZero )
         endif
      endif
c
      return
      end
      subroutine j3xxxx(fi,fo,gaf,gzf,zmass,zwidth , j3)
c
c This subroutine computes the sum of photon and Z currents with the
c suitable weights ( j(W3) = cos(theta_W) j(Z) + sin(theta_W) j(A) ).
c The output j3 is useful as an input of vvvxxx, jvvxxx or w3w3xx.
c The photon propagator is given in Feynman gauge, and the Z propagator
c is given in unitary gauge.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gaf(2)         : fi couplings with A                 gaf
c       complex gzf(2)         : fi couplings with Z                 gzf
c       real    zmass          : mass  of Z
c       real    zwidth         : width of Z
c
c output:
c       complex j3(6)          : W3 current             j^mu(<fo|w3|fi>)
c     
      implicit none
      double complex fi(6),fo(6),j3(6),gaf(2),gzf(2)
      double complex c0l,c1l,c2l,c3l,csl,c0r,c1r,c2r,c3r,csr,dz,ddif
      double complex gn,gz3l,ga3l
      double complex cm2  ! mass**2- I Gamma mass (Fabio)
      double precision q(0:3),zmass,zwidth,zm2,zmw
      double precision q2,da,ww,cw,sw

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      j3(5) = fo(5)-fi(5)
      j3(6) = fo(6)-fi(6)

      q(0) = -dble( j3(5))
      q(1) = -dble( j3(6))
      q(2) = -dimag(j3(6))
      q(3) = -dimag(j3(5))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      zm2 = zmass**2
      zmw = zmass*zwidth


      da = rOne/q2
C      ww = max(dsign(zmw,q2), rZero)
      dz = rOne/dcmplx( q2-zm2, zmw )
      ddif = dcmplx( -zm2, zmw )*da*dz

c ddif is the difference : ddif=da-dz
c  For the running width, use below instead of the above ww,dz and ddif.
c      ww = max( zwidth*q2/zmass, rZero )
c      dz = rOne/dcmplx( q2-zm2, zmw )
c      ddif = dcmplx( -zm2, zmw )*da*dz



      cw = rOne/sqrt(rOne+(gzf(2)/gaf(2))**2)
      sw = sqrt((rOne-cw)*(rOne+cw))
      gn = gaf(2)*sw
      gz3l = gzf(1)*cw
      ga3l = gaf(1)*sw
      c0l =   fo(3)*fi(1)+fo(4)*fi(2)
      c0r =   fo(1)*fi(3)+fo(2)*fi(4)
      c1l = -(fo(3)*fi(2)+fo(4)*fi(1))
      c1r =   fo(1)*fi(4)+fo(2)*fi(3)
      c2l =  (fo(3)*fi(2)-fo(4)*fi(1))*cImag
      c2r = (-fo(1)*fi(4)+fo(2)*fi(3))*cImag
      c3l =  -fo(3)*fi(1)+fo(4)*fi(2)
      c3r =   fo(1)*fi(3)-fo(2)*fi(4)

c     Fabio's implementation of the fixed width
      cm2=dcmplx( zm2, -zmw )
c     csl = (q(0)*c0l-q(1)*c1l-q(2)*c2l-q(3)*c3l)/zm2
c     csr = (q(0)*c0r-q(1)*c1r-q(2)*c2r-q(3)*c3r)/zm2
      csl = (q(0)*c0l-q(1)*c1l-q(2)*c2l-q(3)*c3l)/cm2
      csr = (q(0)*c0r-q(1)*c1r-q(2)*c2r-q(3)*c3r)/cm2
      
      j3(1) =  gz3l*dz*(c0l-csl*q(0))+ga3l*c0l*da
     &       + gn*(c0r*ddif+csr*q(0)*dz)
      j3(2) =  gz3l*dz*(c1l-csl*q(1))+ga3l*c1l*da
     &       + gn*(c1r*ddif+csr*q(1)*dz)
      j3(3) =  gz3l*dz*(c2l-csl*q(2))+ga3l*c2l*da
     &       + gn*(c2r*ddif+csr*q(2)*dz)
      j3(4) =  gz3l*dz*(c3l-csl*q(3))+ga3l*c3l*da
     &       + gn*(c3r*ddif+csr*q(3)*dz)
c
      return
      end
      subroutine jeexxx(eb,ef,shlf,chlf,phi,nhb,nhf,nsf , jee)
c
c This subroutine computes an off-shell photon wavefunction emitted from
c the electron or positron beam, with a special care for the small angle
c region.  The momenta are measured in the laboratory frame, where the
c e- (e+) beam is along the positive (negative) z axis.
c
c input:
c       real    eb             : energy (gev)    of beam  e-/e+
c       real    ef             : energy (gev)    of final e-/e+
c       real    shlf           : sin(theta/2)    of final e-/e+
c       real    chlf           : cos(theta/2)    of final e-/e+
c       real    phi            : azimuthal angle of final e-/e+
c       integer nhb  = -1 or 1 : helicity        of beam  e-/e+
c       integer nhf  = -1 or 1 : helicity        of final e-/e+
c       integer nsf  = -1 or 1 : +1 for electron, -1 for positron
c
c output:
c       complex jee(6)         : off-shell photon          j^mu(<e|a|e>)
c     
      implicit none
      double complex jee(6),coeff
      double precision cs(2),eb,ef,shlf,chlf,phi,alpha,gal,hi,sf,sfh
      double precision x,me2,q2,rfp,rfm,snp,csp,rxc,c,s
      integer nhb,nhf,nsf

      double precision rZero, rHalf, rOne, rTwo, rFour, rOte
      double precision rPi, rIalph
      parameter( rZero = 0.0d0, rHalf = 0.5d0, rOne = 1.0d0 )
      parameter( rTwo = 2.0d0, rFour = 4.0d0, rOte = 128.9d0 )
      parameter( rPi = 3.14159265358979323846d0 )
      parameter( rIalph = 137.0359895d0 )

      double precision me
      parameter( me = 0.51099906d-3 )

c

      alpha = rOne/rOte
      gal = sqrt(alpha*rFour*rPi)

      hi = nhb
      sf = nsf
      sfh = nhb*nsf
      cs((3+nsf)/2) = shlf
      cs((3-nsf)/2) = chlf
c cs(1)=chlf and cs(2)=shlf for electron
c cs(1)=shlf and cs(2)=chlf for positron
      x = ef/eb
      me2 = me**2
      q2 = - rFour*cs(2)**2*(ef*eb-me2)
     &     + sf*(rOne-x)**2/x*(shlf+chlf)*(shlf-chlf)*me2
      rfp = (1+nsf)
      rfm = (1-nsf)
      snp = sin(phi)
      csp = cos(phi)

      if ( nhb.eq.nhf ) then
         rxc = rTwo*x/(rOne-x)*cs(1)**2
         coeff = gal*rTwo*eb*sqrt(x)*cs(2)/q2
     &          *(dcmplx( rfp )-rfm*dcmplx( csp, -snp*hi ))*rHalf
         jee(1) = dcmplx( rZero )
         jee(2) = coeff*dcmplx( (rOne+rxc)*csp, -sfh*snp )
         jee(3) = coeff*dcmplx( (rOne+rxc)*snp,  sfh*csp )
         jee(4) = coeff*(-sf*rxc/cs(1)*cs(2))
      else
         coeff = gal*me/q2/sqrt(x)
     &          *(dcmplx( rfp )+rfm*dcmplx( csp, snp*hi ))*rHalf*hi
         jee(1) = -coeff*(rOne+x)*cs(2)*dcmplx( csp , sfh*snp )
         jee(2) =  coeff*(rOne-x)*cs(1)
         jee(3) =  jee(2)*dcmplx( rZero, sfh )
         jee(4) =  jee(1)*sf*(rOne-x)/(rOne+x)
      endif

      c = (chlf+shlf)*(chlf-shlf)
      s = rTwo*chlf*shlf

      jee(5) = -eb*dcmplx( rOne-x, sf-x*c )
      jee(6) =  eb*x*s*dcmplx( csp, snp )
c
      return
      end
      subroutine jgggtx(va,vb,vc,tc,gc,gt , jgggt)
c
c This subroutine computes an off-shell vector boson current from
c the 4 gluons-tensor boson coupling, corresponding 
c to the color structure f^{x,a,e} f{b,c,e}. 
c
c To obtain the complete amplitude, this subroutine must be called three
c times (once for each color structure) with the following permutations:
c     call jgggtx(va,vb,vc,tc,gc,gt , jgggt1)
c     call jgggtx(vb,vc,va,tc,gc,gt , jgggt2)
c     call jgggtx(vc,va,vb,tc,gc,gt , jgggt3)
c corresponding to 
c	f^{x,a,e} f^{b,c,e}
c	f^{x,b,e} f^{c,a,e}
c	f^{x,c,e} f^{a,b,e}
c
c input:
c       complex va(6)          : boson with adjoint color index a     va
c       complex vb(6)          : boson with adjoint color index b     vb
c       complex vc(6)          : boson with adjoint color index c     vc
c       complex tc(18)         : input tensor                          T
c       real    gc             : coupling constant                    gs
c       complex gt             : coupling constant         gtv=-1/Lambda
c
c output:
c       complex jgggt(6)       : gluon current        j^mu(v:va,vb,vc,T)
c
c- by Q.Li - JAN. 2008
c     
      implicit none
      double complex va(6), vb(6), vc(6), tc(18), gt,  jgggt(6)
      double precision gc

      double complex ft(6,4)
      double complex TVBVC, TVBVD, TVDM(4), TVCM(4)
      double complex T00,T12,T13,T14,T23,T24,T34,EBEC,EBED

      integer a,b,i,j

      double precision pb(4),pc(4),pd(4),pT(4),px(4)

      double precision MET(4,4)
      double complex cZero, d
      double precision rZero, rTwo,px2,r3
      parameter( rZero = 0.0d0, rTwo = 2.0d0, r3=3.0d0)
      parameter(cZero=(0.0d0,0.0d0))


      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      jgggt(5)=va(5)+vb(5)+vc(5)+ft(5,1)
      jgggt(6)=ft(6,1)+va(6)+vb(6)+vc(6)

      px(1) = dreal(jgggt(5))
      px(2) = dreal(jgggt(6))
      px(3) = dimag(jgggt(6))
      px(4) = dimag(jgggt(5))

      pb(1) = dreal(va(5))
      pb(2) = dreal(va(6))
      pb(3) = dimag(va(6))
      pb(4) = dimag(va(5))

      pc(1) = dreal(vb(5))
      pc(2) = dreal(vb(6))
      pc(3) = dimag(vb(6))
      pc(4) = dimag(vb(5))

      pd(1) = dreal(vc(5))
      pd(2) = dreal(vc(6))
      pd(3) = dimag(vc(6))
      pd(4) = dimag(vc(5))

     
      do i=1,4
         do j=1,4
            MET(i,j)=0.0d0
         enddo 
      enddo

      MET(1,1)=1.0d0
      MET(2,2)=-1.0d0
      MET(3,3)=-1.0d0
      MET(4,4)=-1.0d0

      px2 = px(1)**2-px(2)**2-px(3)**2-px(4)**2
      EBEC = va(1)*vb(1)-va(2)*vb(2)-va(3)*vb(3)-va(4)*vb(4)
      EBED = va(1)*vc(1)-va(2)*vc(2)-va(3)*vc(3)-va(4)*vc(4)
      
      d = 1.0d0/dcmplx( px2, 0.0d0 )
      
      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      TVBVC =2.0d0*(ft(1,1)*va(1)*vb(1)
     &+ft(2,2)*va(2)*vb(2)+ft(3,3)*va(3)*vb(3)
     &+ft(4,4)*va(4)*vb(4)
     &)
     & - T12*(va(1)*vb(2) + va(2)*vb(1))
     &          - T13*(va(1)*vb(3) + va(3)*vb(1))
     &          - T14*(va(1)*vb(4) + va(4)*vb(1))
     &          + T23*(va(2)*vb(3) + va(3)*vb(2))
     &          + T24*(va(2)*vb(4) + va(4)*vb(2))
     &          + T34*(va(3)*vb(4) + va(4)*vb(3))

      TVBVD =2.0d0*(ft(1,1)*va(1)*vc(1)
     &+ft(2,2)*va(2)*vc(2)+ft(3,3)*va(3)*vc(3)
     &+ft(4,4)*va(4)*vc(4)
     &)
     & - T12*(va(1)*vc(2) + va(2)*vc(1))
     &          - T13*(va(1)*vc(3) + va(3)*vc(1))
     &          - T14*(va(1)*vc(4) + va(4)*vc(1))
     &          + T23*(va(2)*vc(3) + va(3)*vc(2))
     &          + T24*(va(2)*vc(4) + va(4)*vc(2))
     &          + T34*(va(3)*vc(4) + va(4)*vc(3))

      do j=1,4
         TVDM(j)=
     &MET(j,1)*(ft(1,1)*vc(1)-ft(2,1)*vc(2)
     &-ft(3,1)*vc(3)-ft(4,1)*vc(4))
     &-MET(j,2)*(ft(1,2)*vc(1)-ft(2,2)*vc(2)
     &-ft(3,2)*vc(3)-ft(4,2)*vc(4))
     &-MET(j,3)*(ft(1,3)*vc(1)-ft(2,3)*vc(2)
     &-ft(3,3)*vc(3)-ft(4,3)*vc(4))
     &-MET(j,4)*(ft(1,4)*vc(1)-ft(2,4)*vc(2)
     &-ft(3,4)*vc(3)-ft(4,4)*vc(4))
     &+
     &MET(j,1)*(ft(1,1)*vc(1)-ft(1,2)*vc(2)
     &-ft(1,3)*vc(3)-ft(1,4)*vc(4))
     &-MET(j,2)*(ft(2,1)*vc(1)-ft(2,2)*vc(2)
     &-ft(2,3)*vc(3)-ft(2,4)*vc(4))
     &-MET(j,3)*(ft(3,1)*vc(1)-ft(3,2)*vc(2)
     &-ft(3,3)*vc(3)-ft(3,4)*vc(4))
     &-MET(j,4)*(ft(4,1)*vc(1)-ft(4,2)*vc(2)
     &-ft(4,3)*vc(3)-ft(4,4)*vc(4))
	
	TVCM(j)=
     &MET(j,1)*(ft(1,1)*vb(1)-ft(2,1)*vb(2)
     &-ft(3,1)*vb(3)-ft(4,1)*vb(4))
     &-MET(j,2)*(ft(1,2)*vb(1)-ft(2,2)*vb(2)
     &-ft(3,2)*vb(3)-ft(4,2)*vb(4))
     &-MET(j,3)*(ft(1,3)*vb(1)-ft(2,3)*vb(2)
     &-ft(3,3)*vb(3)-ft(4,3)*vb(4))
     &-MET(j,4)*(ft(1,4)*vb(1)-ft(2,4)*vb(2)
     &-ft(3,4)*vb(3)-ft(4,4)*vb(4))
     &+
     &MET(j,1)*(ft(1,1)*vb(1)-ft(1,2)*vb(2)
     &-ft(1,3)*vb(3)-ft(1,4)*vb(4))
     &-MET(j,2)*(ft(2,1)*vb(1)-ft(2,2)*vb(2)
     &-ft(2,3)*vb(3)-ft(2,4)*vb(4))
     &-MET(j,3)*(ft(3,1)*vb(1)-ft(3,2)*vb(2)
     &-ft(3,3)*vb(3)-ft(3,4)*vb(4))
     &-MET(j,4)*(ft(4,1)*vb(1)-ft(4,2)*vb(2)
     &-ft(4,3)*vb(3)-ft(4,4)*vb(4))

      enddo
	    
      do a=1,4
	
         jgggt(a)=vc(a)*TVBVC-vb(a)*TVBVD
     &        -T00*vc(a)*EBEC+T00*vb(a)*EBED
     &        +EBEC*TVDM(a)-EBED*TVCM(a)
	 
         jgggt(a)=-jgggt(a)*d*gc*gc*gt
         
      enddo 
 
      return
      end
      subroutine jgggxx(w1,w2,w3,g, jw3w)
c
c This subroutine computes an off-shell W+, W-, W3, Z or photon current
c from the four-point gauge boson coupling, including the contributions
c of W exchange diagrams.  The vector propagator is given in Feynman
c gauge for a photon and in unitary gauge for W and Z bosons.  If one
c sets wmass=0.0, then the ggg-->g current is given (see sect 2.9.1 of
c the manual).
c
c input:
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       complex w3(6)          : third  vector                        w3
c       real    g             : first  coupling constant
c                                                  (see the table below)
c
c output:
c       complex jw3w(6)        : W current             j^mu(w':w1,w2,w3)
c
      implicit none
      double complex w1(6),w2(6),w3(6),jw3w(6)
      double complex dw1(0:3),dw2(0:3),dw3(0:3)
      double complex jj(0:3),dv,w32,w13
      double precision p1(0:3),p2(0:3),p3(0:3),q(0:3),g,dg2,q2

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )

c

      jw3w(5) = w1(5)+w2(5)+w3(5)
      jw3w(6) = w1(6)+w2(6)+w3(6)


      dw1(0) = dcmplx(w1(1))
      dw1(1) = dcmplx(w1(2))
      dw1(2) = dcmplx(w1(3))
      dw1(3) = dcmplx(w1(4))
      dw2(0) = dcmplx(w2(1))
      dw2(1) = dcmplx(w2(2))
      dw2(2) = dcmplx(w2(3))
      dw2(3) = dcmplx(w2(4))
      dw3(0) = dcmplx(w3(1))
      dw3(1) = dcmplx(w3(2))
      dw3(2) = dcmplx(w3(3))
      dw3(3) = dcmplx(w3(4))
      p1(0) = dble(      w1(5))
      p1(1) = dble(      w1(6))
      p1(2) = dble(dimag(w1(6)))
      p1(3) = dble(dimag(w1(5)))
      p2(0) = dble(      w2(5))
      p2(1) = dble(      w2(6))
      p2(2) = dble(dimag(w2(6)))
      p2(3) = dble(dimag(w2(5)))
      p3(0) = dble(      w3(5))
      p3(1) = dble(      w3(6))
      p3(2) = dble(dimag(w3(6)))
      p3(3) = dble(dimag(w3(5)))
      q(0) = -(p1(0)+p2(0)+p3(0))
      q(1) = -(p1(1)+p2(1)+p3(1))
      q(2) = -(p1(2)+p2(2)+p3(2))
      q(3) = -(p1(3)+p2(3)+p3(3))

      q2 = q(0)**2 -(q(1)**2 +q(2)**2 +q(3)**2)

      dg2 = dble(g)*dble(g)

      dv = rOne/dcmplx( q2 )

      w32 = dw3(0)*dw2(0)-dw3(1)*dw2(1)-dw3(2)*dw2(2)-dw3(3)*dw2(3)

      w13 = dw1(0)*dw3(0)-dw1(1)*dw3(1)-dw1(2)*dw3(2)-dw1(3)*dw3(3)

      jj(0) = dg2*( dw1(0)*w32 - dw2(0)*w13 )
      jj(1) = dg2*( dw1(1)*w32 - dw2(1)*w13 )
      jj(2) = dg2*( dw1(2)*w32 - dw2(2)*w13 )
      jj(3) = dg2*( dw1(3)*w32 - dw2(3)*w13 )

      jw3w(1) = dcmplx( jj(0)*dv )
      jw3w(2) = dcmplx( jj(1)*dv )
      jw3w(3) = dcmplx( jj(2)*dv )
      jw3w(4) = dcmplx( jj(3)*dv )
c
      return
      end
      subroutine jggxxx(v1,v2,g, jvv)
c
c This subroutine computes an off-shell vector current from the three-
c point gauge boson coupling.  The vector propagator is given in Feynman
c gauge for a massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       real    g              : coupling constant (see the table below)
c
c output:
c       complex jvv(6)         : vector current            j^mu(v:v1,v2)
c
      implicit none
      double complex v1(6),v2(6),jvv(6),j12(0:3)
      double complex sv1,sv2,v12
      double precision p1(0:3),p2(0:3),q(0:3),g,gs,s

c

      jvv(5) = v1(5) + v2(5)
      jvv(6) = v1(6) + v2(6)

      p1(0) = dble( v1(5))
      p1(1) = dble( v1(6))
      p1(2) = dimag(v1(6))
      p1(3) = dimag(v1(5))
      p2(0) = dble( v2(5))
      p2(1) = dble( v2(6))
      p2(2) = dimag(v2(6))
      p2(3) = dimag(v2(5))
      q(0) = -dble( jvv(5))
      q(1) = -dble( jvv(6))
      q(2) = -dimag(jvv(6))
      q(3) = -dimag(jvv(5))
      s = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)

      v12 = v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4)
      sv1 =   (p2(0)-q(0))*v1(1) -(p2(1)-q(1))*v1(2)
     &      - (p2(2)-q(2))*v1(3) -(p2(3)-q(3))*v1(4)
      sv2 = - (p1(0)-q(0))*v2(1) +(p1(1)-q(1))*v2(2)
     &      + (p1(2)-q(2))*v2(3) +(p1(3)-q(3))*v2(4)
      j12(0) = (p1(0)-p2(0))*v12 +sv1*v2(1) +sv2*v1(1)
      j12(1) = (p1(1)-p2(1))*v12 +sv1*v2(2) +sv2*v1(2)
      j12(2) = (p1(2)-p2(2))*v12 +sv1*v2(3) +sv2*v1(3)
      j12(3) = (p1(3)-p2(3))*v12 +sv1*v2(4) +sv2*v1(4)

      gs = -g/s

      jvv(1) = gs*j12(0)
      jvv(2) = gs*j12(1)
      jvv(3) = gs*j12(2)
      jvv(4) = gs*j12(3)
c
      return
      end
      subroutine jiodmx(fi,fo,gc,vmass,vwidth , jio)
c
c This subroutine computes an off-shell vector dipole moment
c (non-renormalizable) current from an external
c fermion pair.  The vector boson propagator is given in Feynman gauge
c for a massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gc(2,2)        : coupling constants                 gvf
c                              : first index is L,R as normal
c                              : second index is EDM,-MDM
c       real    vmass          : mass  of output vector v
c       real    vwidth         : width of output vector v
c
c output:
c       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
c
      implicit none
      double complex fi(6), fo(6), jio(6), c0, c1, c2, c3, d
      double complex gc(2,2), gL, gR
      double precision  q(0:3), vmass, vwidth, q2, vm2, dd

      double complex f1122, f12, f21, f12p21, f12m21
      double complex f3344, f34, f43, f34p43, f34m43
      double complex dumL1, dumL2, dumL3, dumL4
      double complex dumR1, dumR2, dumR3, dumR4
      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      gL = -gc(1,1) + cImag*gc(1,2)
      gR =  gc(2,1) + cImag*gc(2,2)

      jio(5) = fo(5) - fi(5)
      jio(6) = fo(6) - fi(6)

      q(0) = dble( jio(5))
      q(1) = dble( jio(6))
      q(2) = dimag(jio(6))
      q(3) = dimag(jio(5))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2


      f1122  = fo(1)*fi(1) - fo(2)*fi(2)
      f12    = fo(1)*fi(2)
      f21    = fo(2)*fi(1)
      f12p21 = f12 + f21
      f12m21 = f12 - f21
      if ( gc(2,1).ne.cZero .or. gc(2,2).ne.cZero ) then
         f3344  = fo(3)*fi(3) - fo(4)*fi(4)
         f34    = fo(3)*fi(4)
         f43    = fo(4)*fi(3)
         f34p43 = f34 + f43
         f34m43 = f34 - f43
      end if

c note overall (-), since k in vertex is -q above
      dumL1 = -q(2)*f12m21 - cImag*( q(1)*f12p21 + q(3)*f1122 )
      dumL2 =  q(2)*f1122  - cImag*( q(0)*f12p21 - q(3)*f12m21 )
      dumL3 = -q(0)*f12m21 - q(1)*f1122 - q(3)*f12p21
      dumL4 = -q(2)*f12p21 - cImag*( q(0)*f1122  + q(1)*f12m21 )
      if ( gc(2,1).ne.cZero .or. gc(2,2).ne.cZero ) then
         dumR1 =  q(2)*f34m43 + cImag*( q(1)*f34p43 + q(3)*f3344 )
         dumR2 =  q(2)*f3344  + cImag*( q(0)*f34p43 + q(3)*f34m43 )
         dumR3 =  q(0)*f34m43 - q(1)*f3344 - q(3)*f34p43
         dumR4 = -q(2)*f34p43 + cImag*( q(0)*f3344  - q(1)*f34m43 )
      end if

      if ( vmass.ne.rZero ) then

         d = rOne/dcmplx( q2-vm2, vmass*vwidth )

         c0 = gL*dumL1
         c1 = gL*dumL2
         c2 = gL*dumL3
         c3 = gL*dumL4

         if ( gc(2,1).ne.cZero .or.
     &        gc(2,2).ne.cZero      ) then
            c0 = c0 + gR*dumR1
            c1 = c1 + gR*dumR2
            c2 = c2 + gR*dumR3
            c3 = c3 + gR*dumR4
         end if

         jio(1) = c0*d
         jio(2) = c1*d
         jio(3) = c2*d
         jio(4) = c3*d

      else

         dd = rOne/q2

         jio(1) = gL*dumL1
         jio(2) = gL*dumL2
         jio(3) = gL*dumL3
         jio(4) = gL*dumL4

         if ( gc(2,1).ne.cZero .or.
     &        gc(2,2).ne.cZero      ) then
            jio(1) = jio(1) + gR*dumR1
            jio(2) = jio(2) + gR*dumR2
            jio(3) = jio(3) + gR*dumR3
            jio(4) = jio(4) + gR*dumR4
         end if

         jio(1) = jio(1)*dd
         jio(2) = jio(2)*dd
         jio(3) = jio(3)*dd
         jio(4) = jio(4)*dd

      end if
c
      return
      end
      subroutine jiokxx(fi,fo,g,tmass,twidth , uio)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex g(1)           : coupling constant                 -kappa/8
c       real    g(2)           : fermion mass                        m_f
c
c output:
c       complex uio(18)        : outgoing tensor    T_{mu,vu} = T(mu+4*(vu-1))
c     
      implicit none
      double complex fi(18), fo(18), uio(18),g(2)
      double precision fmass,tmass,twidth

      double complex tc(4,4),tc1(4,4),fgamf(4),ffm
      double complex fkslaf,fgamfk(4,4),prop,denom
      double precision k(4),pp,m2,p(4),eta(4,4)
      integer i,j,ii,jj

      double complex ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )

      fmass = dreal(g(2))
      do i=1,16
         uio(i)=(0d0,0d0)
      enddo

      m2 = tmass**2
      uio(17) = fi(5)-fo(5)
      uio(18) = fi(6)-fo(6)

      p(1) = dble( uio(17))/tmass
      p(2) = dble( uio(18))/tmass
      p(3) = dimag(uio(18))/tmass
      p(4) = dimag(uio(17))/tmass

      k(1) = dble( fi(5)+fo(5))
      k(2) = dble( fi(6)+fo(6))
      k(3) = dimag(fi(6)+fo(6))
      k(4) = dimag(fi(5)+fo(5))

      pp = p(1)**2 - p(2)**2 - p(3)**2 - p(4)**2

      denom = dcmplx(pp*m2 - m2,tmass*twidth)

      fgamf(1) =     fi(3)*fo(1)+fi(4)*fo(2)+fi(1)*fo(3)+fi(2)*fo(4)
      fgamf(2) =     fi(4)*fo(1)+fi(3)*fo(2)-fi(2)*fo(3)-fi(1)*fo(4)
      fgamf(3) =ci*(-fi(4)*fo(1)+fi(3)*fo(2)+fi(2)*fo(3)-fi(1)*fo(4))
      fgamf(4) =     fi(3)*fo(1)-fi(4)*fo(2)-fi(1)*fo(3)+fi(2)*fo(4)
      ffm = 2d0*fmass*(fo(1)*fi(1)+fo(2)*fi(2)+fo(3)*fi(3)+fo(4)*fi(4))
      fkslaf  = k(1)* fgamf(1)-k(2)* fgamf(2)-
     &          k(3)* fgamf(3)-k(4)* fgamf(4)
      do i=1,3
         do j=1+1,4
            eta(i,j)=0d0
            eta(j,i)=eta(i,j)
         enddo
      enddo
      eta(1,1)= 1d0
      eta(2,2)=-1d0
      eta(3,3)=-1d0
      eta(4,4)=-1d0

c vertex
c upper triangle:
      do i=2,4
      do j=1,i-1
      tc(i,j)=fgamf(i)*k(j)+fgamf(j)*k(i)
c lower triangle:
      tc(j,i)=tc(i,j)
      enddo
c diagonal terms:
      tc(i,i)=2d0*(fgamf(i)*k(i)+(fkslaf-ffm))
      enddo
      tc(1,1)=2d0*(fgamf(1)*k(1)-(fkslaf-ffm))
      
c make indices upper indices before contracting with propagator
      do i=2,4
         tc(1,i)=-tc(1,i)
         tc(i,1)=-tc(i,1)
      enddo



c multiply by propagator
      do i=1,4
         do j=1,4
            do ii=1,4
               do jj=1,4
                  uio(i+4*(j-1))=uio(i+4*(j-1))+
     &            (
     &                    (eta(i,ii)-p(i)*p(ii))*(eta(j,jj)-p(j)*p(jj))+
     &                    (eta(i,jj)-p(i)*p(jj))*(eta(j,ii)-p(j)*p(ii))-
     &            2d0/3d0*(eta(i,j)-p(i)*p(j))*(eta(ii,jj)-p(ii)*p(jj))
     &            )
     &            *tc(ii,jj)*g(1)/denom
               enddo
            enddo
         enddo
      enddo

      return
      end
      subroutine jiotxx(fi,fo,tc,gc,gt,vmass,vwidth , jiot)
c
c This subroutine computes an off-shell vector boson wavefunction from a
c flowing-out fermion, a flowing-in fermion and a tensor boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex tc(18)         : input    tensor                       T
c       complex gc(2)          : coupling constants                  gvf
c       complex gt             : coupling constant      gtfv=-1/Lambda/2
c       real    vmass          : mass  of output vector v
c       real    vwidth         : width of output vector v
c
c output:
c       complex jiot(6)        : vector boson          j^mu(<fo|v,T|fi>)
c
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex fi(6), fo(6), tc(18), gc(2), gt, jiot(6)
      double precision vmass, vwidth

      double complex ft(6,4)
      double complex d, T00, T12, T13, T14, T23, T24, T34
      double precision pv(4), pv2
      integer i
      
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      double complex cone
      parameter( cone = ( 0.0d0, 1.0d0 ))


      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      jiot(5) = fo(5) + ft(5,1) -fi(5)
      jiot(6) = fo(6) + ft(6,1) -fi(6)

      pv(1) = dreal(jiot(5))
      pv(2) = dreal(jiot(6))
      pv(3) = dimag(jiot(6))
      pv(4) = dimag(jiot(5))
 
      pv2=pv(1)**2-pv(2)**2-pv(3)**2-pv(4)**2
     

      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)
      
      if ( vmass.gt.rZero ) then
         d =  1.0d0/dcmplx( pv2-vmass**2, vmass*vwidth )
      else
         d =  1.0d0/dcmplx( pv2, rZero )
      end if

      if (vmass.gt.rzero) then
         jiot(1) = fi(3)*gc(2)*(-((rtwo*T00*fo(2)*pv(1)*(-pv(2)
     &	 - cone*pv(3)))/vmass**2) + 
     &       fo(1)*((rtwo*T00) 
     &- (rtwo*T00*pv(1)*(pv(1) - pv(4)))/vmass**2)) + 
     &    fi(2)*gc(1)*(-((rtwo*T00*fo(3)*pv(1)
     &*(pv(2) - cone*pv(3)))/vmass**2) + 
     &       fo(4)*((rtwo*T00) - (rtwo*T00*pv(1)
     &*(pv(1) - pv(4)))/vmass**2)) + 
     &    fi(4)*gc(2)*(-((rtwo*T00*fo(1)*pv(1)
     &*(-pv(2) + cone*pv(3)))/vmass**2) + 
     &       fo(2)*((rtwo*T00) 
     &- (rtwo*T00*pv(1)*(pv(1) + pv(4)))/vmass**2)) + 
     &    fi(1)*gc(1)*(-((rtwo*T00*fo(4)*pv(1)*(pv(2)
     & + cone*pv(3)))/vmass**2) + 
     &       fo(3)*((rtwo*T00) 
     &- (rtwo*T00*pv(1)*(pv(1) + pv(4)))/vmass**2)) + 
     &    fi(4)*gc(2)*(fo(1)*(ft(1,2) - cone*ft(1,3) 
     &+ ft(2,1) - cone*ft(3,1))
     & + fo(2)*(-2*ft(1,1) - ft(1,4) - ft(4,1))) + 
     &    fi(1)*gc(1)*(fo(4)*(-ft(1,2) - cone*ft(1,3) - ft(2,1) 
     &- cone*ft(3,1)) + fo(3)*(-2*ft(1,1) - ft(1,4) - ft(4,1))) + 
     &    fi(3)*gc(2)*(fo(2)*(ft(1,2) + cone*ft(1,3) + ft(2,1) 
     &+ cone*ft(3,1)) + fo(1)*(-2*ft(1,1) + ft(1,4) + ft(4,1))) + 
     &    fi(2)*gc(1)*(fo(3)*(-ft(1,2) + cone*ft(1,3) - ft(2,1)
     & + cone*ft(3,1)) + fo(4)*(-2*ft(1,1) + ft(1,4) + ft(4,1))) + 
     &    (pv(1)*(fi(4)*gc(2)*(fo(1)*(-(T12*pv(1)) + cone*T13*pv(1)
     & + T23*(-(cone*pv(2)) + pv(3)) + T24*pv(4) - cone*T34*pv(4) + 
     &               rtwo*pv(2)*ft(2,2) - cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(2)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3) - 
     &T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) - 
     &               rtwo*pv(4)*ft(4,4))) + fi(1)*gc(1)*
     &          (fo(4)*(T12*pv(1) + cone*T13*pv(1) 
     &+ T23*(-(cone*pv(2)) - pv(3)) - T24*pv(4) - cone*T34*pv(4) - 
     &               rtwo*pv(2)*ft(2,2) - cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(3)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3) 
     &- T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) - 
     &               rtwo*pv(4)*ft(4,4))) + fi(3)*gc(2)*
     &          (fo(2)*(-(T12*pv(1)) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) + pv(3)) + T24*pv(4) + cone*T34*pv(4) + 
     &               rtwo*pv(2)*ft(2,2) + cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(1)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) + 
     &               rtwo*pv(4)*ft(4,4))) + fi(2)*gc(1)*
     &          (fo(3)*(T12*pv(1) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) - pv(3)) - T24*pv(4) + cone*T34*pv(4) - 
     &               rtwo*pv(2)*ft(2,2) + cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(4)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3) 
     &+ T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) + 
     &               rtwo*pv(4)*ft(4,4)))))/vmass**2
      
         jiot(2) = fi(3)*gc(2)*(fo(2)*((rtwo*T00) 
     &- (rtwo*T00*pv(2)*(-pv(2) - cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(1)*pv(2)*(pv(1) - pv(4)))/vmass**2) + 
     &    fi(2)*gc(1)*(fo(3)*(-rtwo*T00 - (rtwo*T00*pv(2)*(pv(2)
     & - cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(4)*pv(2)*(pv(1) - pv(4)))/vmass**2) + 
     &    fi(4)*gc(2)*(fo(1)*((rtwo*T00)
     & - (rtwo*T00*pv(2)*(-pv(2) + cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(2)*pv(2)*(pv(1) + pv(4)))/vmass**2) + 
     &    fi(1)*gc(1)*(fo(4)*(-rtwo*T00 - (rtwo*T00*pv(2)*(pv(2) 
     &+ cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(3)*pv(2)*(pv(1) + pv(4)))/vmass**2) + 
     &    fi(4)*gc(2)*(fo(1)*(2*ft(2,2) - cone*ft(2,3) - cone*ft(3,2))
     & + fo(2)*(-ft(1,2) - ft(2,1) - ft(2,4) - ft(4,2))) + 
     &    fi(1)*gc(1)*(fo(4)*(-2*ft(2,2) - cone*ft(2,3) - cone*ft(3,2)) 
     &+ fo(3)*(-ft(1,2) - ft(2,1) - ft(2,4) - ft(4,2))) + 
     &    fi(3)*gc(2)*(fo(2)*(2*ft(2,2) + cone*ft(2,3) + cone*ft(3,2)) 
     &+ fo(1)*(-ft(1,2) - ft(2,1) + ft(2,4) + ft(4,2))) + 
     &    fi(2)*gc(1)*(fo(3)*(-2*ft(2,2) + cone*ft(2,3) + cone*ft(3,2))
     & + fo(4)*(-ft(1,2) - ft(2,1) + ft(2,4) + ft(4,2))) + 
     &    (pv(2)*(fi(4)*gc(2)*(fo(1)*(-(T12*pv(1)) + cone*T13*pv(1)
     & + T23*(-(cone*pv(2)) + pv(3)) + T24*pv(4) - cone*T34*pv(4) + 
     &               rtwo*pv(2)*ft(2,2) - cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(2)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) - 
     &               rtwo*pv(4)*ft(4,4))) + fi(1)*gc(1)*
     &          (fo(4)*(T12*pv(1) + cone*T13*pv(1) + T23*(-(cone*pv(2))
     & - pv(3)) - T24*pv(4) - cone*T34*pv(4) - 
     &               rtwo*pv(2)*ft(2,2) - cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(3)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) - 
     &               rtwo*pv(4)*ft(4,4))) + fi(3)*gc(2)*
     &          (fo(2)*(-(T12*pv(1)) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) + pv(3)) + T24*pv(4) + cone*T34*pv(4) + 
     &               rtwo*pv(2)*ft(2,2) + cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(1)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) + 
     &               rtwo*pv(4)*ft(4,4))) + fi(2)*gc(1)*
     &          (fo(3)*(T12*pv(1) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) - pv(3)) - T24*pv(4) + cone*T34*pv(4) - 
     &               rtwo*pv(2)*ft(2,2) + cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(4)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3) 
     &+ T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) + 
     &               rtwo*pv(4)*ft(4,4)))))/vmass**2

         jiot(3) =fi(3)*gc(2)*(fo(2)*((cone*rtwo*T00)
     &	 - (rtwo*T00*pv(3)*(-pv(2) - cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(1)*pv(3)*(pv(1) - pv(4)))/vmass**2) + 
     &    fi(2)*gc(1)*(fo(3)*((cone*rtwo*T00) 
     &- (rtwo*T00*pv(3)*(pv(2) - cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(4)*pv(3)*(pv(1) - pv(4)))/vmass**2) + 
     &    fi(4)*gc(2)*(fo(1)*(-cone*rtwo*T00 
     &- (rtwo*T00*pv(3)*(-pv(2) + cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(2)*pv(3)*(pv(1) + pv(4)))/vmass**2) + 
     &    fi(1)*gc(1)*(fo(4)*(-cone*rtwo*T00
     & - (rtwo*T00*pv(3)*(pv(2) + cone*pv(3)))/vmass**2) - 
     &       (rtwo*T00*fo(3)*pv(3)*(pv(1) + pv(4)))/vmass**2) + 
     &    fi(4)*gc(2)*(fo(1)*(ft(2,3) + ft(3,2) 
     &- 2*cone*ft(3,3)) + fo(2)*(-ft(1,3) - ft(3,1) 
     &- ft(3,4) - ft(4,3))) + 
     &    fi(1)*gc(1)*(fo(4)*(-ft(2,3) - ft(3,2) 
     &- 2*cone*ft(3,3)) + fo(3)*(-ft(1,3) - ft(3,1) 
     &- ft(3,4) - ft(4,3))) + 
     &    fi(3)*gc(2)*(fo(2)*(ft(2,3) + ft(3,2) + 2*cone*ft(3,3))
     & + fo(1)*(-ft(1,3) - ft(3,1) + ft(3,4) + ft(4,3))) + 
     &    fi(2)*gc(1)*(fo(3)*(-ft(2,3) - ft(3,2) + 2*cone*ft(3,3)) 
     &+ fo(4)*(-ft(1,3) - ft(3,1) + ft(3,4) + ft(4,3))) + 
     &    (pv(3)*(fi(4)*gc(2)*(fo(1)*(-(T12*pv(1)) + cone*T13*pv(1) 
     &+ T23*(-(cone*pv(2)) + pv(3)) + T24*pv(4) - cone*T34*pv(4) + 
     &               rtwo*pv(2)*ft(2,2) - cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(2)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     &- T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) - 
     &               rtwo*pv(4)*ft(4,4))) + fi(1)*gc(1)*
     &          (fo(4)*(T12*pv(1) + cone*T13*pv(1) 
     &+ T23*(-(cone*pv(2)) - pv(3)) - T24*pv(4) - cone*T34*pv(4) - 
     &               rtwo*pv(2)*ft(2,2) - cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(3)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) - 
     &               rtwo*pv(4)*ft(4,4))) + fi(3)*gc(2)*
     &          (fo(2)*(-(T12*pv(1)) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) + pv(3)) + T24*pv(4) + cone*T34*pv(4) + 
     &               rtwo*pv(2)*ft(2,2) + cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(1)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) + 
     &               rtwo*pv(4)*ft(4,4))) + fi(2)*gc(1)*
     &          (fo(3)*(T12*pv(1) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) - pv(3)) - T24*pv(4) + cone*T34*pv(4) - 
     &               rtwo*pv(2)*ft(2,2) + cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(4)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) + 
     &               rtwo*pv(4)*ft(4,4)))))/vmass**2

         jiot(4) =fi(3)*gc(2)*(-((rtwo*T00*fo(2)*(-pv(2) 
     &	- cone*pv(3))*pv(4))/vmass**2) + 
     &       fo(1)*((rtwo*T00) - (rtwo*T00*(pv(1)
     & - pv(4))*pv(4))/vmass**2)) + 
     &    fi(2)*gc(1)*(-((rtwo*T00*fo(3)*(pv(2)
     & - cone*pv(3))*pv(4))/vmass**2) + 
     &       fo(4)*((rtwo*T00) - (rtwo*T00*(pv(1)
     & - pv(4))*pv(4))/vmass**2)) + 
     &    fi(4)*gc(2)*(-((rtwo*T00*fo(1)*(-pv(2) 
     &+ cone*pv(3))*pv(4))/vmass**2) + 
     &       fo(2)*(-rtwo*T00
     & - (rtwo*T00*pv(4)*(pv(1) + pv(4)))/vmass**2)) + 
     &    fi(1)*gc(1)*(-((rtwo*T00*fo(4)*(pv(2)
     & + cone*pv(3))*pv(4))/vmass**2) + 
     &       fo(3)*(-rtwo*T00 - (rtwo*T00*pv(4)*(pv(1)
     & + pv(4)))/vmass**2)) + 
     &    fi(4)*gc(2)*(fo(1)*(ft(2,4) - cone*ft(3,4)
     & + ft(4,2) - cone*ft(4,3)) + fo(2)*(-ft(1,4) 
     &- ft(4,1) - 2*ft(4,4))) + 
     &    fi(1)*gc(1)*(fo(4)*(-ft(2,4) - cone*ft(3,4)
     & - ft(4,2) - cone*ft(4,3)) + fo(3)*(-ft(1,4)
     & - ft(4,1) - 2*ft(4,4))) + 
     &    fi(3)*gc(2)*(fo(2)*(ft(2,4) + cone*ft(3,4) 
     &+ ft(4,2) + cone*ft(4,3)) + fo(1)*(-ft(1,4) 
     &- ft(4,1) + 2*ft(4,4))) + 
     &    fi(2)*gc(1)*(fo(3)*(-ft(2,4) + cone*ft(3,4)
     & - ft(4,2) + cone*ft(4,3)) + fo(4)*(-ft(1,4)
     & - ft(4,1) + 2*ft(4,4))) + 
     &    (pv(4)*(fi(4)*gc(2)*(fo(1)*(-(T12*pv(1)) 
     &+ cone*T13*pv(1) + T23*(-(cone*pv(2)) + pv(3))
     & + T24*pv(4) - cone*T34*pv(4) + 
     &               rtwo*pv(2)*ft(2,2) - cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(2)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) - 
     &               rtwo*pv(4)*ft(4,4))) + fi(1)*gc(1)*
     &          (fo(4)*(T12*pv(1) + cone*T13*pv(1) 
     &+ T23*(-(cone*pv(2)) - pv(3)) - T24*pv(4) - cone*T34*pv(4) - 
     &               rtwo*pv(2)*ft(2,2) - cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(3)*(-(T12*pv(2)) - T24*pv(2) - T13*pv(3)
     & - T34*pv(3) - T14*(-pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) - 
     &               rtwo*pv(4)*ft(4,4))) + fi(3)*gc(2)*
     &          (fo(2)*(-(T12*pv(1)) - cone*T13*pv(1) 
     &+ T23*(cone*pv(2) + pv(3)) + T24*pv(4) + cone*T34*pv(4) + 
     &               rtwo*pv(2)*ft(2,2) + cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(1)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3) 
     &+ T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) + 
     &               rtwo*pv(4)*ft(4,4))) + fi(2)*gc(1)*
     &          (fo(3)*(T12*pv(1) - cone*T13*pv(1) + T23*(cone*pv(2) 
     &- pv(3)) - T24*pv(4) + cone*T34*pv(4) - 
     &               rtwo*pv(2)*ft(2,2) + cone*rtwo*pv(3)*ft(3,3)) + 
     &            fo(4)*(-(T12*pv(2)) + T24*pv(2) - T13*pv(3)
     & + T34*pv(3) - T14*(pv(1) + pv(4)) + rtwo*pv(1)*ft(1,1) + 
     &               rtwo*pv(4)*ft(4,4)))))/vmass**2

      else

         jiot(1) = rtwo*T00*fi(1)*fo(3)*gc(1)
     & + rtwo*T00*fi(2)*fo(4)*gc(1) + rtwo*T00*fi(3)*fo(1)*gc(2) + 
     &    rtwo*T00*fi(4)*fo(2)*gc(2) + fi(4)*gc(2)*(fo(1)*(ft(1,2) 
     &- cone*ft(1,3) + ft(2,1) - cone*ft(3,1)) + 
     &       fo(2)*(-2*ft(1,1) - ft(1,4) - ft(4,1))) + 
     &    fi(1)*gc(1)*(fo(4)*(-ft(1,2) - cone*ft(1,3) - ft(2,1) 
     &- cone*ft(3,1)) + fo(3)*(-2*ft(1,1) - ft(1,4) - ft(4,1))) + 
     &    fi(3)*gc(2)*(fo(2)*(ft(1,2) + cone*ft(1,3) + ft(2,1) 
     &+ cone*ft(3,1)) + fo(1)*(-2*ft(1,1) + ft(1,4) + ft(4,1))) + 
     &    fi(2)*gc(1)*(fo(3)*(-ft(1,2) + cone*ft(1,3) - ft(2,1)
     & + cone*ft(3,1)) + fo(4)*(-2*ft(1,1) + ft(1,4) + ft(4,1)))

         jiot(2) = -(rtwo*T00*fi(2)*fo(3)*gc(1)) 
     &  - rtwo*T00*fi(1)*fo(4)*gc(1)
     &	 + rtwo*T00*fi(4)*fo(1)*gc(2) + 
     &    rtwo*T00*fi(3)*fo(2)*gc(2) + fi(4)*gc(2)*(fo(1)*(2*ft(2,2) 
     &- cone*ft(2,3) - cone*ft(3,2)) + 
     &       fo(2)*(-ft(1,2) - ft(2,1) - ft(2,4) - ft(4,2))) + 
     &    fi(1)*gc(1)*(fo(4)*(-2*ft(2,2) - cone*ft(2,3) - cone*ft(3,2))
     & + fo(3)*(-ft(1,2) - ft(2,1) - ft(2,4) - ft(4,2))) + 
     &    fi(3)*gc(2)*(fo(2)*(2*ft(2,2) + cone*ft(2,3) + cone*ft(3,2))
     & + fo(1)*(-ft(1,2) - ft(2,1) + ft(2,4) + ft(4,2))) + 
     &    fi(2)*gc(1)*(fo(3)*(-2*ft(2,2) + cone*ft(2,3) + cone*ft(3,2)) 
     &+ fo(4)*(-ft(1,2) - ft(2,1) + ft(2,4) + ft(4,2)))

         jiot(3) = cone*rtwo*T00*fi(2)*fo(3)*gc(1) 
     &	- cone*rtwo*T00*fi(1)*fo(4)*gc(1)
     & - cone*rtwo*T00*fi(4)*fo(1)*gc(2) + 
     &    cone*rtwo*T00*fi(3)*fo(2)*gc(2) + fi(4)*gc(2)*
     &     (fo(1)*(ft(2,3) + ft(3,2) - 2*cone*ft(3,3))
     & + fo(2)*(-ft(1,3) - ft(3,1) - ft(3,4) - ft(4,3))) + 
     &    fi(1)*gc(1)*(fo(4)*(-ft(2,3) - ft(3,2) - 2*cone*ft(3,3))
     & + fo(3)*(-ft(1,3) - ft(3,1) - ft(3,4) - ft(4,3))) + 
     &    fi(3)*gc(2)*(fo(2)*(ft(2,3) + ft(3,2) + 2*cone*ft(3,3)) 
     &+ fo(1)*(-ft(1,3) - ft(3,1) + ft(3,4) + ft(4,3))) + 
     &    fi(2)*gc(1)*(fo(3)*(-ft(2,3) - ft(3,2) + 2*cone*ft(3,3)) 
     &+ fo(4)*(-ft(1,3) - ft(3,1) + ft(3,4) + ft(4,3)))
       
         jiot(4) = -(rtwo*T00*fi(1)*fo(3)*gc(1))
     &	 + rtwo*T00*fi(2)*fo(4)*gc(1) + rtwo*T00*fi(3)*fo(1)*gc(2) - 
     &    rtwo*T00*fi(4)*fo(2)*gc(2) + fi(4)*gc(2)*(fo(1)*(ft(2,4) 
     &- cone*ft(3,4) + ft(4,2) - cone*ft(4,3)) + 
     &       fo(2)*(-ft(1,4) - ft(4,1) - 2*ft(4,4))) + 
     &    fi(1)*gc(1)*(fo(4)*(-ft(2,4) - cone*ft(3,4) - ft(4,2)
     & - cone*ft(4,3)) + fo(3)*(-ft(1,4) - ft(4,1) - 2*ft(4,4))) + 
     &    fi(3)*gc(2)*(fo(2)*(ft(2,4) + cone*ft(3,4) + ft(4,2)
     & + cone*ft(4,3)) + fo(1)*(-ft(1,4) - ft(4,1) + 2*ft(4,4))) + 
     &    fi(2)*gc(1)*(fo(3)*(-ft(2,4) + cone*ft(3,4) - ft(4,2) 
     &+ cone*ft(4,3)) + fo(4)*(-ft(1,4) - ft(4,1) + 2*ft(4,4)))  

      endif

      do i = 1,4
         jiot(i) = -jiot(i)*d*gt
      end do

      return
      end
      subroutine jioxxx(fi,fo,gc,vmass,vwidth , jio)
c
c This subroutine computes an off-shell vector current from an external
c fermion pair.  The vector boson propagator is given in Feynman gauge
c for a massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gc(2)          : coupling constants                  gvf
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
c     
      implicit none
      double complex fi(6),fo(6),gc(2),jio(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)


      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

c

      jio(5) = fo(5)-fi(5)
      jio(6) = fo(6)-fi(6)

      q(0) = dble( jio(5))
      q(1) = dble( jio(6))
      q(2) = dimag(jio(6))
      q(3) = dimag(jio(5))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2


      if ( vmass.ne.rZero ) then

         d = rOne/dcmplx( q2-vm2, vmass*vwidth )
c     For the running width, use below instead of the above d.
c     d = rOne/dcmplx( q2-vm2, max(vwidth*q2/vmass,rZero) )


         if ( gc(2).ne.cZero ) then
            c0 =  gc(1)*( fo(3)*fi(1)+fo(4)*fi(2))
     &           +gc(2)*( fo(1)*fi(3)+fo(2)*fi(4))
            c1 = -gc(1)*( fo(3)*fi(2)+fo(4)*fi(1))
     &           +gc(2)*( fo(1)*fi(4)+fo(2)*fi(3))
            c2 =( gc(1)*( fo(3)*fi(2)-fo(4)*fi(1))
     &           +gc(2)*(-fo(1)*fi(4)+fo(2)*fi(3)))*cImag
            c3 =  gc(1)*(-fo(3)*fi(1)+fo(4)*fi(2))
     &           +gc(2)*( fo(1)*fi(3)-fo(2)*fi(4))
         else
            d = d*gc(1)
            c0 =   fo(3)*fi(1)+fo(4)*fi(2)
            c1 =  -fo(3)*fi(2)-fo(4)*fi(1)
            c2 = ( fo(3)*fi(2)-fo(4)*fi(1))*cImag
            c3 =  -fo(3)*fi(1)+fo(4)*fi(2)
         end if

c     Fabio's implementation of the fixed width
         cm2=dcmplx( vm2, -vmass*vwidth )
c     cs = (q(0)*c0-q(1)*c1-q(2)*c2-q(3)*c3)/vm2
         cs = (q(0)*c0-q(1)*c1-q(2)*c2-q(3)*c3)/cm2
         jio(1) = (c0-cs*q(0))*d
         jio(2) = (c1-cs*q(1))*d
         jio(3) = (c2-cs*q(2))*d
         jio(4) = (c3-cs*q(3))*d

      else

         d = dcmplx( rOne/q2, rZero )
         if ( gc(2).ne.cZero ) then
            jio(1) = ( gc(1)*( fo(3)*fi(1)+fo(4)*fi(2))
     &                +gc(2)*( fo(1)*fi(3)+fo(2)*fi(4)) )*d
            jio(2) = (-gc(1)*( fo(3)*fi(2)+fo(4)*fi(1))
     &                +gc(2)*( fo(1)*fi(4)+fo(2)*fi(3)) )*d
            jio(3) = ( gc(1)*( fo(3)*fi(2)-fo(4)*fi(1))
     &                +gc(2)*(-fo(1)*fi(4)+fo(2)*fi(3)))
     &               *d*cImag
            jio(4) = ( gc(1)*(-fo(3)*fi(1)+fo(4)*fi(2))
     &                +gc(2)*( fo(1)*fi(3)-fo(2)*fi(4)) )*d
         else
            d = d*gc(1)
            jio(1) =  ( fo(3)*fi(1)+fo(4)*fi(2))*d
            jio(2) = -( fo(3)*fi(2)+fo(4)*fi(1))*d
            jio(3) =  ( fo(3)*fi(2)-fo(4)*fi(1))*d*cImag
            jio(4) =  (-fo(3)*fi(1)+fo(4)*fi(2))*d
         end if

      end if
c
      return
      end
      subroutine jssxxx(s1,s2,gc,vmass,vwidth , jss)
c
c This subroutine computes an off-shell vector current from the vector-
c scalar-scalar coupling.  The coupling is absent in the minimal SM in
c unitary gauge.  The propagator is given in Feynman gauge for a
c massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex gc             : coupling constant (s1 charge)
c       real    vmass          : mass  of output vector v
c       real    vwidth         : width of output vector v
c
c examples of the coupling constant g for susy particles are as follows:
c   -----------------------------------------------------------
c   |    s1    | (q,i3) of s1  ||   v=A   |   v=Z   |   v=W   |
c   -----------------------------------------------------------
c   | nu~_l    | (  0  , +1/2) ||   ---   |  gzn(1) |  gwf(1) |
c   | e~_l     | ( -1  , -1/2) ||  gal(1) |  gzl(1) |  gwf(1) |
c   | u~_l     | (+2/3 , +1/2) ||  gau(1) |  gzu(1) |  gwf(1) |
c   | d~_l     | (-1/3 , -1/2) ||  gad(1) |  gzd(1) |  gwf(1) |
c   -----------------------------------------------------------
c   | e~_r-bar | ( +1  ,  0  ) || -gal(2) | -gzl(2) | -gwf(2) |
c   | u~_r-bar | (-2/3 ,  0  ) || -gau(2) | -gzu(2) | -gwf(2) |
c   | d~_r-bar | (+1/3 ,  0  ) || -gad(2) | -gzd(2) | -gwf(2) |
c   -----------------------------------------------------------
c where the s1 charge is defined by the flowing-OUT quantum number.
c
c output:
c       complex jss(6)         : vector current            j^mu(v:s1,s2)
c     
      implicit none
      double complex s1(3),s2(3),gc,jss(6),dg,adg
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision pp(0:3),pa(0:3),q(0:3),vmass,vwidth
      double precision q2,vm2,mp2,ma2,m2d

      double precision rZero
      parameter( rZero = 0.0d0 )

c

      jss(5) = s1(2)+s2(2)
      jss(6) = s1(3)+s2(3)

      q(0) = dble( jss(5))
      q(1) = dble( jss(6))
      q(2) = dimag(jss(6))
      q(3) = dimag(jss(5))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2


      if ( vmass.ne.rZero ) then

         dg = gc/dcmplx( q2-vm2, vmass*vwidth )
c  For the running width, use below instead of the above dg.
c         dg = g/dcmplx( q2-vm2, max(vwidth*q2/vmass,rZero) )


         adg = dg*s1(1)*s2(1)

         pp(0) = dble( s1(2))
         pp(1) = dble( s1(3))
         pp(2) = dimag(s1(3))
         pp(3) = dimag(s1(2))
         pa(0) = dble( s2(2))
         pa(1) = dble( s2(3))
         pa(2) = dimag(s2(3))
         pa(3) = dimag(s2(2))
         mp2 = pp(0)**2-(pp(1)**2+pp(2)**2+pp(3)**2)
         ma2 = pa(0)**2-(pa(1)**2+pa(2)**2+pa(3)**2)
         m2d = mp2-ma2

c     Fabio's implementation of the fixed width
         cm2=dcmplx( vm2, -vmass*vwidth )
c     jss(1) = adg*( (pp(0)-pa(0)) - q(0)*m2d/vm2)
c     jss(2) = adg*( (pp(1)-pa(1)) - q(1)*m2d/vm2)
c     jss(3) = adg*( (pp(2)-pa(2)) - q(2)*m2d/vm2)
c     jss(4) = adg*( (pp(3)-pa(3)) - q(3)*m2d/vm2)
         jss(1) = adg*( (pp(0)-pa(0)) - q(0)*m2d/cm2)
         jss(2) = adg*( (pp(1)-pa(1)) - q(1)*m2d/cm2)
         jss(3) = adg*( (pp(2)-pa(2)) - q(2)*m2d/cm2)
         jss(4) = adg*( (pp(3)-pa(3)) - q(3)*m2d/cm2)

      else

         adg = gc*s1(1)*s2(1)/q2

         jss(1) = adg*dble( s1(2)-s2(2))
         jss(2) = adg*dble( s1(3)-s2(3))
         jss(3) = adg*dimag(s1(3)-s2(3))
         jss(4) = adg*dimag(s1(2)-s2(2))

      endif
c
      return
      end
      subroutine jvshxx(vc,sc,gc,xm,xw , jvs)
c
c- by RF - Mar. 2006
c
c
c This subroutine computes an off-shell vector current from the vector-
c vector-(pseudo-)scalar effective coupling.  
c input:
c       complex vc(6)          : input vector
c       complex sc(3)          : input scalar
c       real xm                : mass outgoing vector
c       real xw                : width outgoing vector
c       complex gc(2)          : coupling constant: gc(1) scalar
c                                                   gc(2) pseudo-scalar
c
c output:
c       complex jvs(6)         : outgoing vector current
c
c
      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex vc(DIM),sc(DIM),jvs(DIM),jvs1(DIM),jvs2(DIM)
      double complex qvc,dg
      double precision qp1,p12,p13,p14,p23,p24,p34
      double precision p1(4),q(4),xm,xw,q2
      double complex gc(2)

      jvs(5) = vc(5)+sc(2)
      jvs(6) = vc(6)+sc(3)

      p1(1) = dble( vc(5))
      p1(2) = dble( vc(6))
      p1(3) = dimag(vc(6))
      p1(4) = dimag(vc(5))

      q(1) = -dble( jvs(5))
      q(2) = -dble( jvs(6))
      q(3) = -dimag(jvs(6))
      q(4) = -dimag(jvs(5))

      q2 = q(1)**2 - q(2)**2 - q(3)**2 - q(4)**2

      jvs1(1) = (0D0,0D0)
      jvs1(2) = (0D0,0D0)
      jvs1(3) = (0D0,0D0)
      jvs1(4) = (0D0,0D0)
      jvs2(1) = (0D0,0D0)
      jvs2(2) = (0D0,0D0)
      jvs2(3) = (0D0,0D0)
      jvs2(4) = (0D0,0D0)

      if(xm.eq.0d0)then
         dg = sc(1) /q2
      else
         dg = sc(1)/ dcmplx( q2-xm**2 , xm*xw )
      endif

      if (gc(1).NE.(0D0,0D0)) then
         qvc = vc(1)*q(1) - vc(2)*q(2) - vc(3)*q(3) - vc(4)*q(4)	
         qp1 = q(1)*p1(1) - q(2)*p1(2) - q(3)*p1(3) - q(4)*p1(4)	
         
         jvs1(1) = -gc(1)* (vc(1)*qp1 - p1(1)*qvc)
         jvs1(2) = -gc(1)* (vc(2)*qp1 - p1(2)*qvc)
         jvs1(3) = -gc(1)* (vc(3)*qp1 - p1(3)*qvc)
         jvs1(4) = -gc(1)* (vc(4)*qp1 - p1(4)*qvc)
      endif
 
      if (gc(2).NE.(0D0,0D0)) then
         p12 = p1(1)*q(2) - p1(2)*q(1)
         p13 = p1(1)*q(3) - p1(3)*q(1)
         p14 = p1(1)*q(4) - p1(4)*q(1)
         p23 = p1(2)*q(3) - p1(3)*q(2)
         p24 = p1(2)*q(4) - p1(4)*q(2)
         p34 = p1(3)*q(4) - p1(4)*q(3)


         jvs2(1)=  gc(2)* (           -vc(2)*p34 +vc(3)*p24 -vc(4)*p23)
         jvs2(2)= -gc(2)* ( vc(1)*p34            -vc(3)*p14 +vc(4)*p13)
         jvs2(3)= -gc(2)* (-vc(1)*p24 +vc(2)*p14            -vc(4)*p12)
         jvs2(4)= -gc(2)* ( vc(1)*p23 -vc(2)*p13 +vc(3)*p12           )
      endif


      jvs(1) = dg * (jvs1(1) + jvs2(1))
      jvs(2) = dg * (jvs1(2) + jvs2(2))
      jvs(3) = dg * (jvs1(3) + jvs2(3))
      jvs(4) = dg * (jvs1(4) + jvs2(4))

      return
      end
      subroutine jvsshx(vc,sc1,sc2,g1,xm,xw , jvs)
c
c- by RF - Mar. 2006
c
c
c This subroutine computes an off-shell vector current from the vector-
c vector-Higgs-Higgs effective coupling.  
c
      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex vc(DIM),sc1(DIM),sc2(DIM),jvs(DIM)
      double complex jvs1(DIM),jvs2(DIM),dg,qvc
      double precision qp1,p12,p13,p14,p23,p24,p34
      double precision p1(4),q(4),xm,xw,q2,vm2
      double complex g1(2)

      jvs(5) = vc(5)+sc1(2)+sc2(2)
      jvs(6) = vc(6)+sc1(3)+sc2(3)

      p1(1) = dble( vc(5))
      p1(2) = dble( vc(6))
      p1(3) = dimag(vc(6))
      p1(4) = dimag(vc(5))

      q(1) = -dble( jvs(5))
      q(2) = -dble( jvs(6))
      q(3) = -dimag(jvs(6))
      q(4) = -dimag(jvs(5))

      q2 = q(1)**2 - q(2)**2 - q(3)**2 - q(4)**2

      jvs1(1) = (0D0,0D0)
      jvs1(2) = (0D0,0D0)
      jvs1(3) = (0D0,0D0)
      jvs1(4) = (0D0,0D0)
      jvs2(1) = (0D0,0D0)
      jvs2(2) = (0D0,0D0)
      jvs2(3) = (0D0,0D0)
      jvs2(4) = (0D0,0D0)

      dg = sc1(1)*sc2(1) /q2

      if (g1(1).NE.(0D0,0D0)) then
         qvc = vc(1)*q(1) - vc(2)*q(2) - vc(3)*q(3) - vc(4)*q(4)	
         qp1 = q(1)*p1(1) - q(2)*p1(2) - q(3)*p1(3) - q(4)*p1(4)	
         
         jvs1(1) = g1(1)* (vc(1)*qp1 - p1(1)*qvc)
         jvs1(2) = g1(1)* (vc(2)*qp1 - p1(2)*qvc)
         jvs1(3) = g1(1)* (vc(3)*qp1 - p1(3)*qvc)
         jvs1(4) = g1(1)* (vc(4)*qp1 - p1(4)*qvc)
      endif
 
      if (g1(2).NE.(0D0,0D0)) then
         p12 = p1(1)*q(2) - p1(2)*q(1)
         p13 = p1(1)*q(3) - p1(3)*q(1)
         p14 = p1(1)*q(4) - p1(4)*q(1)
         p23 = p1(2)*q(3) - p1(3)*q(2)
         p24 = p1(2)*q(4) - p1(4)*q(2)
         p34 = p1(3)*q(4) - p1(4)*q(3)


         jvs2(1)= - g1(2)* (-vc(2)*p34 +vc(3)*p24 -vc(4)*p23)
         jvs2(2)=   g1(2)* ( vc(1)*p34 -vc(3)*p14 +vc(4)*p13)
         jvs2(3)=   g1(2)* (-vc(1)*p24 +vc(2)*p14 -vc(4)*p12)
         jvs2(4)=   g1(2)* ( vc(1)*p23 -vc(2)*p13 +vc(3)*p12)
      endif


      jvs(1) = dg * (jvs1(1) + jvs2(1))
      jvs(2) = dg * (jvs1(2) + jvs2(2))
      jvs(3) = dg * (jvs1(3) + jvs2(3))
      jvs(4) = dg * (jvs1(4) + jvs2(4))

      return
      end
      subroutine jvssxx(vc,s1,s2,gc,vmass,vwidth , jvss)
c
c This subroutine computes an off-shell vector current from the vector-
c vector-scalar-scalar coupling.  The vector propagator is given in
c Feynman gauge for a massless vector and in unitary gauge for a massive
c vector.
c
c input:
c       complex vc(6)          : input  vector                        v
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex gc             : coupling constant                 gvvhh
c       real    vmass          : mass  of output vector v'
c       real    vwidth         : width of output vector v'
c
c output:
c       complex jvss(6)        : vector current         j^mu(v':v,s1,s2)
c     
      implicit none
      double complex vc(6),s1(3),s2(3),gc,jvss(6),dg
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex vk

      double precision rZero
      parameter( rZero = 0.0d0 )

c

      jvss(5) = vc(5)+s1(2)+s2(2)
      jvss(6) = vc(6)+s1(3)+s2(3)

      q(0) = dble( jvss(5))
      q(1) = dble( jvss(6))
      q(2) = dimag(jvss(6))
      q(3) = dimag(jvss(5))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2


      if ( vmass.ne.rZero ) then

         dg = gc*s1(1)*s2(1)/dcmplx( q2-vm2, vmass*vwidth )
c  For the running width, use below instead of the above dg.
c         dg = gc*s1(1)*s2(1)/cmplx( q2-vm2 , max( vwidth*q2/vmass ,rZero))

c     Fabio's implementation of the fixed width
         cm2=dcmplx( vm2, -vmass*vwidth )
c     vk = (q(0)*vc(1)-q(1)*vc(2)-q(2)*vc(3)-q(3)*vc(4))/vm2
         vk = (q(0)*vc(1)-q(1)*vc(2)-q(2)*vc(3)-q(3)*vc(4))/cm2

         jvss(1) = dg*(vc(1)-vk*q(0))
         jvss(2) = dg*(vc(2)-vk*q(1))
         jvss(3) = dg*(vc(3)-vk*q(2))
         jvss(4) = dg*(vc(4)-vk*q(3))

      else

         dg = gc*s1(1)*s2(1)/q2

         jvss(1) = dg*vc(1)
         jvss(2) = dg*vc(2)
         jvss(3) = dg*vc(3)
         jvss(4) = dg*vc(4)

      endif
c
      return
      end
      subroutine jvstxx(ga,tc,g,xm,xw,jw)
c
c- by RF - Feb. 2006
c
c This subroutine computes the off-shell vector current
c in terms of the vector w1 and tensor jt 
c
c     input:
c         complex ga                  : incoming vector boson (gluon)
c         complex tc                  : incoming tensor particle
c         real    g                   : coupling constant
c
c     output:
c         complex jw                  : outgoing vector current (gluon)
c
c     not used:
c         real xm, xw
c

      implicit none

c dimension of the current set to arbitrary length
      integer DIM
      parameter (DIM=18)
c      include "dimension.inc"
      double complex ga(DIM),jw(DIM),tc(DIM)

      double complex gt1(4), gt2(4), dv, q2
      double precision q(0:3),g
      double precision xm,xw
	
      double precision sqrTwo
      parameter( sqrTwo = 1.41421356237309514547462185873882845044d0 )

      jw(5) = ga(5) + tc(17)
      jw(6) = ga(6) + tc(18)

      q(0) = - dble( jw(5))
      q(1) = - dble( jw(6))
      q(2) = - dimag(jw(6))
      q(3) = - dimag(jw(5))

      q2 = q(0)**2 -q(1)**2 -q(2)**2 -q(3)**2


c Gluon propagator:

      dv = - g /sqrTwo /dcmplx( q2 )

cfax--> madgraph bug compensation

      dv = -dv

c First take the inner product of the incoming vector with the
c second index of the tensor

      gt2(1) = ga(1)*tc( 1) - ga(2)*tc( 2) - ga(3)*tc( 3) - ga(4)*tc( 4)
      gt2(2) = ga(1)*tc( 5) - ga(2)*tc( 6) - ga(3)*tc( 7) - ga(4)*tc( 8)
      gt2(3) = ga(1)*tc( 9) - ga(2)*tc(10) - ga(3)*tc(11) - ga(4)*tc(12)
      gt2(4) = ga(1)*tc(13) - ga(2)*tc(14) - ga(3)*tc(15) - ga(4)*tc(16)

c and with the first index of the tensor

      gt1(1) = ga(1)*tc( 1) - ga(2)*tc( 5) - ga(3)*tc( 9) - ga(4)*tc(13)
      gt1(2) = ga(1)*tc( 2) - ga(2)*tc( 6) - ga(3)*tc(10) - ga(4)*tc(14)
      gt1(3) = ga(1)*tc( 3) - ga(2)*tc( 7) - ga(3)*tc(11) - ga(4)*tc(15)
      gt1(4) = ga(1)*tc( 4) - ga(2)*tc( 8) - ga(3)*tc(12) - ga(4)*tc(16)

c The current is the difference of gt1 and gt2 with the remaining
c tensor indices the indices of the outgoing vector current

      jw(1) = dv * (gt1(1) - gt2(1))
      jw(2) = dv * (gt1(2) - gt2(2))
      jw(3) = dv * (gt1(3) - gt2(3))
      jw(4) = dv * (gt1(4) - gt2(4))

      return
      end

      subroutine jvsxxx(vc,sc,gc,vmass,vwidth , jvs)
c
c This subroutine computes an off-shell vector current from the vector-
c vector-scalar coupling.  The vector propagator is given in Feynman
c gauge for a massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex vc(6)          : input vector                          v
c       complex sc(3)          : input scalar                          s
c       complex gc             : coupling constant                  gvvh
c       real    vmass          : mass  of output vector v'
c       real    vwidth         : width of output vector v'
c
c output:
c       complex jvs(6)         : vector current             j^mu(v':v,s)
c     
      implicit none
      double complex vc(6),sc(3),gc,jvs(6),dg,vk
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision q(0:3),vmass,vwidth,q2,vm2

      double precision rZero
      parameter( rZero = 0.0d0 )

c

      jvs(5) = vc(5)+sc(2)
      jvs(6) = vc(6)+sc(3)

      q(0) = dble( jvs(5))
      q(1) = dble( jvs(6))
      q(2) = dimag(jvs(6))
      q(3) = dimag(jvs(5))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2


      if ( vmass.ne.rZero ) then

         dg = gc*sc(1)/dcmplx( q2-vm2, vmass*vwidth )
c  For the running width, use below instead of the above dg.
c         dg = g*sc(1)/dcmplx( q2-vm2, max(vwidth*q2/vmass,rZero) )

c     Fabio's implementation of the fixed width
         cm2=dcmplx( vm2, -vmass*vwidth )
c     vk = (-q(0)*vc(1)+q(1)*vc(2)+q(2)*vc(3)+q(3)*vc(4))/vm2
         vk = (-q(0)*vc(1)+q(1)*vc(2)+q(2)*vc(3)+q(3)*vc(4))/cm2

         jvs(1) = dg*(q(0)*vk+vc(1))
         jvs(2) = dg*(q(1)*vk+vc(2))
         jvs(3) = dg*(q(2)*vk+vc(3))
         jvs(4) = dg*(q(3)*vk+vc(4))

      else

         dg=gc*sc(1)/q2

         jvs(1) = dg*vc(1)
         jvs(2) = dg*vc(2)
         jvs(3) = dg*vc(3)
         jvs(4) = dg*vc(4)

      endif
c
      return
      end
      subroutine jvtaxx(ga,tc,g,xm1,xw,jw)
c
c- by RF - Feb. 2006
c
c This subroutine computes the off-shell vector current
c in terms of the vector w1 and tensor jt 
c
c     input:
c         complex ga                  : incoming vector boson
c         complex tc                  : incoming tensor particle
c         real    g                   : coupling constant
c
c     output:
c         complex jw                  : outgoing vector current
c
c     not used:
c         real xm1, xw
c

      implicit none

c dimension of the current set to arbitrary length
c      integer DIM
c      parameter (DIM=18)
      include "dimension.inc"
      double complex ga(DIM),jw(DIM),tc(DIM)
      double complex gt1(4), gt2(4)
      double precision q(0:3),g,q2,dv
      double precision xm1,xw
      double precision sqrTwo
      parameter( sqrTwo = 1.41421356237309514547462185873882845044d0 )

      jw(5) = ga(5) + tc(17)
      jw(6) = ga(6) + tc(18)

      q(0) = - dble( jw(5))
      q(1) = - dble( jw(6))
      q(2) = - dimag(jw(6))
      q(3) = - dimag(jw(5))

      q2 = q(0)**2 -q(1)**2 -q(2)**2 -q(3)**2



c Gluon propagator:

      dv = - g /sqrTwo / q2


c First take the inner product of the incoming vector with the
c second index of the tensor

      gt2(1) = ga(1)*tc( 1) - ga(2)*tc( 2) - ga(3)*tc( 3) - ga(4)*tc( 4)
      gt2(2) = ga(1)*tc( 5) - ga(2)*tc( 6) - ga(3)*tc( 7) - ga(4)*tc( 8)
      gt2(3) = ga(1)*tc( 9) - ga(2)*tc(10) - ga(3)*tc(11) - ga(4)*tc(12)
      gt2(4) = ga(1)*tc(13) - ga(2)*tc(14) - ga(3)*tc(15) - ga(4)*tc(16)

c and with the first index of the tensor

      gt1(1) = ga(1)*tc( 1) - ga(2)*tc( 5) - ga(3)*tc( 9) - ga(4)*tc(13)
      gt1(2) = ga(1)*tc( 2) - ga(2)*tc( 6) - ga(3)*tc(10) - ga(4)*tc(14)
      gt1(3) = ga(1)*tc( 3) - ga(2)*tc( 7) - ga(3)*tc(11) - ga(4)*tc(15)
      gt1(4) = ga(1)*tc( 4) - ga(2)*tc( 8) - ga(3)*tc(12) - ga(4)*tc(16)

c The current is the difference of gt1 and gt2 with the remaining
c tensor indices the indices of the outgoing vector current

      jw(1) = dv * (gt1(1) - gt2(1))
      jw(2) = dv * (gt1(2) - gt2(2))
      jw(3) = dv * (gt1(3) - gt2(3))
      jw(4) = dv * (gt1(4) - gt2(4))

      return
      end

      subroutine jvtcxx(vc,tc,gt,vmass,vwidth , jvt)
c-------------------CP3  2009.10-----------------
c
c This subroutine computes an off-shell vector current from 
c the coupling of two gauge bosons and a non-propagating tensor boson.
c
c input:
c       complex vc(6)          : input vector                                                     v
c       complex tc(18)         : input non-propagating tensor                          T
c       complex gt            : coupling constant         gt=gs
c       real    vmass         : mass  of output vector  v'
c       real    vwidth         :  width of output vector    v'
c
c output:
c       complex jvt(6)         : vector current             j^mu(v':v,T)
c
      implicit none
      double complex vc(6), tc(18), jvt(6)
      double precision vmass, vwidth,gt

      double complex ft(6,4)
      double precision pv2(4), pp2

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)
 
      jvt(5) = vc(5)+ft(5,1)
      jvt(6) = vc(6)+ft(6,1)
      pv2(1) = dreal(jvt(5))
      pv2(2) = dreal(jvt(6))
      pv2(3) = dimag(jvt(6))
      pv2(4) = dimag(jvt(5))
      pp2 = pv2(1)**2 - pv2(2)**2 - pv2(3)**2 - pv2(4)**2


        do i=1,4
            jvt(i)=-(ft(1,i)*vc(1)) + ft(2,i)*vc(2) + ft(3,i)*vc(3) 
     & + ft(4,i)*vc(4) 
            jvt(i)=jvt(i)*gt/dcmplx(pp2, 0.0d0)
         enddo

         

      return
      end
      subroutine jvtxxx(vc,tc,gt,vmass,vwidth , jvt)
c
c This subroutine computes an off-shell vector current from 
c the coupling of two gauge bosons and a tensor boson.
c
c input:
c       complex vc(6)          : input vector                          v
c       complex tc(18)         : input tensor                          T
c       complex gt             : coupling constant         gtv=-1/Lambda
c       real    vmass          : mass  of output vector v'
c       real    vwidth         : width of output vector v'
c
c output:
c       complex jvt(6)         : vector current             j^mu(v':v,T)
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex vc(6), tc(18), gt, jvt(6)
      double precision vmass, vwidth

      double complex ft(6,4),TVM(4),TKM(4)
      double precision MET(4,4)
      double complex T12, T13, T14, T23, T24, T34, T00
      double complex K2V1,K1V1
      double complex TKK,TK2V1, dum
      double precision pv1(4), pv2(4), F, pp2

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)
 
      jvt(5) = vc(5)+ft(5,1)
      jvt(6) = vc(6)+ft(6,1)

      pv1(1) = -dreal(vc(5))
      pv1(2) = -dreal(vc(6))
      pv1(3) = -dimag(vc(6))
      pv1(4) = -dimag(vc(5))

      pv2(1) = dreal(jvt(5))
      pv2(2) = dreal(jvt(6))
      pv2(3) = dimag(jvt(6))
      pv2(4) = dimag(jvt(5))

      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      K2V1 = pv2(1)*vc(1) - pv2(2)*vc(2) - pv2(3)*vc(3) - pv2(4)*vc(4)
      K1V1 = pv1(1)*vc(1) - pv1(2)*vc(2) - pv1(3)*vc(3) - pv1(4)*vc(4)
      F = pv1(1)*pv2(1) - pv1(2)*pv2(2) - pv1(3)*pv2(3) - pv1(4)*pv2(4)
      pp2 = pv2(1)**2 - pv2(2)**2 - pv2(3)**2 - pv2(4)**2

      TKK   = cZero
      TK2V1 = cZero

      do i = 1,4
         dum   = ft(i,i)*pv2(i)
         TKK   = TKK   + dum*pv1(i)
         dum   = ft(i,i)*vc(i)
         TK2V1 = TK2V1 + dum*pv2(i)
      end do
      
      TKK   = rTwo*TKK
      TK2V1 = rTwo*TK2V1
      
      TKK = TKK - T12*(pv1(1)*pv2(2) + pv1(2)*pv2(1))
     &          - T13*(pv1(1)*pv2(3) + pv1(3)*pv2(1))
     &          - T14*(pv1(1)*pv2(4) + pv1(4)*pv2(1))
     &          + T23*(pv1(2)*pv2(3) + pv1(3)*pv2(2))
     &          + T24*(pv1(2)*pv2(4) + pv1(4)*pv2(2))
     &          + T34*(pv1(3)*pv2(4) + pv1(4)*pv2(3))

      TK2V1 = TK2V1 - T12*(vc(1)*pv2(2) + vc(2)*pv2(1))
     &              - T13*(vc(1)*pv2(3) + vc(3)*pv2(1))
     &              - T14*(vc(1)*pv2(4) + vc(4)*pv2(1))
     &              + T23*(vc(2)*pv2(3) + vc(3)*pv2(2))
     &              + T24*(vc(2)*pv2(4) + vc(4)*pv2(2))
     &              + T34*(vc(3)*pv2(4) + vc(4)*pv2(3))


      do j=1,4

         TVM(j) =
     &MET(j,1)*(ft(1,1)*vc(1)-ft(2,1)*vc(2)
     &-ft(3,1)*vc(3)-ft(4,1)*vc(4))
     &-MET(j,2)*(ft(1,2)*vc(1)-ft(2,2)*vc(2)
     &-ft(3,2)*vc(3)-ft(4,2)*vc(4))
     &-MET(j,3)*(ft(1,3)*vc(1)-ft(2,3)*vc(2)
     &-ft(3,3)*vc(3)-ft(4,3)*vc(4))
     &-MET(j,4)*(ft(1,4)*vc(1)-ft(2,4)*vc(2)
     &-ft(3,4)*vc(3)-ft(4,4)*vc(4))
     &+
     &MET(j,1)*(ft(1,1)*vc(1)-ft(1,2)*vc(2)
     &-ft(1,3)*vc(3)-ft(1,4)*vc(4))
     &-MET(j,2)*(ft(2,1)*vc(1)-ft(2,2)*vc(2)
     &-ft(2,3)*vc(3)-ft(2,4)*vc(4))
     &-MET(j,3)*(ft(3,1)*vc(1)-ft(3,2)*vc(2)
     &-ft(3,3)*vc(3)-ft(3,4)*vc(4))
     &-MET(j,4)*(ft(4,1)*vc(1)-ft(4,2)*vc(2)
     &-ft(4,3)*vc(3)-ft(4,4)*vc(4))

         TKM(j) =
     &MET(j,1)*(ft(1,1)*pv1(1)-ft(2,1)*pv1(2)
     &-ft(3,1)*pv1(3)-ft(4,1)*pv1(4))
     &-MET(j,2)*(ft(1,2)*pv1(1)-ft(2,2)*pv1(2)
     &-ft(3,2)*pv1(3)-ft(4,2)*pv1(4))
     &-MET(j,3)*(ft(1,3)*pv1(1)-ft(2,3)*pv1(2)
     &-ft(3,3)*pv1(3)-ft(4,3)*pv1(4))
     &-MET(j,4)*(ft(1,4)*pv1(1)-ft(2,4)*pv1(2)
     &-ft(3,4)*pv1(3)-ft(4,4)*pv1(4))
     &+
     &MET(j,1)*(ft(1,1)*pv1(1)-ft(1,2)*pv1(2)
     &-ft(1,3)*pv1(3)-ft(1,4)*pv1(4))
     &-MET(j,2)*(ft(2,1)*pv1(1)-ft(2,2)*pv1(2)
     &-ft(2,3)*pv1(3)-ft(2,4)*pv1(4))
     &-MET(j,3)*(ft(3,1)*pv1(1)-ft(3,2)*pv1(2)
     &-ft(3,3)*pv1(3)-ft(3,4)*pv1(4))
     &-MET(j,4)*(ft(4,1)*pv1(1)-ft(4,2)*pv1(2)
     &-ft(4,3)*pv1(3)-ft(4,4)*pv1(4))
     
      enddo

      if ( vmass.ne.rZero ) then

         do i=1,4

            jvt(i) = -(vmass**2+F)*T00*vc(i)
     &+T00*K2V1*(pv1(i)+(1.0d0+F/vmass**2)*pv2(i)
     &-F/vmass**2*pv2(i))
     &+(vmass**2+F)*TVM(i)
     &-TK2V1*pv1(i)
     &-TK2V1*pv2(i)*(1.0d0+F/vmass**2)
     &+TKK*vc(i)
     &-K2V1*TKM(i)
     &+F/vmass**2*TK2V1*pv2(i)
	 
            jvt(i)=jvt(i)*gt/dcmplx(pp2-vmass**2, vmass*vwidth )

         enddo

      else
	
         do i=1,4

            jvt(i) = -F*T00*vc(i)
     &+K1V1*T00*(pv1(i)+pv2(i))
     &+T00*K2V1*(pv1(i)+pv2(i))
     &+F*TVM(i)
     &-TK2V1*(pv1(i)+pv2(i))
     &+TKK*vc(i)
     &-(K2V1+K1V1)*TKM(i)
  
            jvt(i)=jvt(i)*gt/dcmplx(pp2, 0.0d0)
            
         enddo

         
      endif

      return
      end
      subroutine jvvkxx(wm,wp,g,tmass,twidth, tc)
c
c This subroutine computes an off-shell tensor boson from two vector
c gauge bosons.
c
c input:
c       complex wm(6)          : vector               flow-in  V
c       complex wp(6)          : vector               flow-out V~
c       complex g(1)           : coupling constant    -kappa/2
c       real    g(2)           : V boson mass          m_V
c       real    tmass          : t boson mass          m_T
c       real    twidth         : t boson width         W_T
c
c output:
c       complex tc(18)        : tensor               KK mode T
c     
      implicit none
      double complex wm(18), wp(18), tc(18),g(2)
      double precision  tmass,twidth,vmass

      double precision pwm(4),pwp(4),preC,m2,p(4),pp,eta(4,4)
      integer i,j,ii,jj
      double complex C(4,4),D(4,4),wpwm,wpeta(4),wmeta(4)
      double complex pwpwm,wppwm,uio(4,4),denom

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      logical firsttime
      data firsttime/.true./

      vmass = dreal(g(2))
      if (firsttime)then
         write(*,*)'----------------------------------------------'
         write(*,*)'Using the jvvkxx HELAS routine. This routine  '
         write(*,*)'is only tested for gg>tt~ (sub)process.       '
         write(*,*)'----------------------------------------------'
         firsttime=.false.
      endif
c
      do i=1,16
         tc(i)=(0d0,0d0)
      enddo

      pwm(1) = dreal(wm(5))
      pwm(2) = dreal(wm(6))
      pwm(3) = dimag(wm(6))
      pwm(4) = dimag(wm(5))
      pwp(1) = dreal(wp(5))
      pwp(2) = dreal(wp(6))
      pwp(3) = dimag(wp(6))
      pwp(4) = dimag(wp(5))

      tc(17) = wm(5)+wp(5)
      tc(18) = wm(6)+wp(6)

      p(1) = -dble( tc(17))/tmass
      p(2) = -dble( tc(18))/tmass
      p(3) = -dimag(tc(18))/tmass
      p(4) = -dimag(tc(17))/tmass

      pp = p(1)**2 - p(2)**2 - p(3)**2 - p(4)**2

      do i=1,3
         do j=1+1,4
            eta(i,j)=0d0
            eta(j,i)=eta(i,j)
         enddo
      enddo
      eta(1,1)= 1d0
      eta(2,2)=-1d0
      eta(3,3)=-1d0
      eta(4,4)=-1d0

      m2 = tmass**2
      denom = dcmplx(pp*m2 - m2,tmass*twidth)

      preC = pwm(1)*pwp(1)-pwm(2)*pwp(2)-pwm(3)*pwp(3)-pwm(4)*pwp(4)

      if (vmass.ne.0d0)then
         preC=preC+vmass
      endif

      wpwm  =  wm(1)* wp(1)-  wm(2)* wp(2)-  wm(3)* wp(3)-  wm(4)* wp(4)
      pwpwm =  wm(1)*pwp(1)-  wm(2)*pwp(2)-  wm(3)*pwp(3)-  wm(4)*pwp(4)
      wppwm = pwm(1)* wp(1)- pwm(2)* wp(2)- pwm(3)* wp(3)- pwm(4)* wp(4)

      do i=1,4
         wpeta(i) = -wp(i)
         wmeta(i) = -wm(i)
      enddo
      wpeta(1) = -wpeta(1)
      wmeta(1) = -wmeta(1)

      do i=1,4
         do j=1,4
            C(i,j) = wpeta(i)*wmeta(j)+wpeta(j)*wmeta(i)-eta(i,j)*wpwm
            D(i,j) = eta(i,j)*pwpwm*wppwm
     & -wmeta(i)*pwm(j)*pwpwm-wpeta(i)*pwp(j)*wppwm+wpwm*pwm(i)*pwp(j)
     & -wmeta(j)*pwm(i)*pwpwm-wpeta(j)*pwp(i)*wppwm+wpwm*pwm(j)*pwp(i)
         enddo
      enddo

      do i=1,4
         do j=1,4
            uio(i,j)=preC*C(i,j)+D(i,j)
         enddo
      enddo
c make indices upper indices before contracting with propagator
      do i=2,4
         uio(1,i)=-uio(1,i)
         uio(i,1)=-uio(i,1)
      enddo

c multiply by propagator
      do i=1,4
         do j=1,4
            do ii=1,4
               do jj=1,4
                  tc(i+4*(j-1))=tc(i+4*(j-1))+
     &            (
     &                    (eta(i,ii)-p(i)*p(ii))*(eta(j,jj)-p(j)*p(jj))+
     &                    (eta(i,jj)-p(i)*p(jj))*(eta(j,ii)-p(j)*p(ii))-
     &            2d0/3d0*(eta(i,j)-p(i)*p(j))*(eta(ii,jj)-p(ii)*p(jj))
     &            )
     &            *uio(ii,jj)*g(1)/denom
               enddo
            enddo
         enddo
      enddo


      return
      end
      subroutine jvvsxx(ga,gb,sc,g1,g2,xm,xw,jggs)
c
c- by RF - Mar. 2006	
c
c This subroutine computes an off-shell vector current from the coupling
c of three gauge bosons and a scalar particle. The outgoing vector current
c is given in the Feynman gauge.
c
c input:
c       complex ga(6)       : first  incoming vector (gluon)
c       complex gb(6)       : second incoming vector (gluon)
c       complex sc(3)       : incoming scalar        (Higgs)
c       real    g1          : coupling constant      (QCD)
c       complex g2(2)       : coupling constant: g2(1) scalar
c                                                g2(2) pseudo-scalar
c
c output:
c       complex jggs(6)         : vector current
c not used:
c       real    xm, xw
c
 
      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex ga(DIM),gb(DIM),jggs(DIM),sc(DIM)
      double complex jggs1(DIM),jggs2(DIM)
      double complex sva,svb,vab,j12(0:3)
      double complex v12,v13,v14,v23,v24,v34,dv
      double precision p1(0:3),p2(0:3),q(0:3),q2,p12q(4)
      double precision g1,xm,xw
      double complex g2(2)


      jggs(5) = ga(5) + gb(5) + sc(2)
      jggs(6) = ga(6) + gb(6) + sc(3)

      p1(0) = dble( ga(5))
      p1(1) = dble( ga(6))
      p1(2) = dimag(ga(6))
      p1(3) = dimag(ga(5))

      p2(0) = dble( gb(5))
      p2(1) = dble( gb(6))
      p2(2) = dimag(gb(6))
      p2(3) = dimag(gb(5))

      q(0) = -dble( jggs(5))
      q(1) = -dble( jggs(6))
      q(2) = -dimag(jggs(6))
      q(3) = -dimag(jggs(5))

      q2 = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2

      dv = g1 * sc(1) /q2

      jggs1(1) = (0D0,0D0)
      jggs1(2) = (0D0,0D0)
      jggs1(3) = (0D0,0D0)
      jggs1(4) = (0D0,0D0)
      jggs2(1) = (0D0,0D0)
      jggs2(2) = (0D0,0D0)
      jggs2(3) = (0D0,0D0)
      jggs2(4) = (0D0,0D0)


      if (g2(1).NE.(0D0,0D0)) then
      vab = ga(1)*gb(1) - ga(2)*gb(2) - ga(3)*gb(3) - ga(4)*gb(4)
      sva =   (p2(0)-q(0)) *ga(1) - (p2(1)-q(1)) *ga(2)
     &      - (p2(2)-q(2)) *ga(3) - (p2(3)-q(3)) *ga(4)
      svb = - (p1(0)-q(0)) *gb(1) + (p1(1)-q(1)) *gb(2)
     &      + (p1(2)-q(2)) *gb(3) + (p1(3)-q(3)) *gb(4)

      jggs1(1)= g2(1)*((p1(0)-p2(0))*vab + sva*gb(1) + svb*ga(1))
      jggs1(2)= g2(1)*((p1(1)-p2(1))*vab + sva*gb(2) + svb*ga(2))
      jggs1(3)= g2(1)*((p1(2)-p2(2))*vab + sva*gb(3) + svb*ga(3))
      jggs1(4)= g2(1)*((p1(3)-p2(3))*vab + sva*gb(4) + svb*ga(4))
      endif
      
      if (g2(2).NE.(0D0,0D0)) then

      p12q(1) = p1(0) + p2(0) + q(0)
      p12q(2) = p1(1) + p2(1) + q(1)
      p12q(3) = p1(2) + p2(2) + q(2)
      p12q(4) = p1(3) + p2(3) + q(3)

      v12 = ga(1)*gb(2) - ga(2)*gb(1)
      v13 = ga(1)*gb(3) - ga(3)*gb(1)
      v14 = ga(1)*gb(4) - ga(4)*gb(1)
      v23 = ga(2)*gb(3) - ga(3)*gb(2)
      v24 = ga(2)*gb(4) - ga(4)*gb(2)
      v34 = ga(3)*gb(4) - ga(4)*gb(3)

      jggs2(1) =   g2(2)*(   v23*p12q(4) - v24*p12q(3) + v34*p12q(2) )
      jggs2(2) = - g2(2)*( - v13*p12q(4) + v14*p12q(3) - v34*p12q(1) )
      jggs2(3) = - g2(2)*(   v12*p12q(4) - v14*p12q(2) + v24*p12q(1) )
      jggs2(4) = - g2(2)*( - v12*p12q(3) + v13*p12q(2) - v23*p12q(1) )
      endif


      jggs(1) = dv * (jggs1(1) + jggs2(1))
      jggs(2) = dv * (jggs1(2) + jggs2(2))
      jggs(3) = dv * (jggs1(3) + jggs2(3))
      jggs(4) = dv * (jggs1(4) + jggs2(4))


      return
      end
      subroutine jvvtlx(ga,gb,sc,g1,g2,xm,xw,jggs)
c
c- by RF - Feb. 2006	
c
c This subroutine computes an off-shell vector current from the coupling
c of three gauge bosons and a scalar particle. The outgoing vector
c current is given in the Feynman gauge.
c
c input:
c       complex ga(6)       : first  incoming vector (gluon)
c       complex gb(6)       : second incoming vector (gluon)
c       complex sc(3)       : incoming scalar        (Higgs)
c       real    g1          : coupling constant      (QCD)
c       complex g2(2)       : coupling constant      (Higgs Effct. Thr.)
c
c output:
c       complex jggs(6)         : vector current
c not used:
c       real    xm, xw
c
 
      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex ga(DIM),gb(DIM),jggs(DIM),sc(DIM),g2(2)
      double complex jggs1(DIM),jggs2(DIM)
      double complex sva,svb,vab,j12(0:3)
      double complex v12,v13,v14,v23,v24,v34,dv
      double precision p1(0:3),p2(0:3),q(0:3),q2,p12q(4)
      double precision g1,xm,xw


      jggs(5) = ga(5) + gb(5) + sc(2)
      jggs(6) = ga(6) + gb(6) + sc(3)

      p1(0) = dble( ga(5))
      p1(1) = dble( ga(6))
      p1(2) = dimag(ga(6))
      p1(3) = dimag(ga(5))

      p2(0) = dble( gb(5))
      p2(1) = dble( gb(6))
      p2(2) = dimag(gb(6))
      p2(3) = dimag(gb(5))

      q(0) = -dble( jggs(5))
      q(1) = -dble( jggs(6))
      q(2) = -dimag(jggs(6))
      q(3) = -dimag(jggs(5))

      q2 = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2

      dv = g1 * sc(1) /q2

      jggs1(1) = (0D0,0D0)
      jggs1(2) = (0D0,0D0)
      jggs1(3) = (0D0,0D0)
      jggs1(4) = (0D0,0D0)
      jggs2(1) = (0D0,0D0)
      jggs2(2) = (0D0,0D0)
      jggs2(3) = (0D0,0D0)
      jggs2(4) = (0D0,0D0)


      if (g2(1).NE.(0D0,0D0)) then
      vab = ga(1)*gb(1) - ga(2)*gb(2) - ga(3)*gb(3) - ga(4)*gb(4)
      sva =   (p2(0)-q(0)) *ga(1) - (p2(1)-q(1)) *ga(2)
     &      - (p2(2)-q(2)) *ga(3) - (p2(3)-q(3)) *ga(4)
      svb = - (p1(0)-q(0)) *gb(1) + (p1(1)-q(1)) *gb(2)
     &      + (p1(2)-q(2)) *gb(3) + (p1(3)-q(3)) *gb(4)

      jggs1(1)= g2(1)*((p1(0)-p2(0))*vab+sva*gb(1)+svb*ga(1))
      jggs1(2)= g2(1)*((p1(1)-p2(1))*vab+sva*gb(2)+svb*ga(2))
      jggs1(3)= g2(1)*((p1(2)-p2(2))*vab+sva*gb(3)+svb*ga(3))
      jggs1(4)= g2(1)*((p1(3)-p2(3))*vab+sva*gb(4)+svb*ga(4))
      endif
      
      if (g2(2).NE.(0D0,0D0)) then

      p12q(1) = p1(0) + p2(0) + q(0)
      p12q(2) = p1(1) + p2(1) + q(1)
      p12q(3) = p1(2) + p2(2) + q(2)
      p12q(4) = p1(3) + p2(3) + q(3)

      v12 = ga(1)*gb(2) - ga(2)*gb(1)
      v13 = ga(1)*gb(3) - ga(3)*gb(1)
      v14 = ga(1)*gb(4) - ga(4)*gb(1)
      v23 = ga(2)*gb(3) - ga(3)*gb(2)
      v24 = ga(2)*gb(4) - ga(4)*gb(2)
      v34 = ga(3)*gb(4) - ga(4)*gb(3)

      jggs2(1) =   g2(2)*(   v23*p12q(4) - v24*p12q(3) + v34*p12q(2) )
      jggs2(2) = - g2(2)*( - v13*p12q(4) + v14*p12q(3) - v34*p12q(1) )
      jggs2(3) = - g2(2)*(   v12*p12q(4) - v14*p12q(2) + v24*p12q(1) )
      jggs2(4) = - g2(2)*( - v12*p12q(3) + v13*p12q(2) - v23*p12q(1) )
      endif


      jggs(1) = dv * (jggs1(1) + jggs2(1))
      jggs(2) = dv * (jggs1(2) + jggs2(2))
      jggs(3) = dv * (jggs1(3) + jggs2(3))
      jggs(4) = dv * (jggs1(4) + jggs2(4))


      return
      end
      subroutine jvvtxx(va,vb,tc,gc,gt,vmass,vwidth , jvvt)
c
c This subroutine computes an off-shell vector boson wavefunction from
c the four-point coupling of three gauge bosons and a tensor boson.
c
c input:
c       complex va(6)          : first  vector                        va
c       complex vb(6)          : second vector                        vb
c       complex tc(18)         : input  tensor                         T
c       real    gc             : coupling constant       gs (for gluons)
c       complex gt             : coupling constant         gtv=-1/Lambda
c       real    vmass          : mass  of output vector v
c       real    vwidth         : width of output vector v

c
c output:
c       complex jvvt(6)        : vector current          j^mu(v:va,vb,T)
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex va(6), vb(6), tc(18), gt, jvvt(6)
      double precision gc,vmass,vwidth

      double complex ft(6,4)
      double complex TV1M(4),TV2M(4),TK12M(4)
      double precision MET(4,4)
      double complex d,T00, T12, T13, T14, T23, T24, T34
      double complex V1V2,K1V2, K2V1
     &,K3V1,K3V2
      double precision K1K2

      double complex TV12,TKV1,TKV2,
     &TK3V1,TK3V2,TK312,dum
      double precision pva(4), pvb(4), pvc(4),
     &pv32,p31(4),p23(4),p12(4),K1K3,K2K3

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)


      jvvt(5) = va(5)+vb(5)+ft(5,1)
      jvvt(6) = va(6)+vb(6)+ft(6,1)

      pva(1) = dreal(va(5))
      pva(2) = dreal(va(6))
      pva(3) = dimag(va(6))
      pva(4) = dimag(va(5))

      pvb(1) = dreal(vb(5))
      pvb(2) = dreal(vb(6))
      pvb(3) = dimag(vb(6))
      pvb(4) = dimag(vb(5))

      pvc(1) = -dreal(jvvt(5))
      pvc(2) = -dreal(jvvt(6))
      pvc(3) = -dimag(jvvt(6))
      pvc(4) = -dimag(jvvt(5))
	
      pv32=pvc(1)**2-pvc(2)**2-pvc(3)**2-pvc(4)**2
      if ( vmass.gt.rZero ) then
         d =  gc/dcmplx( pv32-vmass**2, vmass*vwidth )
      else
         d =  gc/dcmplx( pv32, rZero )
      end if

      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) = 1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      p31(1) = pvc(1)-pva(1)
      p31(2) = pvc(2)-pva(2)
      p31(3) = pvc(3)-pva(3)
      p31(4) = pvc(4)-pva(4)
      
      p12(1) = pva(1)-pvb(1)
      p12(2) = pva(2)-pvb(2)
      p12(3) = pva(3)-pvb(3)
      p12(4) = pva(4)-pvb(4)
      
      p23(1) = pvb(1)-pvc(1)
      p23(2) = pvb(2)-pvc(2)
      p23(3) = pvb(3)-pvc(3)
      p23(4) = pvb(4)-pvc(4)
      
      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)
      
      V1V2 =  va(1)*vb(1) -  va(2)*vb(2) -  va(3)*vb(3) -  va(4)*vb(4)
      K1V2 = pva(1)*vb(1) - pva(2)*vb(2) - pva(3)*vb(3) - pva(4)*vb(4)
      K2V1 = pvb(1)*va(1) - pvb(2)*va(2) - pvb(3)*va(3) - pvb(4)*va(4)
      K3V1 = pvc(1)*va(1) - pvc(2)*va(2) - pvc(3)*va(3) - pvc(4)*va(4)
      K3V2 = pvc(1)*vb(1) - pvc(2)*vb(2) - pvc(3)*vb(3) - pvc(4)*vb(4)
      
      K1K3 = pva(1)*pvc(1)-pva(2)*pvc(2)-pva(3)*pvc(3)-pva(4)*pvc(4)
      K2K3 = pvb(1)*pvc(1)-pvb(2)*pvc(2)-pvb(3)*pvc(3)-pvb(4)*pvc(4)
      
      TV12 = rtwo*(ft(1,1)*va(1)*vb(1)+ft(2,2)*va(2)*vb(2)
     &+ft(3,3)*va(3)*vb(3)+ft(4,4)*va(4)*vb(4))


      TKV1 = rtwo*(ft(1,1)*p23(1)*va(1)+ft(2,2)*p23(2)*va(2)
     &+ft(3,3)*p23(3)*va(3)+ft(4,4)*p23(4)*va(4))

      TKV2 = rtwo*(ft(1,1)*p31(1)*vb(1)+ft(2,2)*p31(2)*vb(2)
     &+ft(3,3)*p31(3)*vb(3)+ft(4,4)*p31(4)*vb(4))

     
      TK3V1 = rtwo*(ft(1,1)*pvc(1)*va(1)+ft(2,2)*pvc(2)*va(2)
     &+ft(3,3)*pvc(3)*va(3)+ft(4,4)*pvc(4)*va(4))	

      TK3V2 = rtwo*(ft(1,1)*pvc(1)*vb(1)+ft(2,2)*pvc(2)*vb(2)
     &+ft(3,3)*pvc(3)*vb(3)+ft(4,4)*pvc(4)*vb(4))	


      TK312 = rtwo*(ft(1,1)*pvc(1)*p12(1)+ft(2,2)*pvc(2)*p12(2)
     &+ft(3,3)*pvc(3)*p12(3)+ft(4,4)*pvc(4)*p12(4))	


      TV12 = TV12 - T12*(va(1)*vb(2) + va(2)*vb(1))
     &          - T13*(va(1)*vb(3) + va(3)*vb(1))
     &          - T14*(va(1)*vb(4) + va(4)*vb(1))
     &          + T23*(va(2)*vb(3) + va(3)*vb(2))
     &          + T24*(va(2)*vb(4) + va(4)*vb(2))
     &          + T34*(va(3)*vb(4) + va(4)*vb(3))


      TKV1 = TKV1 - T12*(p23(1)*va(2) + p23(2)*va(1))
     &              - T13*(p23(1)*va(3) + p23(3)*va(1))
     &              - T14*(p23(1)*va(4) + p23(4)*va(1))
     &              + T23*(p23(2)*va(3) + p23(3)*va(2))
     &              + T24*(p23(2)*va(4) + p23(4)*va(2))
     &              + T34*(p23(3)*va(4) + p23(4)*va(3))

      TKV2 = TKV2 - T12*(p31(1)*vb(2) + p31(2)*vb(1))
     &              - T13*(p31(1)*vb(3) + p31(3)*vb(1))
     &              - T14*(p31(1)*vb(4) + p31(4)*vb(1))
     &              + T23*(p31(2)*vb(3) + p31(3)*vb(2))
     &              + T24*(p31(2)*vb(4) + p31(4)*vb(2))
     &              + T34*(p31(3)*vb(4) + p31(4)*vb(3))

      TK3V1 = TK3V1 - T12*(pvc(1)*va(2) + pvc(2)*va(1))
     &              - T13*(pvc(1)*va(3) + pvc(3)*va(1))
     &              - T14*(pvc(1)*va(4) + pvc(4)*va(1))
     &              + T23*(pvc(2)*va(3) + pvc(3)*va(2))
     &              + T24*(pvc(2)*va(4) + pvc(4)*va(2))
     &              + T34*(pvc(3)*va(4) + pvc(4)*va(3))

      TK3V2 = TK3V2 - T12*(pvc(1)*vb(2) + pvc(2)*vb(1))
     &              - T13*(pvc(1)*vb(3) + pvc(3)*vb(1))
     &              - T14*(pvc(1)*vb(4) + pvc(4)*vb(1))
     &              + T23*(pvc(2)*vb(3) + pvc(3)*vb(2))
     &              + T24*(pvc(2)*vb(4) + pvc(4)*vb(2))
     &              + T34*(pvc(3)*vb(4) + pvc(4)*vb(3))

      TK312 = TK312 - T12*(pvc(1)*p12(2) + pvc(2)*p12(1))
     &              - T13*(pvc(1)*p12(3) + pvc(3)*p12(1))
     &              - T14*(pvc(1)*p12(4) + pvc(4)*p12(1))
     &              + T23*(pvc(2)*p12(3) + pvc(3)*p12(2))
     &              + T24*(pvc(2)*p12(4) + pvc(4)*p12(2))
     &              + T34*(pvc(3)*p12(4) + pvc(4)*p12(3))

      do j=1,4

         TV1M(j) =
     &MET(j,1)*(ft(1,1)*va(1)-ft(2,1)*va(2)
     &-ft(3,1)*va(3)-ft(4,1)*va(4))
     &-MET(j,2)*(ft(1,2)*va(1)-ft(2,2)*va(2)
     &-ft(3,2)*va(3)-ft(4,2)*va(4))
     &-MET(j,3)*(ft(1,3)*va(1)-ft(2,3)*va(2)
     &-ft(3,3)*va(3)-ft(4,3)*va(4))
     &-MET(j,4)*(ft(1,4)*va(1)-ft(2,4)*va(2)
     &-ft(3,4)*va(3)-ft(4,4)*va(4))
     &+
     &MET(j,1)*(ft(1,1)*va(1)-ft(1,2)*va(2)
     &-ft(1,3)*va(3)-ft(1,4)*va(4))
     &-MET(j,2)*(ft(2,1)*va(1)-ft(2,2)*va(2)
     &-ft(2,3)*va(3)-ft(2,4)*va(4))
     &-MET(j,3)*(ft(3,1)*va(1)-ft(3,2)*va(2)
     &-ft(3,3)*va(3)-ft(3,4)*va(4))
     &-MET(j,4)*(ft(4,1)*va(1)-ft(4,2)*va(2)
     &-ft(4,3)*va(3)-ft(4,4)*va(4))

         TV2M(j) =
     &MET(j,1)*(ft(1,1)*vb(1)-ft(2,1)*vb(2)
     &-ft(3,1)*vb(3)-ft(4,1)*vb(4))
     &-MET(j,2)*(ft(1,2)*vb(1)-ft(2,2)*vb(2)
     &-ft(3,2)*vb(3)-ft(4,2)*vb(4))
     &-MET(j,3)*(ft(1,3)*vb(1)-ft(2,3)*vb(2)
     &-ft(3,3)*vb(3)-ft(4,3)*vb(4))
     &-MET(j,4)*(ft(1,4)*vb(1)-ft(2,4)*vb(2)
     &-ft(3,4)*vb(3)-ft(4,4)*vb(4))
     &+
     &MET(j,1)*(ft(1,1)*vb(1)-ft(1,2)*vb(2)
     &-ft(1,3)*vb(3)-ft(1,4)*vb(4))
     &-MET(j,2)*(ft(2,1)*vb(1)-ft(2,2)*vb(2)
     &-ft(2,3)*vb(3)-ft(2,4)*vb(4))
     &-MET(j,3)*(ft(3,1)*vb(1)-ft(3,2)*vb(2)
     &-ft(3,3)*vb(3)-ft(3,4)*vb(4))
     &-MET(j,4)*(ft(4,1)*vb(1)-ft(4,2)*vb(2)
     &-ft(4,3)*vb(3)-ft(4,4)*vb(4))


         TK12M(j) =
     &MET(j,1)*(ft(1,1)*p12(1)-ft(2,1)*p12(2)
     &-ft(3,1)*p12(3)-ft(4,1)*p12(4))
     &-MET(j,2)*(ft(1,2)*p12(1)-ft(2,2)*p12(2)
     &-ft(3,2)*p12(3)-ft(4,2)*p12(4))
     &-MET(j,3)*(ft(1,3)*p12(1)-ft(2,3)*p12(2)
     &-ft(3,3)*p12(3)-ft(4,3)*p12(4))
     &-MET(j,4)*(ft(1,4)*p12(1)-ft(2,4)*p12(2)
     &-ft(3,4)*p12(3)-ft(4,4)*p12(4))
     &+
     &MET(j,1)*(ft(1,1)*p12(1)-ft(1,2)*p12(2)
     &-ft(1,3)*p12(3)-ft(1,4)*p12(4))
     &-MET(j,2)*(ft(2,1)*p12(1)-ft(2,2)*p12(2)
     &-ft(2,3)*p12(3)-ft(2,4)*p12(4))
     &-MET(j,3)*(ft(3,1)*p12(1)-ft(3,2)*p12(2)
     &-ft(3,3)*p12(3)-ft(3,4)*p12(4))
     &-MET(j,4)*(ft(4,1)*p12(1)-ft(4,2)*p12(2)
     &-ft(4,3)*p12(3)-ft(4,4)*p12(4))
     
      enddo
        
      do i=1,4
         jvvt(i) = TV12*p12(i)+TKV1*vb(i)+TKV2*va(i)
     &+(-p12(i)*V1V2-vb(i)*K2V1+vb(i)*K3V1
     &  +va(i)*K1V2-va(i)*K3V2)*T00
     &-K1V2*TV1M(i)+K3V2*TV1M(i)
     &+V1V2*TK12M(i)
     &+K2V1*TV2M(i)-K3V1*TV2M(i)
      enddo 

      if ( vmass.gt.rZero ) then
         do i=1,4
            jvvt(i) = jvvt(i)
     &+(K1V2*pvc(i)
     &-K3V2*TKV1*pvc(i)-K3V2*TK3V1*pvc(i)
     &-K1K3*TV12*pvc(i)+K2K3*TV12*pvc(i)
     &-V1V2*TK312*pvc(i)
     &-K2V1*TK3V2*pvc(i)+K3V1*TK3V2*pvc(i)
     &-K3V1*TKV2*pvc(i))/vmass**2
     &+(-pvc(i)*K1V2*K3V1+pvc(i)*K2V1*K3V2
     &+pvc(i)*V1V2*K1K3-pvc(i)*V1V2*K2K3
     &)/vmass**2*T00
         enddo
           
      endif	
      
      do i=1,4
         jvvt(i) = -jvvt(i) * d*gt
      enddo
      
      return
      end
      subroutine jvvxxx(v1,v2,g,vmass,vwidth , jvv)
c
c This subroutine computes an off-shell vector current from the three-
c point gauge boson coupling.  The vector propagator is given in Feynman
c gauge for a massless vector and in unitary gauge for a massive vector.
c
c input:
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       real    g              : coupling constant (see the table below)
c       real    vmass          : mass  of output vector v
c       real    vwidth         : width of output vector v
c
c the possible sets of the inputs are as follows:
c    ------------------------------------------------------------------
c    |   v1   |   v2   |  jvv   |      g       |   vmass  |  vwidth   |
c    ------------------------------------------------------------------
c    |   W-   |   W+   |  A/Z   |  gwwa/gwwz   | 0./zmass | 0./zwidth |
c    | W3/A/Z |   W-   |  W+    | gw/gwwa/gwwz |   wmass  |  wwidth   |
c    |   W+   | W3/A/Z |  W-    | gw/gwwa/gwwz |   wmass  |  wwidth   |
c    ------------------------------------------------------------------
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex jvv(6)         : vector current            j^mu(v:v1,v2)
c     
      implicit none
      double complex v1(6),v2(6),jvv(6),j12(0:3),js,dg
      double complex sv1,sv2,s11,s12,s21,s22,v12
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision p1(0:3),p2(0:3),q(0:3),g,vmass,vwidth,gs,s
      double precision vm2,m1,m2

      double precision rZero
      parameter( rZero = 0.0d0 )
 
c

      jvv(5) = v1(5)+v2(5)
      jvv(6) = v1(6)+v2(6)

      p1(0) =  dble( v1(5))
      p1(1) =  dble( v1(6))
      p1(2) =  dimag(v1(6))
      p1(3) =  dimag(v1(5))
      p2(0) =  dble( v2(5))
      p2(1) =  dble( v2(6))
      p2(2) =  dimag(v2(6))
      p2(3) =  dimag(v2(5))
      q(0)  = -dble( jvv(5))
      q(1)  = -dble( jvv(6))
      q(2)  = -dimag(jvv(6))
      q(3)  = -dimag(jvv(5))
      s = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2


      v12 = v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4)
      sv1 =   (p2(0)-q(0))*v1(1) -(p2(1)-q(1))*v1(2)
     &      - (p2(2)-q(2))*v1(3) -(p2(3)-q(3))*v1(4)
      sv2 = - (p1(0)-q(0))*v2(1) +(p1(1)-q(1))*v2(2)
     &      + (p1(2)-q(2))*v2(3) +(p1(3)-q(3))*v2(4)
      j12(0) = (p1(0)-p2(0))*v12 +sv1*v2(1) +sv2*v1(1)
      j12(1) = (p1(1)-p2(1))*v12 +sv1*v2(2) +sv2*v1(2)
      j12(2) = (p1(2)-p2(2))*v12 +sv1*v2(3) +sv2*v1(3)
      j12(3) = (p1(3)-p2(3))*v12 +sv1*v2(4) +sv2*v1(4)

      if ( vmass.ne.rZero ) then

         m1 = p1(0)**2-(p1(1)**2+p1(2)**2+p1(3)**2)
         m2 = p2(0)**2-(p2(1)**2+p2(2)**2+p2(3)**2)
         s11 = p1(0)*v1(1)-p1(1)*v1(2)-p1(2)*v1(3)-p1(3)*v1(4)
         s12 = p1(0)*v2(1)-p1(1)*v2(2)-p1(2)*v2(3)-p1(3)*v2(4)
         s21 = p2(0)*v1(1)-p2(1)*v1(2)-p2(2)*v1(3)-p2(3)*v1(4)
         s22 = p2(0)*v2(1)-p2(1)*v2(2)-p2(2)*v2(3)-p2(3)*v2(4)

c     Fabio's implementation of the fixed width
         cm2=dcmplx( vm2, -vmass*vwidth )
c     js = (v12*(-m1+m2) +s11*s12 -s21*s22)/vm2
         js = (v12*(-m1+m2) +s11*s12 -s21*s22)/cm2
        
         dg = -g/dcmplx( s-vm2, vmass*vwidth )

c  For the running width, use below instead of the above dg.
c         dg = -g/dcmplx( s-vm2, max(vwidth*s/vmass,rZero) )

         jvv(1) = dg*(j12(0)-q(0)*js)
         jvv(2) = dg*(j12(1)-q(1)*js)
         jvv(3) = dg*(j12(2)-q(2)*js)
         jvv(4) = dg*(j12(3)-q(3)*js)

      else

         gs = -g/s

         jvv(1) = gs*j12(0)
         jvv(2) = gs*j12(1)
         jvv(3) = gs*j12(2)
         jvv(4) = gs*j12(3)

      end if
c
      return
      end
      subroutine jw3wnx(w1,w2,w3,g1,g2,vmass,vwidth , jw3w)
c
c This subroutine computes an off-shell W+, W-, W3, Z or photon current
c from the four-point gauge boson coupling.  The vector propagator is
c given in Feynman gauge for a photon and in unitary gauge for W and
c Z bosons.  If one sets wmass=0.0, then the ggg-->g current is given
c (see sect 2.9.1 of the manual).
c
c input:
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       complex w3(6)          : third  vector                        w3
c       real    g1             : first  coupling constant
c       real    g2             : second coupling constant
c                                                  (see the table below)
c       real    wmass          : mass  of internal W
c       real    wwidth         : width of internal W
c       real    vmass          : mass  of output W'
c       real    vwidth         : width of output W'
c
c the possible sets of the inputs are as follows:
c   -------------------------------------------------------------------
c   |  w1  |  w2  |  w3  | g1 | g2 |wmass|wwidth|vmass|vwidth || jw3w |
c   -------------------------------------------------------------------
c   |  W-  |  W3  |  W+  | gw |gwwz|wmass|wwidth|zmass|zwidth ||  Z   |
c   |  W-  |  W3  |  W+  | gw |gwwa|wmass|wwidth|  0. |  0.   ||  A   |
c   |  W-  |  Z   |  W+  |gwwz|gwwz|wmass|wwidth|zmass|zwidth ||  Z   |
c   |  W-  |  Z   |  W+  |gwwz|gwwa|wmass|wwidth|  0. |  0.   ||  A   |
c   |  W-  |  A   |  W+  |gwwa|gwwz|wmass|wwidth|zmass|zwidth ||  Z   |
c   |  W-  |  A   |  W+  |gwwa|gwwa|wmass|wwidth|  0. |  0.   ||  A   |
c   -------------------------------------------------------------------
c   |  W3  |  W-  |  W3  | gw | gw |wmass|wwidth|wmass|wwidth ||  W+  |
c   |  W3  |  W+  |  W3  | gw | gw |wmass|wwidth|wmass|wwidth ||  W-  |
c   |  W3  |  W-  |  Z   | gw |gwwz|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  W3  |  W+  |  Z   | gw |gwwz|wmass|wwidth|wmass|wwidth ||  W-  |
c   |  W3  |  W-  |  A   | gw |gwwa|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  W3  |  W+  |  A   | gw |gwwa|wmass|wwidth|wmass|wwidth ||  W-  |
c   |  Z   |  W-  |  Z   |gwwz|gwwz|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  Z   |  W+  |  Z   |gwwz|gwwz|wmass|wwidth|wmass|wwidth ||  W-  |
c   |  Z   |  W-  |  A   |gwwz|gwwa|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  Z   |  W+  |  A   |gwwz|gwwa|wmass|wwidth|wmass|wwidth ||  W-  |
c   |  A   |  W-  |  A   |gwwa|gwwa|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  A   |  W+  |  A   |gwwa|gwwa|wmass|wwidth|wmass|wwidth ||  W-  |
c   -------------------------------------------------------------------
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex jw3w(6)        : W current             j^mu(w':w1,w2,w3)
c     
      implicit none
      double complex w1(6),w2(6),w3(6),jw3w(6)
      double complex dw1(0:3),dw2(0:3),dw3(0:3)
      double complex jj(0:3),j4(0:3),dv,w12,w32,w13,jq
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision g1,g2,vmass,vwidth
      double precision p1(0:3),p2(0:3),p3(0:3),q(0:3)
      double precision dg2,dmv,dwv,mv2,q2

      double precision rZero, rOne, rTwo
      parameter( rZero = 0.0d0, rOne = 1.0d0, rTwo = 2.0d0 )

c

      jw3w(5) = w1(5)+w2(5)+w3(5)
      jw3w(6) = w1(6)+w2(6)+w3(6)

      dw1(0) = dcmplx(w1(1))
      dw1(1) = dcmplx(w1(2))
      dw1(2) = dcmplx(w1(3))
      dw1(3) = dcmplx(w1(4))
      dw2(0) = dcmplx(w2(1))
      dw2(1) = dcmplx(w2(2))
      dw2(2) = dcmplx(w2(3))
      dw2(3) = dcmplx(w2(4))
      dw3(0) = dcmplx(w3(1))
      dw3(1) = dcmplx(w3(2))
      dw3(2) = dcmplx(w3(3))
      dw3(3) = dcmplx(w3(4))
      p1(0) = dble(      w1(5))
      p1(1) = dble(      w1(6))
      p1(2) = dble(dimag(w1(6)))
      p1(3) = dble(dimag(w1(5)))
      p2(0) = dble(      w2(5))
      p2(1) = dble(      w2(6))
      p2(2) = dble(dimag(w2(6)))
      p2(3) = dble(dimag(w2(5)))
      p3(0) = dble(      w3(5))
      p3(1) = dble(      w3(6))
      p3(2) = dble(dimag(w3(6)))
      p3(3) = dble(dimag(w3(5)))
      q(0) = -(p1(0)+p2(0)+p3(0))
      q(1) = -(p1(1)+p2(1)+p3(1))
      q(2) = -(p1(2)+p2(2)+p3(2))
      q(3) = -(p1(3)+p2(3)+p3(3))

      q2 = q(0)**2 -(q(1)**2 +q(2)**2 +q(3)**2)
c     Neil edited this to allow 3-site couplings.
c     Now only g1 is important.  g2 does nothing.
c      dg2 = dble(g1)*dble(g2)
c      dg2 = dble(g1)
c     End Neil's edit.

c Benj's modif in order to have FR running
      if(g1.eq.rzero) dg2=g2
      if(g2.eq.rzero) dg2=g1
c End of Benj'S modif

      dmv = dble(vmass)
      dwv = dble(vwidth)
      mv2 = dmv**2


      if ( vmass.eq.rZero ) then
         dv = rOne/dcmplx( q2 )
      else
         dv = rOne/dcmplx( q2-mv2, dmv*dwv )
      endif

c  For the running width, use below instead of the above dv.
c      dv = rOne/dcmplx( q2-mv2 , max(dwv*q2/dmv,rZero) )

      w12=dw1(0)*dw2(0)-dw1(1)*dw2(1)-dw1(2)*dw2(2)-dw1(3)*dw2(3)
      w32=dw3(0)*dw2(0)-dw3(1)*dw2(1)-dw3(2)*dw2(2)-dw3(3)*dw2(3)

      w13=dw1(0)*dw3(0)-dw1(1)*dw3(1)-dw1(2)*dw3(2)-dw1(3)*dw3(3)
      
      j4(0) = dg2*( dw1(0)*w32 + dw3(0)*w12 - rTwo*dw2(0)*w13 )
      j4(1) = dg2*( dw1(1)*w32 + dw3(1)*w12 - rTwo*dw2(1)*w13 )
      j4(2) = dg2*( dw1(2)*w32 + dw3(2)*w12 - rTwo*dw2(2)*w13 )
      j4(3) = dg2*( dw1(3)*w32 + dw3(3)*w12 - rTwo*dw2(3)*w13 )

      jj(0) = j4(0)
      jj(1) = j4(1)
      jj(2) = j4(2)
      jj(3) = j4(3)

      if ( vmass.ne.rZero ) then

c     Fabio's implementation of the fixed width
         cm2=dcmplx( mv2, -dmv*dwv)
c     jq = (jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/mv2
         jq = (jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/cm2
         
         jw3w(1) = dcmplx( (jj(0)-jq*q(0))*dv )
         jw3w(2) = dcmplx( (jj(1)-jq*q(1))*dv )
         jw3w(3) = dcmplx( (jj(2)-jq*q(2))*dv )
         jw3w(4) = dcmplx( (jj(3)-jq*q(3))*dv )

      else

         jw3w(1) = dcmplx( jj(0)*dv )
         jw3w(2) = dcmplx( jj(1)*dv )
         jw3w(3) = dcmplx( jj(2)*dv )
         jw3w(4) = dcmplx( jj(3)*dv )
      end if
c
      return
      end
      subroutine jw3wxx(w1,w2,w3,g1,g2,vmass,vwidth , jw3w)
c
c This subroutine computes an off-shell W+, W-, W3, Z or photon current
c from the four-point gauge boson coupling.  The vector propagator is
c given in Feynman gauge for a photon and in unitary gauge for W and
c Z bosons.  If one sets wmass=0.0, then the ggg-->g current is given
c (see sect 2.9.1 of the manual).
c
c input:
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       complex w3(6)          : third  vector                        w3
c       real    g1             : first  coupling constant
c       real    g2             : second coupling constant
c                                                  (see the table below)
c       real    wmass          : mass  of internal W
c       real    wwidth         : width of internal W
c       real    vmass          : mass  of output W'
c       real    vwidth         : width of output W'
c
c the possible sets of the inputs are as follows:
c   -------------------------------------------------------------------
c   |  w1  |  w2  |  w3  | g1 | g2 |wmass|wwidth|vmass|vwidth || jw3w |
c   -------------------------------------------------------------------
c   |  W-  |  W3  |  W+  | gw |gwwz|wmass|wwidth|zmass|zwidth ||  Z   |
c   |  W-  |  W3  |  W+  | gw |gwwa|wmass|wwidth|  0. |  0.   ||  A   |
c   |  W-  |  Z   |  W+  |gwwz|gwwz|wmass|wwidth|zmass|zwidth ||  Z   |
c   |  W-  |  Z   |  W+  |gwwz|gwwa|wmass|wwidth|  0. |  0.   ||  A   |
c   |  W-  |  A   |  W+  |gwwa|gwwz|wmass|wwidth|zmass|zwidth ||  Z   |
c   |  W-  |  A   |  W+  |gwwa|gwwa|wmass|wwidth|  0. |  0.   ||  A   |
c   -------------------------------------------------------------------
c   |  W3  |  W-  |  W3  | gw | gw |wmass|wwidth|wmass|wwidth ||  W+  |
c   |  W3  |  W+  |  W3  | gw | gw |wmass|wwidth|wmass|wwidth ||  W-  |
c   |  W3  |  W-  |  Z   | gw |gwwz|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  W3  |  W+  |  Z   | gw |gwwz|wmass|wwidth|wmass|wwidth ||  W-  |
c   |  W3  |  W-  |  A   | gw |gwwa|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  W3  |  W+  |  A   | gw |gwwa|wmass|wwidth|wmass|wwidth ||  W-  |
c   |  Z   |  W-  |  Z   |gwwz|gwwz|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  Z   |  W+  |  Z   |gwwz|gwwz|wmass|wwidth|wmass|wwidth ||  W-  |
c   |  Z   |  W-  |  A   |gwwz|gwwa|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  Z   |  W+  |  A   |gwwz|gwwa|wmass|wwidth|wmass|wwidth ||  W-  |
c   |  A   |  W-  |  A   |gwwa|gwwa|wmass|wwidth|wmass|wwidth ||  W+  |
c   |  A   |  W+  |  A   |gwwa|gwwa|wmass|wwidth|wmass|wwidth ||  W-  |
c   -------------------------------------------------------------------
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex jw3w(6)        : W current             j^mu(w':w1,w2,w3)
c     
      implicit none
      double complex w1(6),w2(6),w3(6),jw3w(6)
      double complex dw1(0:3),dw2(0:3),dw3(0:3)
      double complex jj(0:3),j4(0:3),dv,w12,w32,w13,jq
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision g1,g2,vmass,vwidth
      double precision p1(0:3),p2(0:3),p3(0:3),q(0:3)
      double precision dg2,dmv,dwv,mv2,q2

      double precision rZero, rOne, rTwo
      parameter( rZero = 0.0d0, rOne = 1.0d0, rTwo = 2.0d0 )

c

      jw3w(5) = w1(5)+w2(5)+w3(5)
      jw3w(6) = w1(6)+w2(6)+w3(6)

      dw1(0) = dcmplx(w1(1))
      dw1(1) = dcmplx(w1(2))
      dw1(2) = dcmplx(w1(3))
      dw1(3) = dcmplx(w1(4))
      dw2(0) = dcmplx(w2(1))
      dw2(1) = dcmplx(w2(2))
      dw2(2) = dcmplx(w2(3))
      dw2(3) = dcmplx(w2(4))
      dw3(0) = dcmplx(w3(1))
      dw3(1) = dcmplx(w3(2))
      dw3(2) = dcmplx(w3(3))
      dw3(3) = dcmplx(w3(4))
      p1(0) = dble(      w1(5))
      p1(1) = dble(      w1(6))
      p1(2) = dble(dimag(w1(6)))
      p1(3) = dble(dimag(w1(5)))
      p2(0) = dble(      w2(5))
      p2(1) = dble(      w2(6))
      p2(2) = dble(dimag(w2(6)))
      p2(3) = dble(dimag(w2(5)))
      p3(0) = dble(      w3(5))
      p3(1) = dble(      w3(6))
      p3(2) = dble(dimag(w3(6)))
      p3(3) = dble(dimag(w3(5)))
      q(0) = -(p1(0)+p2(0)+p3(0))
      q(1) = -(p1(1)+p2(1)+p3(1))
      q(2) = -(p1(2)+p2(2)+p3(2))
      q(3) = -(p1(3)+p2(3)+p3(3))

      q2 = q(0)**2 -(q(1)**2 +q(2)**2 +q(3)**2)
      dg2 = dble(g1)*dble(g2)
      dmv = dble(vmass)
      dwv = dble(vwidth)
      mv2 = dmv**2


      if ( vmass.eq.rZero ) then
         dv = rOne/dcmplx( q2 )
      else
         dv = rOne/dcmplx( q2-mv2, dmv*dwv )
      endif

c  For the running width, use below instead of the above dv.
c      dv = rOne/dcmplx( q2-mv2 , max(dwv*q2/dmv,rZero) )

      w12=dw1(0)*dw2(0)-dw1(1)*dw2(1)-dw1(2)*dw2(2)-dw1(3)*dw2(3)
      w32=dw3(0)*dw2(0)-dw3(1)*dw2(1)-dw3(2)*dw2(2)-dw3(3)*dw2(3)

      w13=dw1(0)*dw3(0)-dw1(1)*dw3(1)-dw1(2)*dw3(2)-dw1(3)*dw3(3)
      
      j4(0) = dg2*( dw1(0)*w32 + dw3(0)*w12 - rTwo*dw2(0)*w13 )
      j4(1) = dg2*( dw1(1)*w32 + dw3(1)*w12 - rTwo*dw2(1)*w13 )
      j4(2) = dg2*( dw1(2)*w32 + dw3(2)*w12 - rTwo*dw2(2)*w13 )
      j4(3) = dg2*( dw1(3)*w32 + dw3(3)*w12 - rTwo*dw2(3)*w13 )

      jj(0) = j4(0)
      jj(1) = j4(1)
      jj(2) = j4(2)
      jj(3) = j4(3)

      if ( vmass.ne.rZero ) then

c     Fabio's implementation of the fixed width
         cm2=dcmplx( mv2, -dmv*dwv)
c     jq = (jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/mv2
         jq = (jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/cm2
         
         jw3w(1) = dcmplx( (jj(0)-jq*q(0))*dv )
         jw3w(2) = dcmplx( (jj(1)-jq*q(1))*dv )
         jw3w(3) = dcmplx( (jj(2)-jq*q(2))*dv )
         jw3w(4) = dcmplx( (jj(3)-jq*q(3))*dv )

      else

         jw3w(1) = dcmplx( jj(0)*dv )
         jw3w(2) = dcmplx( jj(1)*dv )
         jw3w(3) = dcmplx( jj(2)*dv )
         jw3w(4) = dcmplx( jj(3)*dv )
      end if
c
      return
      end
      subroutine jwwwnx(w1,w2,w3,gwwa,gwwz,wmass,wwidth , jwww)
c
c This subroutine computes an off-shell W+/W- current from the four-
c point gauge boson coupling.  The vector propagators for the output
c W and the internal Z bosons are given in unitary gauge, and that of
c the internal photon is given in Feynman gauge.
c
c input:
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       complex w3(6)          : third  vector                        w3
c       real    gwwa           : coupling constant of W and A       gwwa
c       real    gwwz           : coupling constant of W and Z       gwwz
c       real    zmass          : mass  of internal Z
c       real    zwidth         : width of internal Z
c       real    wmass          : mass  of output W
c       real    wwidth         : width of output W
c
c the possible sets of the inputs are as follows:
c   -------------------------------------------------------------------
c   |  w1  |  w2  |  w3  |gwwa|gwwz|zmass|zwidth|wmass|wwidth || jwww |
c   -------------------------------------------------------------------
c   |  W-  |  W+  |  W-  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  W+  |
c   |  W+  |  W-  |  W+  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  W-  |
c   -------------------------------------------------------------------
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex jwww(6)        : W current             j^mu(w':w1,w2,w3)
c     
      implicit none
      double complex w1(6),w2(6),w3(6),jwww(6)
      double complex dw1(0:3),dw2(0:3),dw3(0:3),jj(0:3)
      double complex dw,w12,w32,w13,jq
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision gwwa,gwwz,wmass,wwidth
      double precision p1(0:3),p2(0:3),p3(0:3),q(0:3)
      double precision dgwwa2,dgwwz2,dgw2,dmw,dww,mw2,q2

      double precision rZero, rOne, rTwo
      parameter( rZero = 0.0d0, rOne = 1.0d0, rTwo = 2.0d0 )

c



      
      jwww(5) = w1(5)+w2(5)+w3(5)
      jwww(6) = w1(6)+w2(6)+w3(6)

      dw1(0) = dcmplx(w1(1))
      dw1(1) = dcmplx(w1(2))
      dw1(2) = dcmplx(w1(3))
      dw1(3) = dcmplx(w1(4))
      dw2(0) = dcmplx(w2(1))
      dw2(1) = dcmplx(w2(2))
      dw2(2) = dcmplx(w2(3))
      dw2(3) = dcmplx(w2(4))
      dw3(0) = dcmplx(w3(1))
      dw3(1) = dcmplx(w3(2))
      dw3(2) = dcmplx(w3(3))
      dw3(3) = dcmplx(w3(4))
      p1(0) = dble(      w1(5))
      p1(1) = dble(      w1(6))
      p1(2) = dble(dimag(w1(6)))
      p1(3) = dble(dimag(w1(5)))
      p2(0) = dble(      w2(5))
      p2(1) = dble(      w2(6))
      p2(2) = dble(dimag(w2(6)))
      p2(3) = dble(dimag(w2(5)))
      p3(0) = dble(      w3(5))
      p3(1) = dble(      w3(6))
      p3(2) = dble(dimag(w3(6)))
      p3(3) = dble(dimag(w3(5)))
      q(0) = -(p1(0)+p2(0)+p3(0))
      q(1) = -(p1(1)+p2(1)+p3(1))
      q(2) = -(p1(2)+p2(2)+p3(2))
      q(3) = -(p1(3)+p2(3)+p3(3))
      q2 = q(0)**2 -(q(1)**2 +q(2)**2 +q(3)**2)
c     Neil edited this file to allow implementation of 3-site model.
c     Now, only gwwa matters and is the full coupling squared.
c      dgwwa2 = dble(gwwa)**2
c      dgwwz2 = dble(gwwz)**2
c      dgw2 = dgwwa2+dgwwz2
c      dgw2 = gwwa
c     End Neil's edit

c Benj's modif in order to have FR running
      dgw2=rZero
      if(gwwa.eq.rZero) dgw2=gwwz
      if(gwwz.eq.rZero) dgw2=gwwa
c End of Benj'S modif



      dmw = dble(wmass)
      dww = dble(wwidth)
      mw2 = dmw**2


c Start of Claude'S modif (Removed minus Sign)
      dw = rOne/dcmplx( q2-mw2, dmw*dww )
c End of Claude'S modif

c  For the running width, use below instead of the above dw.
c      dw = -rOne/dcmplx( q2-mw2 , max(dww*q2/dmw,rZero) )

      w12=dw1(0)*dw2(0)-dw1(1)*dw2(1)-dw1(2)*dw2(2)-dw1(3)*dw2(3)
      w32=dw3(0)*dw2(0)-dw3(1)*dw2(1)-dw3(2)*dw2(2)-dw3(3)*dw2(3)

      w13=dw1(0)*dw3(0)-dw1(1)*dw3(1)-dw1(2)*dw3(2)-dw1(3)*dw3(3)

      jj(0) = dgw2*( dw1(0)*w32 + dw3(0)*w12 - rTwo*dw2(0)*w13 )
      jj(1) = dgw2*( dw1(1)*w32 + dw3(1)*w12 - rTwo*dw2(1)*w13 )
      jj(2) = dgw2*( dw1(2)*w32 + dw3(2)*w12 - rTwo*dw2(2)*w13 )
      jj(3) = dgw2*( dw1(3)*w32 + dw3(3)*w12 - rTwo*dw2(3)*w13 )

c     Fabio's implementation of the fixed width
      cm2=dcmplx( mw2, -dmw*dww )
c     jq = (jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/mw2
      jq = (jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/cm2

      jwww(1) = dcmplx( (jj(0)-jq*q(0))*dw )
      jwww(2) = dcmplx( (jj(1)-jq*q(1))*dw )
      jwww(3) = dcmplx( (jj(2)-jq*q(2))*dw )
      jwww(4) = dcmplx( (jj(3)-jq*q(3))*dw )
c
      return
      end
      subroutine jwwwxx(w1,w2,w3,gwwa,gwwz,wmass,wwidth , jwww)
c
c This subroutine computes an off-shell W+/W- current from the four-
c point gauge boson coupling.  The vector propagators for the output
c W and the internal Z bosons are given in unitary gauge, and that of
c the internal photon is given in Feynman gauge.
c
c input:
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       complex w3(6)          : third  vector                        w3
c       real    gwwa           : coupling constant of W and A       gwwa
c       real    gwwz           : coupling constant of W and Z       gwwz
c       real    zmass          : mass  of internal Z
c       real    zwidth         : width of internal Z
c       real    wmass          : mass  of output W
c       real    wwidth         : width of output W
c
c the possible sets of the inputs are as follows:
c   -------------------------------------------------------------------
c   |  w1  |  w2  |  w3  |gwwa|gwwz|zmass|zwidth|wmass|wwidth || jwww |
c   -------------------------------------------------------------------
c   |  W-  |  W+  |  W-  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  W+  |
c   |  W+  |  W-  |  W+  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  W-  |
c   -------------------------------------------------------------------
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex jwww(6)        : W current             j^mu(w':w1,w2,w3)
c     
      implicit none
      double complex w1(6),w2(6),w3(6),jwww(6)
      double complex dw1(0:3),dw2(0:3),dw3(0:3),jj(0:3)
      double complex dw,w12,w32,w13,jq
      double complex cm2        ! mass**2- I Gamma mass (Fabio)
      double precision gwwa,gwwz,wmass,wwidth
      double precision p1(0:3),p2(0:3),p3(0:3),q(0:3)
      double precision dgwwa2,dgwwz2,dgw2,dmw,dww,mw2,q2

      double precision rZero, rOne, rTwo
      parameter( rZero = 0.0d0, rOne = 1.0d0, rTwo = 2.0d0 )

c

      jwww(5) = w1(5)+w2(5)+w3(5)
      jwww(6) = w1(6)+w2(6)+w3(6)

      dw1(0) = dcmplx(w1(1))
      dw1(1) = dcmplx(w1(2))
      dw1(2) = dcmplx(w1(3))
      dw1(3) = dcmplx(w1(4))
      dw2(0) = dcmplx(w2(1))
      dw2(1) = dcmplx(w2(2))
      dw2(2) = dcmplx(w2(3))
      dw2(3) = dcmplx(w2(4))
      dw3(0) = dcmplx(w3(1))
      dw3(1) = dcmplx(w3(2))
      dw3(2) = dcmplx(w3(3))
      dw3(3) = dcmplx(w3(4))
      p1(0) = dble(      w1(5))
      p1(1) = dble(      w1(6))
      p1(2) = dble(dimag(w1(6)))
      p1(3) = dble(dimag(w1(5)))
      p2(0) = dble(      w2(5))
      p2(1) = dble(      w2(6))
      p2(2) = dble(dimag(w2(6)))
      p2(3) = dble(dimag(w2(5)))
      p3(0) = dble(      w3(5))
      p3(1) = dble(      w3(6))
      p3(2) = dble(dimag(w3(6)))
      p3(3) = dble(dimag(w3(5)))
      q(0) = -(p1(0)+p2(0)+p3(0))
      q(1) = -(p1(1)+p2(1)+p3(1))
      q(2) = -(p1(2)+p2(2)+p3(2))
      q(3) = -(p1(3)+p2(3)+p3(3))
      q2 = q(0)**2 -(q(1)**2 +q(2)**2 +q(3)**2)
      dgwwa2 = dble(gwwa)**2
      dgwwz2 = dble(gwwz)**2
      dgw2 = dgwwa2+dgwwz2
      dmw = dble(wmass)
      dww = dble(wwidth)
      mw2 = dmw**2


      dw = -rOne/dcmplx( q2-mw2, dmw*dww )
c  For the running width, use below instead of the above dw.
c      dw = -rOne/dcmplx( q2-mw2 , max(dww*q2/dmw,rZero) )

      w12=dw1(0)*dw2(0)-dw1(1)*dw2(1)-dw1(2)*dw2(2)-dw1(3)*dw2(3)
      w32=dw3(0)*dw2(0)-dw3(1)*dw2(1)-dw3(2)*dw2(2)-dw3(3)*dw2(3)

      w13=dw1(0)*dw3(0)-dw1(1)*dw3(1)-dw1(2)*dw3(2)-dw1(3)*dw3(3)

      jj(0) = dgw2*( dw1(0)*w32 + dw3(0)*w12 - rTwo*dw2(0)*w13 )
      jj(1) = dgw2*( dw1(1)*w32 + dw3(1)*w12 - rTwo*dw2(1)*w13 )
      jj(2) = dgw2*( dw1(2)*w32 + dw3(2)*w12 - rTwo*dw2(2)*w13 )
      jj(3) = dgw2*( dw1(3)*w32 + dw3(3)*w12 - rTwo*dw2(3)*w13 )

c     Fabio's implementation of the fixed width
      cm2=dcmplx( mw2, -dmw*dww )
c     jq = (jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/mw2
      jq = (jj(0)*q(0)-jj(1)*q(1)-jj(2)*q(2)-jj(3)*q(3))/cm2

      jwww(1) = dcmplx( (jj(0)-jq*q(0))*dw )
      jwww(2) = dcmplx( (jj(1)-jq*q(1))*dw )
      jwww(3) = dcmplx( (jj(2)-jq*q(2))*dw )
      jwww(4) = dcmplx( (jj(3)-jq*q(3))*dw )
c
      return
      end
      subroutine mom2cx(esum,mass1,mass2,costh1,phi1 , p1,p2)
c
c This subroutine sets up two four-momenta in the two particle rest
c frame.
c
c input:
c       real    esum           : energy sum of particle 1 and 2
c       real    mass1          : mass            of particle 1
c       real    mass2          : mass            of particle 2
c       real    costh1         : cos(theta)      of particle 1
c       real    phi1           : azimuthal angle of particle 1
c
c output:
c       real    p1(0:3)        : four-momentum of particle 1
c       real    p2(0:3)        : four-momentum of particle 2
c     
      implicit none
      double precision p1(0:3),p2(0:3),
     &     esum,mass1,mass2,costh1,phi1,md2,ed,pp,sinth1

      double precision rZero, rHalf, rOne, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0 )

c

      md2 = (mass1-mass2)*(mass1+mass2)
      ed = md2/esum
      if ( mass1*mass2.eq.rZero ) then
         pp = (esum-abs(ed))*rHalf
      else
         pp = sqrt((md2/esum)**2-rTwo*(mass1**2+mass2**2)+esum**2)*rHalf
      endif
      sinth1 = sqrt((rOne-costh1)*(rOne+costh1))

      p1(0) = max((esum+ed)*rHalf,rZero)
      p1(1) = pp*sinth1*cos(phi1)
      p1(2) = pp*sinth1*sin(phi1)
      p1(3) = pp*costh1

      p2(0) = max((esum-ed)*rHalf,rZero)
      p2(1) = -p1(1)
      p2(2) = -p1(2)
      p2(3) = -p1(3)
c
      return
      end
      subroutine momntx(energy,mass,costh,phi , p)
c
c This subroutine sets up a four-momentum from the four inputs.
c
c input:
c       real    energy         : energy
c       real    mass           : mass
c       real    costh          : cos(theta)
c       real    phi            : azimuthal angle
c
c output:
c       real    p(0:3)         : four-momentum
c     
      implicit none
      double precision p(0:3),energy,mass,costh,phi,pp,sinth

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )

c

      p(0) = energy

      if ( energy.eq.mass ) then

         p(1) = rZero
         p(2) = rZero
         p(3) = rZero

      else

         pp = sqrt((energy-mass)*(energy+mass))
         sinth = sqrt((rOne-costh)*(rOne+costh))
         p(3) = pp*costh
         if ( phi.eq.rZero ) then
            p(1) = pp*sinth
            p(2) = rZero
         else
            p(1) = pp*sinth*cos(phi)
            p(2) = pp*sinth*sin(phi)
         endif

      endif
c
      return
      end
      subroutine oxxxxx(p,fmass,nhel,nsf , fo)
c
c This subroutine computes a fermion wavefunction with the flowing-OUT
c fermion number.
c
c input:
c       real    p(0:3)         : four-momentum of fermion
c       real    fmass          : mass          of fermion
c       integer nhel = -1 or 1 : helicity      of fermion
c       integer nsf  = -1 or 1 : +1 for particle, -1 for anti-particle
c
c output:
c       complex fo(6)          : fermion wavefunction               <fo|
c     
      implicit none
      double complex fo(6),chi(2)
      double precision p(0:3),sf(2),sfomeg(2),omega(2),fmass,
     &     pp,pp3,sqp0p3,sqm(0:1)
      integer nhel,nsf,nh,ip,im

      double precision rZero, rHalf, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0, rTwo = 2.0d0 )

c

      fo(5) = dcmplx(p(0),p(3))*nsf
      fo(6) = dcmplx(p(1),p(2))*nsf

      nh = nhel*nsf

      if ( fmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))

         if ( pp.eq.rZero ) then
            
            sqm(0) = dsqrt(abs(fmass)) ! possibility of negative fermion masses
            sqm(1) = sign(sqm(0),fmass) ! possibility of negative fermion masses
            ip = -((1+nh)/2)
            im =  (1-nh)/2
            
            fo(1) = im     * sqm(im)
            fo(2) = ip*nsf * sqm(im)
            fo(3) = im*nsf * sqm(-ip)
            fo(4) = ip     * sqm(-ip)
            
         else
            
            pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))
            sf(1) = dble(1+nsf+(1-nsf)*nh)*rHalf
            sf(2) = dble(1+nsf-(1-nsf)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = fmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , -p(2) )/dsqrt(rTwo*pp*pp3)
            endif
            
            fo(1) = sfomeg(2)*chi(im)
            fo(2) = sfomeg(2)*chi(ip)
            fo(3) = sfomeg(1)*chi(im)
            fo(4) = sfomeg(1)*chi(ip)

         endif
         
      else
         
         if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
            sqp0p3 = 0d0
         else
            sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsf
         end if
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), -p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fo(1) = chi(1)
            fo(2) = chi(2)
            fo(3) = dcmplx( rZero )
            fo(4) = dcmplx( rZero )
         else
            fo(1) = dcmplx( rZero )
            fo(2) = dcmplx( rZero )
            fo(3) = chi(2)
            fo(4) = chi(1)
         endif
         
      endif
c
      return
      end
      subroutine pxxxxx(p,tmass,nhel,nst , tc)

c    CP3 2009.NOV

c This subroutine computes a PSEUDOR wavefunction.
c
c input:
c       real    p(0:3)         : four-momentum of tensor boson
c       real    tmass          : mass          of tensor boson
c       integer nhel           : helicity      of tensor boson
c                = -2,-1,0,1,2 : (0 is forbidden if tmass=0.0)
c       integer nst  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex tc(18)         : PSEUDOR  wavefunction    epsilon^mu^nu(t)
c     
      implicit none
      double precision p(0:3), tmass
      integer nhel, nst
      double complex tc(18)

      double complex ft(6,4), ep(4), em(4), e0(4)
      double precision pt, pt2, pp, pzpt, emp, sqh, sqs
      integer i, j

      double precision rZero, rHalf, rOne, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0 )


      tc(1)=NHEL
      tc(17) = dcmplx(p(0),p(3))*nst
      tc(18) = dcmplx(p(1),p(2))*nst

      return
      end
      subroutine rotxxx(p,q , prot)
c
c This subroutine performs the spacial rotation of a four-momentum.
c the momentum p is assumed to be given in the frame where the spacial
c component of q points the positive z-axis.  prot is the momentum p
c rotated to the frame where q is given.
c
c input:
c       real    p(0:3)         : four-momentum p in q(1)=q(2)=0 frame
c       real    q(0:3)         : four-momentum q in the rotated frame
c
c output:
c       real    prot(0:3)      : four-momentum p in the rotated frame
c     
      implicit none
      double precision p(0:3),q(0:3),prot(0:3),qt2,qt,psgn,qq,p1

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )

c
      prot(0) = p(0)

      qt2 = q(1)**2 + q(2)**2


      if ( qt2.eq.rZero ) then
          if ( q(3).eq.rZero ) then
             prot(1) = p(1)
             prot(2) = p(2)
             prot(3) = p(3)
          else
             psgn = dsign(rOne,q(3))
             prot(1) = p(1)*psgn
             prot(2) = p(2)*psgn
             prot(3) = p(3)*psgn
          endif
      else
          qq = sqrt(qt2+q(3)**2)
          qt = sqrt(qt2)
          p1 = p(1)
          prot(1) = q(1)*q(3)/qq/qt*p1 -q(2)/qt*p(2) +q(1)/qq*p(3)
          prot(2) = q(2)*q(3)/qq/qt*p1 +q(1)/qt*p(2) +q(2)/qq*p(3)
          prot(3) =          -qt/qq*p1               +q(3)/qq*p(3)
      endif
c
      return
      end
      subroutine ssssxx(s1,s2,s3,s4,gc , vertex)
c
c This subroutine computes an amplitude of the four-scalar coupling.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex s3(3)          : third  scalar                        s3
c       complex s4(3)          : fourth scalar                        s4
c       complex gc             : coupling constant                 ghhhh
c
c output:
c       complex vertex         : amplitude            gamma(s1,s2,s3,s4)
c     
      implicit none
      double complex s1(3),s2(3),s3(3),s4(3),gc,vertex

c

      vertex = gc*s1(1)*s2(1)*s3(1)*s4(1)
c
      return
      end
      subroutine ssstxx(tc1,tc2,sc,gt,vertex)
c
c- by RF - Feb. 2006 
c
c This subroutine computes an amplitude of the tts coupling.
c
c     input:
c          complex tc1               : Incoming tensor particle
c          complex tc2               : Incoming tensor particle
c          complex sc                : Incoming scalar particle (Higgs)
c          real    gt                : coupling constant for the tts vertex
c
c     output:
c          complex vertex            : amplitude for a tts vertex
c

      implicit none
c--   dimension of the current set to arbitrary length
      INTEGER DIM
      PARAMETER(DIM=18)
c      include "dimension.inc"
      double complex tc1(DIM),tc2(DIM),sc(DIM)

      double complex vertex
      double precision gt

c Take the inner product between the tensor particles
c and multiply it with the scalar particle and the coupling constant.
c Note that the tensor particle is antisymmetric, thus all diagonal terms
c are zero.

      vertex = gt * sc(1) * (
c     &                       + tc1( 1) * tc2( 1)
     &                       - tc1( 2) * tc2( 2)
     &                       - tc1( 3) * tc2( 3)
     &                       - tc1( 4) * tc2( 4)

     &                       - tc1( 5) * tc2( 5)
c     &                       + tc1( 6) * tc2( 6)
     &                       + tc1( 7) * tc2( 7)
     &                       + tc1( 8) * tc2( 8)

     &                       - tc1( 9) * tc2( 9)
     &                       + tc1(10) * tc2(10)
c     &                       + tc1(11) * tc2(11)
     &                       + tc1(12) * tc2(12)

     &                       - tc1(13) * tc2(13)
     &                       + tc1(14) * tc2(14)
     &                       + tc1(15) * tc2(15)
c     &                       + tc1(16) * tc2(16)
     &                                           )

      return
      end
      subroutine sssxxx(s1,s2,s3,gc , vertex)
c
c This subroutine computes an amplitude of the three-scalar coupling.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex s3(3)          : third  scalar                        s3
c       complex gc             : coupling constant                  ghhh
c
c output:
c       complex vertex         : amplitude               gamma(s1,s2,s3)
c     
      implicit none
      double complex s1(3),s2(3),s3(3),gc,vertex

c

      vertex = gc*s1(1)*s2(1)*s3(1)
c
      return
      end
      subroutine sstlxx(s1,s2,t3,gc , vertex)
c- by RF - Mar. 2006
c
c This subroutine computes an amplitude of the three-scalar coupling.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex t3(3)          : internal particle                    s3
c       real    gc             : coupling constant                    gv
c
c output:
c       complex vertex         : amplitude               gamma(s1,s2,s3)
c     
      implicit none

      include "dimension.inc"
      double complex s1(DIM),s2(DIM),t3(DIM),vertex
      double precision gc

      vertex = - gc*s1(1)*s2(1)*t3(1)

      return
      end
      subroutine sstxxx(s1,s2,tc,gt,smass , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two scalar and a tensor boson.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex tc(18)         : input  tensor                         T
c       complex gt             : coupling constant         gts=-1/Lambda
c       real    smass          : scalar mass                         m_s
c
c output:
c       complex vertex         : amplitude                gamma(s1,s2,T)
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex s1(3), s2(3), tc(18), vertex
      double complex gt
      double precision smass

      double complex ft(6,4)
      double complex T12, T13, T14, T23, T24, T34
      double complex TKK
      double precision ps1(4), ps2(4)
	integer i
      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      ps1(1) = dreal(s1(2))
      ps1(2) = dreal(s1(3))
      ps1(3) = dimag(s1(3))
      ps1(4) = dimag(s1(2))

      ps2(1) = -dreal(s2(2))
      ps2(2) = -dreal(s2(3))
      ps2(3) = -dimag(s2(3))
      ps2(4) = -dimag(s2(2))

      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)


      TKK   = cZero
    
      do i = 1,4
         TKK=TKK+ft(i,i)*ps1(i)*ps2(i)
      end do

      TKK   = rTwo*TKK

      TKK = TKK - T12*(ps1(1)*ps2(2) + ps1(2)*ps2(1))
     &          - T13*(ps1(1)*ps2(3) + ps1(3)*ps2(1))
     &          - T14*(ps1(1)*ps2(4) + ps1(4)*ps2(1))
     &          + T23*(ps1(2)*ps2(3) + ps1(3)*ps2(2))
     &          + T24*(ps1(2)*ps2(4) + ps1(4)*ps2(2))
     &          + T34*(ps1(3)*ps2(4) + ps1(4)*ps2(3))


      vertex = TKK+(ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4))
     &	*(smass**2-ps1(1)*ps2(1)+ps1(2)*ps2(2)
     &      +ps1(3)*ps2(3)+ps1(4)*ps2(4))

      vertex = vertex * gt*s1(1)*s2(1)

      return
      end
      subroutine sxxxxx(p,nss , sc)
c
c This subroutine computes a complex SCALAR wavefunction.
c
c input:
c       real    p(0:3)         : four-momentum of scalar boson
c       integer nss  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex sc(3)          : scalar wavefunction                   s
c     
      implicit none
      double complex sc(3)
      double precision p(0:3)
      integer nss

      double precision rOne
      parameter( rOne = 1.0d0 )

c

      sc(1) = dcmplx( rOne )
      sc(2) = dcmplx(p(0),p(3))*nss
      sc(3) = dcmplx(p(1),p(2))*nss
c
      return
      end
      subroutine tpsxxx(t1,t2,sc,gt,xm,xw,vertex)

c  Subroutines for graviton phase space integration
c  KEK 2009.11
c
      implicit none
      double complex t1(18), t2(18), sc(3), vertex, tc2(18)
      double complex gt(2)
      double precision smass,xmass

      double complex ft(6,4),ft2(6,4)
      double precision ps1(4), pt2(4),PT(4),PTT2,PG(0:3)
      integer i
      double complex cZero
      double precision rZero, rTwo,Pi
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )
      double precision PADD
      double precision L_ADD,NADD,MGLOW,MGUP
      double precision YMASS,YWIDTH,xm,xw
      external txxxxx

      Pi=dacos(-1.0d0)
      L_ADD=dimag(gt(1))
      NADD=dreal(gt(1))
      MGUP=dimag(gt(2))
      MGLOW=dreal(gt(2))
      YMASS=xm
      YWIDTH=xw


      ft(1,1) = t1(1)
      ft(1,2) = t1(2)
      ft(1,3) = t1(3)
      ft(1,4) = t1(4)
      ft(2,1) = t1(5)
      ft(2,2) = t1(6)
      ft(2,3) = t1(7)
      ft(2,4) = t1(8)
      ft(3,1) = t1(9)
      ft(3,2) = t1(10)
      ft(3,3) = t1(11)
      ft(3,4) = t1(12)
      ft(4,1) = t1(13)
      ft(4,2) = t1(14)
      ft(4,3) = t1(15)
      ft(4,4) = t1(16)
      ft(5,1) = t1(17)
      ft(6,1) = t1(18)

      ps1(1) = dreal(sc(2))
      ps1(2) = dreal(sc(3))
      ps1(3) = dimag(sc(3))
      ps1(4) = dimag(sc(2))

      pt2(1) = -dreal(t2(17))
      pt2(2) = -dreal(t2(18))
      pt2(3) = -dimag(t2(18))
      pt2(4) = -dimag(t2(17))

      PG(0)=ps1(1)-pt2(1)
      PG(1)=ps1(2)-pt2(2)
      PG(2)=ps1(3)-pt2(3)
      PG(3)=ps1(4)-pt2(4)

      PTT2=PG(0)**2-PG(1)**2-PG(2)**2-PG(3)**2
      xmass=dsqrt(PTT2)

      if(xmass.lt.MGLOW.or.xmass.gt.MGUP) then
      vertex=dcmplx(0.0d0,0.0d0)
      return
      endif

      CALL txxxxx(PG,xmass,INT(t2(1)),+1 , tc2)


       if(INT(NADD).eq.2) then
         PADD=2.0d0*Pi
        elseif(INT(NADD).eq.3) then
         PADD=4.0d0*Pi
        elseif(INT(NADD).eq.4) then
         PADD=2.0d0*Pi**2
        elseif(INT(NADD).eq.5) then
          PADD=8.0d0/3.0d0*Pi**2
        elseif(INT(NADD).eq.6) then
           PADD=Pi**3
        else
        print *, "OUT CASE"
        stop
        endif 

       vertex =dcmplx( PTT2-YMASS**2, YMASS*YWIDTH )* 
     & ( t1(1)*tc2(1)+t1(6)*tc2(6)+t1(11)*tc2(11)+t1(16)*tc2(16)
     & - t1(2)*tc2(2)-t1(3)*tc2(3)-t1(4)*tc2(4)
     & - t1(5)*tc2(5)-t1(9)*tc2(9)-t1(13)*tc2(13)
     & +t1(7)*tc2(7)+t1(8)*tc2(8)+t1(10)*tc2(10)
     & +t1(12)*tc2(12)+t1(14)*tc2(14)+t1(15)*tc2(15)
     & )
     & *dsqrt( 
     &2.0d0*Pi*8.0d0*Pi  ! to compensate the decay phase factor
     &* PADD/L_ADD**NADD*xmass**(NADD-1)  ! density factor for d=4 case
     &/2.0d0/xmass)   ! dm=dm^2/2/m    

      return
      end
      subroutine ttsaxx(tc1,tc2,sc,gc,vertex)
c
c- by RF - Feb. 2006 
c
c This subroutine computes an amplitude of the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex tc2(18)           : Incoming tensor particle
c          complex sc(3)             : Incoming scalar particle (Higgs)
c          complex gc(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c
c     output:
c          complex vertex            : amplitude for a tts vertex
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),tc2(DIM),sc(DIM)
      double complex vertex, gc(2)


c Take the inner product between the tensor particles
c and multiply it with the scalar particle and the coupling constant.
c Note that the tensor particle is antisymmetric, thus all diagonal terms
c are zero.

      if (gc(1).NE.(0D0,0D0)) then

      vertex = - gc(1)*sc(1)* (
     &                 !      + tc1( 1) * tc2( 1)
     &                       - tc1( 2) * tc2( 2)
     &                       - tc1( 3) * tc2( 3)
     &                       - tc1( 4) * tc2( 4)

     &                       - tc1( 5) * tc2( 5)
     &                 !      + tc1( 6) * tc2( 6)
     &                       + tc1( 7) * tc2( 7)
     &                       + tc1( 8) * tc2( 8)

     &                       - tc1( 9) * tc2( 9)
     &                       + tc1(10) * tc2(10)
     &                 !      + tc1(11) * tc2(11)
     &                       + tc1(12) * tc2(12)

     &                       - tc1(13) * tc2(13)
     &                       + tc1(14) * tc2(14)
     &                       + tc1(15) * tc2(15)
     &                 !      + tc1(16) * tc2(16)
     &                                           )


      else
      vertex = (0D0,0D0)
      endif      


      return
      end
      subroutine ttscxx(tc1,tc2,sc,gc,vertex)
c
c- by RF - Feb. 2006 
c  CP3  Modified Nov. 2009 
c
c This subroutine computes an amplitude of the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex tc2(18)           : Incoming tensor particle
c          complex sc(3)             : Incoming scalar particle (Higgs)
c          complex gc(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c
c     output:
c          complex vertex            : amplitude for a tts vertex
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),tc2(DIM),sc(DIM),t1(6,4),t2(6,4)
      double complex vertex, gc(2)


c Take the inner product between the tensor particles
c and multiply it with the scalar particle and the coupling constant.
c Note that the tensor particle is antisymmetric, thus all diagonal terms
c are zero.

      T1(1,1) = tc1(1)
      T1(1,2) = tc1(2)
      T1(1,3) = tc1(3)
      T1(1,4) = tc1(4)
      T1(2,1) = tc1(5)
      T1(2,2) = tc1(6)
      T1(2,3) = tc1(7)
      T1(2,4) = tc1(8)
      T1(3,1) = tc1(9)
      T1(3,2) = tc1(10)
      T1(3,3) = tc1(11)
      T1(3,4) = tc1(12)
      T1(4,1) = tc1(13)
      T1(4,2) = tc1(14)
      T1(4,3) = tc1(15)
      T1(4,4) = tc1(16)
      T1(5,1) = tc1(17)
      T1(6,1) = tc1(18)

      T2(1,1) = tc2(1)
      T2(1,2) = tc2(2)
      T2(1,3) = tc2(3)
      T2(1,4) = tc2(4)
      T2(2,1) = tc2(5)
      T2(2,2) = tc2(6)
      T2(2,3) = tc2(7)
      T2(2,4) = tc2(8)
      T2(3,1) = tc2(9)
      T2(3,2) = tc2(10)
      T2(3,3) = tc2(11)
      T2(3,4) = tc2(12)
      T2(4,1) = tc2(13)
      T2(4,2) = tc2(14)
      T2(4,3) = tc2(15)
      T2(4,4) = tc2(16)
      T2(5,1) = tc2(17)
      T2(6,1) = tc2(18)

      if (gc(1).NE.(0D0,0D0)) then

      vertex =  gc(1)*(SC(1)*T1(1,2)*T2(1,2) - SC(1)*T1(2,1)*T2(1,2) + 
     -  SC(1)*T1(1,3)*T2(1,3) - SC(1)*T1(3,1)*T2(1,3) + 
     -  SC(1)*T1(1,4)*T2(1,4) - SC(1)*T1(4,1)*T2(1,4) - 
     -  SC(1)*T1(1,2)*T2(2,1) + SC(1)*T1(2,1)*T2(2,1) - 
     -  SC(1)*T1(2,3)*T2(2,3) + SC(1)*T1(3,2)*T2(2,3) - 
     -  SC(1)*T1(2,4)*T2(2,4) + SC(1)*T1(4,2)*T2(2,4) - 
     -  SC(1)*T1(1,3)*T2(3,1) + SC(1)*T1(3,1)*T2(3,1) + 
     -  SC(1)*T1(2,3)*T2(3,2) - SC(1)*T1(3,2)*T2(3,2) - 
     -  SC(1)*T1(3,4)*T2(3,4) + SC(1)*T1(4,3)*T2(3,4) - 
     -  SC(1)*T1(1,4)*T2(4,1) + SC(1)*T1(4,1)*T2(4,1) + 
     -  SC(1)*T1(2,4)*T2(4,2) - SC(1)*T1(4,2)*T2(4,2) + 
     -  SC(1)*T1(3,4)*T2(4,3) - SC(1)*T1(4,3)*T2(4,3))


      else
      vertex = (0D0,0D0)
      endif      


      return
      end
      subroutine ttssxx(tc1,tc2,sc1,sc2,g1,g2,vertex)
c
c- by RF - Mar. 2006 
c
c This subroutine computes an amplitude of the tts coupling.
c
c     input:
c          complex tc1               : Incoming tensor particle
c          complex tc2               : Incoming tensor particle
c          complex sc1               : Incoming scalar particle (Higgs)
c          complex sc2               : Incoming scalar particle (Higgs)
c          complex g1(2)             : coupling constant (Higgs Eff. Thr.)
c          real    g2                : coupling constant
c
c     output:
c          complex vertex            : amplitude for a tts vertex
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),tc2(DIM),sc1(DIM),sc2(DIM)

      double complex vertex, g1(2)
      double precision g2


c Take the inner product between the tensor particles
c and multiply it with the scalar particles and the coupling constants.
c Note that the tensor particles are antisymmetric, thus all diagonal terms
c are zero.

      if (g1(1).NE.(0D0,0D0)) then

      vertex = g1(1)*g2*sc1(1)*sc2(1)* (
c     &                       + tc1( 1) * tc2( 1)
     &                       - tc1( 2) * tc2( 2)
     &                       - tc1( 3) * tc2( 3)
     &                       - tc1( 4) * tc2( 4)

     &                       - tc1( 5) * tc2( 5)
c     &                       + tc1( 6) * tc2( 6)
     &                       + tc1( 7) * tc2( 7)
     &                       + tc1( 8) * tc2( 8)

     &                       - tc1( 9) * tc2( 9)
     &                       + tc1(10) * tc2(10)
c     &                       + tc1(11) * tc2(11)
     &                       + tc1(12) * tc2(12)

     &                       - tc1(13) * tc2(13)
     &                       + tc1(14) * tc2(14)
     &                       + tc1(15) * tc2(15)
c     &                       + tc1(16) * tc2(16)
     &                                           )


      else
      vertex = (0D0,0D0)
      endif      


      return
      end
      subroutine ttsxxx(tc1,tc2,sc,gc,vertex)
c
c- by RF - Feb. 2006 
c
c This subroutine computes an amplitude of the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex tc2(18)           : Incoming tensor particle
c          complex sc(3)             : Incoming scalar particle (Higgs)
c          complex gc(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c
c     output:
c          complex vertex            : amplitude for a tts vertex
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),tc2(DIM),sc(DIM)
      double complex vertex, gc(2)


c Take the inner product between the tensor particles
c and multiply it with the scalar particle and the coupling constant.
c Note that the tensor particle is antisymmetric, thus all diagonal terms
c are zero.

      if (gc(1).NE.(0D0,0D0)) then

      vertex = - gc(1)*sc(1)* (
     &                 !      + tc1( 1) * tc2( 1)
     &                       - tc1( 2) * tc2( 2)
     &                       - tc1( 3) * tc2( 3)
     &                       - tc1( 4) * tc2( 4)

     &                       - tc1( 5) * tc2( 5)
     &                 !      + tc1( 6) * tc2( 6)
     &                       + tc1( 7) * tc2( 7)
     &                       + tc1( 8) * tc2( 8)

     &                       - tc1( 9) * tc2( 9)
     &                       + tc1(10) * tc2(10)
     &                 !      + tc1(11) * tc2(11)
     &                       + tc1(12) * tc2(12)

     &                       - tc1(13) * tc2(13)
     &                       + tc1(14) * tc2(14)
     &                       + tc1(15) * tc2(15)
     &                 !      + tc1(16) * tc2(16)
     &                                           )


      else
      vertex = (0D0,0D0)
      endif      


      return
      end
      subroutine tttxxx(tt,ta,tb,gt,vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two gauge bosons and a tensor boson.
c
c input:
c       complex ta(18)          : first  tensor                           ta
c       complex tb(18)          : second tensor                       tb
c       complex tt(18)           : input  tensor                         tt
c       complex gt             : coupling constant         gtv=-1/Lambda
c output:
c       complex vertex         : amplitude                gamma(tt,ta,tb)
c
      implicit none
      double complex ta(18), tb(18), tt(18), gt, vertex
      double complex T1(6,4),T2(6,4),T(6,4)

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      T(1,1) = tt(1)
      T(1,2) = tt(2)
      T(1,3) = tt(3)
      T(1,4) = tt(4)
      T(2,1) = tt(5)
      T(2,2) = tt(6)
      T(2,3) = tt(7)
      T(2,4) = tt(8)
      T(3,1) = tt(9)
      T(3,2) = tt(10)
      T(3,3) = tt(11)
      T(3,4) = tt(12)
      T(4,1) = tt(13)
      T(4,2) = tt(14)
      T(4,3) = tt(15)
      T(4,4) = tt(16)
      T(5,1) = tt(17)
      T(6,1) = tt(18)
 
      T1(1,1) = ta(1)
      T1(1,2) = ta(2)
      T1(1,3) = ta(3)
      T1(1,4) = ta(4)
      T1(2,1) = ta(5)
      T1(2,2) = ta(6)
      T1(2,3) = ta(7)
      T1(2,4) = ta(8)
      T1(3,1) = ta(9)
      T1(3,2) = ta(10)
      T1(3,3) = ta(11)
      T1(3,4) = ta(12)
      T1(4,1) = ta(13)
      T1(4,2) = ta(14)
      T1(4,3) = ta(15)
      T1(4,4) = ta(16)
      T1(5,1) = ta(17)
      T1(6,1) = ta(18)

      T2(1,1) = tb(1)
      T2(1,2) = tb(2)
      T2(1,3) = tb(3)
      T2(1,4) = tb(4)
      T2(2,1) = tb(5)
      T2(2,2) = tb(6)
      T2(2,3) = tb(7)
      T2(2,4) = tb(8)
      T2(3,1) = tb(9)
      T2(3,2) = tb(10)
      T2(3,3) = tb(11)
      T2(3,4) = tb(12)
      T2(4,1) = tb(13)
      T2(4,2) = tb(14)
      T2(4,3) = tb(15)
      T2(4,4) = tb(16)
      T2(5,1) = tb(17)
      T2(6,1) = tb(18)


      vertex =-gt*(T(1,1)*T1(1,2)*T2(1,2) - T(2,2)*T1(1,2)*T2(1,2) + 
     -  T(3,3)*T1(1,2)*T2(1,2) + T(4,4)*T1(1,2)*T2(1,2) - 
     -  T(2,3)*T1(1,3)*T2(1,2) - T(3,2)*T1(1,3)*T2(1,2) - 
     -  T(2,4)*T1(1,4)*T2(1,2) - T(4,2)*T1(1,4)*T2(1,2) - 
     -  T(1,1)*T1(2,1)*T2(1,2) + T(2,2)*T1(2,1)*T2(1,2) - 
     -  T(3,3)*T1(2,1)*T2(1,2) - T(4,4)*T1(2,1)*T2(1,2) + 
     -  T(1,3)*T1(2,3)*T2(1,2) + T(3,1)*T1(2,3)*T2(1,2) + 
     -  T(1,4)*T1(2,4)*T2(1,2) + T(4,1)*T1(2,4)*T2(1,2) + 
     -  T(2,3)*T1(3,1)*T2(1,2) + T(3,2)*T1(3,1)*T2(1,2) - 
     -  T(1,3)*T1(3,2)*T2(1,2) - T(3,1)*T1(3,2)*T2(1,2) + 
     -  T(2,4)*T1(4,1)*T2(1,2) + T(4,2)*T1(4,1)*T2(1,2) - 
     -  T(1,4)*T1(4,2)*T2(1,2) - T(4,1)*T1(4,2)*T2(1,2) - 
     -  T(2,3)*T1(1,2)*T2(1,3) - T(3,2)*T1(1,2)*T2(1,3) + 
     -  T(1,1)*T1(1,3)*T2(1,3) + T(2,2)*T1(1,3)*T2(1,3) - 
     -  T(3,3)*T1(1,3)*T2(1,3) + T(4,4)*T1(1,3)*T2(1,3) - 
     -  T(3,4)*T1(1,4)*T2(1,3) - T(4,3)*T1(1,4)*T2(1,3) + 
     -  T(2,3)*T1(2,1)*T2(1,3) + T(3,2)*T1(2,1)*T2(1,3) - 
     -  T(1,2)*T1(2,3)*T2(1,3) - T(2,1)*T1(2,3)*T2(1,3) - 
     -  T(1,1)*T1(3,1)*T2(1,3) - T(2,2)*T1(3,1)*T2(1,3) + 
     -  T(3,3)*T1(3,1)*T2(1,3) - T(4,4)*T1(3,1)*T2(1,3) + 
     -  T(1,2)*T1(3,2)*T2(1,3) + T(2,1)*T1(3,2)*T2(1,3) + 
     -  T(1,4)*T1(3,4)*T2(1,3) + T(4,1)*T1(3,4)*T2(1,3) + 
     -  T(3,4)*T1(4,1)*T2(1,3) + T(4,3)*T1(4,1)*T2(1,3) - 
     -  T(1,4)*T1(4,3)*T2(1,3) - T(4,1)*T1(4,3)*T2(1,3) - 
     -  T(2,4)*T1(1,2)*T2(1,4) - T(4,2)*T1(1,2)*T2(1,4) - 
     -  T(3,4)*T1(1,3)*T2(1,4) - T(4,3)*T1(1,3)*T2(1,4) + 
     -  T(1,1)*T1(1,4)*T2(1,4) + T(2,2)*T1(1,4)*T2(1,4) + 
     -  T(3,3)*T1(1,4)*T2(1,4) - T(4,4)*T1(1,4)*T2(1,4) + 
     -  T(2,4)*T1(2,1)*T2(1,4) + T(4,2)*T1(2,1)*T2(1,4) - 
     -  T(1,2)*T1(2,4)*T2(1,4) - T(2,1)*T1(2,4)*T2(1,4) + 
     -  T(3,4)*T1(3,1)*T2(1,4) + T(4,3)*T1(3,1)*T2(1,4) - 
     -  T(1,3)*T1(3,4)*T2(1,4) - T(3,1)*T1(3,4)*T2(1,4) - 
     -  T(1,1)*T1(4,1)*T2(1,4) - T(2,2)*T1(4,1)*T2(1,4) - 
     -  T(3,3)*T1(4,1)*T2(1,4) + T(4,4)*T1(4,1)*T2(1,4) + 
     -  T(1,2)*T1(4,2)*T2(1,4) + T(2,1)*T1(4,2)*T2(1,4) + 
     -  T(1,3)*T1(4,3)*T2(1,4) + T(3,1)*T1(4,3)*T2(1,4) - 
     -  T(1,1)*T1(1,2)*T2(2,1) + T(2,2)*T1(1,2)*T2(2,1) - 
     -  T(3,3)*T1(1,2)*T2(2,1) - T(4,4)*T1(1,2)*T2(2,1) + 
     -  T(2,3)*T1(1,3)*T2(2,1) + T(3,2)*T1(1,3)*T2(2,1) + 
     -  T(2,4)*T1(1,4)*T2(2,1) + T(4,2)*T1(1,4)*T2(2,1) + 
     -  T(1,1)*T1(2,1)*T2(2,1) - T(2,2)*T1(2,1)*T2(2,1) + 
     -  T(3,3)*T1(2,1)*T2(2,1) + T(4,4)*T1(2,1)*T2(2,1) - 
     -  T(1,3)*T1(2,3)*T2(2,1) - T(3,1)*T1(2,3)*T2(2,1) - 
     -  T(1,4)*T1(2,4)*T2(2,1) - T(4,1)*T1(2,4)*T2(2,1) - 
     -  T(2,3)*T1(3,1)*T2(2,1) - T(3,2)*T1(3,1)*T2(2,1) + 
     -  T(1,3)*T1(3,2)*T2(2,1) + T(3,1)*T1(3,2)*T2(2,1) - 
     -  T(2,4)*T1(4,1)*T2(2,1) - T(4,2)*T1(4,1)*T2(2,1) + 
     -  T(1,4)*T1(4,2)*T2(2,1) + T(4,1)*T1(4,2)*T2(2,1) + 
     -  T(1,3)*T1(1,2)*T2(2,3) + T(3,1)*T1(1,2)*T2(2,3) - 
     -  T(1,2)*T1(1,3)*T2(2,3) - T(2,1)*T1(1,3)*T2(2,3) - 
     -  T(1,3)*T1(2,1)*T2(2,3) - T(3,1)*T1(2,1)*T2(2,3) + 
     -  T(1,1)*T1(2,3)*T2(2,3) + T(2,2)*T1(2,3)*T2(2,3) + 
     -  T(3,3)*T1(2,3)*T2(2,3) - T(4,4)*T1(2,3)*T2(2,3) + 
     -  T(3,4)*T1(2,4)*T2(2,3) + T(4,3)*T1(2,4)*T2(2,3) + 
     -  T(1,2)*T1(3,1)*T2(2,3) + T(2,1)*T1(3,1)*T2(2,3) - 
     -  T(1,1)*T1(3,2)*T2(2,3) - T(2,2)*T1(3,2)*T2(2,3) - 
     -  T(3,3)*T1(3,2)*T2(2,3) + T(4,4)*T1(3,2)*T2(2,3) - 
     -  T(2,4)*T1(3,4)*T2(2,3) - T(4,2)*T1(3,4)*T2(2,3) - 
     -  T(3,4)*T1(4,2)*T2(2,3) - T(4,3)*T1(4,2)*T2(2,3) + 
     -  T(2,4)*T1(4,3)*T2(2,3) + T(4,2)*T1(4,3)*T2(2,3) + 
     -  T(1,4)*T1(1,2)*T2(2,4) + T(4,1)*T1(1,2)*T2(2,4) - 
     -  T(1,2)*T1(1,4)*T2(2,4) - T(2,1)*T1(1,4)*T2(2,4) - 
     -  T(1,4)*T1(2,1)*T2(2,4) - T(4,1)*T1(2,1)*T2(2,4) + 
     -  T(3,4)*T1(2,3)*T2(2,4) + T(4,3)*T1(2,3)*T2(2,4) + 
     -  T(1,1)*T1(2,4)*T2(2,4) + T(2,2)*T1(2,4)*T2(2,4) - 
     -  T(3,3)*T1(2,4)*T2(2,4) + T(4,4)*T1(2,4)*T2(2,4) - 
     -  T(3,4)*T1(3,2)*T2(2,4) - T(4,3)*T1(3,2)*T2(2,4) + 
     -  T(2,3)*T1(3,4)*T2(2,4) + T(3,2)*T1(3,4)*T2(2,4) + 
     -  T(1,2)*T1(4,1)*T2(2,4) + T(2,1)*T1(4,1)*T2(2,4) - 
     -  T(1,1)*T1(4,2)*T2(2,4) - T(2,2)*T1(4,2)*T2(2,4) + 
     -  T(3,3)*T1(4,2)*T2(2,4) - T(4,4)*T1(4,2)*T2(2,4) - 
     -  T(2,3)*T1(4,3)*T2(2,4) - T(3,2)*T1(4,3)*T2(2,4) + 
     -  T(2,3)*T1(1,2)*T2(3,1) + T(3,2)*T1(1,2)*T2(3,1) - 
     -  T(1,1)*T1(1,3)*T2(3,1) - T(2,2)*T1(1,3)*T2(3,1) + 
     -  T(3,3)*T1(1,3)*T2(3,1) - T(4,4)*T1(1,3)*T2(3,1) + 
     -  T(3,4)*T1(1,4)*T2(3,1) + T(4,3)*T1(1,4)*T2(3,1) - 
     -  T(2,3)*T1(2,1)*T2(3,1) - T(3,2)*T1(2,1)*T2(3,1) + 
     -  T(1,2)*T1(2,3)*T2(3,1) + T(2,1)*T1(2,3)*T2(3,1) + 
     -  T(1,1)*T1(3,1)*T2(3,1) + T(2,2)*T1(3,1)*T2(3,1) - 
     -  T(3,3)*T1(3,1)*T2(3,1) + T(4,4)*T1(3,1)*T2(3,1) - 
     -  T(1,2)*T1(3,2)*T2(3,1) - T(2,1)*T1(3,2)*T2(3,1) - 
     -  T(1,4)*T1(3,4)*T2(3,1) - T(4,1)*T1(3,4)*T2(3,1) - 
     -  T(3,4)*T1(4,1)*T2(3,1) - T(4,3)*T1(4,1)*T2(3,1) + 
     -  T(1,4)*T1(4,3)*T2(3,1) + T(4,1)*T1(4,3)*T2(3,1) - 
     -  T(1,3)*T1(1,2)*T2(3,2) - T(3,1)*T1(1,2)*T2(3,2) + 
     -  T(1,2)*T1(1,3)*T2(3,2) + T(2,1)*T1(1,3)*T2(3,2) + 
     -  T(1,3)*T1(2,1)*T2(3,2) + T(3,1)*T1(2,1)*T2(3,2) - 
     -  T(1,1)*T1(2,3)*T2(3,2) - T(2,2)*T1(2,3)*T2(3,2) - 
     -  T(3,3)*T1(2,3)*T2(3,2) + T(4,4)*T1(2,3)*T2(3,2) - 
     -  T(3,4)*T1(2,4)*T2(3,2) - T(4,3)*T1(2,4)*T2(3,2) - 
     -  T(1,2)*T1(3,1)*T2(3,2) - T(2,1)*T1(3,1)*T2(3,2) + 
     -  T(1,1)*T1(3,2)*T2(3,2) + T(2,2)*T1(3,2)*T2(3,2) + 
     -  T(3,3)*T1(3,2)*T2(3,2) - T(4,4)*T1(3,2)*T2(3,2) + 
     -  T(2,4)*T1(3,4)*T2(3,2) + T(4,2)*T1(3,4)*T2(3,2) + 
     -  T(3,4)*T1(4,2)*T2(3,2) + T(4,3)*T1(4,2)*T2(3,2) - 
     -  T(2,4)*T1(4,3)*T2(3,2) - T(4,2)*T1(4,3)*T2(3,2) + 
     -  T(1,4)*T1(1,3)*T2(3,4) + T(4,1)*T1(1,3)*T2(3,4) - 
     -  T(1,3)*T1(1,4)*T2(3,4) - T(3,1)*T1(1,4)*T2(3,4) - 
     -  T(2,4)*T1(2,3)*T2(3,4) - T(4,2)*T1(2,3)*T2(3,4) + 
     -  T(2,3)*T1(2,4)*T2(3,4) + T(3,2)*T1(2,4)*T2(3,4) - 
     -  T(1,4)*T1(3,1)*T2(3,4) - T(4,1)*T1(3,1)*T2(3,4) + 
     -  T(2,4)*T1(3,2)*T2(3,4) + T(4,2)*T1(3,2)*T2(3,4) + 
     -  T(1,1)*T1(3,4)*T2(3,4) - T(2,2)*T1(3,4)*T2(3,4) + 
     -  T(3,3)*T1(3,4)*T2(3,4) + T(4,4)*T1(3,4)*T2(3,4) + 
     -  T(1,3)*T1(4,1)*T2(3,4) + T(3,1)*T1(4,1)*T2(3,4) - 
     -  T(2,3)*T1(4,2)*T2(3,4) - T(3,2)*T1(4,2)*T2(3,4) - 
     -  T(1,1)*T1(4,3)*T2(3,4) + T(2,2)*T1(4,3)*T2(3,4) - 
     -  T(3,3)*T1(4,3)*T2(3,4) - T(4,4)*T1(4,3)*T2(3,4) + 
     -  T(2,4)*T1(1,2)*T2(4,1) + T(4,2)*T1(1,2)*T2(4,1) + 
     -  T(3,4)*T1(1,3)*T2(4,1) + T(4,3)*T1(1,3)*T2(4,1) - 
     -  T(1,1)*T1(1,4)*T2(4,1) - T(2,2)*T1(1,4)*T2(4,1) - 
     -  T(3,3)*T1(1,4)*T2(4,1) + T(4,4)*T1(1,4)*T2(4,1) - 
     -  T(2,4)*T1(2,1)*T2(4,1) - T(4,2)*T1(2,1)*T2(4,1) + 
     -  T(1,2)*T1(2,4)*T2(4,1) + T(2,1)*T1(2,4)*T2(4,1) - 
     -  T(3,4)*T1(3,1)*T2(4,1) - T(4,3)*T1(3,1)*T2(4,1) + 
     -  T(1,3)*T1(3,4)*T2(4,1) + T(3,1)*T1(3,4)*T2(4,1) + 
     -  T(1,1)*T1(4,1)*T2(4,1) + T(2,2)*T1(4,1)*T2(4,1) + 
     -  T(3,3)*T1(4,1)*T2(4,1) - T(4,4)*T1(4,1)*T2(4,1) - 
     -  T(1,2)*T1(4,2)*T2(4,1) - T(2,1)*T1(4,2)*T2(4,1) - 
     -  T(1,3)*T1(4,3)*T2(4,1) - T(3,1)*T1(4,3)*T2(4,1) - 
     -  T(1,4)*T1(1,2)*T2(4,2) - T(4,1)*T1(1,2)*T2(4,2) + 
     -  T(1,2)*T1(1,4)*T2(4,2) + T(2,1)*T1(1,4)*T2(4,2) + 
     -  T(1,4)*T1(2,1)*T2(4,2) + T(4,1)*T1(2,1)*T2(4,2) - 
     -  T(3,4)*T1(2,3)*T2(4,2) - T(4,3)*T1(2,3)*T2(4,2) - 
     -  T(1,1)*T1(2,4)*T2(4,2) - T(2,2)*T1(2,4)*T2(4,2) + 
     -  T(3,3)*T1(2,4)*T2(4,2) - T(4,4)*T1(2,4)*T2(4,2) + 
     -  T(3,4)*T1(3,2)*T2(4,2) + T(4,3)*T1(3,2)*T2(4,2) - 
     -  T(2,3)*T1(3,4)*T2(4,2) - T(3,2)*T1(3,4)*T2(4,2) - 
     -  T(1,2)*T1(4,1)*T2(4,2) - T(2,1)*T1(4,1)*T2(4,2) + 
     -  T(1,1)*T1(4,2)*T2(4,2) + T(2,2)*T1(4,2)*T2(4,2) - 
     -  T(3,3)*T1(4,2)*T2(4,2) + T(4,4)*T1(4,2)*T2(4,2) + 
     -  T(2,3)*T1(4,3)*T2(4,2) + T(3,2)*T1(4,3)*T2(4,2) - 
     -  T(1,4)*T1(1,3)*T2(4,3) - T(4,1)*T1(1,3)*T2(4,3) + 
     -  T(1,3)*T1(1,4)*T2(4,3) + T(3,1)*T1(1,4)*T2(4,3) + 
     -  T(2,4)*T1(2,3)*T2(4,3) + T(4,2)*T1(2,3)*T2(4,3) - 
     -  T(2,3)*T1(2,4)*T2(4,3) - T(3,2)*T1(2,4)*T2(4,3) + 
     -  T(1,4)*T1(3,1)*T2(4,3) + T(4,1)*T1(3,1)*T2(4,3) - 
     -  T(2,4)*T1(3,2)*T2(4,3) - T(4,2)*T1(3,2)*T2(4,3) - 
     -  T(1,1)*T1(3,4)*T2(4,3) + T(2,2)*T1(3,4)*T2(4,3) - 
     -  T(3,3)*T1(3,4)*T2(4,3) - T(4,4)*T1(3,4)*T2(4,3) - 
     -  T(1,3)*T1(4,1)*T2(4,3) - T(3,1)*T1(4,1)*T2(4,3) + 
     -  T(2,3)*T1(4,2)*T2(4,3) + T(3,2)*T1(4,2)*T2(4,3) + 
     -  T(1,1)*T1(4,3)*T2(4,3) - T(2,2)*T1(4,3)*T2(4,3) + 
     -  T(3,3)*T1(4,3)*T2(4,3) + T(4,4)*T1(4,3)*T2(4,3))

       return
       end

      subroutine uggggx(va,vb,vc,vd,gc,gt,tmass,twidth , ugggg)
c
c This subroutine computes an off-shell tensor boson 
c current from the 4 gluons-tensor boson coupling, 
c corresponding to the color structure f^{a,b,e} f{c,d,e}. 
c
c To obtain the complete amplitude, this subroutine must be called three
c times (once for each color structure) with the following permutations:
c     call uggggx(ga,gb,gc,gd,g,gt,tmass,twidth , ugggg1)
c     call uggggx(ga,gc,gd,gb,g,gt,tmass,twidth , ugggg2)
c     call uggggx(ga,gd,gb,gv,g,gt,tmass,twidth , ugggg3)
c corresponding to 
c	f^{a,b,e} f^{c,d,e}
c	f^{a,c,e} f^{d,b,e}
c	f^{a,d,e} f^{b,c,e}
c
c input:
c       complex va(6)          : boson with adjoint color index a     va
c       complex vb(6)          : boson with adjoint color index b     vb
c       complex vc(6)          : boson with adjoint color index c     vc
c       complex vd(6)          : boson with adjoint color index d     vd
c       real    gc             : coupling constant                    gs
c       complex gt             : coupling constant         gtv=-1/Lambda
c       real    tmass          : mass  of output tensor T
c       real    twidth         : width of output tensor T
c
c output:
c       complex ugggg(18)      : tensor current   j^mu^nu(T:va,vb,vc,vd)
c     
c- by Q.Li - JAN. 2008
c
      implicit none
      double complex va(6), vb(6), vc(6), vd(6), gt, ugggg(18)
      double precision gc, tmass, twidth

      double complex tc(6,4)
      double complex E1E3,E1E4,E2E3,E2E4,E1PT,E2PT,E3PT,E4PT

      integer a,b,i,j

      double precision pa(4), pb(4),pc(4),pd(4),pT(4)
      double precision MET(4,4)
      double complex cZero, d
      double precision rZero, rTwo,PT2,r3
      parameter( rZero = 0.0d0, rTwo = 2.0d0, r3=3.0d0)
      parameter(cZero=(0.0d0,0.0d0))

      
      tc(5,1) = va(5)+vb(5)+vc(5)+vd(5)
      tc(6,1) = va(6)+vb(6)+vc(6)+vd(6)

      pa(1) = dreal(va(5))
      pa(2) = dreal(va(6))
      pa(3) = dimag(va(6))
      pa(4) = dimag(va(5))

      pb(1) = dreal(vb(5))
      pb(2) = dreal(vb(6))
      pb(3) = dimag(vb(6))
      pb(4) = dimag(vb(5))

      pc(1) = dreal(vc(5))
      pc(2) = dreal(vc(6))
      pc(3) = dimag(vc(6))
      pc(4) = dimag(vc(5))

      pd(1) = dreal(vd(5))
      pd(2) = dreal(vd(6))
      pd(3) = dimag(vd(6))
      pd(4) = dimag(vd(5))

     
      pT(1) = dreal(tc(5,1))
      pT(2) = dreal(tc(6,1))
      pT(3) = dimag(tc(6,1))
      pT(4) = dimag(tc(5,1))

      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0

      PT2 = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2
      
      E1PT = pT(1)*va(1)-pT(2)*va(2)-pT(3)*va(3)-pT(4)*va(4)
      E2PT = pT(1)*vb(1)-pT(2)*vb(2)-pT(3)*vb(3)-pT(4)*vb(4)
      E3PT = pT(1)*vc(1)-pT(2)*vc(2)-pT(3)*vc(3)-pT(4)*vc(4)
      E4PT = pT(1)*vd(1)-pT(2)*vd(2)-pT(3)*vd(3)-pT(4)*vd(4)
      
      E1E3 = va(1)*vc(1)-va(2)*vc(2)-va(3)*vc(3)-va(4)*vc(4)
      E1E4 = va(1)*vd(1)-va(2)*vd(2)-va(3)*vd(3)-va(4)*vd(4)
      E2E3 = vb(1)*vc(1)-vb(2)*vc(2)-vb(3)*vc(3)-vb(4)*vc(4)
      E2E4 = vb(1)*vd(1)-vb(2)*vd(2)-vb(3)*vd(3)-vb(4)*vd(4)
   
      if ( tmass.gt.rZero ) then
         d = - 1.0d0/dcmplx( PT2-tmass**2, tmass*twidth )
      else
         d = - 1.0d0/dcmplx( PT2, rZero )
      end if
      
	    
      do a = 1,4
         do b=1,4
            
            tc(a,b) = -2*E2E4*va(b)*vc(a) + 2*E1E4*vb(b)*vc(a)
     &- 2*E2E4*va(a)*vc(b) + 2*E1E4*vb(a)*vc(b) + 2*E2E3*va(b)*vd(a)
     & - 2*E1E3*vb(b)*vd(a) + 2*E2E3*va(a)*vd(b) - 
     &  2*E1E3*vb(a)*vd(b) - 2*E1E4*E2E3*MET(a,b) + 2*E1E3*E2E4*MET(a,b) 
     &- (2*E1PT*E2E4*E3PT*rtwo*MET(a,b))/(r3*tmass**2) + 
     &  (2*E1E4*E2PT*E3PT*rtwo*MET(a,b))/(r3*tmass**2) 
     &+ (2*E1PT*E2E3*E4PT*rtwo*MET(a,b))/(r3*tmass**2)
     & - (2*E1E3*E2PT*E4PT*rtwo*MET(a,b))/(r3*tmass**2) - 
     &  (E1E4*E2E3*PT2*rtwo*MET(a,b))/(r3*tmass**2) 
     &+ (E1E3*E2E4*PT2*rtwo*MET(a,b))/(r3*tmass**2)
     & + (2*E2E4*E3PT*va(b)*PT(a))/tmass**2 - 
     &  (2*E2E3*E4PT*va(b)*PT(a))/tmass**2
     & - (2*E1E4*E3PT*vb(b)*PT(a))/tmass**2
     & + (2*E1E3*E4PT*vb(b)*PT(a))/tmass**2
     & + (2*E1PT*E2E4*vc(b)*PT(a))/tmass**2 - 
     &  (2*E1E4*E2PT*vc(b)*PT(a))/tmass**2
     & - (2*E1PT*E2E3*vd(b)*PT(a))/tmass**2
     & + (2*E1E3*E2PT*vd(b)*PT(a))/tmass**2
     & + (2*E2E4*E3PT*va(a)*PT(b))/tmass**2 - 
     &  (2*E2E3*E4PT*va(a)*PT(b))/tmass**2 
     &- (2*E1E4*E3PT*vb(a)*PT(b))/tmass**2
     & + (2*E1E3*E4PT*vb(a)*PT(b))/tmass**2 
     &+ (2*E1PT*E2E4*vc(a)*PT(b))/tmass**2 - 
     &  (2*E1E4*E2PT*vc(a)*PT(b))/tmass**2 
     &- (2*E1PT*E2E3*vd(a)*PT(b))/tmass**2 
     &+ (2*E1E3*E2PT*vd(a)*PT(b))/tmass**2
     & - (4*E1PT*E2E4*E3PT*PT(a)*PT(b))/tmass**4 + 
     &  (4*E1E4*E2PT*E3PT*PT(a)*PT(b))/tmass**4 
     &+ (4*E1PT*E2E3*E4PT*PT(a)*PT(b))/tmass**4
     & - (4*E1E3*E2PT*E4PT*PT(a)*PT(b))/tmass**4 - 
     &  (2*E1E4*E2E3*PT2*PT(a)*PT(b))/tmass**4
     & + (2*E1E3*E2E4*PT2*PT(a)*PT(b))/tmass**4
     & + (2*E1PT*E2E4*E3PT*rtwo*PT(a)*PT(b))/(r3*tmass**4) - 
     &  (2*E1E4*E2PT*E3PT*rtwo*PT(a)*PT(b))/(r3*tmass**4)
     & - (2*E1PT*E2E3*E4PT*rtwo*PT(a)*PT(b))/(r3*tmass**4)
     & + (2*E1E3*E2PT*E4PT*rtwo*PT(a)*PT(b))/(r3*tmass**4) + 
     &  (E1E4*E2E3*PT2*rtwo*PT(a)*PT(b))/(r3*tmass**4)
     & - (E1E3*E2E4*PT2*rtwo*PT(a)*PT(b))/(r3*tmass**4) 
     &+ (4*E1E4*E2E3*PT(a)*PT(b))/tmass**2 - 
     &  (4*E1E3*E2E4*PT(a)*PT(b))/tmass**2

            tc(a,b) = -tc(a,b)*d/2.0d0*gc*gc*gt
c     2.0 factor from propagator convention

         enddo
      enddo

      ugggg(1) = tc(1,1)
      ugggg(2) = tc(1,2)
      ugggg(3) = tc(1,3)
      ugggg(4) = tc(1,4)
      ugggg(5) = tc(2,1)
      ugggg(6) = tc(2,2)
      ugggg(7) = tc(2,3)
      ugggg(8) = tc(2,4)
      ugggg(9) = tc(3,1)
      ugggg(10) = tc(3,2)
      ugggg(11) = tc(3,3)
      ugggg(12) = tc(3,4)
      ugggg(13) = tc(4,1)
      ugggg(14) = tc(4,2)
      ugggg(15) = tc(4,3)
      ugggg(16) = tc(4,4)
      ugggg(17) = tc(5,1)
      ugggg(18) = tc(6,1)

      return
      end
      subroutine uiovxx(fi,fo,vc,gc,gt,tmass,twidth , uiov)
c
c This subroutine computes an off-shell tensor boson 
c wavefunction from a flowing-out fermion, a flowing-in fermion and 
c a gauge vector boson.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                       v
c       complex gc(2)          : coupling constants                  gvf
c       complex gt             : coupling constant      gtfv=-1/Lambda/2
c       real    tmass          : mass  of output tensor T
c       real    twidth         : width of output tensor T
c
c output:
c       complex uiov(18)       : tensor boson       j^mu^nu(<fo|v,T|fi>)
c
c- by Q.Li - OCT. 2006
c- Added massless tensor - P. de Aquino - Oct. 2009 
c
      implicit none
      double complex fi(6), fo(6), vc(6), gc(2), gt, uiov(18)
      double precision  tmass, twidth

      double complex yiov(6,4)
      double complex d
      double precision pt(4), pt2, KTVC
      double precision MET(4,4)
      integer i,j
      
      double precision rZero, rTwo,r4,r8,r3
      parameter(rZero=0.0d0, rTwo=2.0d0, r4=4.0d0, r8=8.0d0, r3=3.0d0)
      double complex cone
      parameter( cone = ( 0.0d0, 1.0d0 ))


      yiov(5,1) = fo(5) + vc(5) -fi(5)
      yiov(6,1) = fo(6) + vc(6) -fi(6)

      pt(1) = dreal(yiov(5,1))
      pt(2) = dreal(yiov(6,1))
      pt(3) = dimag(yiov(6,1))
      pt(4) = dimag(yiov(5,1))
      
      pt2 = pt(1)**2-pt(2)**2-pt(3)**2-pt(4)**2
      KTVC = pt(1)*vc(1)-pt(2)*vc(2)-pt(3)*vc(3)-pt(4)*vc(4)
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      if ( tmass.eq.rZero ) then
         d =  -1.0d0/dcmplx( pt2, rZero )

         do i = 1,4
	    do j=1,4
            yiov(i,j) = fi(4)*gc(2)*(fo(1)*MET(i,j)*(vc(2) - cone*vc(3)) 
     &+ fo(2)*MET(i,j)*(-vc(1) - vc(4))) 
     &+ fi(1)*gc(1)*(fo(4)*MET(i,j)*(-vc(2) 
     &- cone*vc(3)) + fo(3)*MET(i,j)*(-vc(1) - vc(4))) 
     &+ fi(3)*gc(2)*(fo(2)*MET(i,j)*(vc(2) 
     &+ cone*vc(3)) + fo(1)*MET(i,j)*(-vc(1) + vc(4))) 
     &+ fi(2)*gc(1)*(fo(3)*MET(i,j)*(-vc(2) + cone*vc(3)) 
     &+ fo(4)*MET(i,j)*(-vc(1) + vc(4)))
   
         end do
        end do 
   
         yiov(1,1) = yiov(1,1)
     &- 2*fi(1)*fo(3)*gc(1)*vc(1) 
     &- 2*fi(2)*fo(4)*gc(1)*vc(1) 
     &- 2*fi(3)*fo(1)*gc(2)*vc(1) 
     &- 2*fi(4)*fo(2)*gc(2)*vc(1)
   
         yiov(2,2) = yiov(2,2)
     &+	2*fi(2)*fo(3)*gc(1)*vc(2) 
     &+ 2*fi(1)*fo(4)*gc(1)*vc(2) 
     &- 2*fi(4)*fo(1)*gc(2)*vc(2) 
     &- 2*fi(3)*fo(2)*gc(2)*vc(2)
   
         yiov(3,3) = yiov(3,3)
     &- 2*cone*fi(2)*fo(3)*gc(1)*vc(3) 
     &+ 2*cone*fi(1)*fo(4)*gc(1)*vc(3) 
     &+ 2*cone*fi(4)*fo(1)*gc(2)*vc(3) 
     &- 2*cone*fi(3)*fo(2)*gc(2)*vc(3)
   
         yiov(4,4) = yiov(4,4)
     &+ 2*fi(1)*fo(3)*gc(1)*vc(4) 
     &- 2*fi(2)*fo(4)*gc(1)*vc(4) 
     &- 2*fi(3)*fo(1)*gc(2)*vc(4) 
     &+ 2*fi(4)*fo(2)*gc(2)*vc(4)
   
         yiov(1,2) = yiov(1,2)
     &+ fi(3)*gc(2)*(-(fo(2)*vc(1)) - fo(1)*vc(2)) 
     &+ fi(4)*gc(2)*(-(fo(1)*vc(1)) - fo(2)*vc(2)) 
     &+ fi(1)*gc(1)*(fo(4)*vc(1) - fo(3)*vc(2)) 
     &+ fi(2)*gc(1)*(fo(3)*vc(1) - fo(4)*vc(2))
     
         yiov(1,3) = yiov(1,3)
     &+	fi(3)*gc(2)*(-(cone*fo(2)*vc(1)) - fo(1)*vc(3)) 
     &+ fi(4)*gc(2)*(cone*fo(1)*vc(1) - fo(2)*vc(3)) 
     &+ fi(1)*gc(1)*(cone*fo(4)*vc(1) - fo(3)*vc(3)) 
     &+ fi(2)*gc(1)*(-(cone*fo(3)*vc(1)) - fo(4)*vc(3))
     
         yiov(1,4) = yiov(1,4)
     &+ fi(2)*fo(4)*gc(1)*(-vc(1) - vc(4)) 
     &+ fi(3)*fo(1)*gc(2)*(-vc(1) - vc(4)) 
     &+ fi(1)*fo(3)*gc(1)*(vc(1) - vc(4)) 
     &+ fi(4)*fo(2)*gc(2)*(vc(1) - vc(4))
   
         yiov(2,3) = yiov(2,3)
     &+ fi(3)*fo(2)*gc(2)*(-(cone*vc(2)) - vc(3)) 
     &+ fi(4)*fo(1)*gc(2)*(cone*vc(2) - vc(3)) 
     &+ fi(2)*fo(3)*gc(1)*(-(cone*vc(2)) + vc(3)) 
     &+ fi(1)*fo(4)*gc(1)*(cone*vc(2) + vc(3))
   
         yiov(2,4) = yiov(2,4)
     &+ fi(4)*gc(2)*(fo(2)*vc(2) - fo(1)*vc(4)) 
     &+ fi(3)*gc(2)*(-(fo(1)*vc(2)) - fo(2)*vc(4)) 
     &+ fi(2)*gc(1)*(-(fo(4)*vc(2)) + fo(3)*vc(4)) 
     &+ fi(1)*gc(1)*(fo(3)*vc(2) + fo(4)*vc(4))

         yiov(3,4) = yiov(3,4)
     &+fi(4)*gc(2)*(fo(2)*vc(3) + cone*fo(1)*vc(4)) 
     &+ fi(3)*gc(2)*(-(fo(1)*vc(3)) - cone*fo(2)*vc(4)) 
     &+ fi(2)*gc(1)*(-(fo(4)*vc(3)) - cone*fo(3)*vc(4)) 
     &+ fi(1)*gc(1)*(fo(3)*vc(3) + cone*fo(4)*vc(4))
   
         yiov(2,1) = yiov(1,2)
         yiov(3,1) = yiov(1,3)
         yiov(4,1) = yiov(1,4)
         yiov(3,2) = yiov(2,3)
         yiov(4,2) = yiov(2,4)
         yiov(4,3) = yiov(3,4)
   			
         do i = 1,4
   	    do j=1,4
            yiov(i,j) = -yiov(i,j)*d*gt
         end do
        end do
      
      else
         if ( tmass.gt.rZero ) then
            d =  -1.0d0/dcmplx( pt2-tmass**2, tmass*twidth )
   
   
            do i = 1,4
   	       do j=1,4
                  yiov(i,j) = fi(3)*gc(2)*(fo(2)*(-((KTVC*r4*MET(i,j)*
     &      (-pt(2) - cone*pt(3)))/(r3*Tmass**2)) - 
     &      (KTVC*r8*(-pt(2) - cone*pt(3))*pt(i)*pt(j))/(r3*Tmass**4) +
     &          (pt2*r4*MET(i,j)*(-vc(2) - cone*vc(3)))/(r3*Tmass**2) + 
     &          (pt2*r8*pt(i)*pt(j)*(-vc(2) - cone*vc(3)))/(r3*Tmass**4)
     & - (r4*pt(i)*pt(j)*(-vc(2) - cone*vc(3)))/Tmass**2 + 
     &    (rtwo*(-pt(2) - cone*pt(3))*(pt(j)*vc(i) + pt(i)*vc(j)))
     &/Tmass**2) + 
     &       fo(1)*(-((KTVC*r4*MET(i,j)*(pt(1) - pt(4)))/(r3*Tmass**2))
     & - (KTVC*r8*(pt(1) - pt(4))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(1) - vc(4)))/(r3*Tmass**2)
     & + (pt2*r8*pt(i)*pt(j)*(vc(1) - vc(4)))/(r3*Tmass**4) - 
     &          (r4*pt(i)*pt(j)*(vc(1) - vc(4)))/Tmass**2
     & + (rtwo*(pt(1) - pt(4))*(pt(j)*vc(i) + pt(i)*vc(j)))/Tmass**2))  
     &   + fi(2)*gc(1)*(fo(3)*(-((KTVC*r4*MET(i,j)*(pt(2) 
     &- cone*pt(3)))/(r3*Tmass**2)) - 
     &     (KTVC*r8*(pt(2) - cone*pt(3))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(2) - cone*vc(3)))/(r3*Tmass**2) + 
     &          (pt2*r8*pt(i)*pt(j)*(vc(2) - cone*vc(3)))/(r3*Tmass**4)
     & - (r4*pt(i)*pt(j)*(vc(2) - cone*vc(3)))/Tmass**2 + 
     &    (rtwo*(pt(2) - cone*pt(3))*(pt(j)*vc(i) 
     &+ pt(i)*vc(j)))/Tmass**2) + 
     &       fo(4)*(-((KTVC*r4*MET(i,j)*(pt(1) - pt(4)))/(r3*Tmass**2)) 
     &- (KTVC*r8*(pt(1) - pt(4))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(1) - vc(4)))/(r3*Tmass**2) 
     &+ (pt2*r8*pt(i)*pt(j)*(vc(1) - vc(4)))/(r3*Tmass**4) - 
     &          (r4*pt(i)*pt(j)*(vc(1) - vc(4)))/Tmass**2 
     &+ (rtwo*(pt(1) - pt(4))*(pt(j)*vc(i) + pt(i)*vc(j)))/Tmass**2)) + 
     &    fi(4)*gc(2)*(fo(1)*(-((KTVC*r4*MET(i,j)*(-pt(2) + cone*pt(3)))
     &/(r3*Tmass**2)) - 
     &      (KTVC*r8*(-pt(2) + cone*pt(3))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(-vc(2) + cone*vc(3)))/(r3*Tmass**2) + 
     &          (pt2*r8*pt(i)*pt(j)*(-vc(2) + cone*vc(3)))/(r3*Tmass**4) 
     &- (r4*pt(i)*pt(j)*(-vc(2) + cone*vc(3)))/Tmass**2 + 
     &          (rtwo*(-pt(2) + cone*pt(3))*(pt(j)*vc(i) 
     &+ pt(i)*vc(j)))/Tmass**2) + 
     &       fo(2)*(-((KTVC*r4*MET(i,j)*(pt(1) + pt(4)))/(r3*Tmass**2))
     & - (KTVC*r8*(pt(1) + pt(4))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(1) + vc(4)))/(r3*Tmass**2)
     & + (pt2*r8*pt(i)*pt(j)*(vc(1) + vc(4)))/(r3*Tmass**4) - 
     &          (r4*pt(i)*pt(j)*(vc(1) + vc(4)))/Tmass**2
     & + (rtwo*(pt(1) + pt(4))*(pt(j)*vc(i) + pt(i)*vc(j)))/Tmass**2)) + 
     &    fi(1)*gc(1)*(fo(4)*(-((KTVC*r4*MET(i,j)*(pt(2) + cone*pt(3)))
     &/(r3*Tmass**2)) - 
     &     (KTVC*r8*(pt(2) + cone*pt(3))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(2) + cone*vc(3)))/(r3*Tmass**2) + 
     &          (pt2*r8*pt(i)*pt(j)*(vc(2) + cone*vc(3)))/(r3*Tmass**4)
     & - (r4*pt(i)*pt(j)*(vc(2) + cone*vc(3)))/Tmass**2 + 
     &          (rtwo*(pt(2) + cone*pt(3))*(pt(j)*vc(i)
     & + pt(i)*vc(j)))/Tmass**2) + 
     &       fo(3)*(-((KTVC*r4*MET(i,j)*(pt(1) + pt(4)))/(r3*Tmass**2))
     & - (KTVC*r8*(pt(1) + pt(4))*pt(i)*pt(j))/(r3*Tmass**4) + 
     &          (pt2*r4*MET(i,j)*(vc(1) + vc(4)))/(r3*Tmass**2) 
     &+ (pt2*r8*pt(i)*pt(j)*(vc(1) + vc(4)))/(r3*Tmass**4) - 
     &          (r4*pt(i)*pt(j)*(vc(1) + vc(4)))/Tmass**2 
     &+ (rtwo*(pt(1) + pt(4))*(pt(j)*vc(i) + pt(i)*vc(j)))/Tmass**2))
   
               end do
            end do 
   
            yiov(1,1) = yiov(1,1)
     & +fi(1)*fo(3)*gc(1)*((4*KTVC*pt(1))/Tmass**2 - 2*rtwo*vc(1)) + 
     &  fi(2)*fo(4)*gc(1)*((4*KTVC*pt(1))/Tmass**2 - 2*rtwo*vc(1)) + 
     &  fi(3)*fo(1)*gc(2)*((4*KTVC*pt(1))/Tmass**2 - 2*rtwo*vc(1)) 
     &+ fi(4)*fo(2)*gc(2)*((4*KTVC*pt(1))/Tmass**2 - 2*rtwo*vc(1))
   
            yiov(2,2) = yiov(2,2)
     &	+fi(4)*fo(1)*gc(2)*((4*KTVC*pt(2))/Tmass**2 - 2*rtwo*vc(2)) + 
     &    fi(3)*fo(2)*gc(2)*((4*KTVC*pt(2))/Tmass**2 - 2*rtwo*vc(2)) + 
     &    fi(2)*fo(3)*gc(1)*((-4*KTVC*pt(2))/Tmass**2 + 2*rtwo*vc(2)) + 
     &    fi(1)*fo(4)*gc(1)*((-4*KTVC*pt(2))/Tmass**2 + 2*rtwo*vc(2))

            yiov(3,3) = yiov(3,3)
     &+fi(2)*fo(3)*gc(1)*((4*cone*KTVC*pt(3))/Tmass**2
     & - 2*cone*rtwo*vc(3)) + 
     &    fi(3)*fo(2)*gc(2)*((4*cone*KTVC*pt(3))/Tmass**2 
     &- 2*cone*rtwo*vc(3)) + 
     &    fi(1)*fo(4)*gc(1)*((-4*cone*KTVC*pt(3))/Tmass**2
     & + 2*cone*rtwo*vc(3)) + 
     &    fi(4)*fo(1)*gc(2)*((-4*cone*KTVC*pt(3))/Tmass**2
     & + 2*cone*rtwo*vc(3))
   
            yiov(4,4) = yiov(4,4)
     &+fi(2)*fo(4)*gc(1)*((4*KTVC*pt(4))/Tmass**2 - 2*rtwo*vc(4)) + 
     &    fi(3)*fo(1)*gc(2)*((4*KTVC*pt(4))/Tmass**2 - 2*rtwo*vc(4)) + 
     &    fi(1)*fo(3)*gc(1)*((-4*KTVC*pt(4))/Tmass**2 + 2*rtwo*vc(4)) + 
     &    fi(4)*fo(2)*gc(2)*((-4*KTVC*pt(4))/Tmass**2 + 2*rtwo*vc(4))
   
            yiov(1,2) = yiov(1,2)
     &	+fi(3)*gc(2)*(fo(2)*((2*KTVC*pt(1))/Tmass**2 - rtwo*vc(1))
     & + fo(1)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2))) + 
     &    fi(4)*gc(2)*(fo(1)*((2*KTVC*pt(1))/Tmass**2 - rtwo*vc(1)) 
     &+ fo(2)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2))) + 
     &    fi(1)*gc(1)*(fo(4)*((-2*KTVC*pt(1))/Tmass**2 + rtwo*vc(1))
     & + fo(3)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2))) + 
     &    fi(2)*gc(1)*(fo(3)*((-2*KTVC*pt(1))/Tmass**2 + rtwo*vc(1))
     & + fo(4)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2)))
   
            yiov(1,3) = yiov(1,3)
     &	+fi(3)*gc(2)*(fo(2)*((2*cone*KTVC*pt(1))/Tmass**2 
     &- cone*rtwo*vc(1)) + 
     &       fo(1)*((2*KTVC*pt(3))/Tmass**2 - rtwo*vc(3))) + 
     &    fi(4)*gc(2)*(fo(1)*((-2*cone*KTVC*pt(1))/Tmass**2 
     &+ cone*rtwo*vc(1)) + fo(2)*((2*KTVC*pt(3))/Tmass**2 
     &- rtwo*vc(3))) + 
     &    fi(1)*gc(1)*(fo(4)*((-2*cone*KTVC*pt(1))/Tmass**2
     & + cone*rtwo*vc(1)) + fo(3)*((2*KTVC*pt(3))/Tmass**2 
     &- rtwo*vc(3))) + 
     &    fi(2)*gc(1)*(fo(3)*((2*cone*KTVC*pt(1))/Tmass**2
     & -cone*rtwo*vc(1))+ fo(4)*((2*KTVC*pt(3))/Tmass**2 - rtwo*vc(3)))
   
            yiov(1,4) = yiov(1,4)
     &+fi(1)*fo(3)*gc(1)*((2*KTVC*(-pt(1) + pt(4)))/Tmass**2
     & - rtwo*(-vc(1) + vc(4))) + 
     &     fi(4)*fo(2)*gc(2)*((2*KTVC*(-pt(1) + pt(4)))/Tmass**2
     & - rtwo*(-vc(1) + vc(4))) + 
     &    fi(2)*fo(4)*gc(1)*((2*KTVC*(pt(1) + pt(4)))/Tmass**2 
     &- rtwo*(vc(1) + vc(4))) + 
     &    fi(3)*fo(1)*gc(2)*((2*KTVC*(pt(1) + pt(4)))/Tmass**2 
     &- rtwo*(vc(1) + vc(4)))
   
            yiov(2,3) = yiov(2,3)
     &+fi(1)*fo(4)*gc(1)*((2*KTVC*(-(cone*pt(2)) - pt(3)))/Tmass**2
     & - rtwo*(-(cone*vc(2)) - vc(3))) + 
     &    fi(2)*fo(3)*gc(1)*((2*KTVC*(cone*pt(2) - pt(3)))/Tmass**2
     & - rtwo*(cone*vc(2) - vc(3))) + 
     &    fi(4)*fo(1)*gc(2)*((2*KTVC*(-(cone*pt(2)) + pt(3)))/Tmass**2
     & - rtwo*(-(cone*vc(2)) + vc(3))) + 
     &    fi(3)*fo(2)*gc(2)*((2*KTVC*(cone*pt(2) + pt(3)))/Tmass**2 
     &- rtwo*(cone*vc(2) + vc(3)))
   
            yiov(2,4) = yiov(2,4)
     &+fi(4)*gc(2)*(fo(2)*((-2*KTVC*pt(2))/Tmass**2 + rtwo*vc(2))
     & + fo(1)*((2*KTVC*pt(4))/Tmass**2 - rtwo*vc(4))) + 
     &    fi(3)*gc(2)*(fo(1)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2))
     & + fo(2)*((2*KTVC*pt(4))/Tmass**2 - rtwo*vc(4))) + 
     &    fi(2)*gc(1)*(fo(4)*((2*KTVC*pt(2))/Tmass**2 - rtwo*vc(2)) 
     &+ fo(3)*((-2*KTVC*pt(4))/Tmass**2 + rtwo*vc(4))) + 
     &    fi(1)*gc(1)*(fo(3)*((-2*KTVC*pt(2))/Tmass**2 + rtwo*vc(2))
     & + fo(4)*((-2*KTVC*pt(4))/Tmass**2 + rtwo*vc(4)))

            yiov(3,4) = yiov(3,4)
     &	+fi(3)*gc(2)*(fo(1)*((2*KTVC*pt(3))/Tmass**2 - rtwo*vc(3)) + 
     &       fo(2)*((2*cone*KTVC*pt(4))/Tmass**2 - cone*rtwo*vc(4))) + 
     &    fi(2)*gc(1)*(fo(4)*((2*KTVC*pt(3))/Tmass**2 - rtwo*vc(3)) 
     &+ fo(3)*((2*cone*KTVC*pt(4))/Tmass**2 - cone*rtwo*vc(4))) + 
     &    fi(4)*gc(2)*(fo(2)*((-2*KTVC*pt(3))/Tmass**2 + rtwo*vc(3)) + 
     &       fo(1)*((-2*cone*KTVC*pt(4))/Tmass**2 + cone*rtwo*vc(4))) + 
     &    fi(1)*gc(1)*(fo(3)*((-2*KTVC*pt(3))/Tmass**2 + rtwo*vc(3)) 
     &+ fo(4)*((-2*cone*KTVC*pt(4))/Tmass**2 + cone*rtwo*vc(4)))

            yiov(2,1) = yiov(1,2)
            yiov(3,1) = yiov(1,3)
            yiov(4,1) = yiov(1,4)
            yiov(3,2) = yiov(2,3)
            yiov(4,2) = yiov(2,4)
            yiov(4,3) = yiov(3,4)
			
             do i = 1,4
            do j=1,4
               yiov(i,j) = -yiov(i,j)*d/2.0d0*gt
             end do
            end do
         
         else
            write(*,*) 'invalid tensor mass'
         end if
      end if
  
      uiov(1) = yiov(1,1)
      uiov(2) = yiov(1,2)
      uiov(3) = yiov(1,3)
      uiov(4) = yiov(1,4)
      uiov(5) = yiov(2,1)
      uiov(6) = yiov(2,2)
      uiov(7) = yiov(2,3)
      uiov(8) = yiov(2,4)
      uiov(9) = yiov(3,1)
      uiov(10) = yiov(3,2)
      uiov(11) = yiov(3,3)
      uiov(12) = yiov(3,4)
      uiov(13) = yiov(4,1)
      uiov(14) = yiov(4,2)
      uiov(15) = yiov(4,3)
      uiov(16) = yiov(4,4)
      uiov(17) = yiov(5,1)
      uiov(18) = yiov(6,1)

      return
      end
      subroutine uioxxx(fi,fo,gt,fmass,tmass,twidth , uio)
c
c This subroutine computes an off-shell tensor current from
c the fermion-anti-fermion-tensor boson coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gt             : coupling constant       gtf=-1/Lambda/4
c       real    fmass          : mass  of input fermion
c       real    tmass          : mass  of output tensor T
c       real    twidth         : width of output tensor T
c
c output:
c       complex uio(18)        : tensor current       j^mu^nu(<fo|T|fi>)
c
c- by Q.Li - OCT. 2006
c- Added massless tensor - P. de Aquino - Oct. 2009 
c
      implicit none
      double complex fi(6), fo(6), gt, uio(18)
      double precision fmass, tmass, twidth

      double complex yio(6,4)
      integer i,j
      double precision pi(4),po(4),km(4),kp(4)
      double precision MET(4,4)
      double complex cone,cZero, d, tt1,tt2,tt3
      double precision rZero, rTwo
      double precision KT2,K1KT,K2KT
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ), cone=(0.0d0,1.0d0) )

      
      yio(5,1) = -fi(5)+fo(5)
      yio(6,1) = -fi(6)+fo(6)
      
      pi(1) = dreal(fi(5))
      pi(2) = dreal(fi(6))
      pi(3) = dimag(fi(6))
      pi(4) = dimag(fi(5))
      
      po(1) = dreal(fo(5))
      po(2) = dreal(fo(6))
      po(3) = dimag(fo(6))
      po(4) = dimag(fo(5))
      
      km(1) = dreal(yio(5,1))
      km(2) = dreal(yio(6,1))
      km(3) = dimag(yio(6,1))
      km(4) = dimag(yio(5,1))

      kp(1) = po(1)+pi(1)
      kp(2) = po(2)+pi(2)
      kp(3) = po(3)+pi(3)
      kp(4) = po(4)+pi(4)
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      K1KT = pi(1)*km(1)-pi(2)*km(2)-pi(3)*km(3)-pi(4)*km(4)
      K2KT = po(1)*km(1)-po(2)*km(2)-po(3)*km(3)-po(4)*km(4)
      KT2 = km(1)**2-km(2)**2-km(3)**2-km(4)**2
      
      if ( tmass.eq.rZero ) then
         d = - gt/dcmplx( KT2, rZero )

         tt1 = fi(3)*(-4*fmass*fo(3) + fo(2)*(-kp(2) - cone*kp(3)) 
     &+ fo(1)*(kp(1) - kp(4)))  
     &+ fi(2)*(-4*fmass*fo(2) + fo(3)*(kp(2) - cone*kp(3)) 
     &+ fo(4)*(kp(1) - kp(4)))  
     &+ fi(4)*(-4*fmass*fo(4) + fo(1)*(-kp(2) + cone*kp(3)) 
     &+ fo(2)*(kp(1) + kp(4)))  
     &+ fi(1)*(-4*fmass*fo(1) + fo(4)*(kp(2) + cone*kp(3)) 
     &+ fo(3)*(kp(1) + kp(4)))
      
         do i = 1,4
            do j=1,4
               yio(i,j) = tt1*MET(i,j)
            end do
         enddo
	
         yio(1,1) = yio(1,1)
     &+ 2*fi(3)*fo(1)*kp(1) + 2*fi(4)*fo(2)*kp(1)  
     &+ 2*fi(1)*fo(3)*kp(1) + 2*fi(2)*fo(4)*kp(1)
	
         yio(2,2) = yio(2,2)
     &+ 2*fi(4)*fo(1)*kp(2) + 2*fi(3)*fo(2)*kp(2)  
     &- 2*fi(2)*fo(3)*kp(2) - 2*fi(1)*fo(4)*kp(2)
   
         yio(3,3) = yio(3,3)
     &- 2*cone*fi(4)*fo(1)*kp(3) + 2*cone*fi(3)*fo(2)*kp(3) 
     &+ 2*cone*fi(2)*fo(3)*kp(3) - 2*cone*fi(1)*fo(4)*kp(3)
        
         yio(4,4) = yio(4,4)
     &+ 2*fi(3)*fo(1)*kp(4) - 2*fi(4)*fo(2)*kp(4) 
     &- 2*fi(1)*fo(3)*kp(4) + 2*fi(2)*fo(4)*kp(4)	
   
         yio(1,2) = yio(1,2)
     &+ fi(3)*(fo(2)*kp(1) + fo(1)*kp(2)) 
     &+ fi(4)*(fo(1)*kp(1) + fo(2)*kp(2)) 
     &+ fi(1)*(-(fo(4)*kp(1))+ fo(3)*kp(2)) 
     &+ fi(2)*(-(fo(3)*kp(1)) + fo(4)*kp(2))
   
         yio(1,3) = yio(1,3)
     &+ fi(3)*(cone*fo(2)*kp(1) + fo(1)*kp(3)) 
     &+ fi(4)*(-(cone*fo(1)*kp(1)) + fo(2)*kp(3)) 
     &+ fi(1)*(-(cone*fo(4)*kp(1)) + fo(3)*kp(3)) 
     &+ fi(2)*(cone*fo(3)*kp(1) + fo(4)*kp(3))
     
         yio(1,4) = yio(1,4)
     &+ fi(4)*fo(2)*(-kp(1) + kp(4)) 
     &+ fi(1)*fo(3)*(-kp(1) + kp(4)) 
     &+ fi(3)*fo(1)*(kp(1) + kp(4))
     &+ fi(2)*fo(4)*(kp(1) + kp(4))
     
         yio(2,3) = yio(2,3)
     &+ fi(1)*fo(4)*(-(cone*kp(2)) - kp(3)) 
     &+ fi(2)*fo(3)*(cone*kp(2) - kp(3))  
     &+ fi(4)*fo(1)*(-(cone*kp(2)) + kp(3)) 
     &+ fi(3)*fo(2)*(cone*kp(2) + kp(3))
   
         yio(2,4) = yio(2,4)
     &+ fi(4)*(-(fo(2)*kp(2)) + fo(1)*kp(4)) 
     &+ fi(3)*(fo(1)*kp(2) + fo(2)*kp(4))  
     &+ fi(2)*(fo(4)*kp(2) - fo(3)*kp(4)) 
     &+ fi(1)*(-(fo(3)*kp(2)) - fo(4)*kp(4))
     
         yio(3,4) = yio(3,4)
     &+ fi(4)*(-(fo(2)*kp(3)) - cone*fo(1)*kp(4)) 
     &+ fi(3)*(fo(1)*kp(3) + cone*fo(2)*kp(4)) 
     &+ fi(2)*(fo(4)*kp(3) + cone*fo(3)*kp(4)) 
     &+ fi(1)*(-(fo(3)*kp(3)) - cone*fo(4)*kp(4))
   
         yio(2,1) = yio(1,2)
         yio(3,1) = yio(1,3)
         yio(4,1) = yio(1,4)         
         yio(3,2) = yio(2,3)
         yio(4,2) = yio(2,4)
         yio(4,3) = yio(3,4)
        
         do i = 1,4
            do j=1,4
               yio(i,j) = yio(i,j)*d
            end do
         enddo
   
      else if (tmass.gt.rZero) then
           d = - gt/dcmplx( KT2-tmass**2, tmass*twidth )
            
           tt1 = fi(3)*(((-8.0d0*fmass)/3.d0 
     &	+ (8.0d0*fmass*KT2)/(3.d0*tmass**2))*fo(3)
     & +  fo(2)*((4.0d0*(K1KT + K2KT)
     &            *(-km(2) - cone*km(3)))/(3.d0*tmass**2)
     & - (4.0d0*KT2*(-kp(2) - cone*kp(3)))/(3.d0*tmass**2)) + 
     &    fo(1)*((4.0d0*(K1KT + K2KT)*(km(1) - km(4)))/(3.d0*tmass**2) 
     & - (4.0d0*KT2*(kp(1) - kp(4)))/(3.d0*tmass**2))) + 
     &    fi(2)*(((-8.0d0*fmass)/3.d0 
     &        + (8.0d0*fmass*KT2)/(3.d0*tmass**2))*fo(2)  
     & +  fo(3)*((4.0d0*(K1KT + K2KT)*(km(2) 
     &          - cone*km(3)))/(3.d0*tmass**2) 
     &- (4.0d0*KT2*(kp(2) - cone*kp(3)))/(3.d0*tmass**2)) + 
     &    fo(4)*((4.0d0*(K1KT + K2KT)*(km(1) - km(4)))/(3.d0*tmass**2)
     & - (4.0d0*KT2*(kp(1) - kp(4)))/(3.d0*tmass**2))) + 
     &    fi(4)*(((-8.0d0*fmass)/3.d0 
     &+ (8.0d0*fmass*KT2)/(3.d0*tmass**2))*fo(4)  
     & +  fo(1)*((4.0d0*(K1KT + K2KT)
     &       *(-km(2) + cone*km(3)))/(3.d0*tmass**2)
     & - (4.0d0*KT2*(-kp(2) + cone*kp(3)))/(3.d0*tmass**2)) + 
     &    fo(2)*((4.0d0*(K1KT + K2KT)*(km(1) + km(4)))/(3.d0*tmass**2)
     & - (4.0d0*KT2*(kp(1) + kp(4)))/(3.d0*tmass**2))) + 
     &    fi(1)*(((-8.0d0*fmass)/3.d0 
     &+ (8.0d0*fmass*KT2)/(3.d0*tmass**2))*fo(1) 
     & +  fo(4)*((4.0d0*(K1KT + K2KT)*(km(2)
     &      + cone*km(3)))/(3.d0*tmass**2) 
     &- (4.0d0*KT2*(kp(2) + cone*kp(3)))/(3.d0*tmass**2)) + 
     &    fo(3)*((4.0d0*(K1KT + K2KT)
     &*(km(1) + km(4)))/(3.d0*tmass**2)
     & - (4.0d0*KT2*(kp(1) + kp(4)))/(3.d0*tmass**2)))
 
      
           tt2 = (fi(3)*(((-16.0d0*KT2*fmass)/(3.0D0*tmass**4)
     &	 - (16.d0*fmass)/(3.d0*tmass**2))*fo(3) + 
     &       fo(2)*((8.d0*(K1KT + K2KT)*(-km(2)
     & - cone*km(3)))/(3.d0*tmass**4) 
     &+ (8.d0*KT2*(-kp(2) - cone*kp(3)))/(3.0D0*tmass**4) + 
     &          (4.d0*(-kp(2) - cone*kp(3)))/tmass**2) + 
     &       fo(1)*((8.d0*(K1KT + K2KT)*(km(1)
     & - km(4)))/(3.d0*tmass**4) 
     &+ (8.d0*KT2*(kp(1) - kp(4)))/(3.d0*tmass**4) + 
     &          (4.d0*(kp(1) - kp(4)))/tmass**2)) + fi(2)*
     &     (((-16.d0*KT2*fmass)/(3.d0*tmass**4) 
     &- (16.d0*fmass)/(3.d0*tmass**2))*fo(2) + 
     &       fo(3)*((8.d0*(K1KT + K2KT)*(km(2)
     & - cone*km(3)))/(3.d0*tmass**4) 
     &+ (8.d0*KT2*(kp(2) - cone*kp(3)))/(3.d0*tmass**4) + 
     &          (4.d0*(kp(2) - cone*kp(3)))/tmass**2) + 
     &       fo(4)*((8.d0*(K1KT + K2KT)*(km(1) - km(4)))/(3.d0*tmass**4)
     & + (8.d0*KT2*(kp(1) - kp(4)))/(3.d0*tmass**4) + 
     &          (4.d0*(kp(1) - kp(4)))/tmass**2)) + fi(4)*
     &     (((-16.d0*KT2*fmass)/(3.d0*tmass**4) 
     &- (16.d0*fmass)/(3.d0*tmass**2))*fo(4) + 
     &       fo(1)*((8.d0*(K1KT + K2KT)*(-km(2) 
     &+ cone*km(3)))/(3.d0*tmass**4)
     & + (8.d0*KT2*(-kp(2) + cone*kp(3)))/(3.d0*tmass**4) + 
     &          (4.d0*(-kp(2) + cone*kp(3)))/tmass**2) + 
     &       fo(2)*((8.d0*(K1KT + K2KT)*(km(1) 
     &+ km(4)))/(3.d0*tmass**4)
     & + (8.d0*KT2*(kp(1) + kp(4)))/(3.d0*tmass**4) + 
     &          (4.d0*(kp(1) + kp(4)))/tmass**2)) + fi(1)*
     &     (((-16.d0*KT2*fmass)/(3.d0*tmass**4) 
     &- (16.d0*fmass)/(3.d0*tmass**2))*fo(1) + 
     &       fo(4)*((8.d0*(K1KT + K2KT)*(km(2)
     & + cone*km(3)))/(3.d0*tmass**4)
     & + (8.d0*KT2*(kp(2) + cone*kp(3)))/(3.d0*tmass**4) + 
     &          (4.d0*(kp(2) + cone*kp(3)))/tmass**2) + 
     &       fo(3)*((8.d0*(K1KT + K2KT)
     &*(km(1) + km(4)))/(3.d0*tmass**4)
     & + (8.d0*KT2*(kp(1) + kp(4)))/(3.d0*tmass**4) + 
     &          (4.d0*(kp(1) + kp(4)))/tmass**2)))
   
           tt3 = fi(3)*((-2.0d0*fo(2)*(-km(2) - cone*km(3)))/tmass**2
     &	 - (2.0d0*fo(1)*(km(1) - km(4)))/tmass**2) + 
     &    fi(2)*((-2.0d0*fo(3)*(km(2) - cone*km(3)))/tmass**2 
     &- (2.0d0*fo(4)*(km(1) - km(4)))/tmass**2) + 
     &    fi(4)*((-2.0d0*fo(1)*(-km(2) + cone*km(3)))/tmass**2
     & - (2.0d0*fo(2)*(km(1) + km(4)))/tmass**2) + 
     &    fi(1)*((-2.0d0*fo(4)*(km(2) + cone*km(3)))/tmass**2 
     &- (2.0d0*fo(3)*(km(1) + km(4)))/tmass**2) 
         
           do i = 1,4
              do j=1,4
                 yio(i,j) = tt1*MET(i,j)+tt2*km(i)*km(j)
     &+tt3*(km(i)*kp(j)+km(j)*kp(i))
              end do
           enddo
   	
           yio(1,1) = yio(1,1)
     &+fi(3)*fo(1)*((-4.0d0*(K1KT + K2KT)*km(1))/tmass**2 + 4.0d0*kp(1))  
     & +fi(4)*fo(2)*((-4.d0*(K1KT + K2KT)*km(1))/tmass**2 + 4.0d0*kp(1)) 
     & +fi(1)*fo(3)*((-4.0d0*(K1KT + K2KT)*km(1))/tmass**2+ 4.0d0*kp(1)) 
     &+ fi(2)*fo(4)*((-4.0d0*(K1KT + K2KT)*km(1))/tmass**2+ 4.0d0*kp(1))
   	
           yio(2,2) = yio(2,2)
     &+fi(2)*fo(3)*((4.0d0*(K1KT + K2KT)*km(2))/tmass**2 - 4.d0*kp(2))
     &+fi(1)*fo(4)*((4.0d0*(K1KT + K2KT)*km(2))/tmass**2 - 4.d0*kp(2))  
     &+fi(4)*fo(1)*((-4.d0*(K1KT + K2KT)*km(2))/tmass**2 + 4.d0*kp(2))
     &+fi(3)*fo(2)*((-4.d0*(K1KT + K2KT)*km(2))/tmass**2 + 4.d0*kp(2))
      
           yio(3,3) = yio(3,3)
     &+fi(4)*fo(1)*((4.0d0*cone*(K1KT + K2KT)*km(3))/tmass**2
     & - 4.0d0*cone*kp(3)) + 
     &    fi(1)*fo(4)*((4.0d0*cone*(K1KT + K2KT)*km(3))/tmass**2
     & - 4.0d0*cone*kp(3)) + 
     &    fi(3)*fo(2)*((-4.0d0*cone*(K1KT + K2KT)*km(3))/tmass**2
     & + 4.d0*cone*kp(3)) + 
     &    fi(2)*fo(3)*((-4.d0*cone*(K1KT + K2KT)*km(3))/tmass**2
     & + 4.d0*cone*kp(3))	
      	
           yio(4,4) = yio(4,4)
     &+fi(4)*fo(2)*((4.d0*(K1KT + K2KT)*km(4))/tmass**2 - 4.d0*kp(4)) 
     &+fi(1)*fo(3)*((4.d0*(K1KT + K2KT)*km(4))/tmass**2 - 4.d0*kp(4)) 
     &+fi(3)*fo(1)*((-4.d0*(K1KT + K2KT)*km(4))/tmass**2+ 4.d0*kp(4))
     &+fi(2)*fo(4)*((-4.d0*(K1KT + K2KT)*km(4))/tmass**2+ 4.d0*kp(4))	
         
           yio(1,2) = yio(1,2)
     &+fi(3)*(fo(2)*((-2.d0*(K1KT + K2KT)*km(1))/tmass**2 + 2.d0*kp(1)) 
     &+ fo(1)*((-2.d0*(K1KT + K2KT)*km(2))/tmass**2 + 2.d0*kp(2))) + 
     &  fi(4)*(fo(1)*((-2.d0*(K1KT + K2KT)*km(1))/tmass**2 + 2.d0*kp(1))
     & + fo(2)*((-2.d0*(K1KT + K2KT)*km(2))/tmass**2 + 2.d0*kp(2))) + 
     &  fi(1)*(fo(4)*((2.d0*(K1KT + K2KT)*km(1))/tmass**2 - 2.d0*kp(1))
     & + fo(3)*((-2.d0*(K1KT + K2KT)*km(2))/tmass**2 + 2.d0*kp(2))) + 
     & fi(2)*(fo(3)*((2.d0*(K1KT + K2KT)*km(1))/tmass**2 - 2.d0*kp(1))
     & + fo(4)*((-2.d0*(K1KT + K2KT)*km(2))/tmass**2 + 2.d0*kp(2)))
      
           yio(1,3) = yio(1,3)
     &+fi(3)*(fo(2)*((-2.0d0*cone*(K1KT + K2KT)*km(1))/tmass**2
     & + 2.d0*cone*kp(1)) + 
     &    fo(1)*((-2.d0*(K1KT + K2KT)*km(3))/tmass**2 +2.d0*kp(3))) + 
     &    fi(4)*(fo(1)*((2.d0*cone*(K1KT + K2KT)*km(1))/tmass**2
     & - 2.d0*cone*kp(1)) + 
     &    fo(2)*((-2.d0*(K1KT + K2KT)*km(3))/tmass**2 +2.d0*kp(3))) + 
     &    fi(1)*(fo(4)*((2.d0*cone*(K1KT + K2KT)*km(1))/tmass**2 
     &- 2.d0*cone*kp(1)) + 
     &    fo(3)*((-2.d0*(K1KT + K2KT)*km(3))/tmass**2 +2.d0*kp(3))) + 
     &    fi(2)*(fo(3)*((-2.d0*cone*(K1KT + K2KT)*km(1))/tmass**2 
     &+ 2.d0*cone*kp(1)) + 
     &    fo(4)*((-2.d0*(K1KT + K2KT)*km(3))/tmass**2 + 2.d0*kp(3)))
        
           yio(1,4) = yio(1,4)
     &+fi(4)*fo(2)*((-2.d0*(K1KT + K2KT)*(-km(1) + km(4)))/tmass**2
     & - 2.d0*kp(1) + 2.d0*kp(4)) + 
     & fi(1)*fo(3)*((-2.d0*(K1KT + K2KT)*(-km(1) + km(4)))/tmass**2
     & - 2.d0*kp(1) + 2.d0*kp(4)) + 
     & fi(3)*fo(1)*((-2.d0*(K1KT + K2KT)*(km(1) + km(4)))/tmass**2 
     &+ 2.d0*kp(1) + 2.d0*kp(4)) + 
     & fi(2)*fo(4)*((-2.d0*(K1KT + K2KT)*(km(1) + km(4)))/tmass**2
     & + 2.d0*kp(1) + 2.d0*kp(4))
      
           yio(2,3) = yio(2,3)
     &+fi(1)*fo(4)*((-2.0d0*(K1KT + K2KT)*(-(cone*km(2))
     & - km(3)))/tmass**2 - 2.0d0*cone*kp(2) - 2.d0*kp(3)) + 
     &    fi(2)*fo(3)*((-2.d0*(K1KT + K2KT)*(cone*km(2)
     & - km(3)))/tmass**2 + 2.d0*cone*kp(2) - 2.d0*kp(3)) + 
     &    fi(4)*fo(1)*((-2.d0*(K1KT + K2KT)*(-(cone*km(2))
     & + km(3)))/tmass**2 - 2.d0*cone*kp(2) + 2.d0*kp(3)) + 
     &    fi(3)*fo(2)*((-2.d0*(K1KT + K2KT)*(cone*km(2)
     & + km(3)))/tmass**2 + 2.d0*cone*kp(2) + 2.d0*kp(3))
      
      
           yio(2,4) = yio(2,4)
     &+fi(2)*(fo(4)*((-2.d0*(K1KT + K2KT)*km(2))/tmass**2+2.d0*kp(2)) + 
     &       fo(3)*((2.d0*(K1KT + K2KT)*km(4))/tmass**2 - 2.d0*kp(4))) + 
     &  fi(1)*(fo(3)*((2.d0*(K1KT + K2KT)*km(2))/tmass**2 - 2.d0*kp(2))
     & + fo(4)*((2.d0*(K1KT + K2KT)*km(4))/tmass**2 - 2.d0*kp(4))) + 
     &  fi(4)*(fo(2)*((2.d0*(K1KT + K2KT)*km(2))/tmass**2 - 2.d0*kp(2))
     & + fo(1)*((-2.d0*(K1KT + K2KT)*km(4))/tmass**2 + 2.d0*kp(4))) + 
     & fi(3)*(fo(1)*((-2.d0*(K1KT + K2KT)*km(2))/tmass**2 + 2.d0*kp(2))
     & + fo(2)*((-2.d0*(K1KT + K2KT)*km(4))/tmass**2 + 2.d0*kp(4)))
         
           yio(3,4) = yio(3,4)
     &+fi(4)*(fo(2)*((2.d0*(K1KT + K2KT)*km(3))/tmass**2- 2.d0*kp(3)) + 
     &       fo(1)*((2.d0*cone*(K1KT + K2KT)*km(4))/tmass**2
     & - 2.d0*cone*kp(4))) + 
     &    fi(1)*(fo(3)*((2.d0*(K1KT + K2KT)*km(3))/tmass**2
     & - 2.d0*kp(3)) + 
     &       fo(4)*((2.d0*cone*(K1KT + K2KT)*km(4))/tmass**2
     & - 2.d0*cone*kp(4))) + 
     &    fi(3)*(fo(1)*((-2.d0*(K1KT + K2KT)*km(3))/tmass**2
     & + 2.d0*kp(3)) + 
     &       fo(2)*((-2.d0*cone*(K1KT + K2KT)*km(4))/tmass**2
     & + 2.d0*cone*kp(4))) + 
     &    fi(2)*(fo(4)*((-2.d0*(K1KT + K2KT)*km(3))/tmass**2
     & + 2.d0*kp(3)) + 
     &       fo(3)*((-2.d0*cone*(K1KT + K2KT)*km(4))/tmass**2
     &+ 2.d0*cone*kp(4)))
      
           yio(2,1) = yio(1,2)
           yio(3,1) = yio(1,3)
           yio(4,1) = yio(1,4)
            
           yio(3,2) = yio(2,3)
           yio(4,2) = yio(2,4)
           yio(4,3) = yio(3,4)
           
           do i = 1,4
              do j=1,4
                 yio(i,j) = yio(i,j)*d/2.0d0
              end do
           enddo            
                        
         else 
            write(*,*) 'invalid tensor mass'
            stop
      end if

      uio(1) = yio(1,1)
      uio(2) = yio(1,2)
      uio(3) = yio(1,3)
      uio(4) = yio(1,4)
      uio(5) = yio(2,1)
      uio(6) = yio(2,2)
      uio(7) = yio(2,3)
      uio(8) = yio(2,4)
      uio(9) = yio(3,1)
      uio(10) = yio(3,2)
      uio(11) = yio(3,3)
      uio(12) = yio(3,4)
      uio(13) = yio(4,1)
      uio(14) = yio(4,2)
      uio(15) = yio(4,3)
      uio(16) = yio(4,4)
      uio(17) = yio(5,1)
      uio(18) = yio(6,1)

      return
      end
      subroutine upsxxx(t2,sc,gt, xm,xw,t1)

c  Subroutines for graviton phase space integration
c  KEK 2009.11
c
      implicit none
      double complex t1(18), t2(18), sc(3), vertex, tc2(18)
      double complex gt(2)
      double precision xm,xw,xmass

      double complex ft(6,4),ft2(6,4)
      double precision ps1(4), pt2(4),PT(4),PTT2,PG(0:3)
      integer i
      double complex cZero
      double precision rZero, rTwo,Pi
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )
      double precision PADD
      double precision L_ADD,NADD,MGLOW,MGUP
      external txxxxx

      Pi=dacos(-1.0d0)
      L_ADD=dimag(gt(1))
      NADD=dreal(gt(1))
      MGUP=dimag(gt(2))
      MGLOW=dreal(gt(2))


      ps1(1) = dreal(sc(2))
      ps1(2) = dreal(sc(3))
      ps1(3) = dimag(sc(3))
      ps1(4) = dimag(sc(2))

      pt2(1) = -dreal(t2(17))
      pt2(2) = -dreal(t2(18))
      pt2(3) = -dimag(t2(18))
      pt2(4) = -dimag(t2(17))

      PG(0)=ps1(1)-pt2(1)
      PG(1)=ps1(2)-pt2(2)
      PG(2)=ps1(3)-pt2(3)
      PG(3)=ps1(4)-pt2(4)

      PTT2=PG(0)**2-PG(1)**2-PG(2)**2-PG(3)**2
      xmass=dsqrt(PTT2)

       t1(17) = dcmplx(PG(0),PG(3))
       t1(18) = dcmplx(PG(1), PG(2))

      if(xmass.lt.MGLOW.or.xmass.gt.MGUP) then
      do i=1,16
      t1(i)=dcmplx(0.0d0,0.0d0)
      enddo
      return
      endif

      CALL txxxxx(PG,xmass,INT(t2(1)),+1 , t1)


       if(INT(NADD).eq.2) then
         PADD=2.0d0*Pi
        elseif(INT(NADD).eq.3) then
         PADD=4.0d0*Pi
        elseif(INT(NADD).eq.4) then
         PADD=2.0d0*Pi**2
        elseif(INT(NADD).eq.5) then
          PADD=8.0d0/3.0d0*Pi**2
        elseif(INT(NADD).eq.6) then
           PADD=Pi**3
        else
        print *, "OUT CASE"
        stop
        endif 

 
        do i=1,16
        t1(i)=-1.0d0*t1(i)*
     & dsqrt( 
     & 2.0d0*Pi*8.0d0*Pi  ! to compensate the decay phase factor
     &* PADD/L_ADD**NADD*xmass**(NADD-1)  ! density factor for d=4 case
     &/2.0d0/xmass)   ! dm=dm^2/2/m    
        enddo

      return
      end
      subroutine usslxx(s1,s2,gc,xm,xw , uss)
c- by RF - Mar. 2006
c
c This subroutine computes an internal particle current from the three-
c scalar coupling.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex gc             : coupling constant                  ghhh
c
c output:
c       complex uss(3)         : internal particle current (scalar)
c    
c not used:
c       xm,xw
c

      implicit none
      include "dimension.inc"
      double complex s1(DIM),s2(DIM),uss(DIM)
      double precision xm,xw,gc

      uss(2) = s1(2)+s2(2)
      uss(3) = s1(3)+s2(3)

c the internal particle does not propagate, so no multiplication
c with a propagator necessary.

      uss(1) = - gc*s1(1)*s2(1)
c
      return
      end
      subroutine ussxxx(s1,s2,gt,smass,tmass,twidth , uss)
c
c This subroutine computes an off-shell tensor current 
c from the scalar-scalar-tensor boson coupling.
c
c input:
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex gt             : coupling constant         gts=-1/Lambda
c       real    smass          : scalar mass                         m_s
c       real    tmass          : mass  of output tensor T 
c       real    twidth         : width of output tensor T
c
c output:
c       complex uss(18)        : tensor current         j^mu^nu(T:s1,s2)
c     
c- by Q.Li - OCT. 2006
c
      implicit none
      double complex s1(3), s2(3), gt, uss(18)
      double precision smass, tmass, twidth

      integer i,j
      double complex yss(6,4)
      double precision ps1(4), ps2(4), pT(4)
      double precision MET(4,4)
      double complex cZero, d
      double precision rZero, rTwo
      double precision p1p2,pT2,p1pT,p2pT
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      yss(5,1) = s1(2)+s2(2)
      yss(6,1) = s1(3)+s2(3)
      
      ps1(1) = dreal(s1(2))
      ps1(2) = dreal(s1(3))
      ps1(3) = dimag(s1(3))
      ps1(4) = dimag(s1(2))
      
      ps2(1) = -dreal(s2(2))
      ps2(2) = -dreal(s2(3))
      ps2(3) = -dimag(s2(3))
      ps2(4) = -dimag(s2(2))
      
      pT(1) = dreal(yss(5,1))
      pT(2) = dreal(yss(6,1))
      pT(3) = dimag(yss(6,1))
      pT(4) = dimag(yss(5,1))
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      p1p2 = ps1(1)*ps2(1)-ps1(2)*ps2(2)-ps1(3)*ps2(3)-ps1(4)*ps2(4)
      p1pT = ps1(1)*pT(1)-ps1(2)*pT(2)-ps1(3)*pT(3)-ps1(4)*pT(4)
      p2pT = pT(1)*ps2(1)-pT(2)*ps2(2)-pT(3)*ps2(3)-pT(4)*ps2(4)
      pT2  = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2
      
      if ( tmass.gt.rZero ) then
         d = - gt/dcmplx( pT2-tmass**2, tmass*twidth )
      else
         d = - gt/dcmplx( pT2, rZero )
      end if
    
      do i = 1,4
         do j=1,4
            yss(i,j) = -2.0d0/3.0d0*MET(i,j)*smass**2
     &	-4.0d0/3.0d0*smass**2/tmass**2*pT(i)*pT(j)
     &    +2.0d0/3.0d0*smass**2/tmass**2*PT2*MET(i,j)
     &    +4.0d0/3.0d0*smass**2/tmass**4*PT2*pT(i)*pT(j)
     &    +2.0d0*(ps1(i)*ps2(j)+ps1(j)*ps2(i))
     &    -2.0d0/tmass**2*p1pT*(ps2(i)*pT(j)+ps2(j)*pT(i)) 
     &    -2.0d0/tmass**2*p2pT*(ps1(i)*pT(j)+ps1(j)*pT(i)) 
     &    +4.0d0/3.0d0/tmass**2*MET(i,j)*p1pT*p2pT
     &    +8.0d0/3.0d0/tmass**4*p1pT*p2pT*pT(i)*pT(j)
     &    -2.0d0/3.0d0*p1p2*MET(i,j)
     &    +8.0d0/3.0d0/tmass**2*p1p2*pT(i)*pT(j)
     &    -2.0d0/3.0d0/tmass**2*p1p2*pT2*MET(i,j)
     &    -4.0d0/3.0d0/tmass**4*PT2*p1p2*pT(i)*pT(j)

            yss(i,j) = yss(i,j)*d*s1(1)*s2(1)/2.0d0

         end do
      enddo

      uss(1) = yss(1,1)
      uss(2) = yss(1,2)
      uss(3) = yss(1,3)
      uss(4) = yss(1,4)
      uss(5) = yss(2,1)
      uss(6) = yss(2,2)
      uss(7) = yss(2,3)
      uss(8) = yss(2,4)
      uss(9) = yss(3,1)
      uss(10) = yss(3,2)
      uss(11) = yss(3,3)
      uss(12) = yss(3,4)
      uss(13) = yss(4,1)
      uss(14) = yss(4,2)
      uss(15) = yss(4,3)
      uss(16) = yss(4,4)
      uss(17) = yss(5,1)
      uss(18) = yss(6,1)

      return
      end
      subroutine utsaxx(tc1,sc,gt,xm,xw,jts)
c
c- by RF - Feb. 2006 
c
c This subroutine computes an off-shell tensor current from the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex sc(3)             : Incoming scalar particle
c          complex gt(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c
c     output:
c          complex jts               : off-shell tensor current
c
c     not used:
c          xm, xw
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),jts(DIM),sc(DIM), gt(2)
      double precision xm, xw

c The outgoing tensor current is the same as the incoming multiplied by the
c coupling constant and the scalar particle. Note that the diagonal tensor
c terms are always zero because the tensor particle is anti-symmetric. The 
c tensor particle does not propagate, thus no multiplication with the tensor 
c propagator.

      jts(17) = sc(2) + tc1(17)
      jts(18) = sc(3) + tc1(18)

 
      if (gt(1).NE.(0D0,0D0)) then

         jts( 1) = (0D0,0D0)        !-gt(1) * sc(1) * tc1( 1)
         jts( 2) =  -gt(1) * sc(1) * tc1( 2)
         jts( 3) =  -gt(1) * sc(1) * tc1( 3)
         jts( 4) =  -gt(1) * sc(1) * tc1( 4)

         jts( 5) =  -gt(1) * sc(1) * tc1( 5)
         jts( 6) = (0D0,0D0)        !-gt(1) * sc(1) * tc1( 6)
         jts( 7) =  -gt(1) * sc(1) * tc1( 7)
         jts( 8) =  -gt(1) * sc(1) * tc1( 8)

         jts( 9) =  -gt(1) * sc(1) * tc1( 9)
         jts(10) =  -gt(1) * sc(1) * tc1(10)
         jts(11) = (0D0,0D0)        !-gt(1) * sc(1) * tc1(11)
         jts(12) =  -gt(1) * sc(1) * tc1(12)
         
         jts(13) =  -gt(1) * sc(1) * tc1(13)
         jts(14) =  -gt(1) * sc(1) * tc1(14)
         jts(15) =  -gt(1) * sc(1) * tc1(15)
         jts(16) = (0D0,0D0)        !-gt(1) * sc(1) * tc1(16)

      else
         jts( 1)=(0D0,0D0)
         jts( 2)=(0D0,0D0)
         jts( 3)=(0D0,0D0)
         jts( 4)=(0D0,0D0)
         jts( 5)=(0D0,0D0)
         jts( 6)=(0D0,0D0)
         jts( 7)=(0D0,0D0)
         jts( 8)=(0D0,0D0)
         jts( 9)=(0D0,0D0)
         jts(10)=(0D0,0D0)
         jts(11)=(0D0,0D0)
         jts(12)=(0D0,0D0)
         jts(13)=(0D0,0D0)
         jts(14)=(0D0,0D0)
         jts(15)=(0D0,0D0)
         jts(16)=(0D0,0D0)
      endif


      return
      end
      subroutine utscxx(tc1,sc,gt,xm,xw,jts)
c
c- by RF - Feb. 2006 
c  CP3  Modified Nov. 2009 
c
c This subroutine computes an off-shell tensor current from the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex sc(3)             : Incoming scalar particle
c          complex gt(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c
c     output:
c          complex jts               : off-shell tensor current
c
c     not used:
c          xm, xw
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),jts(DIM),sc(DIM), gt(2),t1(6,4),t2(6,4)
      double precision xm, xw
      integer m1,m3  
      double precision MT(4,4)

c The outgoing tensor current is the same as the incoming multiplied by the
c coupling constant and the scalar particle. Note that the diagonal tensor
c terms are always zero because the tensor particle is anti-symmetric. The 
c tensor particle does not propagate, thus no multiplication with the tensor 
c propagator.

      jts(17) = sc(2) + tc1(17)
      jts(18) = sc(3) + tc1(18)

      T2(1,1) = tc1(1)
      T2(1,2) = tc1(2)
      T2(1,3) = tc1(3)
      T2(1,4) = tc1(4)
      T2(2,1) = tc1(5)
      T2(2,2) = tc1(6)
      T2(2,3) = tc1(7)
      T2(2,4) = tc1(8)
      T2(3,1) = tc1(9)
      T2(3,2) = tc1(10)
      T2(3,3) = tc1(11)
      T2(3,4) = tc1(12)
      T2(4,1) = tc1(13)
      T2(4,2) = tc1(14)
      T2(4,3) = tc1(15)
      T2(4,4) = tc1(16)
      T2(5,1) = tc1(17)
      T2(6,1) = tc1(18)

      do m1=1,4
         do m3=1,4
            MT(m1,m3) = 0.0d0
         enddo 
      enddo
      MT(1,1) =  1.0d0
      MT(2,2) = -1.0d0
      MT(3,3) = -1.0d0
      MT(4,4) = -1.0d0

 
      if (gt(1).NE.(0D0,0D0)) then

       do  m1=1,4
         do m3=1,4
        T1(m1,m3)=gt(1)*(-(MT(1,m3)*MT(2,m1)*SC(1)*T2(1,2)) + 
     -  MT(1,m1)*MT(2,m3)*SC(1)*T2(1,2) - 
     -  MT(1,m3)*MT(3,m1)*SC(1)*T2(1,3) + 
     -  MT(1,m1)*MT(3,m3)*SC(1)*T2(1,3) - 
     -  MT(1,m3)*MT(4,m1)*SC(1)*T2(1,4) + 
     -  MT(1,m1)*MT(4,m3)*SC(1)*T2(1,4) + 
     -  MT(1,m3)*MT(2,m1)*SC(1)*T2(2,1) - 
     -  MT(1,m1)*MT(2,m3)*SC(1)*T2(2,1) + 
     -  MT(2,m3)*MT(3,m1)*SC(1)*T2(2,3) - 
     -  MT(2,m1)*MT(3,m3)*SC(1)*T2(2,3) + 
     -  MT(2,m3)*MT(4,m1)*SC(1)*T2(2,4) - 
     -  MT(2,m1)*MT(4,m3)*SC(1)*T2(2,4) + 
     -  MT(1,m3)*MT(3,m1)*SC(1)*T2(3,1) - 
     -  MT(1,m1)*MT(3,m3)*SC(1)*T2(3,1) - 
     -  MT(2,m3)*MT(3,m1)*SC(1)*T2(3,2) + 
     -  MT(2,m1)*MT(3,m3)*SC(1)*T2(3,2) + 
     -  MT(3,m3)*MT(4,m1)*SC(1)*T2(3,4) - 
     -  MT(3,m1)*MT(4,m3)*SC(1)*T2(3,4) + 
     -  MT(1,m3)*MT(4,m1)*SC(1)*T2(4,1) - 
     -  MT(1,m1)*MT(4,m3)*SC(1)*T2(4,1) - 
     -  MT(2,m3)*MT(4,m1)*SC(1)*T2(4,2) + 
     -  MT(2,m1)*MT(4,m3)*SC(1)*T2(4,2) - 
     -  MT(3,m3)*MT(4,m1)*SC(1)*T2(4,3) + 
     -  MT(3,m1)*MT(4,m3)*SC(1)*T2(4,3))
       enddo
       enddo  

       jts(1) = T1(1,1)
       jts(2) = T1(1,2)
       jts(3) = T1(1,3)
       jts(4) = T1(1,4)
       jts(5) = T1(2,1)
       jts(6) = T1(2,2)
       jts(7) = T1(2,3)
       jts(8) = T1(2,4)
       jts(9) = T1(3,1)
       jts(10) = T1(3,2)
       jts(11) = T1(3,3)
       jts(12) = T1(3,4)
       jts(13) = T1(4,1)
       jts(14) = T1(4,2)
       jts(15) = T1(4,3)
       jts(16) = T1(4,4)


      else
         jts( 1)=(0D0,0D0)
         jts( 2)=(0D0,0D0)
         jts( 3)=(0D0,0D0)
         jts( 4)=(0D0,0D0)
         jts( 5)=(0D0,0D0)
         jts( 6)=(0D0,0D0)
         jts( 7)=(0D0,0D0)
         jts( 8)=(0D0,0D0)
         jts( 9)=(0D0,0D0)
         jts(10)=(0D0,0D0)
         jts(11)=(0D0,0D0)
         jts(12)=(0D0,0D0)
         jts(13)=(0D0,0D0)
         jts(14)=(0D0,0D0)
         jts(15)=(0D0,0D0)
         jts(16)=(0D0,0D0)
      endif


      return
      end
      subroutine utssxx(tc1,sc1,sc2,g1,g2,xm,xw,jts)
c
c- by RF - Mar. 2006 
c
c This subroutine computes an off-shell tensor current from the ttss coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex sc1(3)            : Incoming scalar particle (Higgs)
c          complex sc2(3)            : second incoming scalar particle (Higgs)
c          complex g1(2)             : coupling constant (Higgs effc. theor)
c          real    g2                : coupling constant (include extra Higgs)
c
c     output:
c          complex jts               : off-shell tensor current
c
c     not used:
c          xm, xw
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),jts(DIM),sc1(DIM),sc2(DIM),g1(2)
      double precision g2, xm, xw

c The outgoing tensor current is the same as the incoming multiplied by the
c coupling constants and the scalar particles.
c Note that the diagonal tensor terms are always zero because
c the tensor particle is anti-symmetric.


      jts(17) = sc1(2) + sc2(2) + tc1(17)
      jts(18) = sc1(3) + sc2(3) + tc1(18)

 
      if (g1(1).NE.(0D0,0D0)) then

         jts( 1) = (0D0,0D0)   ! g1(1)* g2 * sc1(1) * sc2(1) * tc1( 1)
         jts( 2) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1( 2)
         jts( 3) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1( 3)
         jts( 4) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1( 4)

         jts( 5) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1( 5)
         jts( 6) = (0D0,0D0)   ! g1(1)* g2 * sc1(1) * sc2(1) * tc1( 6)
         jts( 7) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1( 7)
         jts( 8) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1( 8)

         jts( 9) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1( 9)
         jts(10) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1(10)
         jts(11) = (0D0,0D0)   ! g1(1)* g2 * sc1(1) * sc2(1) * tc1(11)
         jts(12) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1(12)
         
         jts(13) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1(13)
         jts(14) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1(14)
         jts(15) =  g1(1)* g2 * sc1(1) * sc2(1) * tc1(15)
         jts(16) = (0D0,0D0)   ! g1(1)* g2 * sc1(1) * sc2(1) * tc1(16)

      else
         jts( 1)=(0D0,0D0)
         jts( 2)=(0D0,0D0)
         jts( 3)=(0D0,0D0)
         jts( 4)=(0D0,0D0)
         jts( 5)=(0D0,0D0)
         jts( 6)=(0D0,0D0)
         jts( 7)=(0D0,0D0)
         jts( 8)=(0D0,0D0)
         jts( 9)=(0D0,0D0)
         jts(10)=(0D0,0D0)
         jts(11)=(0D0,0D0)
         jts(12)=(0D0,0D0)
         jts(13)=(0D0,0D0)
         jts(14)=(0D0,0D0)
         jts(15)=(0D0,0D0)
         jts(16)=(0D0,0D0)
      endif


      return
      end
      subroutine utsxxx(tc1,sc,gt,xm,xw,jts)
c
c- by RF - Feb. 2006 
c
c This subroutine computes an off-shell tensor current from the tts coupling.
c
c     input:
c          complex tc1(18)           : Incoming tensor particle
c          complex sc(3)             : Incoming scalar particle
c          complex gt(2)             : coupling constant: gt(1) scalar
c                                                         gt(2) not used
c
c     output:
c          complex jts               : off-shell tensor current
c
c     not used:
c          xm, xw
c

      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex tc1(DIM),jts(DIM),sc(DIM), gt(2)
      double precision xm, xw

c The outgoing tensor current is the same as the incoming multiplied by the
c coupling constant and the scalar particle. Note that the diagonal tensor
c terms are always zero because the tensor particle is anti-symmetric. The 
c tensor particle does not propagate, thus no multiplication with the tensor 
c propagator.

      jts(17) = sc(2) + tc1(17)
      jts(18) = sc(3) + tc1(18)

 
      if (gt(1).NE.(0D0,0D0)) then

         jts( 1) = (0D0,0D0)        !-gt(1) * sc(1) * tc1( 1)
         jts( 2) =  -gt(1) * sc(1) * tc1( 2)
         jts( 3) =  -gt(1) * sc(1) * tc1( 3)
         jts( 4) =  -gt(1) * sc(1) * tc1( 4)

         jts( 5) =  -gt(1) * sc(1) * tc1( 5)
         jts( 6) = (0D0,0D0)        !-gt(1) * sc(1) * tc1( 6)
         jts( 7) =  -gt(1) * sc(1) * tc1( 7)
         jts( 8) =  -gt(1) * sc(1) * tc1( 8)

         jts( 9) =  -gt(1) * sc(1) * tc1( 9)
         jts(10) =  -gt(1) * sc(1) * tc1(10)
         jts(11) = (0D0,0D0)        !-gt(1) * sc(1) * tc1(11)
         jts(12) =  -gt(1) * sc(1) * tc1(12)
         
         jts(13) =  -gt(1) * sc(1) * tc1(13)
         jts(14) =  -gt(1) * sc(1) * tc1(14)
         jts(15) =  -gt(1) * sc(1) * tc1(15)
         jts(16) = (0D0,0D0)        !-gt(1) * sc(1) * tc1(16)

      else
         jts( 1)=(0D0,0D0)
         jts( 2)=(0D0,0D0)
         jts( 3)=(0D0,0D0)
         jts( 4)=(0D0,0D0)
         jts( 5)=(0D0,0D0)
         jts( 6)=(0D0,0D0)
         jts( 7)=(0D0,0D0)
         jts( 8)=(0D0,0D0)
         jts( 9)=(0D0,0D0)
         jts(10)=(0D0,0D0)
         jts(11)=(0D0,0D0)
         jts(12)=(0D0,0D0)
         jts(13)=(0D0,0D0)
         jts(14)=(0D0,0D0)
         jts(15)=(0D0,0D0)
         jts(16)=(0D0,0D0)
      endif


      return
      end
      subroutine uttaxx(wt,wt1,gt,wt2)
c-------------------CP3  2009.10-----------------
c
c This subroutine computes an off-shell non-propagating tensor current from 
c the coupling of TTT
c
c input:
c       complex wt(18)           : input tensor                                                       t
c       complex wt1(18)         : input non-propagating tensor                          T
c       complex gt             : coupling constant         gt=-1/Lambda
c output:
c       complex wt2(18)         : non-propagating tensor current  

      implicit none
      double complex wt(18), wt1(18), gt, wt2(18)

      double complex T(6,4),T1(6,4),T2(6,4)
      double precision MT(4,4)
      integer i, j,s,d
      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      T(1,1) = wt(1)
      T(1,2) = wt(2)
      T(1,3) = wt(3)
      T(1,4) = wt(4)
      T(2,1) = wt(5)
      T(2,2) = wt(6)
      T(2,3) = wt(7)
      T(2,4) = wt(8)
      T(3,1) = wt(9)
      T(3,2) = wt(10)
      T(3,3) = wt(11)
      T(3,4) = wt(12)
      T(4,1) = wt(13)
      T(4,2) = wt(14)
      T(4,3) = wt(15)
      T(4,4) = wt(16)
      T(5,1) = wt(17)
      T(6,1) = wt(18)
 
      T1(1,1) = wt1(1)
      T1(1,2) = wt1(2)
      T1(1,3) = wt1(3)
      T1(1,4) = wt1(4)
      T1(2,1) = wt1(5)
      T1(2,2) = wt1(6)
      T1(2,3) = wt1(7)
      T1(2,4) = wt1(8)
      T1(3,1) = wt1(9)
      T1(3,2) = wt1(10)
      T1(3,3) = wt1(11)
      T1(3,4) = wt1(12)
      T1(4,1) = wt1(13)
      T1(4,2) = wt1(14)
      T1(4,3) = wt1(15)
      T1(4,4) = wt1(16)
      T1(5,1) = wt1(17)
      T1(6,1) = wt1(18)

      T2(5,1) = T(5,1)+T1(5,1)
      T2(6,1) = T(6,1)+T1(6,1)


      do i=1,4
         do j=1,4
            MT(i,j) = 0.0d0
         enddo 
      enddo
      MT(1,1) =  1.0d0
      MT(2,2) = -1.0d0
      MT(3,3) = -1.0d0
      MT(4,4) = -1.0d0
      
    
	
         do s=1,4
           do d=1,4
          T2(s,d)=MT(1,s)*MT(2,d)*T(1,1)*T1(1,2) - 
     -  MT(1,d)*MT(2,s)*T(1,1)*T1(1,2) + 
     -  MT(2,s)*MT(3,d)*T(1,3)*T1(1,2) - 
     -  MT(2,d)*MT(3,s)*T(1,3)*T1(1,2) + 
     -  MT(2,s)*MT(4,d)*T(1,4)*T1(1,2) - 
     -  MT(2,d)*MT(4,s)*T(1,4)*T1(1,2) - 
     -  MT(1,s)*MT(2,d)*T(2,2)*T1(1,2) + 
     -  MT(1,d)*MT(2,s)*T(2,2)*T1(1,2) - 
     -  MT(1,s)*MT(3,d)*T(2,3)*T1(1,2) + 
     -  MT(1,d)*MT(3,s)*T(2,3)*T1(1,2) - 
     -  MT(1,s)*MT(4,d)*T(2,4)*T1(1,2) + 
     -  MT(1,d)*MT(4,s)*T(2,4)*T1(1,2) + 
     -  MT(2,s)*MT(3,d)*T(3,1)*T1(1,2) - 
     -  MT(2,d)*MT(3,s)*T(3,1)*T1(1,2) - 
     -  MT(1,s)*MT(3,d)*T(3,2)*T1(1,2) + 
     -  MT(1,d)*MT(3,s)*T(3,2)*T1(1,2) + 
     -  MT(1,s)*MT(2,d)*T(3,3)*T1(1,2) - 
     -  MT(1,d)*MT(2,s)*T(3,3)*T1(1,2) + 
     -  MT(2,s)*MT(4,d)*T(4,1)*T1(1,2) - 
     -  MT(2,d)*MT(4,s)*T(4,1)*T1(1,2) - 
     -  MT(1,s)*MT(4,d)*T(4,2)*T1(1,2) + 
     -  MT(1,d)*MT(4,s)*T(4,2)*T1(1,2) + 
     -  MT(1,s)*MT(2,d)*T(4,4)*T1(1,2) - 
     -  MT(1,d)*MT(2,s)*T(4,4)*T1(1,2) + 
     -  MT(1,s)*MT(3,d)*T(1,1)*T1(1,3) - 
     -  MT(1,d)*MT(3,s)*T(1,1)*T1(1,3) - 
     -  MT(2,s)*MT(3,d)*T(1,2)*T1(1,3) + 
     -  MT(2,d)*MT(3,s)*T(1,2)*T1(1,3) + 
     -  MT(3,s)*MT(4,d)*T(1,4)*T1(1,3) - 
     -  MT(3,d)*MT(4,s)*T(1,4)*T1(1,3) - 
     -  MT(2,s)*MT(3,d)*T(2,1)*T1(1,3) + 
     -  MT(2,d)*MT(3,s)*T(2,1)*T1(1,3) + 
     -  MT(1,s)*MT(3,d)*T(2,2)*T1(1,3) - 
     -  MT(1,d)*MT(3,s)*T(2,2)*T1(1,3) - 
     -  MT(1,s)*MT(2,d)*T(2,3)*T1(1,3) + 
     -  MT(1,d)*MT(2,s)*T(2,3)*T1(1,3) - 
     -  MT(1,s)*MT(2,d)*T(3,2)*T1(1,3) + 
     -  MT(1,d)*MT(2,s)*T(3,2)*T1(1,3) - 
     -  MT(1,s)*MT(3,d)*T(3,3)*T1(1,3) + 
     -  MT(1,d)*MT(3,s)*T(3,3)*T1(1,3) - 
     -  MT(1,s)*MT(4,d)*T(3,4)*T1(1,3) + 
     -  MT(1,d)*MT(4,s)*T(3,4)*T1(1,3) + 
     -  MT(3,s)*MT(4,d)*T(4,1)*T1(1,3) - 
     -  MT(3,d)*MT(4,s)*T(4,1)*T1(1,3) - 
     -  MT(1,s)*MT(4,d)*T(4,3)*T1(1,3) + 
     -  MT(1,d)*MT(4,s)*T(4,3)*T1(1,3) + 
     -  MT(1,s)*MT(3,d)*T(4,4)*T1(1,3) - 
     -  MT(1,d)*MT(3,s)*T(4,4)*T1(1,3) + 
     -  MT(1,s)*MT(4,d)*T(1,1)*T1(1,4) - 
     -  MT(1,d)*MT(4,s)*T(1,1)*T1(1,4) - 
     -  MT(2,s)*MT(4,d)*T(1,2)*T1(1,4) + 
     -  MT(2,d)*MT(4,s)*T(1,2)*T1(1,4) - 
     -  MT(3,s)*MT(4,d)*T(1,3)*T1(1,4) + 
     -  MT(3,d)*MT(4,s)*T(1,3)*T1(1,4) - 
     -  MT(2,s)*MT(4,d)*T(2,1)*T1(1,4) + 
     -  MT(2,d)*MT(4,s)*T(2,1)*T1(1,4) + 
     -  MT(1,s)*MT(4,d)*T(2,2)*T1(1,4) - 
     -  MT(1,d)*MT(4,s)*T(2,2)*T1(1,4) - 
     -  MT(1,s)*MT(2,d)*T(2,4)*T1(1,4) + 
     -  MT(1,d)*MT(2,s)*T(2,4)*T1(1,4) - 
     -  MT(3,s)*MT(4,d)*T(3,1)*T1(1,4) + 
     -  MT(3,d)*MT(4,s)*T(3,1)*T1(1,4) + 
     -  MT(1,s)*MT(4,d)*T(3,3)*T1(1,4) - 
     -  MT(1,d)*MT(4,s)*T(3,3)*T1(1,4) - 
     -  MT(1,s)*MT(3,d)*T(3,4)*T1(1,4) + 
     -  MT(1,d)*MT(3,s)*T(3,4)*T1(1,4) - 
     -  MT(1,s)*MT(2,d)*T(4,2)*T1(1,4) + 
     -  MT(1,d)*MT(2,s)*T(4,2)*T1(1,4) - 
     -  MT(1,s)*MT(3,d)*T(4,3)*T1(1,4) + 
     -  MT(1,d)*MT(3,s)*T(4,3)*T1(1,4) - 
     -  MT(1,s)*MT(4,d)*T(4,4)*T1(1,4) + 
     -  MT(1,d)*MT(4,s)*T(4,4)*T1(1,4) - 
     -  MT(1,s)*MT(2,d)*T(1,1)*T1(2,1) + 
     -  MT(1,d)*MT(2,s)*T(1,1)*T1(2,1) - 
     -  MT(2,s)*MT(3,d)*T(1,3)*T1(2,1) + 
     -  MT(2,d)*MT(3,s)*T(1,3)*T1(2,1) - 
     -  MT(2,s)*MT(4,d)*T(1,4)*T1(2,1) + 
     -  MT(2,d)*MT(4,s)*T(1,4)*T1(2,1) + 
     -  MT(1,s)*MT(2,d)*T(2,2)*T1(2,1) - 
     -  MT(1,d)*MT(2,s)*T(2,2)*T1(2,1) + 
     -  MT(1,s)*MT(3,d)*T(2,3)*T1(2,1) - 
     -  MT(1,d)*MT(3,s)*T(2,3)*T1(2,1) + 
     -  MT(1,s)*MT(4,d)*T(2,4)*T1(2,1) - 
     -  MT(1,d)*MT(4,s)*T(2,4)*T1(2,1) - 
     -  MT(2,s)*MT(3,d)*T(3,1)*T1(2,1) + 
     -  MT(2,d)*MT(3,s)*T(3,1)*T1(2,1) + 
     -  MT(1,s)*MT(3,d)*T(3,2)*T1(2,1) - 
     -  MT(1,d)*MT(3,s)*T(3,2)*T1(2,1) - 
     -  MT(1,s)*MT(2,d)*T(3,3)*T1(2,1) + 
     -  MT(1,d)*MT(2,s)*T(3,3)*T1(2,1) - 
     -  MT(2,s)*MT(4,d)*T(4,1)*T1(2,1) + 
     -  MT(2,d)*MT(4,s)*T(4,1)*T1(2,1) + 
     -  MT(1,s)*MT(4,d)*T(4,2)*T1(2,1) - 
     -  MT(1,d)*MT(4,s)*T(4,2)*T1(2,1) - 
     -  MT(1,s)*MT(2,d)*T(4,4)*T1(2,1) + 
     -  MT(1,d)*MT(2,s)*T(4,4)*T1(2,1) + 
     -  MT(2,s)*MT(3,d)*T(1,1)*T1(2,3) - 
     -  MT(2,d)*MT(3,s)*T(1,1)*T1(2,3) - 
     -  MT(1,s)*MT(3,d)*T(1,2)*T1(2,3) + 
     -  MT(1,d)*MT(3,s)*T(1,2)*T1(2,3) + 
     -  MT(1,s)*MT(2,d)*T(1,3)*T1(2,3) - 
     -  MT(1,d)*MT(2,s)*T(1,3)*T1(2,3) - 
     -  MT(1,s)*MT(3,d)*T(2,1)*T1(2,3) + 
     -  MT(1,d)*MT(3,s)*T(2,1)*T1(2,3) + 
     -  MT(2,s)*MT(3,d)*T(2,2)*T1(2,3) - 
     -  MT(2,d)*MT(3,s)*T(2,2)*T1(2,3) - 
     -  MT(3,s)*MT(4,d)*T(2,4)*T1(2,3) + 
     -  MT(3,d)*MT(4,s)*T(2,4)*T1(2,3) + 
     -  MT(1,s)*MT(2,d)*T(3,1)*T1(2,3) - 
     -  MT(1,d)*MT(2,s)*T(3,1)*T1(2,3) + 
     -  MT(2,s)*MT(3,d)*T(3,3)*T1(2,3) - 
     -  MT(2,d)*MT(3,s)*T(3,3)*T1(2,3) + 
     -  MT(2,s)*MT(4,d)*T(3,4)*T1(2,3) - 
     -  MT(2,d)*MT(4,s)*T(3,4)*T1(2,3) - 
     -  MT(3,s)*MT(4,d)*T(4,2)*T1(2,3) + 
     -  MT(3,d)*MT(4,s)*T(4,2)*T1(2,3) + 
     -  MT(2,s)*MT(4,d)*T(4,3)*T1(2,3) - 
     -  MT(2,d)*MT(4,s)*T(4,3)*T1(2,3) - 
     -  MT(2,s)*MT(3,d)*T(4,4)*T1(2,3) + 
     -  MT(2,d)*MT(3,s)*T(4,4)*T1(2,3) + 
     -  MT(2,s)*MT(4,d)*T(1,1)*T1(2,4) - 
     -  MT(2,d)*MT(4,s)*T(1,1)*T1(2,4) - 
     -  MT(1,s)*MT(4,d)*T(1,2)*T1(2,4) + 
     -  MT(1,d)*MT(4,s)*T(1,2)*T1(2,4) + 
     -  MT(1,s)*MT(2,d)*T(1,4)*T1(2,4) - 
     -  MT(1,d)*MT(2,s)*T(1,4)*T1(2,4) - 
     -  MT(1,s)*MT(4,d)*T(2,1)*T1(2,4) + 
     -  MT(1,d)*MT(4,s)*T(2,1)*T1(2,4) + 
     -  MT(2,s)*MT(4,d)*T(2,2)*T1(2,4) - 
     -  MT(2,d)*MT(4,s)*T(2,2)*T1(2,4) + 
     -  MT(3,s)*MT(4,d)*T(2,3)*T1(2,4) - 
     -  MT(3,d)*MT(4,s)*T(2,3)*T1(2,4) + 
     -  MT(3,s)*MT(4,d)*T(3,2)*T1(2,4) - 
     -  MT(3,d)*MT(4,s)*T(3,2)*T1(2,4) - 
     -  MT(2,s)*MT(4,d)*T(3,3)*T1(2,4) + 
     -  MT(2,d)*MT(4,s)*T(3,3)*T1(2,4) + 
     -  MT(2,s)*MT(3,d)*T(3,4)*T1(2,4) - 
     -  MT(2,d)*MT(3,s)*T(3,4)*T1(2,4) + 
     -  MT(1,s)*MT(2,d)*T(4,1)*T1(2,4) - 
     -  MT(1,d)*MT(2,s)*T(4,1)*T1(2,4) + 
     -  MT(2,s)*MT(3,d)*T(4,3)*T1(2,4) - 
     -  MT(2,d)*MT(3,s)*T(4,3)*T1(2,4) + 
     -  MT(2,s)*MT(4,d)*T(4,4)*T1(2,4) - 
     -  MT(2,d)*MT(4,s)*T(4,4)*T1(2,4) - 
     -  MT(1,s)*MT(3,d)*T(1,1)*T1(3,1) + 
     -  MT(1,d)*MT(3,s)*T(1,1)*T1(3,1) + 
     -  MT(2,s)*MT(3,d)*T(1,2)*T1(3,1) - 
     -  MT(2,d)*MT(3,s)*T(1,2)*T1(3,1) - 
     -  MT(3,s)*MT(4,d)*T(1,4)*T1(3,1) + 
     -  MT(3,d)*MT(4,s)*T(1,4)*T1(3,1) + 
     -  MT(2,s)*MT(3,d)*T(2,1)*T1(3,1) - 
     -  MT(2,d)*MT(3,s)*T(2,1)*T1(3,1) - 
     -  MT(1,s)*MT(3,d)*T(2,2)*T1(3,1) + 
     -  MT(1,d)*MT(3,s)*T(2,2)*T1(3,1) + 
     -  MT(1,s)*MT(2,d)*T(2,3)*T1(3,1) - 
     -  MT(1,d)*MT(2,s)*T(2,3)*T1(3,1) + 
     -  MT(1,s)*MT(2,d)*T(3,2)*T1(3,1) - 
     -  MT(1,d)*MT(2,s)*T(3,2)*T1(3,1) + 
     -  MT(1,s)*MT(3,d)*T(3,3)*T1(3,1) - 
     -  MT(1,d)*MT(3,s)*T(3,3)*T1(3,1) + 
     -  MT(1,s)*MT(4,d)*T(3,4)*T1(3,1) - 
     -  MT(1,d)*MT(4,s)*T(3,4)*T1(3,1) - 
     -  MT(3,s)*MT(4,d)*T(4,1)*T1(3,1) + 
     -  MT(3,d)*MT(4,s)*T(4,1)*T1(3,1) + 
     -  MT(1,s)*MT(4,d)*T(4,3)*T1(3,1) - 
     -  MT(1,d)*MT(4,s)*T(4,3)*T1(3,1) - 
     -  MT(1,s)*MT(3,d)*T(4,4)*T1(3,1) + 
     -  MT(1,d)*MT(3,s)*T(4,4)*T1(3,1) - 
     -  MT(2,s)*MT(3,d)*T(1,1)*T1(3,2) + 
     -  MT(2,d)*MT(3,s)*T(1,1)*T1(3,2) + 
     -  MT(1,s)*MT(3,d)*T(1,2)*T1(3,2) - 
     -  MT(1,d)*MT(3,s)*T(1,2)*T1(3,2) - 
     -  MT(1,s)*MT(2,d)*T(1,3)*T1(3,2) + 
     -  MT(1,d)*MT(2,s)*T(1,3)*T1(3,2) + 
     -  MT(1,s)*MT(3,d)*T(2,1)*T1(3,2) - 
     -  MT(1,d)*MT(3,s)*T(2,1)*T1(3,2) - 
     -  MT(2,s)*MT(3,d)*T(2,2)*T1(3,2) + 
     -  MT(2,d)*MT(3,s)*T(2,2)*T1(3,2) + 
     -  MT(3,s)*MT(4,d)*T(2,4)*T1(3,2) - 
     -  MT(3,d)*MT(4,s)*T(2,4)*T1(3,2) - 
     -  MT(1,s)*MT(2,d)*T(3,1)*T1(3,2) + 
     -  MT(1,d)*MT(2,s)*T(3,1)*T1(3,2) - 
     -  MT(2,s)*MT(3,d)*T(3,3)*T1(3,2) + 
     -  MT(2,d)*MT(3,s)*T(3,3)*T1(3,2) - 
     -  MT(2,s)*MT(4,d)*T(3,4)*T1(3,2) + 
     -  MT(2,d)*MT(4,s)*T(3,4)*T1(3,2) + 
     -  MT(3,s)*MT(4,d)*T(4,2)*T1(3,2) - 
     -  MT(3,d)*MT(4,s)*T(4,2)*T1(3,2) - 
     -  MT(2,s)*MT(4,d)*T(4,3)*T1(3,2) + 
     -  MT(2,d)*MT(4,s)*T(4,3)*T1(3,2) + 
     -  MT(2,s)*MT(3,d)*T(4,4)*T1(3,2) - 
     -  MT(2,d)*MT(3,s)*T(4,4)*T1(3,2) + 
     -  MT(3,s)*MT(4,d)*T(1,1)*T1(3,4) - 
     -  MT(3,d)*MT(4,s)*T(1,1)*T1(3,4) - 
     -  MT(1,s)*MT(4,d)*T(1,3)*T1(3,4) + 
     -  MT(1,d)*MT(4,s)*T(1,3)*T1(3,4) + 
     -  MT(1,s)*MT(3,d)*T(1,4)*T1(3,4) - 
     -  MT(1,d)*MT(3,s)*T(1,4)*T1(3,4) - 
     -  MT(3,s)*MT(4,d)*T(2,2)*T1(3,4) + 
     -  MT(3,d)*MT(4,s)*T(2,2)*T1(3,4) + 
     -  MT(2,s)*MT(4,d)*T(2,3)*T1(3,4) - 
     -  MT(2,d)*MT(4,s)*T(2,3)*T1(3,4) - 
     -  MT(2,s)*MT(3,d)*T(2,4)*T1(3,4) + 
     -  MT(2,d)*MT(3,s)*T(2,4)*T1(3,4) - 
     -  MT(1,s)*MT(4,d)*T(3,1)*T1(3,4) + 
     -  MT(1,d)*MT(4,s)*T(3,1)*T1(3,4) + 
     -  MT(2,s)*MT(4,d)*T(3,2)*T1(3,4) - 
     -  MT(2,d)*MT(4,s)*T(3,2)*T1(3,4) + 
     -  MT(3,s)*MT(4,d)*T(3,3)*T1(3,4) - 
     -  MT(3,d)*MT(4,s)*T(3,3)*T1(3,4) + 
     -  MT(1,s)*MT(3,d)*T(4,1)*T1(3,4) - 
     -  MT(1,d)*MT(3,s)*T(4,1)*T1(3,4) - 
     -  MT(2,s)*MT(3,d)*T(4,2)*T1(3,4) + 
     -  MT(2,d)*MT(3,s)*T(4,2)*T1(3,4) + 
     -  MT(3,s)*MT(4,d)*T(4,4)*T1(3,4) - 
     -  MT(3,d)*MT(4,s)*T(4,4)*T1(3,4) - 
     -  MT(1,s)*MT(4,d)*T(1,1)*T1(4,1) + 
     -  MT(1,d)*MT(4,s)*T(1,1)*T1(4,1) + 
     -  MT(2,s)*MT(4,d)*T(1,2)*T1(4,1) - 
     -  MT(2,d)*MT(4,s)*T(1,2)*T1(4,1) + 
     -  MT(3,s)*MT(4,d)*T(1,3)*T1(4,1) - 
     -  MT(3,d)*MT(4,s)*T(1,3)*T1(4,1) + 
     -  MT(2,s)*MT(4,d)*T(2,1)*T1(4,1) - 
     -  MT(2,d)*MT(4,s)*T(2,1)*T1(4,1) - 
     -  MT(1,s)*MT(4,d)*T(2,2)*T1(4,1) + 
     -  MT(1,d)*MT(4,s)*T(2,2)*T1(4,1) + 
     -  MT(1,s)*MT(2,d)*T(2,4)*T1(4,1) - 
     -  MT(1,d)*MT(2,s)*T(2,4)*T1(4,1) + 
     -  MT(3,s)*MT(4,d)*T(3,1)*T1(4,1) - 
     -  MT(3,d)*MT(4,s)*T(3,1)*T1(4,1) - 
     -  MT(1,s)*MT(4,d)*T(3,3)*T1(4,1) + 
     -  MT(1,d)*MT(4,s)*T(3,3)*T1(4,1) + 
     -  MT(1,s)*MT(3,d)*T(3,4)*T1(4,1) - 
     -  MT(1,d)*MT(3,s)*T(3,4)*T1(4,1) + 
     -  MT(1,s)*MT(2,d)*T(4,2)*T1(4,1) - 
     -  MT(1,d)*MT(2,s)*T(4,2)*T1(4,1) + 
     -  MT(1,s)*MT(3,d)*T(4,3)*T1(4,1) - 
     -  MT(1,d)*MT(3,s)*T(4,3)*T1(4,1) + 
     -  MT(1,s)*MT(4,d)*T(4,4)*T1(4,1) - 
     -  MT(1,d)*MT(4,s)*T(4,4)*T1(4,1) - 
     -  MT(2,s)*MT(4,d)*T(1,1)*T1(4,2) + 
     -  MT(2,d)*MT(4,s)*T(1,1)*T1(4,2) + 
     -  MT(1,s)*MT(4,d)*T(1,2)*T1(4,2) - 
     -  MT(1,d)*MT(4,s)*T(1,2)*T1(4,2) - 
     -  MT(1,s)*MT(2,d)*T(1,4)*T1(4,2) + 
     -  MT(1,d)*MT(2,s)*T(1,4)*T1(4,2) + 
     -  MT(1,s)*MT(4,d)*T(2,1)*T1(4,2) - 
     -  MT(1,d)*MT(4,s)*T(2,1)*T1(4,2) - 
     -  MT(2,s)*MT(4,d)*T(2,2)*T1(4,2) + 
     -  MT(2,d)*MT(4,s)*T(2,2)*T1(4,2) - 
     -  MT(3,s)*MT(4,d)*T(2,3)*T1(4,2) + 
     -  MT(3,d)*MT(4,s)*T(2,3)*T1(4,2) - 
     -  MT(3,s)*MT(4,d)*T(3,2)*T1(4,2) + 
     -  MT(3,d)*MT(4,s)*T(3,2)*T1(4,2) + 
     -  MT(2,s)*MT(4,d)*T(3,3)*T1(4,2) - 
     -  MT(2,d)*MT(4,s)*T(3,3)*T1(4,2) - 
     -  MT(2,s)*MT(3,d)*T(3,4)*T1(4,2) + 
     -  MT(2,d)*MT(3,s)*T(3,4)*T1(4,2) - 
     -  MT(1,s)*MT(2,d)*T(4,1)*T1(4,2) + 
     -  MT(1,d)*MT(2,s)*T(4,1)*T1(4,2) - 
     -  MT(2,s)*MT(3,d)*T(4,3)*T1(4,2) + 
     -  MT(2,d)*MT(3,s)*T(4,3)*T1(4,2) - 
     -  MT(2,s)*MT(4,d)*T(4,4)*T1(4,2) + 
     -  MT(2,d)*MT(4,s)*T(4,4)*T1(4,2) - 
     -  MT(3,s)*MT(4,d)*T(1,1)*T1(4,3) + 
     -  MT(3,d)*MT(4,s)*T(1,1)*T1(4,3) + 
     -  MT(1,s)*MT(4,d)*T(1,3)*T1(4,3) - 
     -  MT(1,d)*MT(4,s)*T(1,3)*T1(4,3) - 
     -  MT(1,s)*MT(3,d)*T(1,4)*T1(4,3) + 
     -  MT(1,d)*MT(3,s)*T(1,4)*T1(4,3) + 
     -  MT(3,s)*MT(4,d)*T(2,2)*T1(4,3) - 
     -  MT(3,d)*MT(4,s)*T(2,2)*T1(4,3) - 
     -  MT(2,s)*MT(4,d)*T(2,3)*T1(4,3) + 
     -  MT(2,d)*MT(4,s)*T(2,3)*T1(4,3) + 
     -  MT(2,s)*MT(3,d)*T(2,4)*T1(4,3) - 
     -  MT(2,d)*MT(3,s)*T(2,4)*T1(4,3) + 
     -  MT(1,s)*MT(4,d)*T(3,1)*T1(4,3) - 
     -  MT(1,d)*MT(4,s)*T(3,1)*T1(4,3) - 
     -  MT(2,s)*MT(4,d)*T(3,2)*T1(4,3) + 
     -  MT(2,d)*MT(4,s)*T(3,2)*T1(4,3) - 
     -  MT(3,s)*MT(4,d)*T(3,3)*T1(4,3) + 
     -  MT(3,d)*MT(4,s)*T(3,3)*T1(4,3) - 
     -  MT(1,s)*MT(3,d)*T(4,1)*T1(4,3) + 
     -  MT(1,d)*MT(3,s)*T(4,1)*T1(4,3) + 
     -  MT(2,s)*MT(3,d)*T(4,2)*T1(4,3) - 
     -  MT(2,d)*MT(3,s)*T(4,2)*T1(4,3) - 
     -  MT(3,s)*MT(4,d)*T(4,4)*T1(4,3) + 
     -  MT(3,d)*MT(4,s)*T(4,4)*T1(4,3)
           enddo
         enddo
       

         do s=1,4
           do d=1,4
          T2(s,d)=-T2(s,d)*gt
           enddo
         enddo
         
      wt2(1) = T2(1,1)
      wt2(2) = T2(1,2)
      wt2(3) = T2(1,3)
      wt2(4) = T2(1,4)
      wt2(5) = T2(2,1)
      wt2(6) = T2(2,2)
      wt2(7) = T2(2,3)
      wt2(8) = T2(2,4)
      wt2(9) = T2(3,1)
      wt2(10) = T2(3,2)
      wt2(11) = T2(3,3)
      wt2(12) = T2(3,4)
      wt2(13) = T2(4,1)
      wt2(14) = T2(4,2)
      wt2(15) = T2(4,3)
      wt2(16) = T2(4,4)
      wt2(17) = T2(5,1)
      wt2(18) = T2(6,1)

      return
      end
      subroutine uttbxx(wt1,wt2,gt,tmass,twidth,wt)
c-------------------CP3  2009.10-----------------
c
c This subroutine computes an off-shell tensor current from 
c the coupling of two non-propagating tensor bosons.
c
c input:
c       complex wt1(18)           : input non-propagating tensor                                                       t
c       complex wt2(18)         :  input non-propagating tensor                      T
c       complex gt             : coupling constant         gt=-1/Lambda
c output:
c       complex wt(18)         :  tensor 

      implicit none
      double complex wt(18), wt1(18), gt, wt2(18)
      double precision pT(4),pT2,tmass,twidth
      double complex T(6,4),T1(6,4),T2(6,4),d
      double precision MT(4,4)
      integer i, j,m,n
      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      T2(1,1) = wt2(1)
      T2(1,2) = wt2(2)
      T2(1,3) = wt2(3)
      T2(1,4) = wt2(4)
      T2(2,1) = wt2(5)
      T2(2,2) = wt2(6)
      T2(2,3) = wt2(7)
      T2(2,4) = wt2(8)
      T2(3,1) = wt2(9)
      T2(3,2) = wt2(10)
      T2(3,3) = wt2(11)
      T2(3,4) = wt2(12)
      T2(4,1) = wt2(13)
      T2(4,2) = wt2(14)
      T2(4,3) = wt2(15)
      T2(4,4) = wt2(16)
      T2(5,1) = wt2(17)
      T2(6,1) = wt2(18)
 
      T1(1,1) = wt1(1)
      T1(1,2) = wt1(2)
      T1(1,3) = wt1(3)
      T1(1,4) = wt1(4)
      T1(2,1) = wt1(5)
      T1(2,2) = wt1(6)
      T1(2,3) = wt1(7)
      T1(2,4) = wt1(8)
      T1(3,1) = wt1(9)
      T1(3,2) = wt1(10)
      T1(3,3) = wt1(11)
      T1(3,4) = wt1(12)
      T1(4,1) = wt1(13)
      T1(4,2) = wt1(14)
      T1(4,3) = wt1(15)
      T1(4,4) = wt1(16)
      T1(5,1) = wt1(17)
      T1(6,1) = wt1(18)

      T(5,1) = T2(5,1)+T1(5,1)
      T(6,1) = T2(6,1)+T1(6,1)


      do i=1,4
         do j=1,4
            MT(i,j) = 0.0d0
         enddo 
      enddo
      MT(1,1) =  1.0d0
      MT(2,2) = -1.0d0
      MT(3,3) = -1.0d0
      MT(4,4) = -1.0d0
      
    
	
         do m=1,4
           do n=1,4

          T(m,n)= -((2*MT(1,m)*MT(1,n)
     -   - (2*MT(m,n))/3.)*T1(1,2)*T2(1,2)) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(1,2)*T2(1,2) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(1,2)*T2(1,2) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(1,2)*T2(1,2) + 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(1,3)*T2(1,2) + 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(1,4)*T2(1,2) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(2,1)*T2(1,2) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(2,1)*T2(1,2) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(2,1)*T2(1,2) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(2,1)*T2(1,2) - 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(2,3)*T2(1,2) - 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(2,4)*T2(1,2) - 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(3,1)*T2(1,2) + 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(3,2)*T2(1,2) - 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(4,1)*T2(1,2) + 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(4,2)*T2(1,2) + 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(1,2)*T2(1,3) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(1,3)*T2(1,3) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(1,3)*T2(1,3) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(1,3)*T2(1,3) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(1,3)*T2(1,3) + 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(1,4)*T2(1,3) - 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(2,1)*T2(1,3) + 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(2,3)*T2(1,3) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(3,1)*T2(1,3) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(3,1)*T2(1,3) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(3,1)*T2(1,3) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(3,1)*T2(1,3) - 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(3,2)*T2(1,3) - 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(3,4)*T2(1,3) - 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(4,1)*T2(1,3) + 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(4,3)*T2(1,3) + 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(1,2)*T2(1,4) + 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(1,3)*T2(1,4) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(1,4)*T2(1,4) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(1,4)*T2(1,4) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(1,4)*T2(1,4) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(1,4)*T2(1,4) - 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(2,1)*T2(1,4) + 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(2,4)*T2(1,4) - 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(3,1)*T2(1,4) + 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(3,4)*T2(1,4) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(4,1)*T2(1,4) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(4,1)*T2(1,4) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(4,1)*T2(1,4) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(4,1)*T2(1,4) - 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(4,2)*T2(1,4) - 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(4,3)*T2(1,4) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(1,2)*T2(2,1) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(1,2)*T2(2,1) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(1,2)*T2(2,1) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(1,2)*T2(2,1) - 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(1,3)*T2(2,1) - 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(1,4)*T2(2,1) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(2,1)*T2(2,1) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(2,1)*T2(2,1) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(2,1)*T2(2,1) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(2,1)*T2(2,1) + 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(2,3)*T2(2,1) + 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(2,4)*T2(2,1) + 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(3,1)*T2(2,1) - 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(3,2)*T2(2,1) + 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(4,1)*T2(2,1) - 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(4,2)*T2(2,1) - 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(1,2)*T2(2,3) + 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(1,3)*T2(2,3) + 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(2,1)*T2(2,3) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(2,3)*T2(2,3) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(2,3)*T2(2,3) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(2,3)*T2(2,3) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(2,3)*T2(2,3) - 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(2,4)*T2(2,3) - 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(3,1)*T2(2,3) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(3,2)*T2(2,3) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(3,2)*T2(2,3) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(3,2)*T2(2,3) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(3,2)*T2(2,3) + 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(3,4)*T2(2,3) + 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(4,2)*T2(2,3) - 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(4,3)*T2(2,3) - 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(1,2)*T2(2,4) + 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(1,4)*T2(2,4) + 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(2,1)*T2(2,4) - 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(2,3)*T2(2,4) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(2,4)*T2(2,4) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(2,4)*T2(2,4) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(2,4)*T2(2,4) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(2,4)*T2(2,4) + 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(3,2)*T2(2,4) - 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(3,4)*T2(2,4) - 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(4,1)*T2(2,4) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(4,2)*T2(2,4) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(4,2)*T2(2,4) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(4,2)*T2(2,4) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(4,2)*T2(2,4) + 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(4,3)*T2(2,4) - 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(1,2)*T2(3,1) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(1,3)*T2(3,1) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(1,3)*T2(3,1) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(1,3)*T2(3,1) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(1,3)*T2(3,1) - 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(1,4)*T2(3,1) + 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(2,1)*T2(3,1) - 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(2,3)*T2(3,1) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(3,1)*T2(3,1) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(3,1)*T2(3,1) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(3,1)*T2(3,1) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(3,1)*T2(3,1) + 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(3,2)*T2(3,1) + 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(3,4)*T2(3,1) + 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(4,1)*T2(3,1) - 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(4,3)*T2(3,1) + 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(1,2)*T2(3,2) - 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(1,3)*T2(3,2) - 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(2,1)*T2(3,2) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(2,3)*T2(3,2) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(2,3)*T2(3,2) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(2,3)*T2(3,2) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(2,3)*T2(3,2) + 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(2,4)*T2(3,2) + 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(3,1)*T2(3,2) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(3,2)*T2(3,2) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(3,2)*T2(3,2) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(3,2)*T2(3,2) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(3,2)*T2(3,2) - 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(3,4)*T2(3,2) - 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(4,2)*T2(3,2) + 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(4,3)*T2(3,2) - 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(1,3)*T2(3,4) + 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(1,4)*T2(3,4) + 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(2,3)*T2(3,4) - 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(2,4)*T2(3,4) + 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(3,1)*T2(3,4) - 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(3,2)*T2(3,4) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(3,4)*T2(3,4) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(3,4)*T2(3,4) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(3,4)*T2(3,4) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(3,4)*T2(3,4) - 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(4,1)*T2(3,4) + 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(4,2)*T2(3,4) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(4,3)*T2(3,4) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(4,3)*T2(3,4) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(4,3)*T2(3,4) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(4,3)*T2(3,4) - 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(1,2)*T2(4,1) - 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(1,3)*T2(4,1) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(1,4)*T2(4,1) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(1,4)*T2(4,1) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(1,4)*T2(4,1) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(1,4)*T2(4,1) + 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(2,1)*T2(4,1) - 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(2,4)*T2(4,1) + 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(3,1)*T2(4,1) - 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(3,4)*T2(4,1) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(4,1)*T2(4,1) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(4,1)*T2(4,1) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(4,1)*T2(4,1) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(4,1)*T2(4,1) + 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(4,2)*T2(4,1) + 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(4,3)*T2(4,1) + 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(1,2)*T2(4,2) - 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(1,4)*T2(4,2) - 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(2,1)*T2(4,2) + 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(2,3)*T2(4,2) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(2,4)*T2(4,2) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(2,4)*T2(4,2) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(2,4)*T2(4,2) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(2,4)*T2(4,2) - 
     -  2*(MT(3,n)*MT(4,m) + MT(3,m)*MT(4,n))*T1(3,2)*T2(4,2) + 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(3,4)*T2(4,2) + 
     -  2*(MT(1,n)*MT(2,m) + MT(1,m)*MT(2,n))*T1(4,1)*T2(4,2) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(4,2)*T2(4,2) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(4,2)*T2(4,2) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(4,2)*T2(4,2) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(4,2)*T2(4,2) - 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(4,3)*T2(4,2) + 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(1,3)*T2(4,3) - 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(1,4)*T2(4,3) - 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(2,3)*T2(4,3) + 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(2,4)*T2(4,3) - 
     -  2*(MT(1,n)*MT(4,m) + MT(1,m)*MT(4,n))*T1(3,1)*T2(4,3) + 
     -  2*(MT(2,n)*MT(4,m) + MT(2,m)*MT(4,n))*T1(3,2)*T2(4,3) + 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(3,4)*T2(4,3) - 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(3,4)*T2(4,3) + 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(3,4)*T2(4,3) + 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(3,4)*T2(4,3) + 
     -  2*(MT(1,n)*MT(3,m) + MT(1,m)*MT(3,n))*T1(4,1)*T2(4,3) - 
     -  2*(MT(2,n)*MT(3,m) + MT(2,m)*MT(3,n))*T1(4,2)*T2(4,3) - 
     -  (2*MT(1,m)*MT(1,n) - (2*MT(m,n))/3.)*T1(4,3)*T2(4,3) + 
     -  (2*MT(2,m)*MT(2,n) + (2*MT(m,n))/3.)*T1(4,3)*T2(4,3) - 
     -  (2*MT(3,m)*MT(3,n) + (2*MT(m,n))/3.)*T1(4,3)*T2(4,3) - 
     -  (2*MT(4,m)*MT(4,n) + (2*MT(m,n))/3.)*T1(4,3)*T2(4,3)
           enddo
         enddo
       
       pT(1) = dreal(T(5,1))
       pT(2) = dreal(T(6,1))
       pT(3) = dimag(T(6,1))
       pT(4) = dimag(T(5,1))
       pT2 = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2


         do m=1,4
           do n=1,4
          T(m,n)=T(m,n)
     - +((-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(1,2)*
     -   T2(1,2) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(1,2)*T2(1,2)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(1,2)*T2(1,2)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(1,2)*T2(1,2)
     -   + 2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(1,3)*T2(1,2) + 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(1,4)*T2(1,2) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(2,1)*
     -   T2(1,2) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(2,1)*T2(1,2)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(2,1)*T2(1,2)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(2,1)*T2(1,2)
     -   - 2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(2,3)*T2(1,2) - 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(2,4)*T2(1,2) - 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(3,1)*T2(1,2) + 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(3,2)*T2(1,2) - 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(4,1)*T2(1,2) + 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(4,2)*T2(1,2) + 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(1,2)*T2(1,3) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(1,3)*
     -   T2(1,3) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(1,3)*T2(1,3)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(1,3)*T2(1,3)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(1,3)*T2(1,3)
     -   + 2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(1,4)*T2(1,3) - 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(2,1)*T2(1,3) + 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(2,3)*T2(1,3) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(3,1)*
     -   T2(1,3) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(3,1)*T2(1,3)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(3,1)*T2(1,3)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(3,1)*T2(1,3)
     -   - 2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(3,2)*T2(1,3) - 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(3,4)*T2(1,3) - 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(4,1)*T2(1,3) + 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(4,3)*T2(1,3) + 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(1,2)*T2(1,4) + 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(1,3)*T2(1,4) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(1,4)*
     -   T2(1,4) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(1,4)*T2(1,4)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(1,4)*T2(1,4)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(1,4)*T2(1,4)
     -   - 2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(2,1)*T2(1,4) + 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(2,4)*T2(1,4) - 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(3,1)*T2(1,4) + 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(3,4)*T2(1,4) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(4,1)*
     -   T2(1,4) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(4,1)*T2(1,4)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(4,1)*T2(1,4)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(4,1)*T2(1,4)
     -   - 2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(4,2)*T2(1,4) - 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(4,3)*T2(1,4) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(1,2)*
     -   T2(2,1) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(1,2)*T2(2,1)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(1,2)*T2(2,1)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(1,2)*T2(2,1)
     -   - 2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(1,3)*T2(2,1) - 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(1,4)*T2(2,1) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(2,1)*
     -   T2(2,1) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(2,1)*T2(2,1)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(2,1)*T2(2,1)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(2,1)*T2(2,1)
     -   + 2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(2,3)*T2(2,1) + 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(2,4)*T2(2,1) + 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(3,1)*T2(2,1) - 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(3,2)*T2(2,1) + 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(4,1)*T2(2,1) - 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(4,2)*T2(2,1) - 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(1,2)*T2(2,3) + 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(1,3)*T2(2,3) + 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(2,1)*T2(2,3) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(2,3)*
     -   T2(2,3) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(2,3)*T2(2,3)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(2,3)*T2(2,3)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(2,3)*T2(2,3)
     -   - 2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(2,4)*T2(2,3) - 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(3,1)*T2(2,3) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(3,2)*
     -   T2(2,3) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(3,2)*T2(2,3)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(3,2)*T2(2,3)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(3,2)*T2(2,3)
     -   + 2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(3,4)*T2(2,3) + 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(4,2)*T2(2,3) - 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(4,3)*T2(2,3) - 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(1,2)*T2(2,4) + 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(1,4)*T2(2,4) + 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(2,1)*T2(2,4) - 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(2,3)*T2(2,4) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(2,4)*
     -   T2(2,4) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(2,4)*T2(2,4)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(2,4)*T2(2,4)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(2,4)*T2(2,4)
     -   + 2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(3,2)*T2(2,4) - 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(3,4)*T2(2,4) - 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(4,1)*T2(2,4) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(4,2)*
     -   T2(2,4) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(4,2)*T2(2,4)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(4,2)*T2(2,4)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(4,2)*T2(2,4)
     -   + 2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(4,3)*T2(2,4) - 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(1,2)*T2(3,1) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(1,3)*
     -   T2(3,1) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(1,3)*T2(3,1)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(1,3)*T2(3,1)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(1,3)*T2(3,1)
     -   - 2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(1,4)*T2(3,1) + 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(2,1)*T2(3,1) - 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(2,3)*T2(3,1) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(3,1)*
     -   T2(3,1) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(3,1)*T2(3,1)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(3,1)*T2(3,1)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(3,1)*T2(3,1)
     -   + 2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(3,2)*T2(3,1) + 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(3,4)*T2(3,1) + 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(4,1)*T2(3,1) - 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(4,3)*T2(3,1) + 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(1,2)*T2(3,2) - 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(1,3)*T2(3,2) - 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(2,1)*T2(3,2) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(2,3)*
     -   T2(3,2) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(2,3)*T2(3,2)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(2,3)*T2(3,2)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(2,3)*T2(3,2)
     -   + 2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(2,4)*T2(3,2) + 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(3,1)*T2(3,2) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(3,2)*
     -   T2(3,2) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(3,2)*T2(3,2)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(3,2)*T2(3,2)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(3,2)*T2(3,2)
     -   - 2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(3,4)*T2(3,2) - 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(4,2)*T2(3,2) + 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(4,3)*T2(3,2) - 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(1,3)*T2(3,4) + 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(1,4)*T2(3,4) + 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(2,3)*T2(3,4) - 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(2,4)*T2(3,4) + 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(3,1)*T2(3,4) - 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(3,2)*T2(3,4) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(3,4)*
     -   T2(3,4) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(3,4)*T2(3,4)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(3,4)*T2(3,4)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(3,4)*T2(3,4)
     -   - 2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(4,1)*T2(3,4) + 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(4,2)*T2(3,4) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(4,3)*
     -   T2(3,4) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(4,3)*T2(3,4)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(4,3)*T2(3,4)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(4,3)*T2(3,4)
     -   - 2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(1,2)*T2(4,1) - 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(1,3)*T2(4,1) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(1,4)*
     -   T2(4,1) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(1,4)*T2(4,1)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(1,4)*T2(4,1)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(1,4)*T2(4,1)
     -   + 2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(2,1)*T2(4,1) - 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(2,4)*T2(4,1) + 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(3,1)*T2(4,1) - 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(3,4)*T2(4,1) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(4,1)*
     -   T2(4,1) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(4,1)*T2(4,1)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(4,1)*T2(4,1)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(4,1)*T2(4,1)
     -   + 2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(4,2)*T2(4,1) + 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(4,3)*T2(4,1) + 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(1,2)*T2(4,2) - 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(1,4)*T2(4,2) - 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(2,1)*T2(4,2) + 
     -  2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(2,3)*T2(4,2) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(2,4)*
     -   T2(4,2) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(2,4)*T2(4,2)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(2,4)*T2(4,2)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(2,4)*T2(4,2)
     -   - 2*(-(pT(4)*pT(n)*MT(3,m)) - pT(4)*pT(m)*MT(3,n) - 
     -     pT(3)*pT(n)*MT(4,m) - pT(3)*pT(m)*MT(4,n) + 
     -     (2*pT(3)*pT(4)*MT(m,n))/3.)*T1(3,2)*T2(4,2) + 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(3,4)*T2(4,2) + 
     -  2*(-(pT(2)*pT(n)*MT(1,m)) - pT(2)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(2,m) - pT(1)*pT(m)*MT(2,n) + 
     -     (2*pT(1)*pT(2)*MT(m,n))/3.)*T1(4,1)*T2(4,2) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(4,2)*
     -   T2(4,2) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(4,2)*T2(4,2)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(4,2)*T2(4,2)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(4,2)*T2(4,2)
     -   - 2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(4,3)*T2(4,2) + 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(1,3)*T2(4,3) - 
     -  2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(1,4)*T2(4,3) - 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(2,3)*T2(4,3) + 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(2,4)*T2(4,3) - 
     -  2*(-(pT(4)*pT(n)*MT(1,m)) - pT(4)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(4,m) - pT(1)*pT(m)*MT(4,n) + 
     -     (2*pT(1)*pT(4)*MT(m,n))/3.)*T1(3,1)*T2(4,3) + 
     -  2*(-(pT(4)*pT(n)*MT(2,m)) - pT(4)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(4,m) - pT(2)*pT(m)*MT(4,n) + 
     -     (2*pT(2)*pT(4)*MT(m,n))/3.)*T1(3,2)*T2(4,3) + 
     -  (2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) - 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(3,4)*
     -   T2(4,3) + (-2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) + 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(3,4)*T2(4,3)
     -   + (2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) - 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(3,4)*T2(4,3)
     -   + (2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) - 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)*T1(3,4)*T2(4,3)
     -   + 2*(-(pT(3)*pT(n)*MT(1,m)) - pT(3)*pT(m)*MT(1,n) - 
     -     pT(1)*pT(n)*MT(3,m) - pT(1)*pT(m)*MT(3,n) + 
     -     (2*pT(1)*pT(3)*MT(m,n))/3.)*T1(4,1)*T2(4,3) - 
     -  2*(-(pT(3)*pT(n)*MT(2,m)) - pT(3)*pT(m)*MT(2,n) - 
     -     pT(2)*pT(n)*MT(3,m) - pT(2)*pT(m)*MT(3,n) + 
     -     (2*pT(2)*pT(3)*MT(m,n))/3.)*T1(4,2)*T2(4,3) + 
     -  (-2*(-(pT(1)*pT(n)*MT(1,m)) - pT(1)*pT(m)*MT(1,n)) + 
     -     (2*(-(pT(m)*pT(n)) - pT(1)**2*MT(m,n)))/3.)*T1(4,3)*
     -   T2(4,3) + (2*(-(pT(2)*pT(n)*MT(2,m)) - 
     -        pT(2)*pT(m)*MT(2,n)) - 
     -     (2*(pT(m)*pT(n) - pT(2)**2*MT(m,n)))/3.)*T1(4,3)*T2(4,3)
     -   + (-2*(-(pT(3)*pT(n)*MT(3,m)) - pT(3)*pT(m)*MT(3,n)) + 
     -     (2*(pT(m)*pT(n) - pT(3)**2*MT(m,n)))/3.)*T1(4,3)*T2(4,3)
     -   + (-2*(-(pT(4)*pT(n)*MT(4,m)) - pT(4)*pT(m)*MT(4,n)) + 
     -     (2*(pT(m)*pT(n) - pT(4)**2*MT(m,n)))/3.)
     -    *T1(4,3)*T2(4,3))/tmass*2

           enddo
         enddo



        if ( tmass.gt.rZero ) then
         d = - 1.0d0/dcmplx( pT2-tmass**2, tmass*twidth )
        else
         d = - 1.0d0/dcmplx( pT2, rZero )
        end if

         do m=1,4
           do n=1,4
          T(m,n)=T(m,n)*gt*d/2.0d0
           enddo
         enddo
c     2.0 factor from propagator convention

         
      wt(1) = T(1,1)
      wt(2) = T(1,2)
      wt(3) = T(1,3)
      wt(4) = T(1,4)
      wt(5) = T(2,1)
      wt(6) = T(2,2)
      wt(7) = T(2,3)
      wt(8) = T(2,4)
      wt(9) = T(3,1)
      wt(10) = T(3,2)
      wt(11) = T(3,3)
      wt(12) = T(3,4)
      wt(13) = T(4,1)
      wt(14) = T(4,2)
      wt(15) = T(4,3)
      wt(16) = T(4,4)
      wt(17) = T(5,1)
      wt(18) = T(6,1)

      return
      end
      subroutine uvvaxx(w1,w2,g,xm1,xm2,xw,jt)
c
c- by RF - Feb. 2006
c
c This subroutine computes the portion of the off-shell current
c for the color-octect tensor jt in terms of w1 and w2 
c
c input:
c       complex w1(6)          : first  vector                        w1
c       complex w2(6)          : second vector                        w2
c       real    g              : first  coupling constant
c       real    xm2            : not used
c       real    xm1            : not used
c       real    xw             : not used
c
c output:
c       complex jt(18)        : tensor current  j^(mu,nu)(w':w1,w2,w3)
c
      implicit none

c dimension of the current set to arbitrary length
c      integer DIM
c      parameter (DIM=18)
      include "dimension.inc"
      double complex w1(DIM),w2(DIM),jt(DIM)
      double precision xm1,xm2,xw,g,s2g
      double precision sqrTwo
      parameter( sqrTwo = 1.41421356237309514547462185873882845044d0 )


c the tensor particle does not propagate, so no propagator needed.

      s2g = g / sqrTwo
      
      jt( 1) = (0D0,0D0) ! s2g * (w1(1)*w2(1)-w1(1)*w2(1))
      jt( 2) =  s2g * (w1(1)*w2(2)-w1(2)*w2(1))
      jt( 3) =  s2g * (w1(1)*w2(3)-w1(3)*w2(1))
      jt( 4) =  s2g * (w1(1)*w2(4)-w1(4)*w2(1))

      jt( 5) =  s2g * (w1(2)*w2(1)-w1(1)*w2(2))
      jt( 6) = (0D0,0D0) ! s2g * (w1(2)*w2(2)-w1(2)*w2(2))
      jt( 7) =  s2g * (w1(2)*w2(3)-w1(3)*w2(2))
      jt( 8) =  s2g * (w1(2)*w2(4)-w1(4)*w2(2))

      jt( 9) =  s2g * (w1(3)*w2(1)-w1(1)*w2(3))
      jt(10) =  s2g * (w1(3)*w2(2)-w1(2)*w2(3))
      jt(11) = (0D0,0D0) ! s2g * (w1(3)*w2(3)-w1(3)*w2(3))
      jt(12) =  s2g * (w1(3)*w2(4)-w1(4)*w2(3))

      jt(13) =  s2g * (w1(4)*w2(1)-w1(1)*w2(4))
      jt(14) =  s2g * (w1(4)*w2(2)-w1(2)*w2(4))
      jt(15) =  s2g * (w1(4)*w2(3)-w1(3)*w2(4))
      jt(16) = (0D0,0D0) ! s2g * (w1(4)*w2(4)-w1(4)*w2(4))

      jt(17) = w1(5) + w2(5)
      jt(18) = w1(6) + w2(6)

      return
      end
      subroutine uvvcxx(v1,v2,gt,vmass,tmass,twidth , uvvh)
c
c This subroutine computes an off-shell tensor current from 
c the two gauge bosons-pseudo tensor boson coupling.
c
c input:
c       complex v1(3)          : first  vector                        v1
c       complex v2(3)          : second vector                        v2
c       real    gt             : coupling constant                gtv= gs
c       real    vmass          : vector boson mass                   m_v
c       real    tmass          : mass  of output tensor T
c       real    twidth         : width of output tensor T
c
c output:
c       complex uvv(18)        : tensor current         j^mu^nu(T:v1,v2)
c
      implicit none
      double complex v1(6), v2(6), uvvh(18)
      double precision vmass, tmass, twidth,gt
      integer i,j
      double complex yvv(6,4)
      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      yvv(5,1) = v1(5)+v2(5)
      yvv(6,1) = v1(6)+v2(6)

      do i=1,4
      do j=1,4
      yvv(i,j)=gt*v1(i)*v2(j)
      enddo
      enddo

      uvvh(1) = yvv(1,1)
      uvvh(2) = yvv(1,2)
      uvvh(3) = yvv(1,3)
      uvvh(4) = yvv(1,4)
      uvvh(5) = yvv(2,1)
      uvvh(6) = yvv(2,2)
      uvvh(7) = yvv(2,3)
      uvvh(8) = yvv(2,4)
      uvvh(9) = yvv(3,1)
      uvvh(10) = yvv(3,2)
      uvvh(11) = yvv(3,3)
      uvvh(12) = yvv(3,4)
      uvvh(13) = yvv(4,1)
      uvvh(14) = yvv(4,2)
      uvvh(15) = yvv(4,3)
      uvvh(16) = yvv(4,4)
      uvvh(17) = yvv(5,1)
      uvvh(18) = yvv(6,1)

      return
      end
      subroutine uvvvlx(ga,gb,gc,g1,g2,xm,xw,jhvvv)
c
c- by RF - Mar. 2006
c
c This subroutine computes an off-shell (non-propagating) scalar particle
c from three incoming vector bosons of the coupling of three gauge bosons.
c
c input:
c       complex ga(6)          : first  incoming vector
c       complex gb(6)          : second incoming vector
c       complex gc(6)          : third  incoming vector
c       real    g1             : coupling constant     (QCD)
c       complex g2(2)          : coupling constant     (Higgs Effct. Thr.)
c
c output:
c       complex jhvvv(3)       : output scalar
c
c not used:
c       xm,xw
c

      implicit none

c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex ga(DIM),gb(DIM),gc(DIM),jhvvv(DIM),g2(2)

      double complex dvertx, vertex, vertex1, vertex2
      double complex vab, vbc, vca, v123, v124, v134, v234
      double complex pgagb, pgagc, pgbga, pgbgc, pgcga, pgcgb
      double precision pga(0:3),pgb(0:3),pgc(0:3),pabc(4)
      double precision g1,xm,xw, q2, q(4)

      pga(0) = dble( ga(5))
      pga(1) = dble( ga(6))
      pga(2) = dimag(ga(6))
      pga(3) = dimag(ga(5))

      pgb(0) = dble( gb(5))
      pgb(1) = dble( gb(6))
      pgb(2) = dimag(gb(6))
      pgb(3) = dimag(gb(5))

      pgc(0) = dble( gc(5))
      pgc(1) = dble( gc(6))
      pgc(2) = dimag(gc(6))
      pgc(3) = dimag(gc(5))

      vertex1 = (0D0,0D0)
      vertex2 = (0D0,0D0)


      jhvvv(2) = ga(5) + gb(5) + gc(5)
      jhvvv(3) = ga(6) + gb(6) + gc(6)


c The internal particle does not propagate, so no multiplication
c with a propagator.

      if (g2(1).NE.(0D0,0D0)) then
      vab = ga(1)*gb(1)-ga(2)*gb(2)-ga(3)*gb(3)-ga(4)*gb(4)
      vbc = gb(1)*gc(1)-gb(2)*gc(2)-gb(3)*gc(3)-gb(4)*gc(4)
      vca = gc(1)*ga(1)-gc(2)*ga(2)-gc(3)*ga(3)-gc(4)*ga(4)

      pgagb = pga(0)*gb(1) - pga(1)*gb(2) - pga(2)*gb(3) - pga(3)*gb(4)
      pgagc = pga(0)*gc(1) - pga(1)*gc(2) - pga(2)*gc(3) - pga(3)*gc(4)
      pgbga = pgb(0)*ga(1) - pgb(1)*ga(2) - pgb(2)*ga(3) - pgb(3)*ga(4)
      pgbgc = pgb(0)*gc(1) - pgb(1)*gc(2) - pgb(2)*gc(3) - pgb(3)*gc(4)
      pgcga = pgc(0)*ga(1) - pgc(1)*ga(2) - pgc(2)*ga(3) - pgc(3)*ga(4)
      pgcgb = pgc(0)*gb(1) - pgc(1)*gb(2) - pgc(2)*gb(3) - pgc(3)*gb(4)

      dvertx = vab*(pgagc-pgbgc) + vbc*(pgbga-pgcga) + vca*(pgcgb-pgagb)
      vertex1= dvertx * g2(1)
      endif

      if (g2(2).NE.(0D0,0D0)) then
      pabc(1) = pga(0) + pgb(0) + pgc(0)
      pabc(2) = pga(1) + pgb(1) + pgc(1)
      pabc(3) = pga(2) + pgb(2) + pgc(2)
      pabc(4) = pga(3) + pgb(3) + pgc(3)

      v123 =   ga(1)*gb(2)*gc(3) - ga(1)*gb(3)*gc(2) - ga(2)*gb(1)*gc(3)
     &       + ga(2)*gb(3)*gc(1) + ga(3)*gb(1)*gc(2) - ga(3)*gb(2)*gc(1)
      v124 = - ga(1)*gb(2)*gc(4) + ga(1)*gb(4)*gc(2) + ga(2)*gb(1)*gc(4)
     &       - ga(2)*gb(4)*gc(1) - ga(4)*gb(1)*gc(2) + ga(4)*gb(2)*gc(1)
      v134 =   ga(1)*gb(3)*gc(4) - ga(1)*gb(4)*gc(3) - ga(3)*gb(1)*gc(4)
     &       + ga(3)*gb(4)*gc(1) + ga(4)*gb(1)*gc(3) - ga(4)*gb(3)*gc(1)
      v234 = - ga(2)*gb(3)*gc(4) + ga(2)*gb(4)*gc(3) + ga(3)*gb(2)*gc(4)
     &       - ga(3)*gb(4)*gc(2) - ga(4)*gb(2)*gc(3) + ga(4)*gb(3)*gc(2)


      vertex2= g2(2) * (  v123*pabc(4) + v124*pabc(3)
     &                  + v134*pabc(2) + v234*pabc(1) )
      endif

      jhvvv(1) = g1 * (vertex1 + vertex2)

      return
      end
      subroutine uvvvxx(va,vb,vc,gc,gt,tmass,twidth , uvvv)
c
c This subroutine computes an off-shell tensor current 
c from the four-point coupling of three gauge bosons and a tensor boson.
c
c input:
c       complex va(6)          : first  vector                        va
c       complex vb(6)          : second vector                        vb
c       complex vc(6)          : third  vector                        vc      
c       real gc                : coupling constant       gs (for gluons)
c       complex gt             : coupling constant         gtv=-1/Lambda
c       real tmass             : mass  of output tensor T
c       real twidth            : width of output tensor T 
c
c output:
c       complex uvvv(18)       : tensor current      j^mu^nu(T:va,vb,vc)
c
c- by Q.Li - OCT. 2006
c- Added massless tensor - P. de Aquino - Oct. 2009 
c     
      implicit none
      double complex va(6), vb(6), vc(6), gt, uvvv(18)
      double precision gc, tmass, twidth

      double complex yvvv(6,4)
      double precision MET(4,4)
      double complex d
      double complex V1V2,V1V3,V2V3,KTV1,KTV2,KTV3,
     &K12V3,K23V1,K31V2
      double precision KTK12,KTK31,KTK23 

      double precision pva(4), pvb(4), pvc(4),pt(4),
     &pt2,p31(4),p23(4),p12(4)

      integer i, j

      double complex cZero
      double precision rZero, r2, r3,r4
      parameter( rZero = 0.0d0, r2 = 2.0d0, r3=3.d0,r4=4.d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )


      yvvv(5,1) = va(5)+vb(5)+vc(5)
      yvvv(6,1) = va(6)+vb(6)+vc(6)

      pva(1) = dreal(va(5))
      pva(2) = dreal(va(6))
      pva(3) = dimag(va(6))
      pva(4) = dimag(va(5))

      pvb(1) = dreal(vb(5))
      pvb(2) = dreal(vb(6))
      pvb(3) = dimag(vb(6))
      pvb(4) = dimag(vb(5))

      pvc(1) = dreal(vc(5))
      pvc(2) = dreal(vc(6))
      pvc(3) = dimag(vc(6))
      pvc(4) = dimag(vc(5))

      pt(1) = dreal(yvvv(5,1))
      pt(2) = dreal(yvvv(6,1))
      pt(3) = dimag(yvvv(6,1))
      pt(4) = dimag(yvvv(5,1))
	
      pt2=pt(1)**2-pt(2)**2-pt(3)**2-pt(4)**2
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      p31(1) = pvc(1)-pva(1)
      p31(2) = pvc(2)-pva(2)
      p31(3) = pvc(3)-pva(3)
      p31(4) = pvc(4)-pva(4)
      
      p12(1) = pva(1)-pvb(1)
      p12(2) = pva(2)-pvb(2)
      p12(3) = pva(3)-pvb(3)
      p12(4) = pva(4)-pvb(4)
      
      p23(1) = pvb(1)-pvc(1)
      p23(2) = pvb(2)-pvc(2)
      p23(3) = pvb(3)-pvc(3)
      p23(4) = pvb(4)-pvc(4)
      
	
      V1V2 =  va(1)*vb(1) -  va(2)*vb(2) -  va(3)*vb(3) -  va(4)*vb(4)
      V1V3 =  va(1)*vc(1) -  va(2)*vc(2) -  va(3)*vc(3) -  va(4)*vc(4)
      V2V3 =  vc(1)*vb(1) -  vc(2)*vb(2) -  vc(3)*vb(3) -  vc(4)*vb(4)
      K31V2 = p31(1)*vb(1) - p31(2)*vb(2) - p31(3)*vb(3) - p31(4)*vb(4)
      K12V3 = p12(1)*vc(1) - p12(2)*vc(2) - p12(3)*vc(3) - p12(4)*vc(4)
      K23V1 = p23(1)*va(1) - p23(2)*va(2) - p23(3)*va(3) - p23(4)*va(4)
      
      KTV1 = pt(1)*va(1) - pt(2)*va(2) - pt(3)*va(3) - pt(4)*va(4)
      KTV2 = pt(1)*vb(1) - pt(2)*vb(2) - pt(3)*vb(3) - pt(4)*vb(4)
      KTV3 = pt(1)*vc(1) - pt(2)*vc(2) - pt(3)*vc(3) - pt(4)*vc(4)
      
      KTK12 =pt(1)*p12(1)-pt(2)*p12(2)-pt(3)*p12(3)-pt(4)*p12(4)
      KTK31 =pt(1)*p31(1)-pt(2)*p31(2)-pt(3)*p31(3)-pt(4)*p31(4)
      KTK23 =pt(1)*p23(1)-pt(2)*p23(2)-pt(3)*p23(3)-pt(4)*p23(4)
      
      if ( tmass.eq.rZero ) then
         d =  -gc/dcmplx( pt2, rZero )  
      
      do i=1,4
         do j=1,4

            yvvv(i,j) =  -(K12V3*V1V2*MET(i,j)) 
     &- K31V2*V1V3*MET(i,j) - K23V1*V2V3*MET(i,j)  
     &+ V2V3*p23(j)*va(i) + V2V3*p23(i)*va(j) 
     &+ V1V3*p31(j)*vb(i) + K12V3*va(j)*vb(i)  
     &+ V1V3*p31(i)*vb(j) + K12V3*va(i)*vb(j) 
     &+ V1V2*p12(j)*vc(i) + K31V2*va(j)*vc(i)  
     &+ K23V1*vb(j)*vc(i) + V1V2*p12(i)*vc(j) 
     &+ K31V2*va(i)*vc(j) + K23V1*vb(i)*vc(j)
	
            yvvv(i,j) = -yvvv(i,j)*d*gt
	
         enddo
      enddo

      else
         if ( tmass.gt.rZero ) then
            d =  -gc/dcmplx( pt2-tmass**2, tmass*twidth )

         do i=1,4
            do j=1,4

              yvvv(i,j) = (r2*K12V3*KTV1*KTV2*r2*MET(i,j))/(r3*tmass**2)
     &	 + (r2*K31V2*KTV1*KTV3*r2*MET(i,j))/(r3*tmass**2) + 
     &  (r2*K23V1*KTV2*KTV3*r2*MET(i,j))/(r3*tmass**2) 
     &- r2*K12V3*V1V2*MET(i,j) + 
     &  (r2*KTK12*KTV3*r2*V1V2*MET(i,j))/(r3*tmass**2) 
     &- (K12V3*pt2*r2*V1V2*MET(i,j))/(r3*tmass**2) - 
     &  r2*K31V2*V1V3*MET(i,j)
     & + (r2*KTK31*KTV2*r2*V1V3*MET(i,j))/(r3*tmass**2) - 
     &  (K31V2*pt2*r2*V1V3*MET(i,j))/(r3*tmass**2) 
     &- r2*K23V1*V2V3*MET(i,j) + 
     &  (r2*KTK23*KTV1*r2*V2V3*MET(i,j))/(r3*tmass**2)
     & - (K23V1*pt2*r2*V2V3*MET(i,j))/(r3*tmass**2) - 
     &  (r2*KTV3*V1V2*p12(j)*pt(i))/tmass**2 
     &- (r2*KTV1*V2V3*p23(j)*pt(i))/tmass**2 
     &- (r2*KTV2*V1V3*p31(j)*pt(i))/tmass**2 - 
     &  (r2*KTV3*V1V2*p12(i)*pt(j))/tmass**2
     & - (r2*KTV1*V2V3*p23(i)*pt(j))/tmass**2 
     &- (r2*KTV2*V1V3*p31(i)*pt(j))/tmass**2 + 
     &  (r4*K12V3*KTV1*KTV2*pt(i)*pt(j))/tmass**4
     & + (r4*K31V2*KTV1*KTV3*pt(i)*pt(j))/tmass**4 + 
     &  (r4*K23V1*KTV2*KTV3*pt(i)*pt(j))/tmass**4 
     &- (r2*K12V3*KTV1*KTV2*r2*pt(i)*pt(j))/(r3*tmass**4) - 
     &  (r2*K31V2*KTV1*KTV3*r2*pt(i)*pt(j))/(r3*tmass**4)
     & - (r2*K23V1*KTV2*KTV3*r2*pt(i)*pt(j))/(r3*tmass**4) + 
     &  (r4*KTK12*KTV3*V1V2*pt(i)*pt(j))/tmass**4 
     &- (r2*K12V3*pt2*V1V2*pt(i)*pt(j))/tmass**4 - 
     &  (r2*KTK12*KTV3*r2*V1V2*pt(i)*pt(j))/(r3*tmass**4)
     & + (K12V3*pt2*r2*V1V2*pt(i)*pt(j))/(r3*tmass**4) + 
     &  (r4*K12V3*V1V2*pt(i)*pt(j))/tmass**2 
     &+ (r4*KTK31*KTV2*V1V3*pt(i)*pt(j))/tmass**4 - 
     &  (r2*K31V2*pt2*V1V3*pt(i)*pt(j))/tmass**4
     & - (r2*KTK31*KTV2*r2*V1V3*pt(i)*pt(j))/(r3*tmass**4) + 
     &  (K31V2*pt2*r2*V1V3*pt(i)*pt(j))/(r3*tmass**4)
     & + (r4*K31V2*V1V3*pt(i)*pt(j))/tmass**2 + 
     &  (r4*KTK23*KTV1*V2V3*pt(i)*pt(j))/tmass**4
     & - (r2*K23V1*pt2*V2V3*pt(i)*pt(j))/tmass**4 - 
     &  (r2*KTK23*KTV1*r2*V2V3*pt(i)*pt(j))/(r3*tmass**4)
     & + (K23V1*pt2*r2*V2V3*pt(i)*pt(j))/(r3*tmass**4) + 
     &  (r4*K23V1*V2V3*pt(i)*pt(j))/tmass**2 
     &+ r2*V2V3*p23(j)*va(i) - (r2*K12V3*KTV2*pt(j)*va(i))/tmass**2 - 
     &  (r2*K31V2*KTV3*pt(j)*va(i))/tmass**2 
     &- (r2*KTK23*V2V3*pt(j)*va(i))/tmass**2 + 2*V2V3*p23(i)*va(j) - 
     &  (r2*K12V3*KTV2*pt(i)*va(j))/tmass**2 
     &- (r2*K31V2*KTV3*pt(i)*va(j))/tmass**2 
     &- (r2*KTK23*V2V3*pt(i)*va(j))/tmass**2 + 
     &  r2*V1V3*p31(j)*vb(i) - (r2*K12V3*KTV1*pt(j)*vb(i))/tmass**2 
     &- (r2*K23V1*KTV3*pt(j)*vb(i))/tmass**2 - 
     &  (r2*KTK31*V1V3*pt(j)*vb(i))/tmass**2 
     &+ r2*K12V3*va(j)*vb(i) + r2*V1V3*p31(i)*vb(j)
     & - (r2*K12V3*KTV1*pt(i)*vb(j))/tmass**2 - 
     &  (r2*K23V1*KTV3*pt(i)*vb(j))/tmass**2 
     &- (r2*KTK31*V1V3*pt(i)*vb(j))/tmass**2 
     &+ r2*K12V3*va(i)*vb(j) + r2*V1V2*p12(j)*vc(i) - 
     &  (r2*K31V2*KTV1*pt(j)*vc(i))/tmass**2 
     &- (r2*K23V1*KTV2*pt(j)*vc(i))/tmass**2 
     &- (r2*KTK12*V1V2*pt(j)*vc(i))/tmass**2 + 
     &  r2*K31V2*va(j)*vc(i) + r2*K23V1*vb(j)*vc(i) 
     &+ r2*V1V2*p12(i)*vc(j) - (r2*K31V2*KTV1*pt(i)*vc(j))/tmass**2 - 
     &  (r2*K23V1*KTV2*pt(i)*vc(j))/tmass**2
     & - (r2*KTK12*V1V2*pt(i)*vc(j))/tmass**2
     & + r2*K31V2*va(i)*vc(j) + r2*K23V1*vb(i)*vc(j)
	
              yvvv(i,j) = -yvvv(i,j)*d/2.0d0*gt
	
            enddo
         enddo
         else
            write(*,*) '’nvalid tensor mass'
         end if
      end if
      
      uvvv(1) = yvvv(1,1)
      uvvv(2) = yvvv(1,2)
      uvvv(3) = yvvv(1,3)
      uvvv(4) = yvvv(1,4)
      uvvv(5) = yvvv(2,1)
      uvvv(6) = yvvv(2,2)
      uvvv(7) = yvvv(2,3)
      uvvv(8) = yvvv(2,4)
      uvvv(9) = yvvv(3,1)
      uvvv(10) = yvvv(3,2)
      uvvv(11) = yvvv(3,3)
      uvvv(12) = yvvv(3,4)
      uvvv(13) = yvvv(4,1)
      uvvv(14) = yvvv(4,2)
      uvvv(15) = yvvv(4,3)
      uvvv(16) = yvvv(4,4)
      uvvv(17) = yvvv(5,1)
      uvvv(18) = yvvv(6,1)

      return
      end
      subroutine uvvxxx(v1,v2,gt,vmass,tmass,twidth , uvv)
c
c This subroutine computes an off-shell tensor current from 
c the two gauge bosons-tensor boson coupling.
c
c input:
c       complex v1(3)          : first  vector                        v1
c       complex v2(3)          : second vector                        v2
c       complex gt             : coupling constant         gtv=-1/Lambda
c       real    vmass          : vector boson mass                   m_v
c       real    tmass          : mass  of output tensor T
c       real    twidth         : width of output tensor T
c
c output:
c       complex uvv(18)        : tensor current         j^mu^nu(T:v1,v2)
c
c- by Q.Li - OCT. 2006
c- Added massless tensor - P. de Aquino - Oct. 2009 
c
      implicit none
      double complex v1(6), v2(6), gt, uvv(18)
      double precision vmass, tmass, twidth

      double complex yvv(6,4)
      double complex KTE1,KTE2,K2E1,K1E2,K1E1,K2E2,E1E2
      integer i,j
      double precision pv1(4), pv2(4), pT(4)
      double precision MET(4,4)
      double complex cZero, d
      double precision rZero, rTwo
      double precision K1K2,KT2,K1KT,K2KT
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      yvv(5,1) = v1(5)+v2(5)
      yvv(6,1) = v1(6)+v2(6)
      
      pv1(1) = dreal(v1(5))
      pv1(2) = dreal(v1(6))
      pv1(3) = dimag(v1(6))
      pv1(4) = dimag(v1(5))
      
      pv2(1) = dreal(v2(5))
      pv2(2) = dreal(v2(6))
      pv2(3) = dimag(v2(6))
      pv2(4) = dimag(v2(5))
      
      pT(1) = dreal(yvv(5,1))
      pT(2) = dreal(yvv(6,1))
      pT(3) = dimag(yvv(6,1))
      pT(4) = dimag(yvv(5,1))
      
      do i=1,4
         do j=1,4
            MET(i,j) = 0.0d0
         enddo 
      enddo
      
      MET(1,1) =  1.0d0
      MET(2,2) = -1.0d0
      MET(3,3) = -1.0d0
      MET(4,4) = -1.0d0
      
      K1K2 = pv1(1)*pv2(1)-pv1(2)*pv2(2)-pv1(3)*pv2(3)-pv1(4)*pv2(4)
      K1KT = pv1(1)*pT(1)-pv1(2)*pT(2)-pv1(3)*pT(3)-pv1(4)*pT(4)
      K2KT = pT(1)*pv2(1)-pT(2)*pv2(2)-pT(3)*pv2(3)-pT(4)*pv2(4)
      KT2 = pT(1)**2-pT(2)**2-pT(3)**2-pT(4)**2
      
      KTE1 = pT(1)*v1(1)-pT(2)*v1(2)-pT(3)*v1(3)-pT(4)*v1(4)
      KTE2 = pT(1)*v2(1)-pT(2)*v2(2)-pT(3)*v2(3)-pT(4)*v2(4)
      K1E1 = pv1(1)*v1(1)-pv1(2)*v1(2)-pv1(3)*v1(3)-pv1(4)*v1(4)
      K1E2 = pv1(1)*v2(1)-pv1(2)*v2(2)-pv1(3)*v2(3)-pv1(4)*v2(4)
      K2E1 = pv2(1)*v1(1)-pv2(2)*v1(2)-pv2(3)*v1(3)-pv2(4)*v1(4)
      K2E2 = pv2(1)*v2(1)-pv2(2)*v2(2)-pv2(3)*v2(3)-pv2(4)*v2(4)
      E1E2 = v2(1)*v1(1)-v2(2)*v1(2)-v2(3)*v1(3)-v2(4)*v1(4)
      
      if ( tmass.eq.rZero ) then
         d = - gt/dcmplx( KT2, rZero )
    
         do i = 1,4
            do j=1,4
               yvv(i,j) = -(E1E2*K1K2*MET(i,j)) + K1E2*K2E1*MET(i,j)  
     &+  E1E2*pv1(j)*pv2(i) + E1E2*pv1(i)*pv2(j)  
     &-  K1E2*pv2(j)*v1(i) - K1E2*pv2(i)*v1(j)  
     &-  K2E1*pv1(j)*v2(i) + k1k2*v1(j)*v2(i)  
     &-  K2E1*pv1(i)*v2(j) + k1k2*v1(i)*v2(j) 
               
            
               if ( vmass.ne.rZero ) then
                  yvv(i,j) = yvv(i,j)
     &              + vmass**2*v1(j)*v2(i) + vmass**2*v1(i)*v2(j)
               else
c     gauge fixing term for zero mass photon/gluon
                  yvv(i,j) =  yvv(i,j)
     &             -(K1E1*K2E2*MET(i,j)) - K2E2*pv2(j)*v1(i) - 
     &  K2E2*pv2(i)*v1(j) - K1E1*pv1(j)*v2(i) - 
     &  K1E1*pv1(i)*v2(j)
               endif
               yvv(i,j) = yvv(i,j)*d
            end do
         enddo
      
      else if ( tmass.gt.rZero ) then
            d = - gt/dcmplx( KT2-tmass**2, tmass*twidth )
    
            do i = 1,4
               do j=1,4
                  yvv(i,j) = 2.0d0*K1K2*(v1(i)*v2(j)+v1(j)*v2(i))
     &-2.0d0*K1K2*KTE2/tmass**2*(PT(i)*v1(j)+PT(j)*v1(i))
     &-2.0d0*K1E2*(pv2(i)*v1(j)+pv2(j)*v1(i))
     &+2.0d0*K1E2*K2KT/tmass**2*(PT(i)*v1(j)+PT(j)*v1(i))
     &-2.0d0/3.0d0*E1E2*K1K2*MET(i,j)
     &+8.0d0/3.0d0*K1K2*E1E2/tmass**2*PT(i)*PT(j)
     &+2.0d0*E1E2*(pv1(i)*pv2(j)+pv1(j)*pv2(i)) 
     &-2.0d0*K1K2*KTE1/tmass**2*(PT(i)*v2(j)+PT(j)*v2(i))
     &-2.0d0*K2E1*(pv1(i)*v2(j)+pv1(j)*v2(i))
     &+4.0d0*K1K2*KTE1*KTE2/3.0d0/tmass**2*MET(i,j)
     &+8.0d0*K1K2*KTE1*KTE2/3.0d0/tmass**4*PT(i)*PT(j)
     &+2.0d0*K2E1*KTE2/tmass**2*(pv1(i)*PT(j)+pv1(j)*PT(i))
     &+2.0d0*KTE1*K1E2/tmass**2*(pv2(i)*PT(j)+pv2(j)*PT(i))
     &+2.0d0*K2E1*K1E2*MET(i,j)
     &-4.0d0*K2E1*K1E2/tmass**2*PT(i)*PT(j)
     &-2.0d0/3.0d0/tmass**2*KT2*E1E2*K1K2*MET(i,j)
     &-4.0d0/3.0d0/tmass**4*KT2*E1E2*K1K2*PT(i)*PT(j)
     &+2.0d0/3.0d0/tmass**2*K2E1*K1E2*KT2*MET(i,j)
     &+4.0d0/3.0d0/tmass**4*K2E1*K1E2*KT2*PT(i)*PT(j)
     &-2.0d0/tmass**2*E1E2*K1KT*(pv2(i)*PT(j)+pv2(j)*PT(i))
     &+2.0d0/tmass**2*K2E1*K1KT*(PT(i)*v2(j)+PT(j)*v2(i))
     &-4.0d0/3.0d0/tmass**2*K2E1*KTE2*K1KT*MET(i,j)
     &-8.0d0/3.0d0/tmass**4*K2E1*KTE2*K1KT*PT(i)*PT(j)
     &-2.0d0/tmass**2*E1E2*K2KT*(PT(i)*pv1(j)+PT(j)*pv1(i))
     &-4.0d0/3.0d0/tmass**2*KTE1*K1E2*K2KT*MET(i,j)
     &-8.0d0/3.0d0/tmass**4*KTE1*K1E2*K2KT*PT(i)*PT(j)
     &+4.0d0/3.0d0/tmass**2*E1E2*K2KT*K1KT*MET(i,j)
     &+8.0d0/3.0d0/tmass**4*E1E2*K2KT*K1KT*PT(i)*PT(j)
     &-4.0d0/3.0d0*E1E2*K1K2*MET(i,j)
     &+4.0d0/3.0d0/tmass**2*E1E2*K1K2
     &*PT(i)*PT(j)
            
                  if ( vmass.ne.rZero ) then
                     yvv(i,j) = 
     &                 yvv(i,j)+2*vmass**2*(v1(i)*v2(j)+v1(j)*v2(i))
     &-2.0d0/3.0d0*vmass**2*E1E2*MET(i,j)
     &+8.0d0/3.0d0/tmass**2*vmass**2*E1E2*PT(i)*PT(j)
     &-2.0d0/tmass**2*vmass**2*KTE1*(PT(i)*v2(j)+PT(j)*v2(i)) 
     &-2.0d0/tmass**2*vmass**2*KTE2*(PT(i)*v1(j)+PT(j)*v1(i)) 
     &+4.0d0/3.0d0/tmass**2*vmass**2*KTE1*KTE2*MET(i,j)
     &+8.0d0/3.0d0/tmass**4*vmass**2*KTE1*KTE2*PT(i)*PT(j)
     &-2.0d0/3.0d0/tmass**2*vmass**2*KT2*E1E2*MET(i,j)
     &-4.0d0/3.0d0/tmass**4*vmass**2*KT2*E1E2*PT(i)*PT(j)
                  else
c     gauge fixing term for zero mass photon/gluon
                     yvv(i,j) = 
     &                yvv(i,j)-2.0d0*K1E1*(pv1(i)*v2(j)+pv1(j)*v2(i))
     &+2.0d0/tmass**2*KTE2*K1E1*(PT(i)*pv1(j)+PT(j)*pv1(i))
     &+2.0d0/3.0d0*K1E2*K1E1*MET(i,j)
     &-8.0d0/3.0d0/tmass**2*K1E2*K1E1*PT(i)*PT(j)
     &-2.0d0/3.0d0*K2E2*K1E1*MET(i,j)
     &-4.0d0/3.0d0/tmass**2*K2E2*K1E1*PT(i)*PT(j)
     &+2.0d0/3.0d0/tmass**2*K1E2*KT2*K1E1*MET(i,j)
     &+4.0d0/3.0d0/tmass**4*K1E2*KT2*K1E1*PT(i)*PT(j)
     &+2.0d0/3.0d0/tmass**2*K2E2*KT2*K1E1*MET(i,j)
     &+4.0d0/3.0d0/tmass**4*K2E2*KT2*K1E1*PT(i)*PT(j)
     &+2.0d0/tmass**2*K1KT*K1E1*(PT(i)*v2(j)+PT(j)*v2(i))
     &-4.0d0/3.0d0/tmass**2*K1KT*KTE2*K1E1*MET(i,j)
     &-8.0d0/3.0d0/tmass**4*K1KT*KTE2*K1E1*PT(i)*PT(j)
     &-2.0d0*K2E2*(pv2(i)*v1(j)+pv2(j)*v1(i))
     &+2.0d0/tmass**2*K2E2*KTE1*(pv2(i)*PT(j)+pv2(j)*PT(i))
     &+2.0d0/3.0d0*K2E1*K2E2*MET(i,j)
     &-8.0d0/3.0d0/tmass**2*K2E1*K2E2*PT(i)*PT(j)
     &+2.0d0/3.0d0/tmass**2*KT2*K2E2*K2E1*MET(i,j)
     &+4.0d0/3.0d0/tmass**4*KT2*K2E2*K2E1*PT(i)*PT(j)
     &+2.0d0/tmass**2*K2E2*K2KT*(PT(i)*v1(j)+PT(j)*v1(i))
     &-4.0d0/3.0d0/tmass**2*K2E2*K2KT*KTE1*MET(i,j)
     &-8.0d0/3.0d0/tmass**4*K2E2*K2KT*KTE1
     &*PT(i)*PT(j)
                  endif

                  yvv(i,j) = yvv(i,j)*d/2.0d0

               end do
            enddo
         else
            write(*,*) 'invalid tensor mass'
            stop
      end if
      
      uvv(1) = yvv(1,1)
      uvv(2) = yvv(1,2)
      uvv(3) = yvv(1,3)
      uvv(4) = yvv(1,4)
      uvv(5) = yvv(2,1)
      uvv(6) = yvv(2,2)
      uvv(7) = yvv(2,3)
      uvv(8) = yvv(2,4)
      uvv(9) = yvv(3,1)
      uvv(10) = yvv(3,2)
      uvv(11) = yvv(3,3)
      uvv(12) = yvv(3,4)
      uvv(13) = yvv(4,1)
      uvv(14) = yvv(4,2)
      uvv(15) = yvv(4,3)
      uvv(16) = yvv(4,4)
      uvv(17) = yvv(5,1)
      uvv(18) = yvv(6,1)

      return
      end
      subroutine vssxxx(vc,s1,s2,gc , vertex)
c
c This subroutine computes an amplitude from the vector-scalar-scalar
c coupling.  The coupling is absent in the minimal SM in unitary gauge.
c
c       complex vc(6)          : input  vector                        v
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex gc             : coupling constant (s1 charge)
c
c examples of the coupling constant gc for SUSY particles are as follows:
c   -----------------------------------------------------------
c   |    s1    | (q,i3) of s1  ||   v=a   |   v=z   |   v=w   |
c   -----------------------------------------------------------
c   | nu~_l    | (  0  , +1/2) ||   ---   |  gzn(1) |  gwf(1) |
c   | e~_l     | ( -1  , -1/2) ||  gal(1) |  gzl(1) |  gwf(1) |
c   | u~_l     | (+2/3 , +1/2) ||  gau(1) |  gzu(1) |  gwf(1) |
c   | d~_l     | (-1/3 , -1/2) ||  gad(1) |  gzd(1) |  gwf(1) |
c   -----------------------------------------------------------
c   | e~_r-bar | ( +1  ,  0  ) || -gal(2) | -gzl(2) | -gwf(2) |
c   | u~_r-bar | (-2/3 ,  0  ) || -gau(2) | -gzu(2) | -gwf(2) |
c   | d~_r-bar | (+1/3 ,  0  ) || -gad(2) | -gzd(2) | -gwf(2) |
c   -----------------------------------------------------------
c where the s1 charge is defined by the flowing-OUT quantum number.
c
c output:
c       complex vertex         : amplitude                gamma(v,s1,s2)
c     
      implicit none
      double complex vc(6),s1(3),s2(3),gc,vertex
      double precision p(0:3)

c

      p(0) = dble( s1(2)-s2(2))
      p(1) = dble( s1(3)-s2(3))
      p(2) = dimag(s1(3)-s2(3))
      p(3) = dimag(s1(2)-s2(2))

      vertex = gc*s1(1)*s2(1)
     &        *(vc(1)*p(0)-vc(2)*p(1)-vc(3)*p(2)-vc(4)*p(3))
c
      return
      end
      subroutine vvshxx(v1,v2,sc,gc , vertex)
c
c- by RF - Mar. 2006
c
c
c This subroutine computes an amplitude of the vector-vector-
c (pseudo-)scalar effective coupling.
c
c input:
c       complex v1(6)          : first  vector
c       complex v2(6)          : second vector
c       complex sc(3)          : input  scalar
c       complex gc(2)          : coupling constant: gc(1) scalar
c                                                   gc(2) pseudo-scalar
c
c output:
c       complex vertex         : amplitude
c     
      implicit none
c--   dimension of the current set to arbitrary length
c     INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex v1(DIM),v2(DIM),sc(DIM),vertex,vertex1,vertex2
      double complex v12,p2v1,p1v2,v13,v14,v23,v24,v34
      double precision p12,p13,p14,p23,p24,p34
      double precision p1(0:3),p2(0:3)
      double complex gc(2)


      p1(0) = dble( v1(5))
      p1(1) = dble( v1(6))
      p1(2) = dimag(v1(6))
      p1(3) = dimag(v1(5))

      p2(0) = dble( v2(5))
      p2(1) = dble( v2(6))
      p2(2) = dimag(v2(6))
      p2(3) = dimag(v2(5))

      vertex1 = (0D0,0D0)
      vertex2 = (0D0,0D0)

      if (gc(1).NE.(0D0,0D0)) then

         v12  = v1(1)*v2(1) - v1(2)*v2(2) - v1(3)*v2(3) - v1(4)*v2(4)
         p12  = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
         p2v1 = v1(1)*p2(0) - v1(2)*p2(1) - v1(3)*p2(2) - v1(4)*p2(3)
         p1v2 = p1(0)*v2(1) - p1(1)*v2(2) - p1(2)*v2(3) - p1(3)*v2(4)	

         vertex1 = - gc(1)*(v12*p12 - p2v1*p1v2)
      endif

      if (gc(2).NE.(0D0,0D0)) then
          p12 = p1(0)*p2(1) - p1(1)*p2(0)
          p13 = p1(0)*p2(2) - p1(2)*p2(0)
          p14 = p1(0)*p2(3) - p1(3)*p2(0)
          p23 = p1(1)*p2(2) - p1(2)*p2(1)
          p24 = p1(1)*p2(3) - p1(3)*p2(1)
          p34 = p1(2)*p2(3) - p1(3)*p2(2)

          v12 = v1(1)*v2(2) - v1(2)*v2(1)
          v13 = v1(1)*v2(3) - v1(3)*v2(1)
          v14 = v1(1)*v2(4) - v1(4)*v2(1)
          v23 = v1(2)*v2(3) - v1(3)*v2(2)
          v24 = v1(2)*v2(4) - v1(4)*v2(2)
          v34 = v1(3)*v2(4) - v1(4)*v2(3)

          vertex2 = gc(2)*( v12*p34 - v13*p24 + v14*p23
     &                     +v23*p14 - v24*p13 + v34*p12 )
      endif
       
      vertex = sc(1)*(vertex1 + vertex2)

      return
      end
      subroutine vvsshx(v1,v2,sc1,sc2,g1, vertex)
c
c- by RF - Mar. 2006
c
c
c This subroutine computes an amplitude of the vector-vector-Higgs-Higgs 
c effective coupling.
c
c input:
c       complex v1(6)          : first  vector                        
c       complex v2(6)          : second vector                        
c       complex sc1(3)         : first  scalar                        
c       complex sc2(3)         : second scalar
c       complex g1(2)          : first coupling constant                 
c
c output:
c       complex vertex         : amplitude                gamma(v1,v2,s,s)
c     
      implicit none
c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex v1(DIM),v2(DIM),sc1(DIM),sc2(DIM),vertex
      double complex vertex1,vertex2,g1(2)
      double complex v12,p2v1,p1v2,v13,v14,v23,v24,v34
      double precision p12,p13,p14,p23,p24,p34
      double precision p1(0:3),p2(0:3)


      p1(0) = dble( v1(5))
      p1(1) = dble( v1(6))
      p1(2) = dimag(v1(6))
      p1(3) = dimag(v1(5))

      p2(0) = dble( v2(5))
      p2(1) = dble( v2(6))
      p2(2) = dimag(v2(6))
      p2(3) = dimag(v2(5))

      vertex1 = (0D0,0D0)
      vertex2 = (0D0,0D0)

      if (g1(1).NE.(0D0,0D0)) then

         v12  = v1(1)*v2(1) - v1(2)*v2(2) - v1(3)*v2(3) - v1(4)*v2(4)
         p12  = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
         p2v1 = v1(1)*p2(0) - v1(2)*p2(1) - v1(3)*p2(2) - v1(4)*p2(3)
         p1v2 = p1(0)*v2(1) - p1(1)*v2(2) - p1(2)*v2(3) - p1(3)*v2(4)	

         vertex1 = g1(1) *(v12*p12 - p2v1*p1v2)
      endif

      if (g1(2).NE.(0D0,0D0)) then
          p12 = p1(0)*p2(1) - p1(1)*p2(0)
          p13 = p1(0)*p2(2) - p1(2)*p2(0)
          p14 = p1(0)*p2(3) - p1(3)*p2(0)
          p23 = p1(1)*p2(2) - p1(2)*p2(1)
          p24 = p1(1)*p2(3) - p1(3)*p2(1)
          p34 = p1(2)*p2(3) - p1(3)*p2(2)

          v12 = v1(1)*v2(2) - v1(2)*v2(1)
          v13 = v1(1)*v2(3) - v1(3)*v2(1)
          v14 = v1(1)*v2(4) - v1(4)*v2(1)
          v23 = v1(2)*v2(3) - v1(3)*v2(2)
          v24 = v1(2)*v2(4) - v1(4)*v2(2)
          v34 = v1(3)*v2(4) - v1(4)*v2(3)

          vertex2 = - g1(2)*( v12*p34 - v13*p24 + v14*p23
     &                       +v23*p14 - v24*p13 + v34*p12 )
      endif
       
      vertex = sc1(1)*sc2(1)*(vertex1 + vertex2)

      return
      end
      subroutine vvssxx(v1,v2,s1,s2,gc , vertex)
c
c This subroutine computes an amplitude of the vector-vector-scalar-
c scalar coupling.
c
c input:
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       complex s1(3)          : first  scalar                        s1
c       complex s2(3)          : second scalar                        s2
c       complex gc             : coupling constant                 gvvhh
c
c output:
c       complex vertex         : amplitude            gamma(v1,v2,s1,s2)
c     
      implicit none
      double complex v1(6),v2(6),s1(3),s2(3),gc,vertex

c

      vertex = gc*s1(1)*s2(1)
     &        *(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
c
      return
      end
      subroutine vvstxx(ga,gb,tc,g, vertex)
c
c- by RF - Feb. 2006
c
c This subroutine computes the portion of the amplitude of the four-point 
c coupling of 2 massless color octet gauge bosons (gluons) 
C with an effective color octect antisymmetric tensor.
c                                                                       
c input:
c     complex ga(6)                       : first vector  (gluon)
c     complex gb(6)                       : second vector (gluon)
c     complex tc(18)                      : tensor current
c     real    g                           : coupling constant
c 
c output:
c     complex vertex                      : amplitude
c

      implicit none

c dimension of the current set to arbitrary length
      integer DIM
      parameter (DIM=18)
c      include "dimension.inc"
      double complex ga(DIM),gb(DIM),tc(DIM)

      double precision xm,xw,g

      double precision sqrTwo
      parameter( sqrTwo=1.41421356237309514547462185873882845044d0 )

      double complex dvertx, vertex


      dvertx = ! + tc( 1)*( ga(1)*gb(1) - ga(1)*gb(1) )  Always zero
     &          - tc( 2)*( ga(1)*gb(2) - ga(2)*gb(1) ) 
     &          - tc( 3)*( ga(1)*gb(3) - ga(3)*gb(1) ) 
     &          - tc( 4)*( ga(1)*gb(4) - ga(4)*gb(1) ) 
     
     &          - tc( 5)*( ga(2)*gb(1) - ga(1)*gb(2) ) 
     &         ! + tc( 6)*( ga(2)*gb(2) - ga(2)*gb(2) )  Always zero
     &          + tc( 7)*( ga(2)*gb(3) - ga(3)*gb(2) ) 
     &          + tc( 8)*( ga(2)*gb(4) - ga(4)*gb(2) ) 
     
     &          - tc( 9)*( ga(3)*gb(1) - ga(1)*gb(3) ) 
     &          + tc(10)*( ga(3)*gb(2) - ga(2)*gb(3) ) 
     &         ! + tc(11)*( ga(3)*gb(3) - ga(3)*gb(3) )  Always zero
     &          + tc(12)*( ga(3)*gb(4) - ga(4)*gb(3) ) 
     
     &          - tc(13)*( ga(4)*gb(1) - ga(1)*gb(4) ) 
     &          + tc(14)*( ga(4)*gb(2) - ga(2)*gb(4) ) 
     &          + tc(15)*( ga(4)*gb(3) - ga(3)*gb(4) ) 
     &         ! + tc(16)*( ga(4)*gb(4) - ga(4)*gb(4) )  Always zero

      vertex = g * dvertx /sqrTwo


      return
      end
      subroutine vvsxxx(v1,v2,sc,gc , vertex)
c
c This subroutine computes an amplitude of the vector-vector-scalar
c coupling.
c
c input:
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       complex sc(3)          : input  scalar                        s
c       complex gc             : coupling constant                  gvvh
c
c output:
c       complex vertex         : amplitude                gamma(v1,v2,s)
c     
      implicit none
      double complex v1(6),v2(6),sc(3),gc,vertex

c

      vertex = gc*sc(1)
     &        *(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
c
      return
      end
      subroutine vvtaxx(ga,gb,tc,g,xm, vertex)
c
c- by RF - Feb. 2006
c
c This subroutine computes the portion of the amplitude of the four-point 
c coupling of 2 massless color octet gauge bosons
C with an effective color octect antisymmetric tensor.
c                                                                       
c input:
c     complex ga(6)                       : first vector  (gluon)
c     complex gb(6)                       : second vector (gluon)
c     complex tc(18)                      : tensor current
c     real    g                           : coupling constant
c 
c output:
c     complex vertex                      : amplitude
c

      implicit none

c dimension of the current set to arbitrary length
c      integer DIM
c      parameter (DIM=18)
      include "dimension.inc"
      double complex ga(DIM),gb(DIM),tc(DIM)
      double precision xm,xw,g
      double precision sqrTwo
      parameter( sqrTwo=1.41421356237309514547462185873882845044d0 )

      double complex dvertx, vertex


      dvertx = ! + tc( 1)*( ga(1)*gb(1) - ga(1)*gb(1) )  Always zero
     &          - tc( 2)*( ga(1)*gb(2) - ga(2)*gb(1) ) 
     &          - tc( 3)*( ga(1)*gb(3) - ga(3)*gb(1) ) 
     &          - tc( 4)*( ga(1)*gb(4) - ga(4)*gb(1) ) 
     
     &          - tc( 5)*( ga(2)*gb(1) - ga(1)*gb(2) ) 
     &         ! + tc( 6)*( ga(2)*gb(2) - ga(2)*gb(2) )  Always zero
     &          + tc( 7)*( ga(2)*gb(3) - ga(3)*gb(2) ) 
     &          + tc( 8)*( ga(2)*gb(4) - ga(4)*gb(2) ) 
     
     &          - tc( 9)*( ga(3)*gb(1) - ga(1)*gb(3) ) 
     &          + tc(10)*( ga(3)*gb(2) - ga(2)*gb(3) ) 
     &         ! + tc(11)*( ga(3)*gb(3) - ga(3)*gb(3) )  Always zero
     &          + tc(12)*( ga(3)*gb(4) - ga(4)*gb(3) ) 
     
     &          - tc(13)*( ga(4)*gb(1) - ga(1)*gb(4) ) 
     &          + tc(14)*( ga(4)*gb(2) - ga(2)*gb(4) ) 
     &          + tc(15)*( ga(4)*gb(3) - ga(3)*gb(4) ) 
     &         ! + tc(16)*( ga(4)*gb(4) - ga(4)*gb(4) )  Always zero

      vertex = g * dvertx /sqrTwo


      return
      end
      subroutine vvtcxx(v1,v2,tc,gt,vmass , vertex)
c-------------------CP3  2009.10-----------------
c This subroutine computes an amplitude of the three-point coupling of
c two gauge bosons and a non-propagating tensor boson.
c
c input:
c       complex v1(6)          : first  vector                            v1
c       complex v2(6)          : second vector                        v2
c       complex tc(18)         : input  tensor                           T
c       complex gt             : coupling constant                   gt=gs
c       real    vmass          : vector boson mass                   m_v
c
c output:
c       complex vertex         : amplitude                gamma(v1,v2,T)
c
c     
      implicit none
      double complex v1(6), v2(6), tc(18), vertex
      double precision vmass, gt

      double complex ft(6,4)
      double complex  dum

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

        vertex=ft(1,1)*v1(1)*v2(1) -ft(2,1)*v1(2)*v2(1) - 
     -  ft(3,1)*v1(3)*v2(1) - ft(4,1)*v1(4)*v2(1) - 
     -  ft(1,2)*v1(1)*v2(2) + ft(2,2)*v1(2)*v2(2) + 
     -  ft(3,2)*v1(3)*v2(2) + ft(4,2)*v1(4)*v2(2) - 
     -  ft(1,3)*v1(1)*v2(3) + ft(2,3)*v1(2)*v2(3) + 
     -  ft(3,3)*v1(3)*v2(3) + ft(4,3)*v1(4)*v2(3) - 
     -  ft(1,4)*v1(1)*v2(4) + ft(2,4)*v1(2)*v2(4) + 
     -  ft(3,4)*v1(3)*v2(4) + ft(4,4)*v1(4)*v2(4)

      vertex = vertex * gt

      return
      end
      subroutine vvtxkk(wm,wp,tc,g,vmass , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two gauge bosons and a Kaluza-Klein tensor boson.
c
c input:
c       complex wm(6)          : vector               flow-in  V
c       complex wp(6)          : vector               flow-out V~
c       complex tc(6,4)        : tensor               KK mode T
c       real    g              : coupling constant    -kappa/2
c       real    vmass          : V boson mass          m_V
c
c output:
c       complex vertex         : amplitude            gamma(wm,wp,tc)
c     
      implicit none
      double complex wm(6), wp(6), tc(6,4), vertex
      double precision g, vmass

      double complex T12, T13, T14, T23, T24, T34
      double complex V1V2, k1V2, k2V1
      double complex Tkk, TVV, Tk1V2, Tk2V1, dum
      double precision pwm(4), pwp(4), F

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )
c
      pwm(1) = dreal(wm(5))
      pwm(2) = dreal(wm(6))
      pwm(3) = dimag(wm(6))
      pwm(4) = dimag(wm(5))
      pwp(1) = dreal(wp(5))
      pwp(2) = dreal(wp(6))
      pwp(3) = dimag(wp(6))
      pwp(4) = dimag(wp(5))

      T12 = tc(1,2) + tc(2,1)
      T13 = tc(1,3) + tc(3,1)
      T14 = tc(1,4) + tc(4,1)
      T23 = tc(2,3) + tc(3,2)
      T24 = tc(2,4) + tc(4,2)
      T34 = tc(3,4) + tc(4,3)

      V1V2 =  wm(1)*wp(1) -  wm(2)*wp(2) -  wm(3)*wp(3) -  wm(4)*wp(4)
      k1V2 = pwm(1)*wp(1) - pwm(2)*wp(2) - pwm(3)*wp(3) - pwm(4)*wp(4)
      k2V1 = pwp(1)*wm(1) - pwp(2)*wm(2) - pwp(3)*wm(3) - pwp(4)*wm(4)

      F = pwm(1)*pwp(1) - pwm(2)*pwp(2) - pwm(3)*pwp(3) - pwm(4)*pwp(4)
      if ( vmass.ne.rZero ) then
         F = F + vmass**2
      end if

      Tkk   = cZero
      TVV   = cZero
      Tk1V2 = cZero
      Tk2V1 = cZero

      do i = 1,4
         dum   = tc(i,i)*pwm(i)
         Tkk   = Tkk   + dum*pwp(i)
         Tk1V2 = Tk1V2 + dum*wp(i)
         dum   = tc(i,i)*wm(i)
         TVV   = TVV   + dum*wp(i)
         Tk2V1 = Tk2V1 + dum*pwp(i)
      end do

      Tkk   = rTwo*Tkk
      TVV   = rTwo*TVV
      Tk1V2 = rTwo*Tk1V2
      Tk2V1 = rTwo*Tk2V1

      Tkk = Tkk - T12*(pwm(1)*pwp(2) + pwm(2)*pwp(1))
     &          - T13*(pwm(1)*pwp(3) + pwm(3)*pwp(1))
     &          - T14*(pwm(1)*pwp(4) + pwm(4)*pwp(1))
     &          + T23*(pwm(2)*pwp(3) + pwm(3)*pwp(2))
     &          + T24*(pwm(2)*pwp(4) + pwm(4)*pwp(2))
     &          + T34*(pwm(3)*pwp(4) + pwm(4)*pwp(3))

      Tk1V2 = Tk1V2 - T12*(pwm(1)*wp(2) + pwm(2)*wp(1))
     &              - T13*(pwm(1)*wp(3) + pwm(3)*wp(1))
     &              - T14*(pwm(1)*wp(4) + pwm(4)*wp(1))
     &              + T23*(pwm(2)*wp(3) + pwm(3)*wp(2))
     &              + T24*(pwm(2)*wp(4) + pwm(4)*wp(2))
     &              + T34*(pwm(3)*wp(4) + pwm(4)*wp(3))

      TVV = TVV - T12*(wm(1)*wp(2) + wm(2)*wp(1))
     &          - T13*(wm(1)*wp(3) + wm(3)*wp(1))
     &          - T14*(wm(1)*wp(4) + wm(4)*wp(1))
     &          + T23*(wm(2)*wp(3) + wm(3)*wp(2))
     &          + T24*(wm(2)*wp(4) + wm(4)*wp(2))
     &          + T34*(wm(3)*wp(4) + wm(4)*wp(3))

      Tk2V1 = Tk2V1 - T12*(wm(1)*pwp(2) + wm(2)*pwp(1))
     &              - T13*(wm(1)*pwp(3) + wm(3)*pwp(1))
     &              - T14*(wm(1)*pwp(4) + wm(4)*pwp(1))
     &              + T23*(wm(2)*pwp(3) + wm(3)*pwp(2))
     &              + T24*(wm(2)*pwp(4) + wm(4)*pwp(2))
     &              + T34*(wm(3)*pwp(4) + wm(4)*pwp(3))

      vertex =  (tc(1,1)-tc(2,2)-tc(3,3)-tc(4,4))*( k1V2*k2V1 - V1V2*F )
     &        + F*TVV + V1V2*Tkk - k2V1*Tk1V2 - k1V2*Tk2V1

C      vertex = F*TVV + V1V2*Tkk - k2V1*Tk1V2 - k1V2*Tk2V1

      vertex = vertex * g
c
      return
      end
      subroutine vvtxxx(v1,v2,tc,gt,vmass , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two gauge bosons and a tensor boson.
c
c input:
c       complex v1(6)          : first  vector                        v1
c       complex v2(6)          : second vector                        v2
c       complex tc(18)         : input  tensor                         T
c       complex gt             : coupling constant         gtv=-1/Lambda
c       real    vmass          : vector boson mass                   m_v
c
c output:
c       complex vertex         : amplitude                gamma(v1,v2,T)
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex v1(6), v2(6), tc(18), gt, vertex
      double precision vmass

      double complex ft(6,4)
      double complex T12, T13, T14, T23, T24, T34
      double complex V1V2, K1V2, K2V1
c     new
     &,K1V1,K2V2
c     new
      double complex TKK, TVV, TK1V2, TK2V1, dum
      double precision pv1(4), pv2(4), F

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      pv1(1) = dreal(v1(5))
      pv1(2) = dreal(v1(6))
      pv1(3) = dimag(v1(6))
      pv1(4) = dimag(v1(5))
      pv2(1) = dreal(v2(5))
      pv2(2) = dreal(v2(6))
      pv2(3) = dimag(v2(6))
      pv2(4) = dimag(v2(5))

      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      V1V2 =  v1(1)*v2(1) -  v1(2)*v2(2) -  v1(3)*v2(3) -  v1(4)*v2(4)
      K1V2 = pv1(1)*v2(1) - pv1(2)*v2(2) - pv1(3)*v2(3) - pv1(4)*v2(4)
      K2V1 = pv2(1)*v1(1) - pv2(2)*v1(2) - pv2(3)*v1(3) - pv2(4)*v1(4)
c     new
      K1V1 = pv1(1)*v1(1) - pv1(2)*v1(2) - pv1(3)*v1(3) - pv1(4)*v1(4)
      K2V2 = pv2(1)*v2(1) - pv2(2)*v2(2) - pv2(3)*v2(3) - pv2(4)*v2(4)
c     new

      F = pv1(1)*pv2(1) - pv1(2)*pv2(2) - pv1(3)*pv2(3) - pv1(4)*pv2(4)
      if ( vmass.ne.rZero ) then
         F = F + vmass**2
      end if

      TKK   = cZero
      TVV   = cZero
      TK1V2 = cZero
      TK2V1 = cZero

      do i = 1,4
         dum   = ft(i,i)*pv1(i)
         TKK   = TKK   + dum*pv2(i)
         TK1V2 = TK1V2 + dum*v2(i)
         dum   = ft(i,i)*v1(i)
         TVV   = TVV   + dum*v2(i)
         TK2V1 = TK2V1 + dum*pv2(i)
      end do

      TKK   = rTwo*TKK
      TVV   = rTwo*TVV
      TK1V2 = rTwo*TK1V2
      TK2V1 = rTwo*TK2V1

      TKK = TKK - T12*(pv1(1)*pv2(2) + pv1(2)*pv2(1))
     &          - T13*(pv1(1)*pv2(3) + pv1(3)*pv2(1))
     &          - T14*(pv1(1)*pv2(4) + pv1(4)*pv2(1))
     &          + T23*(pv1(2)*pv2(3) + pv1(3)*pv2(2))
     &          + T24*(pv1(2)*pv2(4) + pv1(4)*pv2(2))
     &          + T34*(pv1(3)*pv2(4) + pv1(4)*pv2(3))

      TK1V2 = TK1V2 - T12*(pv1(1)*v2(2) + pv1(2)*v2(1))
     &              - T13*(pv1(1)*v2(3) + pv1(3)*v2(1))
     &              - T14*(pv1(1)*v2(4) + pv1(4)*v2(1))
     &              + T23*(pv1(2)*v2(3) + pv1(3)*v2(2))
     &              + T24*(pv1(2)*v2(4) + pv1(4)*v2(2))
     &              + T34*(pv1(3)*v2(4) + pv1(4)*v2(3))

      TVV = TVV - T12*(v1(1)*v2(2) + v1(2)*v2(1))
     &          - T13*(v1(1)*v2(3) + v1(3)*v2(1))
     &          - T14*(v1(1)*v2(4) + v1(4)*v2(1))
     &          + T23*(v1(2)*v2(3) + v1(3)*v2(2))
     &          + T24*(v1(2)*v2(4) + v1(4)*v2(2))
     &          + T34*(v1(3)*v2(4) + v1(4)*v2(3))

      TK2V1 = TK2V1 - T12*(v1(1)*pv2(2) + v1(2)*pv2(1))
     &              - T13*(v1(1)*pv2(3) + v1(3)*pv2(1))
     &              - T14*(v1(1)*pv2(4) + v1(4)*pv2(1))
     &              + T23*(v1(2)*pv2(3) + v1(3)*pv2(2))
     &              + T24*(v1(2)*pv2(4) + v1(4)*pv2(2))
     &              + T34*(v1(3)*pv2(4) + v1(4)*pv2(3))

      vertex =  (ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4))*( K1V2*K2V1 - V1V2*F )
     &        + F*TVV + V1V2*TKK - K2V1*TK1V2 - K1V2*TK2V1

C      vertex = F*TVV + V1V2*TKK - K2V1*TK1V2 - K1V2*TK2V1

c     new, additonal gauge fixing term in Feyman gauge
      if ( vmass.eq.rZero ) then
         vertex = vertex 
     &+(ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4))*(K1V1*K1V2+K2V1*K2V2+K1V1*K2V2)
     &-K1V1*TK1V2-K2V2*TK2V1	   
      endif	
c     new    

      vertex = vertex * gt

      return
      end
      subroutine vvvkxx(wm,wp,tc,g, vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c two gauge bosons and a Kaluza-Klein tensor boson.
c
c input:
c       complex wm(6)          : vector               flow-in  V
c       complex wp(6)          : vector               flow-out V~
c       complex tc(6,4)        : tensor               KK mode T
c       complex g(1)           : coupling constant    -kappa/2
c       real    g(2)           : V boson mass          m_V
c
c output:
c       complex vertex         : amplitude            gamma(wm,wp,tc)
c     
      implicit none
      double complex wm(18), wp(18), tc(18), vertex,g(2)
      double precision vmass

      double complex T12, T13, T14, T23, T24, T34
      double complex V1V2, k1V2, k2V1
      double complex Tkk, TVV, Tk1V2, Tk2V1, dum
      double precision pwm(4), pwp(4), F

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )
c
      vmass = dreal(g(2))
      pwm(1) = dreal(wm(5))
      pwm(2) = dreal(wm(6))
      pwm(3) = dimag(wm(6))
      pwm(4) = dimag(wm(5))
      pwp(1) = dreal(wp(5))
      pwp(2) = dreal(wp(6))
      pwp(3) = dimag(wp(6))
      pwp(4) = dimag(wp(5))

      T12 = tc( 2) + tc( 5)
      T13 = tc( 3) + tc( 9)
      T14 = tc( 4) + tc(13)
      T23 = tc( 7) + tc(10)
      T24 = tc( 8) + tc(14)
      T34 = tc(12) + tc(15)

      V1V2 =  wm(1)*wp(1) -  wm(2)*wp(2) -  wm(3)*wp(3) -  wm(4)*wp(4)
      k1V2 = pwm(1)*wp(1) - pwm(2)*wp(2) - pwm(3)*wp(3) - pwm(4)*wp(4)
      k2V1 = pwp(1)*wm(1) - pwp(2)*wm(2) - pwp(3)*wm(3) - pwp(4)*wm(4)

      F = pwm(1)*pwp(1) - pwm(2)*pwp(2) - pwm(3)*pwp(3) - pwm(4)*pwp(4)
      if ( vmass.ne.rZero ) then
         F = F + vmass**2
      end if

      Tkk   = cZero
      TVV   = cZero
      Tk1V2 = cZero
      Tk2V1 = cZero

      do i = 1,4
         dum   = tc(i+4*(i-1))*pwm(i)
         Tkk   = Tkk   + dum*pwp(i)
         Tk1V2 = Tk1V2 + dum*wp(i)
         dum   = tc(i+4*(i-1))*wm(i)
         TVV   = TVV   + dum*wp(i)
         Tk2V1 = Tk2V1 + dum*pwp(i)
      end do

      Tkk   = rTwo*Tkk
      TVV   = rTwo*TVV
      Tk1V2 = rTwo*Tk1V2
      Tk2V1 = rTwo*Tk2V1

      Tkk = Tkk - T12*(pwm(1)*pwp(2) + pwm(2)*pwp(1))
     &          - T13*(pwm(1)*pwp(3) + pwm(3)*pwp(1))
     &          - T14*(pwm(1)*pwp(4) + pwm(4)*pwp(1))
     &          + T23*(pwm(2)*pwp(3) + pwm(3)*pwp(2))
     &          + T24*(pwm(2)*pwp(4) + pwm(4)*pwp(2))
     &          + T34*(pwm(3)*pwp(4) + pwm(4)*pwp(3))

      Tk1V2 = Tk1V2 - T12*(pwm(1)*wp(2) + pwm(2)*wp(1))
     &              - T13*(pwm(1)*wp(3) + pwm(3)*wp(1))
     &              - T14*(pwm(1)*wp(4) + pwm(4)*wp(1))
     &              + T23*(pwm(2)*wp(3) + pwm(3)*wp(2))
     &              + T24*(pwm(2)*wp(4) + pwm(4)*wp(2))
     &              + T34*(pwm(3)*wp(4) + pwm(4)*wp(3))

      TVV = TVV - T12*(wm(1)*wp(2) + wm(2)*wp(1))
     &          - T13*(wm(1)*wp(3) + wm(3)*wp(1))
     &          - T14*(wm(1)*wp(4) + wm(4)*wp(1))
     &          + T23*(wm(2)*wp(3) + wm(3)*wp(2))
     &          + T24*(wm(2)*wp(4) + wm(4)*wp(2))
     &          + T34*(wm(3)*wp(4) + wm(4)*wp(3))

      Tk2V1 = Tk2V1 - T12*(wm(1)*pwp(2) + wm(2)*pwp(1))
     &              - T13*(wm(1)*pwp(3) + wm(3)*pwp(1))
     &              - T14*(wm(1)*pwp(4) + wm(4)*pwp(1))
     &              + T23*(wm(2)*pwp(3) + wm(3)*pwp(2))
     &              + T24*(wm(2)*pwp(4) + wm(4)*pwp(2))
     &              + T34*(wm(3)*pwp(4) + wm(4)*pwp(3))

      vertex =  (tc(1)-tc(6)-tc(11)-tc(16))*( k1V2*k2V1 - V1V2*F )
     &        + F*TVV + V1V2*Tkk - k2V1*Tk1V2 - k1V2*Tk2V1

C      vertex = F*TVV + V1V2*Tkk - k2V1*Tk1V2 - k1V2*Tk2V1

      vertex = vertex * g(1)
c
      return
      end
      subroutine vvvsxx(ga,gb,gc,sc,g1,g2,vertex)
c
c- by RF - Mar. 2006
c
c This subroutine computes an amplitude of the coupling of three gauge bosons
c and a scalar particle
c
c input:
c       complex ga(6)          : first  incoming vector   (gluon)
c       complex gb(6)          : second incoming vector   (gluon)
c       complex gc(6)          : third  incoming vector   (gluon)
c       complex sc(3)          : incoming scalar particle (Higgs)
c       real    g1             : coupling constant        (QCD)
c       complex g2(2)          : coupling constant: gc(1) scalar
c                                                   gc(2) pseudo-scalar
c
c output:
c       complex vertex         : amplitude  
c

      implicit none

c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex ga(DIM),gb(DIM),gc(DIM),sc(DIM)

      double complex dvertx, vertex, vertex1, vertex2
      double complex vab, vbc, vca, v123, v124, v134, v234
      double complex pgagb, pgagc, pgbga, pgbgc, pgcga, pgcgb
      double precision pga(0:3),pgb(0:3),pgc(0:3),pabc(4)
      double precision g1
      double complex g2(2)

      pga(0) = dble( ga(5))
      pga(1) = dble( ga(6))
      pga(2) = dimag(ga(6))
      pga(3) = dimag(ga(5))

      pgb(0) = dble( gb(5))
      pgb(1) = dble( gb(6))
      pgb(2) = dimag(gb(6))
      pgb(3) = dimag(gb(5))

      pgc(0) = dble( gc(5))
      pgc(1) = dble( gc(6))
      pgc(2) = dimag(gc(6))
      pgc(3) = dimag(gc(5))

      vertex1 = (0D0,0D0)
      vertex2 = (0D0,0D0)

      if (g2(1).NE.(0D0,0D0)) then
      vab = ga(1)*gb(1)-ga(2)*gb(2)-ga(3)*gb(3)-ga(4)*gb(4)
      vbc = gb(1)*gc(1)-gb(2)*gc(2)-gb(3)*gc(3)-gb(4)*gc(4)
      vca = gc(1)*ga(1)-gc(2)*ga(2)-gc(3)*ga(3)-gc(4)*ga(4)

      pgagb = pga(0)*gb(1) - pga(1)*gb(2) - pga(2)*gb(3) - pga(3)*gb(4)
      pgagc = pga(0)*gc(1) - pga(1)*gc(2) - pga(2)*gc(3) - pga(3)*gc(4)
      pgbga = pgb(0)*ga(1) - pgb(1)*ga(2) - pgb(2)*ga(3) - pgb(3)*ga(4)
      pgbgc = pgb(0)*gc(1) - pgb(1)*gc(2) - pgb(2)*gc(3) - pgb(3)*gc(4)
      pgcga = pgc(0)*ga(1) - pgc(1)*ga(2) - pgc(2)*ga(3) - pgc(3)*ga(4)
      pgcgb = pgc(0)*gb(1) - pgc(1)*gb(2) - pgc(2)*gb(3) - pgc(3)*gb(4)

      dvertx = vab*(pgagc-pgbgc) + vbc*(pgbga-pgcga) + vca*(pgcgb-pgagb)
      vertex1= dvertx * g2(1)
      endif

      if (g2(2).NE.(0D0,0D0)) then
      pabc(1) = pga(0) + pgb(0) + pgc(0)
      pabc(2) = pga(1) + pgb(1) + pgc(1)
      pabc(3) = pga(2) + pgb(2) + pgc(2)
      pabc(4) = pga(3) + pgb(3) + pgc(3)

      v123 =   ga(1)*gb(2)*gc(3) - ga(1)*gb(3)*gc(2) - ga(2)*gb(1)*gc(3)
     &       + ga(2)*gb(3)*gc(1) + ga(3)*gb(1)*gc(2) - ga(3)*gb(2)*gc(1)
      v124 = - ga(1)*gb(2)*gc(4) + ga(1)*gb(4)*gc(2) + ga(2)*gb(1)*gc(4)
     &       - ga(2)*gb(4)*gc(1) - ga(4)*gb(1)*gc(2) + ga(4)*gb(2)*gc(1)
      v134 =   ga(1)*gb(3)*gc(4) - ga(1)*gb(4)*gc(3) - ga(3)*gb(1)*gc(4)
     &       + ga(3)*gb(4)*gc(1) + ga(4)*gb(1)*gc(3) - ga(4)*gb(3)*gc(1)
      v234 = - ga(2)*gb(3)*gc(4) + ga(2)*gb(4)*gc(3) + ga(3)*gb(2)*gc(4)
     &       - ga(3)*gb(4)*gc(2) - ga(4)*gb(2)*gc(3) + ga(4)*gb(3)*gc(2)


      vertex2= g2(2) * (  v123*pabc(4) + v124*pabc(3)
     &                  + v134*pabc(2) + v234*pabc(1) )
      endif

      vertex = g1*sc(1) * (vertex1 + vertex2)

      return
      end
      subroutine vvvtlx(ga,gb,gc,sc,g1,g2,vertex)
c
c- by RF - Feb. 2006
c
c This subroutine computes an amplitude of the coupling of three gauge bosons
c and a scalar particle
c
c input:
c       complex ga(6)          : first  incoming vector   (gluon)
c       complex gb(6)          : second incoming vector   (gluon)
c       complex gc(6)          : third  incoming vector   (gluon)
c       complex sc(3)          : incoming scalar particle (Higgs)
c       real    g1             : coupling constant        (QCD)
c       complex g2(2)          : coupling constant        (Higgs Effct. Thr.)
c
c output:
c       complex vertex         : amplitude  
c

      implicit none

c--   dimension of the current set to arbitrary length
c      INTEGER DIM
c      PARAMETER(DIM=18)
      include "dimension.inc"
      double complex ga(DIM),gb(DIM),gc(DIM),sc(DIM)
      double complex dvertx, vertex, vertex1, vertex2
      double complex vab, vbc, vca, v123, v124, v134, v234
      double complex pgagb, pgagc, pgbga, pgbgc, pgcga, pgcgb
      double precision pga(0:3),pgb(0:3),pgc(0:3),pabc(4)
      double precision g1
      double complex   g2(2)

      pga(0) = dble( ga(5))
      pga(1) = dble( ga(6))
      pga(2) = dimag(ga(6))
      pga(3) = dimag(ga(5))

      pgb(0) = dble( gb(5))
      pgb(1) = dble( gb(6))
      pgb(2) = dimag(gb(6))
      pgb(3) = dimag(gb(5))

      pgc(0) = dble( gc(5))
      pgc(1) = dble( gc(6))
      pgc(2) = dimag(gc(6))
      pgc(3) = dimag(gc(5))

      vertex1 = (0D0,0D0)
      vertex2 = (0D0,0D0)

      if (g2(1).NE.(0D0,0D0)) then
      vab = ga(1)*gb(1)-ga(2)*gb(2)-ga(3)*gb(3)-ga(4)*gb(4)
      vbc = gb(1)*gc(1)-gb(2)*gc(2)-gb(3)*gc(3)-gb(4)*gc(4)
      vca = gc(1)*ga(1)-gc(2)*ga(2)-gc(3)*ga(3)-gc(4)*ga(4)

      pgagb = pga(0)*gb(1) - pga(1)*gb(2) - pga(2)*gb(3) - pga(3)*gb(4)
      pgagc = pga(0)*gc(1) - pga(1)*gc(2) - pga(2)*gc(3) - pga(3)*gc(4)
      pgbga = pgb(0)*ga(1) - pgb(1)*ga(2) - pgb(2)*ga(3) - pgb(3)*ga(4)
      pgbgc = pgb(0)*gc(1) - pgb(1)*gc(2) - pgb(2)*gc(3) - pgb(3)*gc(4)
      pgcga = pgc(0)*ga(1) - pgc(1)*ga(2) - pgc(2)*ga(3) - pgc(3)*ga(4)
      pgcgb = pgc(0)*gb(1) - pgc(1)*gb(2) - pgc(2)*gb(3) - pgc(3)*gb(4)

      dvertx = vab*(pgagc-pgbgc) + vbc*(pgbga-pgcga) + vca*(pgcgb-pgagb)
      vertex1= dvertx * g2(1)
      endif

      if (g2(2).NE.(0D0,0D0)) then
      pabc(1) = pga(0) + pgb(0) + pgc(0)
      pabc(2) = pga(1) + pgb(1) + pgc(1)
      pabc(3) = pga(2) + pgb(2) + pgc(2)
      pabc(4) = pga(3) + pgb(3) + pgc(3)

      v123 =   ga(1)*gb(2)*gc(3) - ga(1)*gb(3)*gc(2) - ga(2)*gb(1)*gc(3)
     &       + ga(2)*gb(3)*gc(1) + ga(3)*gb(1)*gc(2) - ga(3)*gb(2)*gc(1)
      v124 = - ga(1)*gb(2)*gc(4) + ga(1)*gb(4)*gc(2) + ga(2)*gb(1)*gc(4)
     &       - ga(2)*gb(4)*gc(1) - ga(4)*gb(1)*gc(2) + ga(4)*gb(2)*gc(1)
      v134 =   ga(1)*gb(3)*gc(4) - ga(1)*gb(4)*gc(3) - ga(3)*gb(1)*gc(4)
     &       + ga(3)*gb(4)*gc(1) + ga(4)*gb(1)*gc(3) - ga(4)*gb(3)*gc(1)
      v234 = - ga(2)*gb(3)*gc(4) + ga(2)*gb(4)*gc(3) + ga(3)*gb(2)*gc(4)
     &       - ga(3)*gb(4)*gc(2) - ga(4)*gb(2)*gc(3) + ga(4)*gb(3)*gc(2)


      vertex2= g2(2) * (  v123*pabc(4) + v124*pabc(3)
     &                    + v134*pabc(2) + v234*pabc(1) )
      endif

      vertex = g1*sc(1) * (vertex1 + vertex2)

      return
      end
      subroutine vvvtxx(va,vb,vc,tc,gc,gt , vertex)
c
c This subroutine computes an amplitude of the four-point coupling of
c three gauge bosons and a tensor boson.
c
c input:
c       complex va(6)          : first  vector                        va
c       complex vb(6)          : second vector                        vb
c       complex vc(6)          : third  vector                        vc
c       complex tc(18)         : input  tensor                         T
c       real    gc             : coupling constant       gs (for gluons)
c       complex gt             : coupling constant         gtv=-1/Lambda
c
c output:
c       complex vertex         : amplitude             gamma(va,vb,vc,T)
c
c- by Q.Li - OCT. 2006
c     
      implicit none
      double complex va(6), vb(6), vc(6),  tc(18), gt, vertex
      double precision gc
 
      double complex ft(6,4)
      double complex T00, T12, T13, T14, T23, T24, T34
      double complex V1V2,V1V3,V2V3, K1V2, K1V3, K2V1
     &,K2V3,K3V1,K3V2

      double complex TV12,TV13,TV23,TKV1,TKV2,TKV3, dum
      double precision pva(4), pvb(4), pvc(4),p31(4),p23(4),p12(4)

      integer i, j

      double complex cZero
      double precision rZero, rTwo
      parameter( rZero = 0.0d0, rTwo = 2.0d0 )
      parameter( cZero = ( 0.0d0, 0.0d0 ) )

      
      ft(1,1) = tc(1)
      ft(1,2) = tc(2)
      ft(1,3) = tc(3)
      ft(1,4) = tc(4)
      ft(2,1) = tc(5)
      ft(2,2) = tc(6)
      ft(2,3) = tc(7)
      ft(2,4) = tc(8)
      ft(3,1) = tc(9)
      ft(3,2) = tc(10)
      ft(3,3) = tc(11)
      ft(3,4) = tc(12)
      ft(4,1) = tc(13)
      ft(4,2) = tc(14)
      ft(4,3) = tc(15)
      ft(4,4) = tc(16)
      ft(5,1) = tc(17)
      ft(6,1) = tc(18)

      pva(1) = dreal(va(5))
      pva(2) = dreal(va(6))
      pva(3) = dimag(va(6))
      pva(4) = dimag(va(5))

      pvb(1) = dreal(vb(5))
      pvb(2) = dreal(vb(6))
      pvb(3) = dimag(vb(6))
      pvb(4) = dimag(vb(5))

      pvc(1) = dreal(vc(5))
      pvc(2) = dreal(vc(6))
      pvc(3) = dimag(vc(6))
      pvc(4) = dimag(vc(5))

      p31(1) = pvc(1)-pva(1)
      p31(2) = pvc(2)-pva(2)
      p31(3) = pvc(3)-pva(3)
      p31(4) = pvc(4)-pva(4)
      
      p12(1) = pva(1)-pvb(1)
      p12(2) = pva(2)-pvb(2)
      p12(3) = pva(3)-pvb(3)
      p12(4) = pva(4)-pvb(4)
      
      p23(1) = pvb(1)-pvc(1)
      p23(2) = pvb(2)-pvc(2)
      p23(3) = pvb(3)-pvc(3)
      p23(4) = pvb(4)-pvc(4)
      
      T00 = ft(1,1)-ft(2,2)-ft(3,3)-ft(4,4)
      T12 = ft(1,2) + ft(2,1)
      T13 = ft(1,3) + ft(3,1)
      T14 = ft(1,4) + ft(4,1)
      T23 = ft(2,3) + ft(3,2)
      T24 = ft(2,4) + ft(4,2)
      T34 = ft(3,4) + ft(4,3)

      V1V2 =  va(1)*vb(1) -  va(2)*vb(2) -  va(3)*vb(3) -  va(4)*vb(4)
      V1V3 =  va(1)*vc(1) -  va(2)*vc(2) -  va(3)*vc(3) -  va(4)*vc(4)
      V2V3 =  vc(1)*vb(1) -  vc(2)*vb(2) -  vc(3)*vb(3) -  vc(4)*vb(4)
      K1V2 = pva(1)*vb(1) - pva(2)*vb(2) - pva(3)*vb(3) - pva(4)*vb(4)
      K1V3 = pva(1)*vc(1) - pva(2)*vc(2) - pva(3)*vc(3) - pva(4)*vc(4)
      K2V1 = pvb(1)*va(1) - pvb(2)*va(2) - pvb(3)*va(3) - pvb(4)*va(4)
      K2V3 = pvb(1)*vc(1) - pvb(2)*vc(2) - pvb(3)*vc(3) - pvb(4)*vc(4)
      K3V1 = pvc(1)*va(1) - pvc(2)*va(2) - pvc(3)*va(3) - pvc(4)*va(4)
      K3V2 = pvc(1)*vb(1) - pvc(2)*vb(2) - pvc(3)*vb(3) - pvc(4)*vb(4)


      TV12   = cZero
      TV13   = cZero
      TV23   = cZero
      TKV1   = cZero
      TKV2   = cZero
      TKV3   = cZero

      TV12 = rtwo*(ft(1,1)*va(1)*vb(1)+ft(2,2)*va(2)*vb(2)
     &+ft(3,3)*va(3)*vb(3)+ft(4,4)*va(4)*vb(4))

      TV13 = rtwo*(ft(1,1)*va(1)*vc(1)+ft(2,2)*va(2)*vc(2)
     &+ft(3,3)*va(3)*vc(3)+ft(4,4)*va(4)*vc(4))

      TV23 = rtwo*(ft(1,1)*vb(1)*vc(1)+ft(2,2)*vb(2)*vc(2)
     &+ft(3,3)*vb(3)*vc(3)+ft(4,4)*vb(4)*vc(4))

      TKV1 = rtwo*(ft(1,1)*p23(1)*va(1)+ft(2,2)*p23(2)*va(2)
     &+ft(3,3)*p23(3)*va(3)+ft(4,4)*p23(4)*va(4))

      TKV2 = rtwo*(ft(1,1)*p31(1)*vb(1)+ft(2,2)*p31(2)*vb(2)
     &+ft(3,3)*p31(3)*vb(3)+ft(4,4)*p31(4)*vb(4))
     
      TKV3 = rtwo*(ft(1,1)*p12(1)*vc(1)+ft(2,2)*p12(2)*vc(2)
     &+ft(3,3)*p12(3)*vc(3)+ft(4,4)*p12(4)*vc(4))	


      TV12 = TV12 - T12*(va(1)*vb(2) + va(2)*vb(1))
     &            - T13*(va(1)*vb(3) + va(3)*vb(1))
     &            - T14*(va(1)*vb(4) + va(4)*vb(1))
     &            + T23*(va(2)*vb(3) + va(3)*vb(2))
     &            + T24*(va(2)*vb(4) + va(4)*vb(2))
     &            + T34*(va(3)*vb(4) + va(4)*vb(3))

      TV13 = TV13 - T12*(va(1)*vc(2) + va(2)*vc(1))
     &            - T13*(va(1)*vc(3) + va(3)*vc(1))
     &            - T14*(va(1)*vc(4) + va(4)*vc(1))
     &            + T23*(va(2)*vc(3) + va(3)*vc(2))
     &            + T24*(va(2)*vc(4) + va(4)*vc(2))
     &            + T34*(va(3)*vc(4) + va(4)*vc(3))

      TV23 = TV23 - T12*(vb(1)*vc(2) + vb(2)*vc(1))
     &            - T13*(vb(1)*vc(3) + vb(3)*vc(1))
     &            - T14*(vb(1)*vc(4) + vb(4)*vc(1))
     &            + T23*(vb(2)*vc(3) + vb(3)*vc(2))
     &            + T24*(vb(2)*vc(4) + vb(4)*vc(2))
     &            + T34*(vb(3)*vc(4) + vb(4)*vc(3))


      TKV1 = TKV1 - T12*(p23(1)*va(2) + p23(2)*va(1))
     &            - T13*(p23(1)*va(3) + p23(3)*va(1))
     &            - T14*(p23(1)*va(4) + p23(4)*va(1))
     &            + T23*(p23(2)*va(3) + p23(3)*va(2))
     &            + T24*(p23(2)*va(4) + p23(4)*va(2))
     &            + T34*(p23(3)*va(4) + p23(4)*va(3))

      TKV2 = TKV2 - T12*(p31(1)*vb(2) + p31(2)*vb(1))
     &            - T13*(p31(1)*vb(3) + p31(3)*vb(1))
     &            - T14*(p31(1)*vb(4) + p31(4)*vb(1))
     &            + T23*(p31(2)*vb(3) + p31(3)*vb(2))
     &            + T24*(p31(2)*vb(4) + p31(4)*vb(2))
     &            + T34*(p31(3)*vb(4) + p31(4)*vb(3))

      TKV3 = TKV3 - T12*(p12(1)*vc(2) + p12(2)*vc(1))
     &            - T13*(p12(1)*vc(3) + p12(3)*vc(1))
     &            - T14*(p12(1)*vc(4) + p12(4)*vc(1))
     &            + T23*(p12(2)*vc(3) + p12(3)*vc(2))
     &            + T24*(p12(2)*vc(4) + p12(4)*vc(2))
     &            + T34*(p12(3)*vc(4) + p12(4)*vc(3))


      vertex = TKV3*V1V2-T00*K1V3*V1V2+T00*K2V3*V1V2+TKV2*V1V3
     &+TV23*(K2V1-K3V1)+TKV1*V2V3-T00*K2V1*V2V3+T00*K3V1*V2V3
     &-TV13*(K1V2-K3V2)+T00*K1V2*V1V3-T00*V1V3*K3V2+TV12*(K1V3-K2V3)

      vertex= -vertex * gc*gt

      return
      end
      subroutine vvvxxx(wm,wp,w3,g , vertex)
c
c This subroutine computes an amplitude of the three-point coupling of
c the gauge bosons.
c
c input:
c       complex wm(6)          : vector               flow-out W-
c       complex wp(6)          : vector               flow-out W+
c       complex w3(6)          : vector               j3 or A    or Z
c       real    g              : coupling constant    gw or gwwa or gwwz
c
c output:
c       complex vertex         : amplitude               gamma(wm,wp,w3)
c     
      implicit none
      double complex wm(6),wp(6),w3(6),vertex,
     &     xv1,xv2,xv3,v12,v23,v31,p12,p13,p21,p23,p31,p32
      double precision pwm(0:3),pwp(0:3),pw3(0:3),g

      double precision rZero, rTenth
      parameter( rZero = 0.0d0, rTenth = 0.1d0 )

c
      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pw3(0) = dble( w3(5))
      pw3(1) = dble( w3(6))
      pw3(2) = dimag(w3(6))
      pw3(3) = dimag(w3(5))


      v12 = wm(1)*wp(1)-wm(2)*wp(2)-wm(3)*wp(3)-wm(4)*wp(4)
      v23 = wp(1)*w3(1)-wp(2)*w3(2)-wp(3)*w3(3)-wp(4)*w3(4)
      v31 = w3(1)*wm(1)-w3(2)*wm(2)-w3(3)*wm(3)-w3(4)*wm(4)
      xv1 = rZero
      xv2 = rZero
      xv3 = rZero

      if ( abs(wm(1)).ne.rZero ) then
         if ( abs(wm(1)).ge.max(abs(wm(2)),abs(wm(3)),abs(wm(4)))
     &        *rTenth )
     &      xv1 = pwm(0)/wm(1)
      endif
      if ( abs(wp(1)).ne.rZero) then
         if ( abs(wp(1)).ge.max(abs(wp(2)),abs(wp(3)),abs(wp(4)))
     &        *rTenth )
     &      xv2 = pwp(0)/wp(1)
      endif
      if ( abs(w3(1)).ne.rZero) then
         if ( abs(w3(1)).ge.max(abs(w3(2)),abs(w3(3)),abs(w3(4)))
     &        *rTenth )
     &      xv3 = pw3(0)/w3(1)
      endif

      p12 = (pwm(0)-xv1*wm(1))*wp(1)-(pwm(1)-xv1*wm(2))*wp(2)
     &     -(pwm(2)-xv1*wm(3))*wp(3)-(pwm(3)-xv1*wm(4))*wp(4)
      p13 = (pwm(0)-xv1*wm(1))*w3(1)-(pwm(1)-xv1*wm(2))*w3(2)
     &     -(pwm(2)-xv1*wm(3))*w3(3)-(pwm(3)-xv1*wm(4))*w3(4)
      p21 = (pwp(0)-xv2*wp(1))*wm(1)-(pwp(1)-xv2*wp(2))*wm(2)
     &     -(pwp(2)-xv2*wp(3))*wm(3)-(pwp(3)-xv2*wp(4))*wm(4)
      p23 = (pwp(0)-xv2*wp(1))*w3(1)-(pwp(1)-xv2*wp(2))*w3(2)
     &     -(pwp(2)-xv2*wp(3))*w3(3)-(pwp(3)-xv2*wp(4))*w3(4)
      p31 = (pw3(0)-xv3*w3(1))*wm(1)-(pw3(1)-xv3*w3(2))*wm(2)
     &     -(pw3(2)-xv3*w3(3))*wm(3)-(pw3(3)-xv3*w3(4))*wm(4)
      p32 = (pw3(0)-xv3*w3(1))*wp(1)-(pw3(1)-xv3*w3(2))*wp(2)
     &     -(pw3(2)-xv3*w3(3))*wp(3)-(pw3(3)-xv3*w3(4))*wp(4)

      vertex = -(v12*(p13-p23)+v23*(p21-p31)+v31*(p32-p12))*g
c
      return
      end
      subroutine vxxxxx(p,vmass,nhel,nsv , vc)
c
c This subroutine computes a VECTOR wavefunction.
c
c input:
c       real    p(0:3)         : four-momentum of vector boson
c       real    vmass          : mass          of vector boson
c       integer nhel = -1, 0, 1: helicity      of vector boson
c                                (0 is forbidden if vmass=0.0)
c       integer nsv  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex vc(6)          : vector wavefunction       epsilon^mu(v)
c     
      implicit none
      double complex vc(6)
      double precision p(0:3),vmass,hel,hel0,pt,pt2,pp,pzpt,emp,sqh
      integer nhel,nsv,nsvahl

      double precision rZero, rHalf, rOne, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0 )
      
c

      sqh = dsqrt(rHalf)
      hel = dble(nhel)
      nsvahl = nsv*dabs(hel)
      pt2 = p(1)**2+p(2)**2
      pp = min(p(0),dsqrt(pt2+p(3)**2))
      pt = min(pp,dsqrt(pt2))

      vc(5) = dcmplx(p(0),p(3))*nsv
      vc(6) = dcmplx(p(1),p(2))*nsv


      if ( vmass.ne.rZero ) then

         hel0 = rOne-dabs(hel)

         if ( pp.eq.rZero ) then

            vc(1) = dcmplx( rZero )
            vc(2) = dcmplx(-hel*sqh )
            vc(3) = dcmplx( rZero , nsvahl*sqh )
            vc(4) = dcmplx( hel0 )

         else

            emp = p(0)/(vmass*pp)
            vc(1) = dcmplx( hel0*pp/vmass )
            vc(4) = dcmplx( hel0*p(3)*emp+hel*pt/pp*sqh )
            if ( pt.ne.rZero ) then
               pzpt = p(3)/(pp*pt)*sqh*hel
               vc(2) = dcmplx( hel0*p(1)*emp-p(1)*pzpt , 
     &                         -nsvahl*p(2)/pt*sqh       )
               vc(3) = dcmplx( hel0*p(2)*emp-p(2)*pzpt ,  
     &                          nsvahl*p(1)/pt*sqh       )
            else
               vc(2) = dcmplx( -hel*sqh )
               vc(3) = dcmplx( rZero , nsvahl*sign(sqh,p(3)) )
            endif

         endif

      else

         pp = p(0)
         pt = sqrt(p(1)**2+p(2)**2)
         vc(1) = dcmplx( rZero )
         vc(4) = dcmplx( hel*pt/pp*sqh )
         if ( pt.ne.rZero ) then
            pzpt = p(3)/(pp*pt)*sqh*hel
            vc(2) = dcmplx( -p(1)*pzpt , -nsv*p(2)/pt*sqh )
            vc(3) = dcmplx( -p(2)*pzpt ,  nsv*p(1)/pt*sqh )
         else
            vc(2) = dcmplx( -hel*sqh )
            vc(3) = dcmplx( rZero , nsv*sign(sqh,p(3)) )
         endif

      endif
c
      return
      end
      subroutine w3w3nx(wm,w31,wp,w32,g31,g32, vertex)
c
c This subroutine computes an amplitude of the four-point coupling of
c the W-, W+ and two W3/Z/A.
c If one sets wmass=0.0, then the gggg vertex is given
c (see sect 2.9.1 of the manual).
c
c input:
c       complex wm(0:3)        : flow-out W-                         wm
c       complex w31(0:3)       : first    W3/Z/A                     w31
c       complex wp(0:3)        : flow-out W+                         wp
c       complex w32(0:3)       : second   W3/Z/A                     w32
c       real    g31            : coupling of w31 with W-/W+
c       real    g32            : coupling of w32 with W-/W+
c                                                  (see the table below)
c       real    wmass          : mass  of W
c       real    wwidth         : width of W
c
c the possible sets of the inputs are as follows:
c   -------------------------------------------
c   |  wm  |  w31 |  wp  |  w32 |  g31 |  g32 |
c   -------------------------------------------
c   |  W-  |  W3  |  W+  |  W3  |  gw  |  gw  |
c   |  W-  |  W3  |  W+  |  Z   |  gw  | gwwz |
c   |  W-  |  W3  |  W+  |  A   |  gw  | gwwa |
c   |  W-  |  Z   |  W+  |  Z   | gwwz | gwwz |
c   |  W-  |  Z   |  W+  |  A   | gwwz | gwwa |
c   |  W-  |  A   |  W+  |  A   | gwwa | gwwa |
c   -------------------------------------------
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex vertex         : amplitude          gamma(wm,w31,wp,w32)
c     
      implicit none
      double complex wm(6),w31(6),wp(6),w32(6),vertex
      double complex dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),dvertx
      double complex v12,v13,v14,v23,v24,v34
      double precision pwm(0:3),pw31(0:3),pwp(0:3),pw32(0:3)
      double precision g31,g32,gtemp

      double precision rZero, rOne, rTwo
      parameter( rZero = 0.0d0, rOne = 1.0d0, rTwo = 2.0d0 )

c


c Benj's modif in order to have FR running
      gtemp=rZero
      if(g31.eq.rzero) gtemp=g32
      if(g32.eq.rzero) gtemp=g31
c End of Benj'S modif


      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pw31(0) = dble( w31(5))
      pw31(1) = dble( w31(6))
      pw31(2) = dimag(w31(6))
      pw31(3) = dimag(w31(5))
      pw32(0) = dble( w32(5))
      pw32(1) = dble( w32(6))
      pw32(2) = dimag(w32(6))
      pw32(3) = dimag(w32(5))


      dv1(0) = dcmplx(wm(1))
      dv1(1) = dcmplx(wm(2))
      dv1(2) = dcmplx(wm(3))
      dv1(3) = dcmplx(wm(4))
      dv2(0) = dcmplx(w31(1))
      dv2(1) = dcmplx(w31(2))
      dv2(2) = dcmplx(w31(3))
      dv2(3) = dcmplx(w31(4))
      dv3(0) = dcmplx(wp(1))
      dv3(1) = dcmplx(wp(2))
      dv3(2) = dcmplx(wp(3))
      dv3(3) = dcmplx(wp(4))
      dv4(0) = dcmplx(w32(1))
      dv4(1) = dcmplx(w32(2))
      dv4(2) = dcmplx(w32(3))
      dv4(3) = dcmplx(w32(4))

      v12 = dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
      v13 = dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
      v14 = dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
      v23 = dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
      v24 = dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
      v34 = dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)

      dvertx = v12*v34 + v14*v23 - rTwo*v13*v24

c     Neil edited this to allow 3 site coupling.
c     Now only g31 is important.  g32 does nothing.
c      vertex = dcmplx( dvertx ) * (g31*g32)
c      vertex = dcmplx( dvertx ) * (g31)
c     End Neil's edit
c

c Start of Benj's modif to have FR running
      vertex = dcmplx( dvertx ) * (gtemp)
c End of Benj's modif
      return
      end
      subroutine w3w3xx(wm,w31,wp,w32,g31,g32, vertex)
c
c This subroutine computes an amplitude of the four-point coupling of
c the W-, W+ and two W3/Z/A.
c If one sets wmass=0.0, then the gggg vertex is given
c (see sect 2.9.1 of the manual).
c
c input:
c       complex wm(0:3)        : flow-out W-                         wm
c       complex w31(0:3)       : first    W3/Z/A                     w31
c       complex wp(0:3)        : flow-out W+                         wp
c       complex w32(0:3)       : second   W3/Z/A                     w32
c       real    g31            : coupling of w31 with W-/W+
c       real    g32            : coupling of w32 with W-/W+
c                                                  (see the table below)
c       real    wmass          : mass  of W
c       real    wwidth         : width of W
c
c the possible sets of the inputs are as follows:
c   -------------------------------------------
c   |  wm  |  w31 |  wp  |  w32 |  g31 |  g32 |
c   -------------------------------------------
c   |  W-  |  W3  |  W+  |  W3  |  gw  |  gw  |
c   |  W-  |  W3  |  W+  |  Z   |  gw  | gwwz |
c   |  W-  |  W3  |  W+  |  A   |  gw  | gwwa |
c   |  W-  |  Z   |  W+  |  Z   | gwwz | gwwz |
c   |  W-  |  Z   |  W+  |  A   | gwwz | gwwa |
c   |  W-  |  A   |  W+  |  A   | gwwa | gwwa |
c   -------------------------------------------
c where all the bosons are defined by the flowing-OUT quantum number.
c
c output:
c       complex vertex         : amplitude          gamma(wm,w31,wp,w32)
c     
      implicit none
      double complex wm(6),w31(6),wp(6),w32(6),vertex
      double complex dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),dvertx
      double complex v12,v13,v14,v23,v24,v34
      double precision pwm(0:3),pw31(0:3),pwp(0:3),pw32(0:3)
      double precision g31,g32

      double precision rZero, rOne, rTwo
      parameter( rZero = 0.0d0, rOne = 1.0d0, rTwo = 2.0d0 )

c
      pwm(0) = dble( wm(5))
      pwm(1) = dble( wm(6))
      pwm(2) = dimag(wm(6))
      pwm(3) = dimag(wm(5))
      pwp(0) = dble( wp(5))
      pwp(1) = dble( wp(6))
      pwp(2) = dimag(wp(6))
      pwp(3) = dimag(wp(5))
      pw31(0) = dble( w31(5))
      pw31(1) = dble( w31(6))
      pw31(2) = dimag(w31(6))
      pw31(3) = dimag(w31(5))
      pw32(0) = dble( w32(5))
      pw32(1) = dble( w32(6))
      pw32(2) = dimag(w32(6))
      pw32(3) = dimag(w32(5))


      dv1(0) = dcmplx(wm(1))
      dv1(1) = dcmplx(wm(2))
      dv1(2) = dcmplx(wm(3))
      dv1(3) = dcmplx(wm(4))
      dv2(0) = dcmplx(w31(1))
      dv2(1) = dcmplx(w31(2))
      dv2(2) = dcmplx(w31(3))
      dv2(3) = dcmplx(w31(4))
      dv3(0) = dcmplx(wp(1))
      dv3(1) = dcmplx(wp(2))
      dv3(2) = dcmplx(wp(3))
      dv3(3) = dcmplx(wp(4))
      dv4(0) = dcmplx(w32(1))
      dv4(1) = dcmplx(w32(2))
      dv4(2) = dcmplx(w32(3))
      dv4(3) = dcmplx(w32(4))

      v12 = dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
      v13 = dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
      v14 = dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
      v23 = dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
      v24 = dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
      v34 = dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)

      dvertx = v12*v34 + v14*v23 - rTwo*v13*v24
      
      vertex = dcmplx( dvertx ) * (g31*g32)
c
      return
      end
      subroutine wwwwnx(wm1,wp1,wm2,wp2,gwwz,gwwa , vertex)

c
c This subroutine computes an amplitude of the four-point W-/W+ coupling.
c
c input:
c       complex wm1(0:3)       : first  flow-out W-                  wm1
c       complex wp1(0:3)       : first  flow-out W+                  wp1
c       complex wm2(0:3)       : second flow-out W-                  wm2
c       complex wp2(0:3)       : second flow-out W+                  wp2
c       real    gwwa           : coupling constant of W and A       gwwa
c       real    gwwz           : coupling constant of W and Z       gwwz
c       real    zmass          : mass  of Z
c       real    zwidth         : width of Z
c
c output:
c       complex vertex         : amplitude        gamma(wm1,wp1,wm2,wp2)
c     
      implicit none
      double complex wm1(6),wp1(6),wm2(6),wp2(6),vertex
      double complex dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),dvertx
      double complex v12,v13,v14,v23,v24,v34
      double precision pwm1(0:3),pwp1(0:3),pwm2(0:3),pwp2(0:3)
      double precision gwwa,gwwz,gtemp

      double precision rZero, rOne, rTwo
      parameter( rZero = 0.0d0, rOne = 1.0d0, rTwo = 2.0d0 )

c

c Benj's modif in order to have FR running
      gtemp=rZero
      if(gwwa.eq.rZero) gtemp=gwwz
      if(gwwz.eq.rZero) gtemp=gwwa
c End of Benj'S modif


      pwm1(0) = dble( wm1(5))
      pwm1(1) = dble( wm1(6))
      pwm1(2) = dimag(wm1(6))
      pwm1(3) = dimag(wm1(5))
      pwp1(0) = dble( wp1(5))
      pwp1(1) = dble( wp1(6))
      pwp1(2) = dimag(wp1(6))
      pwp1(3) = dimag(wp1(5))
      pwm2(0) = dble( wm2(5))
      pwm2(1) = dble( wm2(6))
      pwm2(2) = dimag(wm2(6))
      pwm2(3) = dimag(wm2(5))
      pwp2(0) = dble( wp2(5))
      pwp2(1) = dble( wp2(6))
      pwp2(2) = dimag(wp2(6))
      pwp2(3) = dimag(wp2(5))


      dv1(0) = dcmplx(wm1(1))
      dv1(1) = dcmplx(wm1(2))
      dv1(2) = dcmplx(wm1(3))
      dv1(3) = dcmplx(wm1(4))
      dv2(0) = dcmplx(wp1(1))
      dv2(1) = dcmplx(wp1(2))
      dv2(2) = dcmplx(wp1(3))
      dv2(3) = dcmplx(wp1(4))
      dv3(0) = dcmplx(wm2(1))
      dv3(1) = dcmplx(wm2(2))
      dv3(2) = dcmplx(wm2(3))
      dv3(3) = dcmplx(wm2(4))
      dv4(0) = dcmplx(wp2(1))
      dv4(1) = dcmplx(wp2(2))
      dv4(2) = dcmplx(wp2(3))
      dv4(3) = dcmplx(wp2(4))

      v12 = dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
      v13 = dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
      v14 = dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
      v23 = dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
      v24 = dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
      v34 = dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)

c      dvertx = (v12*v34 + v14*v23 - rTwo*v13*v24)*(gwwa**2+gwwz**2)
c     Neil edited this vertex to allow implementation of 3-site model.
c     Now, the coupling gwwa is the full coupling squared of this vertex.
c      dvertx = (v12*v34 + v14*v23 - rTwo*v13*v24)*(gwwa)
c     End Neil's edit

c Start of Benj'S modif
      dvertx = (v12*v34 + v14*v23 - rTwo*v13*v24)*(gtemp)
c End of Benj'S modif

c Start of Claude'S modif (Removed minus sign)
      vertex = dcmplx( dvertx )
c End of Claude'S modif
c
      return
      end
      subroutine wwwwxx(wm1,wp1,wm2,wp2,gwwa,gwwz , vertex)
c
c This subroutine computes an amplitude of the four-point W-/W+ coupling.
c
c input:
c       complex wm1(0:3)       : first  flow-out W-                  wm1
c       complex wp1(0:3)       : first  flow-out W+                  wp1
c       complex wm2(0:3)       : second flow-out W-                  wm2
c       complex wp2(0:3)       : second flow-out W+                  wp2
c       real    gwwa           : coupling constant of W and A       gwwa
c       real    gwwz           : coupling constant of W and Z       gwwz
c       real    zmass          : mass  of Z
c       real    zwidth         : width of Z
c
c output:
c       complex vertex         : amplitude        gamma(wm1,wp1,wm2,wp2)
c     
      implicit none
      double complex wm1(6),wp1(6),wm2(6),wp2(6),vertex
      double complex dv1(0:3),dv2(0:3),dv3(0:3),dv4(0:3),dvertx
      double complex v12,v13,v14,v23,v24,v34
      double precision pwm1(0:3),pwp1(0:3),pwm2(0:3),pwp2(0:3)
      double precision gwwa,gwwz

      double precision rZero, rOne, rTwo
      parameter( rZero = 0.0d0, rOne = 1.0d0, rTwo = 2.0d0 )

c
      pwm1(0) = dble( wm1(5))
      pwm1(1) = dble( wm1(6))
      pwm1(2) = dimag(wm1(6))
      pwm1(3) = dimag(wm1(5))
      pwp1(0) = dble( wp1(5))
      pwp1(1) = dble( wp1(6))
      pwp1(2) = dimag(wp1(6))
      pwp1(3) = dimag(wp1(5))
      pwm2(0) = dble( wm2(5))
      pwm2(1) = dble( wm2(6))
      pwm2(2) = dimag(wm2(6))
      pwm2(3) = dimag(wm2(5))
      pwp2(0) = dble( wp2(5))
      pwp2(1) = dble( wp2(6))
      pwp2(2) = dimag(wp2(6))
      pwp2(3) = dimag(wp2(5))


      dv1(0) = dcmplx(wm1(1))
      dv1(1) = dcmplx(wm1(2))
      dv1(2) = dcmplx(wm1(3))
      dv1(3) = dcmplx(wm1(4))
      dv2(0) = dcmplx(wp1(1))
      dv2(1) = dcmplx(wp1(2))
      dv2(2) = dcmplx(wp1(3))
      dv2(3) = dcmplx(wp1(4))
      dv3(0) = dcmplx(wm2(1))
      dv3(1) = dcmplx(wm2(2))
      dv3(2) = dcmplx(wm2(3))
      dv3(3) = dcmplx(wm2(4))
      dv4(0) = dcmplx(wp2(1))
      dv4(1) = dcmplx(wp2(2))
      dv4(2) = dcmplx(wp2(3))
      dv4(3) = dcmplx(wp2(4))

      v12 = dv1(0)*dv2(0)-dv1(1)*dv2(1)-dv1(2)*dv2(2)-dv1(3)*dv2(3)
      v13 = dv1(0)*dv3(0)-dv1(1)*dv3(1)-dv1(2)*dv3(2)-dv1(3)*dv3(3)
      v14 = dv1(0)*dv4(0)-dv1(1)*dv4(1)-dv1(2)*dv4(2)-dv1(3)*dv4(3)
      v23 = dv2(0)*dv3(0)-dv2(1)*dv3(1)-dv2(2)*dv3(2)-dv2(3)*dv3(3)
      v24 = dv2(0)*dv4(0)-dv2(1)*dv4(1)-dv2(2)*dv4(2)-dv2(3)*dv4(3)
      v34 = dv3(0)*dv4(0)-dv3(1)*dv4(1)-dv3(2)*dv4(2)-dv3(3)*dv4(3)

      dvertx = (v12*v34 + v14*v23 - rTwo*v13*v24)*(gwwa**2+gwwz**2)

      vertex = -dcmplx( dvertx )
c
      return
      end
