c -*-Fortran-*-
      integer nintervals
      parameter (nintervals=50)
      real * 8 xgrid(0:nintervals,ndiminteg),
     1   xgrid0(0:nintervals,ndiminteg),ymax(nintervals,ndiminteg),
     2   ymaxrat(nintervals,ndiminteg),xgridrm(0:nintervals,ndiminteg),
     3   xgrid0rm(0:nintervals,ndiminteg),ymaxrm(nintervals,ndiminteg),
     4   ymaxratrm(nintervals,ndiminteg),
     5   xmmm(0:nintervals,ndiminteg),xmmmrm(0:nintervals,ndiminteg),
     6   xint,xintrm,xacc(0:nintervals,ndiminteg),
     7   xaccrm(0:nintervals,ndiminteg)
      integer nhits(1:nintervals,ndiminteg),
     1     nhitsrm(1:nintervals,ndiminteg)
      integer ifold(ndiminteg),ifoldrm(ndiminteg)
      common/cgengrids/xgrid,xgrid0,xint,ymax,ymaxrat,xmmm,xgridrm,
     1     xgrid0rm,xintrm,ymaxrm,ymaxratrm,xmmmrm,xacc,xaccrm,
     1     ifold,ifoldrm,nhits,nhitsrm
      save /cgengrids/
