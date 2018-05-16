c -*-Fortran-*-
      real * 8 xgrid(0:50,ndiminteg),ymax(50,ndiminteg),
     1     ymaxrat(50,ndiminteg),xgridrm(0:50,ndiminteg),
     2     ymaxrm(50,ndiminteg),ymaxratrm(50,ndiminteg),
     3     xmmm(0:50,ndiminteg),xmmmrm(0:50,ndiminteg)
      integer ifold(ndiminteg),ifoldrm(ndiminteg)
      common/cgengrids/xgrid,ymax,ymaxrat,xmmm,xgridrm,ymaxrm,ymaxratrm,
     #                xmmmrm,ifold,ifoldrm
