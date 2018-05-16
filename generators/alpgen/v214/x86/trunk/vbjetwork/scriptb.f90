program script
!
  implicit none
!
  integer nh,nw,nz,nj
  character*6 t(2)
  character*12 v(2)
  character*28 cc
  character*10 dif_fil
  integer j1,j2
!
  t=(/'test   ','test90 '/); v=(/'vbjetsgen   ','vbjetsgen90V'/); dif_fil='diffeven  '
  write(1,9)'echo start>',dif_fil
!
  nh=1;nz=2;nw=0;nj=1
!
  do nh=0,8
    do nw=0,8-nh
      do nz=0,8-nw-nh
        do nj=0,8-nw-nh-nz
          if(nw+nz.eq.0.or.nw+nz+nh+nj.ne.1) cycle
          if(nj.gt.3) cycle
          do j1=2,2
            do j2=1,1
              write(1,1)v(j1),' <<fine'
              write(1,2)j2
              write(1,3)t(j1)
              write(1,2)1
              write(1,4)'1000d0 5'
              write(1,5)nz,nw,nh,nj
              write(1,6)'100 '
              if(nj.ne.0) write(1,7)'20 2.5 0.7'
              write(1,6)'1 1.'
              write(1,2)0
              write(1,4)'100    1'
              write(1,4)'1000    '
              write(1,2)0
              write(1,6)'fine'
            enddo
          enddo
          write(1,8)'echo "','nz,nw,nh,nj',nz,nw,nh,nj,' ">>',dif_fil
          write(1,10)'diff test.stat test90.stat >>',dif_fil
          write(1,11)'diff test.wgt test90.wgt >>',dif_fil
          write(1,11)'diff test.unw test90.unw >>',dif_fil
        enddo
      enddo
    enddo
  enddo
!
 1    format(a12,a7) 
 2    format(i2)
 3    format(a6)
 7    format(a10)
 4    format(a8)
 6    format(a4)
 5    format(4(i2,2x))
 8    format(a6,a12,4(2x,i2),a4,a10)
 9    format(a11,a10)
 10   format(a29,a10)
 11   format(a27,a10)
!
end program script
