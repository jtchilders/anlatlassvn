#/bin/sh
# add a warning at the beginning of the file
echo "c     !: This should be as the sigreal.f file in the main directory. The only
c     !: difference is the presence here of the genrad variable and the
c     !: corresponding common block." > tmp.f

# copy sigreal.f from common files
cat ../sigreal.f >> tmp.f

# add genrad variable and correpsponding common blocks in proper places
sed '/subroutine sigreal_rad(/,/subroutine sigreal_btl(/s/call realgr(/genrad=.true. !:\n\t\t call realgr(/g' tmp.f > tmp2.f; mv tmp2.f tmp.f

sed '/subroutine sigreal_rad(/,/include/s/implicit none/implicit none\n      logical genrad !:\
      common\/cgenrad\/genrad !:/g' tmp.f > tmp2.f; mv tmp2.f tmp.f

sed '/subroutine sigreal_btl0(/,/subroutine fillmomenta(/s/call realgr(/genrad=.false. !:\n\t\t call realgr(/g' tmp.f > tmp2.f; mv tmp2.f tmp.f

sed '/subroutine sigreal_btl0(/,/include/s/implicit none/implicit none\n      logical genrad !:\
      common\/cgenrad\/genrad !:/g' tmp.f > tmp2.f; mv tmp2.f tmp.f

# rename file
mv tmp.f sigreal.f

# delete (subsitute with null) every occurrence of abc 
# between lines containing aaa and xxx
# 2 positional expression separated by comma
# sed '/aaa/,/xxx/s/abc//g' file > file.mod
