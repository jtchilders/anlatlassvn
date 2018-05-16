#/bin/sh
cat ../sigremnants.f | perl -pe 's/external valid_emitter/external valid_emitter\n\
      real \*8 finitemtcorr\n\
      external finitemtcorr\n/;
s/hc2/hc2\*\n     # finitemtcorr\(\)/gi' >  sigremnants_ggH.f

mv sigremnants_ggH.f tmpfile
echo "c###################################################### " >sigremnants_ggH.f  
echo "c# THIS IS AN AUTOMATICALLY GENERATED FILE. DO NOT EDIT!" >>sigremnants_ggH.f
echo "c###################################################### " >>sigremnants_ggH.f
cat tmpfile >> sigremnants_ggH.f
rm tmpfile