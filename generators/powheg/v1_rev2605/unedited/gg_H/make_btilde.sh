#/bin/sh
cat ../btilde.f | perl -pe 's/         call allborn/c set the mt dependency in Born cross section\
         call setbornmassdep\
c evaluate Born cross section with the chosen top mass\
         call allborn/;
s/         call btildeborn\(resborn\)/         call btildeborn\(resborn\)\
c re-evaluate Born amplitudes in large top mass limit, if needed\
         call setbornmass2inf/;
s/     www=www0\*hc2/     real \*8 finitemtcorr\
      external finitemtcorr\
      www=www0\*hc2/;
s/         call gen_born_phsp\(xborn\)/         call gen_born_phsp\(xborn\)\
         www=www*finitemtcorr()\
         wwwtot=wwwtot*finitemtcorr()/'>  btilde_ggH.f

mv btilde_ggH.f tmpfile
echo "c###################################################### " >btilde_ggH.f  
echo "c# THIS IS AN AUTOMATICALLY GENERATED FILE. DO NOT EDIT!" >>btilde_ggH.f
echo "c###################################################### " >>btilde_ggH.f
cat tmpfile >> btilde_ggH.f
rm tmpfile