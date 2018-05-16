#/usr/bin/env bash
#set -o xtrace
PDF_URL=http://www.hepforge.org/archive/lhapdf/pdfsets/5.9.1
COPTS_OPTIONAL=-fPIC
LOPTS_OPTIONAL=
INSTALL_PATH=$PWD/local
if [[ "$1" == "debug" ]];
then
   echo "BlueGeneQ Make: Setting debug compiler option"
   COPTS_OPTIONAL=-g
fi
if [[ "$1" == "profile" ]];
then
   echo "BlueGeneQ Make: Setting profile compiler option"
   COPTS_OPTIONAL=-pg
   LOPTS_OPTIONAL=-pg
fi
echo $COPTS_OPTIONAL
echo $LOPTS_OPTIONAL
mkdir -p $INSTALL_PATH
./configure --prefix=$INSTALL_PATH CFLAGS="$COPPS_OPTIONAL" CXXFLAGS="$COPTS_OPTIONAL" LDFLAGS="$LOPTS_OPTIONAL" --enable-shared=no
make
make install
mkdir -p $INSTALL_PATH/share/lhapdf/PDFsets
./bin/lhapdf-getdata --repo=$PDF_URL --dest=$INSTALL_PATH/share/lhapdf/PDFsets cteq6ll.LHpdf
