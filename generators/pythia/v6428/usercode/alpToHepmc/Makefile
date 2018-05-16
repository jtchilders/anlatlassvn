

LHAPDF_PATH=/users/hpcusers/svn/tools/lhapdf-5.9.1/local
PYTHIA_PATH=/users/hpcusers/svn/generators/pythia/v6428/bgq/trunk
ALPGEN_PATH=/users/hpcusers/svn/generators/alpgen/v214/unedited

HEPMCPATH=/users/hpcusers/svn/tools/HepMC-2.06.08
HEPMC_COPTS=-I $(HEPMCPATH)
HEPMC_LOPTS=$(HEPMCPATH)/src/.libs/*.o $(HEPMCPATH)/fio/.libs/*.o

TCLAPPATH=/users/hpcusers/svn/tools/tclap-1.2.1
TCLAP_COPTS=-I $(TCLAPPATH)/include

GCC=g++
CFLAGS=-c -fPIC $(HEPMC_COPTS) -g $(TCLAP_COPTS) -Ilibfgz-0.2
LFLAGS=$(HEPMC_LOPTS) -L$(LHAPDF_PATH)/lib -lLHAPDF -lgfortran -L./libfgz-0.2 -lfgz -lz
GFF=gfortran
GFFCFLAGS=-fPIC -c -g

EXE=alpToHepmc

target: ${EXE} setup.sh

${EXE}: ${EXE}.o alpsho.o atopyt.o pythia-6.4.28.o libfgz-0.2/libfgz.a
	${GCC} -o $@ $^ ${LFLAGS}

${EXE}.o: ${EXE}.cpp
	${GCC} ${CFLAGS} $<

alpsho.o: alpsho.f alpsho.inc
	${GFF} $(GFFCFLAGS) $<

atopyt.o: atopyt.f alpsho.inc
	${GFF} ${GFFCFLAGS} $<

pythia-6.4.28.o: ${PYTHIA_PATH}/pythia-6.4.28.f
	${GFF} -fno-align-commons -cpp -DUSING_ALPGEN $(GFFCFLAGS) $<

setup.sh:
	@echo 'here'
	@echo "export LD_LIBRARY_PATH=$${LD_LIBRARY_PATH}:$(HEPMC_PATH)/lib:$(LHAPDF_PATH)/lib" > setup.sh

libfgz-0.2/libfgz.a: libfgz-0.2/*.F 
	make -C libfgz-0.2/

clean:
	rm -f *.o ${EXE}
