

HEPMC_PATH=/users/hpcusers/svn/tools/HepMC-2.06.08/local
LHAPDF_PATH=/users/hpcusers/svn/tools/lhapdf-5.9.1/local
PYTHIA_PATH=/users/hpcusers/svn/generators/pythia/v6428/bgq/trunk

GCC=g++
CFLAGS=-c -fPIC -I$(HEPMC_PATH)/include -g
LFLAGS=-L$(HEPMC_PATH)/lib -lHepMC -lHepMCfio -L$(LHAPDF_PATH)/lib -lLHAPDF -lgfortran
GFF=gfortran
GFFCFLAGS=-fPIC -c -g

EXE=pythiaToHepmc

target: ${EXE} setup.sh

${EXE}: ${EXE}.o pythia-6.4.28.o
	${GCC} -o $@ $^ ${LFLAGS} ${LFLAGS}

${EXE}.o: ${EXE}.cpp
	${GCC} ${CFLAGS} $<

pythia-6.4.28.o: ${PYTHIA_PATH}/pythia-6.4.28.f
	${GFF} -cpp $(GFFCFLAGS) $<

setup.sh:
	@echo 'here'
	@echo "export LD_LIBRARY_PATH=$${LD_LIBRARY_PATH}:$(HEPMC_PATH)/lib:$(LHAPDF_PATH)/lib" > setup.sh


clean:
	rm -f *.o ${EXE}
