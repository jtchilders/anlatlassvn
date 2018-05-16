# ROOT Stuff
ROOTC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)

#HEPMCPATH=/home/parton/tools/HepMC-2.06.08
HEPMCPATH=/users/hpcusers/svn/tools/HepMC-2.06.08
HEPMC_COPTS=-I $(HEPMCPATH)
HEPMC_LOPTS=$(HEPMCPATH)/src/.libs/*.o

#PYTHIAPATH=/home/parton/generators/pythia8180/unedited
PYTHIAPATH=/users/hpcusers/svn/generators/pythia/v8180/bgq/trunk
PYTHIA_COPTS=-I $(PYTHIAPATH)/include -I $(PYTHIAPATH)/examples
PYTHIA_LOPTS=$(PYTHIAPATH)/tmp/archive/*.o

#TCLAPPATH=/home/parton/tools/tclap-1.2.1
TCLAPPATH=/users/hpcusers/svn/tools/tclap-1.2.1
TCLAP_COPTS=-I $(TCLAPPATH)/include


COPTS=-c -g $(PYTHIA_COPTS) $(HEPMC_COPTS) $(TCLAP_COPTS)
LOPTS=$(PYTHIA_LOPTS) $(HEPMC_LOPTS)

GCC=g++

BINARY=runPythiaOnAlpgen

target: $(BINARY)

# Link binary executables
$(BINARY): $(BINARY).o
	@echo " "
	@echo "--------->>>>>> Linking" $<
	$(GCC) $(LOPTS) $< -o $@
	@echo "----->>> Done."


# Build binary Source
$(BINARY).o: $(BINARY).cpp
	@echo " "
	@echo "--------->>>>>> Compiling" $*.o
	$(GCC) -MD $(COPTS) $< -o $@
	@echo "----->>> Done."

clean:
	@echo " "
	@echo "--------->>>>>> Cleaning"
	rm -f $(BINARY).o $(BINARY)
	@echo "----->>> Done."

