# ROOT Stuff
ROOTC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)

COPTS_OPTIONAL=
LOPTS_OPTIONAL=

HEPMCPATH=/users/hpcusers/svn/tools/HepMC-2.06.08
HEPMC_COPTS=-I $(HEPMCPATH)
HEPMC_LOPTS=-L$(HEPMCPATH)/local/lib -lHepMC -lHepMCfio

PYTHIAPATH=/users/hpcusers/svn/generators/pythia/v8180/unedited
PYTHIA_COPTS=-I $(PYTHIAPATH)/include -I $(PYTHIAPATH)/examples 
PYTHIA_LOPTS=-L$(PYTHIAPATH)/lib/archive -lpythia8 -llhapdfdummy -lpythia8tohepmc

TCLAPPATH=/users/hpcusers/svn/tools/tclap-1.2.1
TCLAP_COPTS=-I $(TCLAPPATH)/include

COPTS=-c -g $(COPTS_OPTIONAL) $(PYTHIA_COPTS) $(TCLAP_COPTS) $(ROOTC)
LOPTS=$(PYTHIA_LOPTS) $(LOPTS_OPTIONAL) $(ROOTLIB)

GCC=g++

BINARY=plotUnw

target: $(BINARY)

# Link binary executables
$(BINARY): $(BINARY).o
	@echo " "
	@echo "--------->>>>>> Linking" $<
	$(GCC) $(LOPTS) $< $(LOPTS) -o $@
	@echo "----->>> Done."


# Build binary Source
$(BINARY).o: $(BINARY).cpp GeneratorInput.h
	@echo " "
	@echo "--------->>>>>> Compiling" $*.o
	$(GCC) -MD $(COPTS) $< -o $@
	@echo "----->>> Done."

clean:
	@echo " "
	@echo "--------->>>>>> Cleaning"
	rm -f $(BINARY).o $(BINARY)
	@echo "----->>> Done."

