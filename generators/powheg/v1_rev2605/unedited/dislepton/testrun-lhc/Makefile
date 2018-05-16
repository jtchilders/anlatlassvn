# Quick start guide: enter
#    make -j4
# to compile the executables and run them for the example.

# In case of problems and for details on the build settings, please see the
# ../Makefile in the project directory.

pythia: nlo
	@echo "\n===> building main-PYTHIA-lhef"
	cd .. && $(MAKE) main-PYTHIA-lhef
	@echo "\n===> running main-PYTHIA-lhef"
	../main-PYTHIA-lhef

nlo: clean
	@echo "\n===> building pwhg_main"
	cd .. && $(MAKE) pwhg_main
	@echo "\n===> running pwhg_main"
	../pwhg_main

clean:
	@echo "\n===> cleaning up old results"
	rm -rf FlavRegList pwgbtlgrid.top pwgstat.dat pwgxgrid.dat pwhg_checklimits\
	       pwgevents.lhe pwgborngrid.top pwgcounters.dat pwggrid.dat \
	       pwghistnorms.top pwgNLO.top pwgpwhgalone-output.top pwgubound.dat \
	       LHEF_analysis.top POWHEG+PYTHIA-output.top \
	       pwgPOWHEG+PYTHIA-output.top \
	       *~

	       