BUILT_SOURCES = done-all-links

AUTOMAKE_OPTIONS = -Wno-portability

DIRLINKS = $(top_srcdir)/Analysis \
	   $(top_srcdir)/Decay \
           $(top_srcdir)/Hadronization \
           $(top_srcdir)/MatrixElement \
           $(top_srcdir)/Models \
           $(top_srcdir)/PDF \
           $(top_srcdir)/PDT \
           $(top_srcdir)/Shower \
           $(top_srcdir)/DipoleShower \
           $(top_srcdir)/Exsample2 \
           $(top_srcdir)/UnderlyingEvent \
           $(top_srcdir)/Utilities

LOOPTOOLHEADERS = $(top_srcdir)/Looptools/include/clooptools.h

CLEANFILES = done-all-links 

done-all-links: $(DIRLINKS) $(LOOPTOOLHEADERS)
	mkdir -p Herwig++/Config Herwig++/Looptools
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) Herwig++
	$(LN_S) -f $(addprefix ../../, $(LOOPTOOLHEADERS)) Herwig++/Looptools
	touch done-all-links

install-data-local:
	find Herwig++ -follow \( -name '*.h' -or -name '*.icc' \
	-or -name '*.tcc' -or -name '*.fh' -or -name '*.xh' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/Herwig++

clean-local:
	rm -rf Herwig++
