
REALVERSION=$(shell tools/getversion.sh)
VERSION=$(shell tools/rversion.sh)

all: splicing_$(VERSION).tar.gz

CFILES=$(wildcard splicing/src/*.[ch]) $(wildcard splicing/src/*.pmt)
RFILES=$(wildcard splicing/R/*.R)
MANFILES=$(wildcard splicing/man/*.Rd)

splicing/DESCRIPTION: tools/DESCRIPTION
	sed 's/^Version: .*$$/Version: '$(VERSION)'/' $<     | \
	sed 's/^Date: .*$$/Date: '`date "+%Y-%m-%d"`'/' > $@

splicing/src/splicing_version.h: tools/splicing_version.h
	sed 's/@VERSION@/'$(REALVERSION)'/g' $< >$@
	touch splicing/src/version.c

.PHONY: splicing/DESCRIPTION splicing/src/splicing_version.h

splicing_$(VERSION).tar.gz: splicing/DESCRIPTION splicing/NAMESPACE \
		$(CFILES) $(RFILES) $(MANFILES) splicing/configure.in \
		splicing/src/splicing_version.h
	rm -f splicing/src/Makevars
	cd splicing && autoconf && autoheader
	R CMD build splicing
