
all: splicing_0.1.tar.gz

CFILES=$(wildcard splicing/src/*.[ch]) $(wildcard splicing/src/*.pmt)
RFILES=$(wildcard splicing/R/*.R)
MANFILES=$(wildcard splicing/man/*.Rd)

splicing_0.1.tar.gz: splicing/DESCRIPTION splicing/NAMESPACE \
		$(CFILES) $(RFILES) $(MANFILES) splicing/configure.in
	cd splicing && autoconf && autoheader
	R CMD build splicing
