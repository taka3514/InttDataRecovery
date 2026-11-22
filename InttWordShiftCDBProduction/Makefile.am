AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  InttWordShiftCDBProduction.h

lib_LTLIBRARIES = \
  libInttWordShiftCDBProduction.la

libInttWordShiftCDBProduction_la_SOURCES = \
  InttWordShiftCDBProduction.cc

libInttWordShiftCDBProduction_la_LIBADD = \
  -lphool \
  -lfun4all \
  -lintt \
  -ltrack_reco \
  -lcalo_io \
  -lg4dst \
  -lSubsysReco

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libInttWordShiftCDBProduction.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
