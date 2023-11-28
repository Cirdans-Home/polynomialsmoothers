# -*- mode: makefile -*-

AMGDIR=../install
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

LINKOPT=
EXEDIR=./runs

psctoolkit: anisotropy.o data_input.o
	$(FLINK) $(LINKOPT) anisotropy.o data_input.o -o anisopsblas $(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv anisopsblas $(EXEDIR)

anisotropy.o: data_input.o

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf runs/anisopsblas

# Many suffixes are covered by implicit rules, but you may need to write custom rules
# such as these if you use suffixes that do not have implicit rules.
# https://www.gnu.org/software/make/manual/html_node/Catalogue-of-Rules.html#Catalogue-of-Rules
% : %.F90
	$(LINK.F) $(FINCLUDES) -o $@ $^ $(LDLIBS)
%.o: %.F90
	$(COMPILE.F) $(FINCLUDES) $(OUTPUT_OPTION) $<
% : %.cxx
	$(LINK.cc) -o $@ $^ $(LDLIBS)
%.o: %.cxx
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
%.o : %.cu
	$(CUDAC) -c $(CPPFLAGS) $(CUDAC_FLAGS) $(CUDA_INCLUDE) -o $@ $<
