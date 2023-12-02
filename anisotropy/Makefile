# -*- mode: makefile -*-

AMGDIR=../install
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec -lpsb_ext
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

LINKOPT=
EXEDIR=./runs

psctoolkit: anisotropy.o data_input.o
	$(FLINK) $(LINKOPT) anisotropy.o data_input.o -o anisopsblas $(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv anisopsblas $(EXEDIR)

psctoolkitgpu: anisotropy_gpu.o data_input.o
	$(FLINK) $(LINKOPT) anisotropy_gpu.o data_input.o -o anisopsblascuda $(AMG_LIBS) -lpsb_cuda $(PSBGPULDLIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv anisopsblascuda $(EXEDIR)

anisotropy.o: data_input.o

anisotropy_gpu.o: data_input.o

anisotropy_gpu.o: anisotropy_gpu.F90
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -DCUDA_MODE -c $< -o $@

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf runs/anisopsblas
	rm -rf runs/anisopsblascuda

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
