AMGDIR=../install
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec -lpsb_ext
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

LINKOPT=
EXEDIR=./runs

LINKOPT=
EXEDIR=./runs
DGEN3D=amg_d_pde3d_base_mod.o amg_d_pde3d_exp_mod.o amg_d_pde3d_gauss_mod.o amg_d_pde3d_box_mod.o 	


all: amg_d_pde3d

amg_d_pde3d: amg_d_pde3d.o amg_d_genpde_mod.o $(DGEN3D) data_input.o
	$(FLINK) $(LINKOPT) amg_d_pde3d.o  amg_d_genpde_mod.o $(DGEN3D) data_input.o -o amg_d_pde3d $(AMG_LIBS) -lpsb_cuda $(PSBGPULDLIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_d_pde3d $(EXEDIR)

amg_d_pde3d.o: data_input.o

amg_d_pde3d.o: amg_d_genpde_mod.o $(DGEN3D)

amg_d_pde3d.o: amg_d_pde3d.F90
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -DCUDA_MODE -c $< -o $@

clean:
	/bin/rm -f data_input.o *.o *$(.mod)\
        $(EXEDIR)/amg_d_pde3d  
        
# Many suffixes are covered by implicit rules, but you may need to write custom rules
# such as these if you use suffixes that do not have implicit rules.
# https://www.gnu.org/software/make/manual/html_node/Catalogue-of-Rules.html#Catalogue-of-Rules
% : %.F90
	$(LINK.F) $(FINCLUDES) -o $@ $^ $(LDLIBS)
%.o: %.F90
	$(COMPILE.F) $(FDEFINES) $(FINCLUDES) $(OUTPUT_OPTION) $<
% : %.cxx
	$(LINK.cc) -o $@ $^ $(LDLIBS)
%.o: %.cxx
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
%.o : %.cu
	$(CUDAC) -c $(CPPFLAGS) $(CUDAC_FLAGS) $(CUDA_INCLUDE) -o $@ $<
