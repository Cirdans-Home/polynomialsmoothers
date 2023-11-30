anisopsblas:
	ln -s /home/durastante/polynomialsmoothers/anisotropy/runs/anisopsblas anisopsblas

# V-Cycle with VBM Aggregation 4 Iteration of l1-jacobi-chebyshev 4 and 30 iteration of l1-Jacobi as coarse solver
genscript_cheby4_soc1:
	./gen_script.sh run-cheby4.sh
run_cheby4_soc1: genscript_cheby4_soc1
	sbatch run-cheby4-1.sh
	sbatch run-cheby4-2.sh
	sbatch run-cheby4-4.sh
	sbatch run-cheby4-8.sh
	sbatch run-cheby4-16.sh
	sbatch run-cheby4-32.sh
	sbatch run-cheby4-64.sh
	sbatch run-cheby4-128.sh
	sbatch run-cheby4-216.sh

# V-Cycle with VBM Aggregation 4 Iteration of l1-jacobi-Optimized chebyshev 4 and 30 iteration of l1-Jacobi as coarse solver
genscript_optcheby4_soc1:
	./gen_script.sh run-optcheby4.sh
run_optcheby4_soc1: genscript_optcheby4_soc1
	sbatch run-optcheby4-1.sh
	sbatch run-optcheby4-2.sh
	sbatch run-optcheby4-4.sh
	sbatch run-optcheby4-8.sh
	sbatch run-optcheby4-16.sh
	sbatch run-optcheby4-32.sh	
	sbatch run-optcheby4-64.sh
	sbatch run-optcheby4-128.sh
	sbatch run-optcheby4-216.sh

clean:
	rm anisopsblas
	rm run-cheby4-*
	rm run-optcheby4-*
	rm *.out
