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
	
# V-Cycle with VBM Aggregation 4 Iteration of l1-jacobi-Optimized chebyshev 4 and 30 iteration of l1-Jacobi as coarse solver
genscript_optcheby1_soc1:
	./gen_script.sh run-optcheby1.sh
run_optcheby1_soc1: genscript_optcheby1_soc1
	sbatch run-optcheby1-1.sh
	sbatch run-optcheby1-2.sh
	sbatch run-optcheby1-4.sh
	sbatch run-optcheby1-8.sh
	sbatch run-optcheby1-16.sh
	sbatch run-optcheby1-32.sh      
	sbatch run-optcheby1-64.sh
	sbatch run-optcheby1-128.sh
	
# V-Cycle with Matching Aggregation 4 Iteration of l1-jacobi-chebyshev 4 and 30 iteration of l1-Jacobi as coarse solver
genscript_cheby4_match:
	./gen_script.sh run-cheby4-match.sh
run_cheby4_match: genscript_cheby4_match
	sbatch run-cheby4-match-1.sh
	sbatch run-cheby4-match-2.sh
	sbatch run-cheby4-match-4.sh
	sbatch run-cheby4-match-8.sh
	sbatch run-cheby4-match-16.sh
	sbatch run-cheby4-match-32.sh
	sbatch run-cheby4-match-64.sh
	sbatch run-cheby4-match-128.sh
        
# V-Cycle with Matching Aggregation 4 Iteration of l1-jacobi-Optimized chebyshev 4 and 30 iteration of l1-Jacobi as coarse solver
genscript_optcheby4_match:
	./gen_script.sh run-optcheby4-match.sh
run_optcheby4_match: genscript_optcheby4_match
	sbatch run-optcheby4-match-1.sh
	sbatch run-optcheby4-match-2.sh
	sbatch run-optcheby4-match-4.sh
	sbatch run-optcheby4-match-8.sh
	sbatch run-optcheby4-match-16.sh
	sbatch run-optcheby4-match-32.sh      
	sbatch run-optcheby4-match-64.sh
	sbatch run-optcheby4-match-128.sh
	
# V-Cycle with Matching Aggregation 4 Iteration of l1-jacobi-Optimized chebyshev 4 and 30 iteration of l1-Jacobi as coarse solver
genscript_optcheby1_match:
	./gen_script.sh run-optcheby1-match.sh
run_optcheby1_match: genscript_optcheby1_match
	sbatch run-optcheby1-match-1.sh
	sbatch run-optcheby1-match-2.sh
	sbatch run-optcheby1-match-4.sh
	sbatch run-optcheby1-match-8.sh
	sbatch run-optcheby1-match-16.sh
	sbatch run-optcheby1-match-32.sh      
	sbatch run-optcheby1-match-64.sh
	sbatch run-optcheby1-match-128.sh

clean:
	rm anisopsblas
	rm *.out
