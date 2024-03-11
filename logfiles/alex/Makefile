anisopsblas:
	ln -s /home/hpc/ihpc/ihpc100h/polynomialsmoothers/anisotropy/runs/anisopsblas anisopsblas
	ln -s /home/hpc/ihpc/ihpc100h/polynomialsmoothers/anisotropy/runs/anisopsblascuda anisopsblascuda

genmatch:
	./gen_script.sh run-match-gpu.sh 4

runmatch:
	sbatch run-match-gpu-1.sh
	sbatch run-match-gpu-2.sh
	sbatch run-match-gpu-4.sh
	sbatch run-match-gpu-8.sh
	sbatch run-match-gpu-16.sh
	sbatch run-match-gpu-32.sh
	sbatch run-match-gpu-64.sh

gensoc:
	./gen_script.sh run-soc1-gpu.sh 4
	
gensocmpi:
	./gen_script_mpi.sh run-soc1-mpi.sh

runsoc:
	sbatch run-soc1-gpu-1.sh
	sbatch run-soc1-gpu-2.sh
	sbatch run-soc1-gpu-4.sh
	sbatch run-soc1-gpu-8.sh
	sbatch run-soc1-gpu-16.sh
	sbatch run-soc1-gpu-32.sh
	sbatch run-soc1-gpu-64.sh

cleanscript:
	rm run-match-gpu-*
	rm run-soc1-gpu-*

clean:
	rm anisopsblas
	rm anisopsblascuda
	rm *.out
