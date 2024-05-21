anisopsblas:
	ln -s /home/durastante/polynomialsmoothers/anisotropy/runs/anisopsblascuda anisopsblascuda
laplacian:
	ln -s /home/durastante/polynomialsmoothers/3dlaplacian/runs/amg_d_pde3d 3dlaplacian
	ln -s /home/durastante/polynomialsmoothers/3dlaplacian/runs/amg_d_pde3d_multi 3dlaplacian_multi

runsoc1:
	sbatch run-soc1-1.sh
	sbatch run-soc1-2.sh
	sbatch run-soc1-4.sh
	sbatch run-soc1-8.sh
	sbatch run-soc1-16.sh

runmatch:
	sbatch run-match-1.sh
	sbatch run-match-2.sh
	sbatch run-match-4.sh
	sbatch run-match-8.sh
	sbatch run-match-16.sh

3dlaprunmatch:
	sbatch run-3dlap-match-gpu-1.sh
	sbatch run-3dlap-match-gpu-2.sh
	sbatch run-3dlap-match-gpu-4.sh
	sbatch run-3dlap-match-gpu-8.sh
	sbatch run-3dlap-match-gpu-16.sh

3dlaprunsoc1:
	sbatch run-3dlap-soc1-gpu-1.sh
	sbatch run-3dlap-soc1-gpu-2.sh
	sbatch run-3dlap-soc1-gpu-4.sh
	sbatch run-3dlap-soc1-gpu-8.sh
	sbatch run-3dlap-soc1-gpu-16.sh

clean:
	rm anisopsblascuda
	rm *.out
