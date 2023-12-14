anisopsblas:
	ln -s /home/it4i-fdurast/polynomialsmoothers/anisotropy/runs/anisopsblas anisopsblas
	ln -s /home/it4i-fdurast/polynomialsmoothers/anisotropy/runs/anisopsblascuda anisopsblascuda

genmatch:
	./gen_script.sh run-match-gpu.sh

gensoc:
	./gen_script.sh run-soc1-gpu.sh

cleanscript:
	rm run-match-gpu-*
	rm run-soc1-gpu-*

clean:
	rm anisopsblas
	rm anisopsblascuda
	rm *.out
