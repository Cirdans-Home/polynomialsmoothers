./configure \
	--with-metisincdir=/opt/metis/include/ \
	--with-metislibdir=/opt/metis/lib/ \
	--with-ipk=4 --with-lpk=8 \
	--with-amddir=/usr/include/suitesparse/ \
	--prefix=/home/cirdan/Documenti/Paper/ParallelSmoothers/polynomialsmoothers/install\
	--enable-openmp \
	--with-cuda=/usr/local/cuda-12.3/ 
	--with-cudacc=50,52,53,60,61,62,70,72,75
