./configure --prefix= \
	--with-ccopt="-O3 -fPIC" \
	--with-fcopt="-O3 -frecursive" \
	--with-cxxopt="-O3" \
	--with-cuda=${CUDA_HOME} \
	--with-cudacc=90

