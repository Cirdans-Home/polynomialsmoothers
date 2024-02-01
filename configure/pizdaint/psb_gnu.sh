./configure --prefix=/users/pdambra/polynomialsmoothers/install \
	--with-metisincdir=/users/pdambra/polynomialsmoothers/install/include \
	--with-metislibdir=/users/pdambra/polynomialsmoothers/install/lib \
	--with-blasdir=/users/pdambra/polynomialsmoothers/install \
	--with-ccopt="-O3 -fPIC" \
	--with-fcopt="-O3 -frecursive" \
	--with-cxxopt="-O3" \
	--with-cuda=${CUDA_HOME} \
	MPICC=cc MPICXX=CC MPIFC=ftn

