./configure --prefix=/home/it4i-fdurast/polynomialsmoothers/install \
	--with-metisincdir=/home/it4i-fdurast/polynomialsmoothers/install/include \
	--with-metislibdir=/home/it4i-fdurast/polynomialsmoothers/install/lib \
	--with-libs="-lGKlib -L/home/it4i-fdurast/polynomialsmoothers/install/lib"\
	--with-ccopt="-O3 -fPIC" \
	--with-fcopt="-O3 -frecursive" \
	--with-cxxopt="-O3" \
	--with-cuda=${CUDA_HOME} \
	--with-cudacc=80

