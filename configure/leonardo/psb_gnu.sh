./configure --prefix=/leonardo/home/userexternal/fdurasta/polynomialsmoothers/install \
	--with-metisincdir=/leonardo/home/userexternal/fdurasta/polynomialsmoothers/install/include \
	--with-metislibdir=/leonardo/home/userexternal/fdurasta/polynomialsmoothers/install/lib \
	--with-ccopt="-O3 -fPIC -march=native" \
	--with-fcopt="-O3 -frecursive -march=native" \
	--with-cxxopt="-O3" \
	--with-cuda=${CUDA_HOME} \
	--with-cudacc=80

