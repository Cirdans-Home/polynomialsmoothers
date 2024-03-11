./configure --prefix=/home/hpc/ihpc/ihpc100h/polynomialsmoothers/install \
	--with-metisincdir=/home/hpc/ihpc/ihpc100h/polynomialsmoothers/install/include \
	--with-metislibdir=/home/hpc/ihpc/ihpc100h/polynomialsmoothers/install/lib \
	--with-ccopt="-Ofast -funroll-all-loops -fPIC -march=znver3" \
	--with-fcopt="-Ofast -funroll-all-loops -frecursive -march=znver3" \
	--with-cxxopt="-Ofast -funroll-all-loops -march=znver3" \
	--with-blasdir=/home/hpc/ihpc/ihpc100h/polynomialsmoothers/install/lib \
	--with-cuda=${CUDA_HOME} \
	--with-cudacc=80

# --with-blas="-m64  -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl" \


