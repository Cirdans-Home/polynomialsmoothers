./configure \
	--with-metisincdir=/software/NUMERICAL/PARMETIS/parmetis-4.0.3/metis/include/ \
	--with-metislibdir=/software/NUMERICAL/PARMETIS/parmetis-4.0.3/metis/build/Linux-x86_64/libmetis/ \
	--with-ipk=4 --with-lpk=4 \
	--with-amddir=/usr/include/suitesparse/ \
	--prefix=/home/fdurastante/polynomialsmoothers/install \
	--with-cuda=/opt/cuda/12.0/ \
	--with-cudacc=50,52,53,60,61,62,70,72,75 \
	--enable-openmp	
