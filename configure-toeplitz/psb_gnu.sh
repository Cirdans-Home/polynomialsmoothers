./configure --prefix=/home/durastante/polynomialsmoothers/install \
	--with-metisincdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/metis-5.1.0-kuao2tv62i5uiifpbtfnx2nqa4zuaaci/include \
	--with-metislibdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/metis-5.1.0-kuao2tv62i5uiifpbtfnx2nqa4zuaaci/lib \
	--with-blasdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/openblas-0.3.20-4ipfq4m2x4ku6edyd6o6sq4zxprh4kmo/lib/ \
	--with-amdlibdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/suite-sparse-5.10.1-klxivyhmvz5ckb3zv4fzjtwk6itursni/lib \
	--with-amdincdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/suite-sparse-5.10.1-klxivyhmvz5ckb3zv4fzjtwk6itursni/include \
	--with-ccopt="-O3 -g -ggdb -fPIC -march=native -mtune=native" \
	--with-fcopt="-O3 -g -ggdb -frecursive -fPIC -march=native -mtune=native" \
	--with-cxxopt="-O3 -g -ggdb -fPIC -march=native -mtune=native" \
        --with-lpk=4 --enable-openmp	

