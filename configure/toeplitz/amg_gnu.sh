./configure --with-psblas="/home/durastante/polynomialsmoothers/install" \
	--prefix="/home/durastante/polynomialsmoothers/install/"\
	--with-umfpacklibdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/suite-sparse-5.10.1-klxivyhmvz5ckb3zv4fzjtwk6itursni/lib \
	--with-umfpackincdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/suite-sparse-5.10.1-klxivyhmvz5ckb3zv4fzjtwk6itursni/include \
	--with-blasdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/openblas-0.3.20-4ipfq4m2x4ku6edyd6o6sq4zxprh4kmo/lib \
	--with-mumpslibdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/mumps-5.4.1-euvufbtdti7ozt65nbcvph5twfrlmqmi/lib \
	--with-mumpsincdir=/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/mumps-5.4.1-euvufbtdti7ozt65nbcvph5twfrlmqmi/include\
	--with-extra-libs="-fopenmp -lopenblas -L/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/openblas-0.3.20-4ipfq4m2x4ku6edyd6o6sq4zxprh4kmo/lib -lscalapack -L/data/software/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-12.2.0/netlib-scalapack-2.2.0-b4njionfcuwwa3heilv7dhvetprj2wmo/lib/"
