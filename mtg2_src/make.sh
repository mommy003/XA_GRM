
  #yum install glibc-static
  #yum install zlib-devel
  #yum install zlib-static

  #libmkl_blacs_ilp64.a -> libmkl_blacs_intelmpi_ilp64.a  in Makefile

  #source /opt/intel/bin/compilervars.sh intel64
  #source /opt/intel/mkl/bin/mklvars.sh intel64
  #source /home/hsc-admhl/intel/bin/compilervars.sh intel64
  #source /home/hsc-admhl/intel/mkl/bin/mklvars.sh intel64
  #source /data/alh-admhl/intel/bin/compilervars.sh intel64
  #source /data/alh-admhl/intel/mkl/bin/mklvars.sh intel64
  source /opt/intel/oneapi/setvars.sh

  # set unlimited stack size
  ulimit -s unlimited

  rm *.o
  make -f Makefile-mtg2

