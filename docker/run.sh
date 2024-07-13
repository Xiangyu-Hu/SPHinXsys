# docker run --rm --cpus 40 --gpus '"device=0,2"' -v ./:/home/SPHinXsys sycl bash run.sh
# docker run --rm --cpus 40 -v ./:/home/SPHinXsys sycl bash run_cpu.sh
source /opt/intel/oneapi/setvars.sh --include-intel-llvm
cd ./build/tests/2d_examples/test_2d_dambreak_sycl/ && ctest -V