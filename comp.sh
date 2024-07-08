#!/bin/bash

while getopts ":s:h" option; do
  case $option in
    s)
      solver_name="$OPTARG"
      ;;
    h)
      echo "Usage: $0 [-s solver_name]"
      exit 1
      ;;
  esac
done
MAIN=""

if [ ! -d "bin" ]; then
  mkdir -p bin
fi
if test -f "solvers/${solver_name}/${solver_name}.cu"; then
    MAIN="solvers/${solver_name}/${solver_name}.cu"
    echo "compiling: ${MAIN}"
    nvcc ${MAIN} solvers/${solver_name}/parse_input.cpp src/field_init.cpp src/field.cpp src/term_init.cpp src/term.cpp src/evolver.cpp src/field_kernels.cu src/term_kernels.cu -lcufft -lcurand -lfftw3f -DWITHCUDA -o ./bin/$solver_name -I /home/paarthgulati/local/fftw3/include -L /home/paarthgulati/local/fftw3/lib
else
    echo "Solver not found; Usage: $0 [-s solver_name]"
fi
