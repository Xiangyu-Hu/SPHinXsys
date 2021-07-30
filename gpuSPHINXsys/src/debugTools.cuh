#ifndef DEBUGTOOLS_CUH
#define DEBUGTOOLS_CUH

#define CUDA_ERROR_CHECK

#define CudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define CudaCheckError() __cudaCheckError(__FILE__, __LINE__)

#include <string>

namespace gpu{

  class cuda_generic_error: public std::runtime_error{
    cudaError_t error_code;
  public:
    cuda_generic_error(std::string msg, cudaError_t err):
      std::runtime_error(msg + ": " + cudaGetErrorString(err) + " - code: " + std::to_string(err)),
      error_code(err){}

    cudaError_t code(){return error_code;}
  };

}

inline void __cudaSafeCall(cudaError err, const char *file, const int line){
  #ifdef CUDA_ERROR_CHECK
  if (cudaSuccess != err){
    cudaGetLastError(); //Reset CUDA error status    
    throw gpu::cuda_generic_error("CudaSafeCall() failed at "+
				    std::string(file) + ":" + std::to_string(line), err);
  }
  #endif
}

inline void __cudaCheckError(const char *file, const int line){
  cudaError err;
  err = cudaGetLastError();
  if(cudaSuccess != err){
    throw gpu::cuda_generic_error("CudaCheckError() failed at "+
				    std::string(file) + ":" + std::to_string(line), err);
  }
}

#endif
