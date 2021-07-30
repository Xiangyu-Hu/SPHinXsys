#ifndef SYSTEM_H
#define SYSTEM_H

#include"defines.h"
#include"debugTools.cuh"
#include"Log.h"
#include<cstring>
#include<ostream>
#include<iostream>
#include<cstdarg>
#include"allocator.h"
#include<memory>

namespace gpu{

  using std::shared_ptr;
  using std::string;

  struct access{
    enum location{cpu, gpu, nodevice};
    enum mode{read, write, readwrite, nomode};
  };

  class insuficient_compute_capability_exception: public std::exception{
    const char* what() const noexcept {
      return "Insuficient compute capability";
    }
  };

  class System{
  public:

    using resource = gpu::device_memory_resource;
    using device_temporary_memory_resource = gpu::pool_memory_resource_adaptor<resource>;

    template<class T>
    using allocator = gpu::polymorphic_allocator<T , device_temporary_memory_resource>;

    enum LogLevel{CRITICAL=0, ERROR, EXCEPTION, WARNING, MESSAGE, STDERR, STDOUT,
          DEBUG, DEBUG1, DEBUG2, DEBUG3, DEBUG4, DEBUG5, DEBUG6, DEBUG7};

  private:

    void initializeCUDA(){
        try{
            CudaSafeCall(cudaFree(0));
            CudaSafeCall(cudaDeviceSynchronize());
            CudaCheckError();
        }
        catch(...){
            log<System::ERROR>("[System] Exception raised at CUDA initialization");
            throw;
        }
        log<System::MESSAGE>("[System] CUDA initialized");
    }

  public:
    System(){
      this->initializeCUDA();
      CudaCheckError();
    }

    void finish(){
      log<DEBUG2>("[System] finish");
      CudaSafeCall(cudaDeviceSynchronize());
      CudaCheckError();
    }

    template<int level, class ...T>
    static inline void log(char const *fmt, T... args){
      Logging::log<level>(fmt, args...);
      if(level == CRITICAL){
          throw std::runtime_error("System encountered an unrecoverable error");
      }
    }
    template<int level>
    static inline void log(const std::string &msg){
      log<level>("%s", msg.c_str());
    }

    template<class T = char>
    static allocator<T> getTemporaryDeviceAllocator(){
      return allocator<T>();
    }

  };

}
#endif
