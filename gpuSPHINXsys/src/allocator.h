#ifndef ALLOCATOR_H
#define ALLOCATOR_H
#include<thrust/device_ptr.h>
#include<thrust/system/cuda/memory.h>
#include"debugTools.cuh"
#include<map>

namespace gpu{

  namespace detail{
    template<class T> using cuda_ptr = thrust::device_ptr<T>;
    template<class T>
    class memory_resource{
    public:
      using pointer = T;
      using max_align_t = long double; //This C++11 alias is not available in std with g++-4.8.5
      pointer allocate(std::size_t bytes, std::size_t alignment = alignof(max_align_t)){
	return do_allocate(bytes, alignment);
      }

      void deallocate(pointer p, std::size_t bytes, std::size_t alignment = alignof(max_align_t)){
	return do_deallocate(p, bytes, alignment);
      }

      bool is_equal(const memory_resource &other) const noexcept{
	return do_is_equal(other);
      }

      virtual pointer do_allocate(std::size_t bytes, std::size_t alignment) = 0;

      virtual void do_deallocate(pointer p, std::size_t bytes, std::size_t alignment) = 0;

      virtual bool do_is_equal(const memory_resource &other) const noexcept{
	return this == &other;
      }

    };

    template<class MR>
    MR* get_default_resource(){
      static MR default_resource;
      return &default_resource;
    }
  }

  //A device memory resource, with cuda raw pointers
  class device_memory_resource : public detail::memory_resource<void*>{
    using super = detail::memory_resource<void*>;
  public:

    using pointer = typename super::pointer;

    virtual pointer do_allocate(std::size_t bytes, std::size_t alignment) override{
      return thrust::raw_pointer_cast(thrust::cuda::malloc<char>(bytes));
    }

    virtual void do_deallocate(pointer p, std::size_t bytes, std::size_t alignment) override{
      thrust::cuda::pointer<void> void_ptr(p);
      thrust::cuda::free(void_ptr);
    }    
    
  };

  //A pool device memory_resource, stores previously allocated blocks in a cache
  // and retrieves them fast when similar ones are allocated again (without calling malloc everytime).
  template<class MR>
  struct pool_memory_resource_adaptor: public detail::memory_resource<typename MR::pointer>{
  private:
    using super = detail::memory_resource<typename MR::pointer>;
    MR* res;
  public:
    using pointer = typename super::pointer;

    ~pool_memory_resource_adaptor(){
      try{
	free_all();
      }
      catch(...){
      }
    }

    pool_memory_resource_adaptor(MR* resource): res(resource){}
    pool_memory_resource_adaptor(): res(detail::get_default_resource<MR>()){}

    using FreeBlocks =  std::multimap<std::ptrdiff_t, void*>;
    using AllocatedBlocks =  std::map<void*, std::ptrdiff_t>;
    FreeBlocks free_blocks;
    AllocatedBlocks allocated_blocks;

    virtual pointer do_allocate( std::size_t bytes, std::size_t alignment) override{
      pointer result;
      std::ptrdiff_t blockSize = 0;
      auto available_blocks = free_blocks.equal_range(bytes);
      auto available_block = available_blocks.first;
      //Look for a block of the same size
      if(available_block == free_blocks.end()){
	available_block = available_blocks.second;
      }
      //Try to find a block greater than requested size
      if(available_block != free_blocks.end() ){
	result = pointer(available_block -> second);
	blockSize = available_block -> first;
	free_blocks.erase(available_block);
      }
      else{
	result = res->do_allocate(bytes, alignment);
	blockSize = bytes;
      }
      allocated_blocks.insert(std::make_pair(thrust::raw_pointer_cast(result), blockSize));
      return result;
    }

    virtual void do_deallocate(pointer p, std::size_t bytes, std::size_t alignment) override{
      auto block = allocated_blocks.find(thrust::raw_pointer_cast(p));
      if(block == allocated_blocks.end()){
	throw  std::system_error(EFAULT, std::generic_category(), "Address is not handled by this instance.");
      }
      std::ptrdiff_t num_bytes = block->second;
      allocated_blocks.erase(block);
      free_blocks.insert(std::make_pair(num_bytes, thrust::raw_pointer_cast(p)));
    }

    virtual bool do_is_equal(const super &other) const noexcept override {
      return res->do_is_equal(other);
    }

    void free_all(){
      for(auto &i: free_blocks) res->do_deallocate(static_cast<pointer>(i.second), i.first, 0);
      for(auto &i: allocated_blocks) res->do_deallocate(static_cast<pointer>(i.first), i.second, 0);
      free_blocks.clear();
      allocated_blocks.clear();
    }

  };

  namespace detail{
    //Takes a pointer type (including smart pointers) and returns a reference to the underlying type
    template<class T> struct pointer_to_lvalue_reference{
    private:
      using element_type = typename std::pointer_traits<T>::element_type;
    public:
      using type = typename std::add_lvalue_reference<element_type>::type;
    };

    //Specialization for special thrust pointer/reference types...
    template<class T> struct pointer_to_lvalue_reference<detail::cuda_ptr<T>>{
      using type = thrust::system::cuda::reference<T>;
    };

    template<class T> struct non_void_value_type{using type = T;};
    template<> struct non_void_value_type<void>{using type = char;};

  }

  //An allocator that can be used for any type using the same underlying memory_resource.
  //pointer type can be specified to work with thrust cuda pointers
  template<class T,
	   class MR = detail::memory_resource<void*>,
	   class void_pointer = T*>
  class polymorphic_allocator{
    MR * res;
  public:
    //C++17 definitions for allocator interface
    using size_type = std::size_t;

    using value_type = T;
    using value_size_type = typename detail::non_void_value_type<value_type>::type;

    using differente_type = std::ptrdiff_t;

    //All of the traits below are deprecated in C++17, but thrust counts on them
    //using void_pointer = T*;
    using pointer = typename std::pointer_traits<void_pointer>::template rebind<value_type>;

    using reference = typename detail::pointer_to_lvalue_reference<pointer>::type;
    using const_reference = typename detail::pointer_to_lvalue_reference<std::add_const<pointer>>::type;

    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;

    template<class Other>
    polymorphic_allocator(Other other) : res(other.resource()){}

    polymorphic_allocator(MR * resource) : res(resource){}

    polymorphic_allocator() : res(detail::get_default_resource<MR>()){}

    MR* resource() const { return this->res;}

    pointer allocate (size_type n) const{
      return static_cast<pointer>(static_cast<value_type*>(this->res->do_allocate(n * sizeof(value_size_type),
										  alignof(value_size_type))));
    }

    void deallocate (pointer p, size_type n = 0) const{
      return this->res->do_deallocate(thrust::raw_pointer_cast(p),
				      n * sizeof(value_size_type),
				      alignof(value_size_type));
    }
  };

  //thrust versions are not reliable atm, std::pmr is C++17 which CUDA 10.1 does not support yet and is not thrust::cuda compatible (thrust expects device allocators to return thurst::cuda::pointers)
  //template<class T> using memory_resource = thrust::mr::memory_resource<T>;
  //template<class T> using polymorphic_resource_adaptor = thrust::mr::polymorphic_adaptor_resource<T>;
  //template<typename T, class MR> using allocator = thrust::mr::allocator<T, MR>;
  //template<class T> using polymorphic_allocator = thrust::mr::polymorphic_allocator<T>;

}

#endif
