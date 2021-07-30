#ifndef PARTICLESORTER_CUH
#define PARTICLESORTER_CUH

#include"Box.cuh"
#include"Grid.cuh"
#include"System.h"
#include"debugTools.cuh"
#include<cub/cub.cuh>

namespace gpu{

  namespace Sorter{

    struct MortonHash{
      Grid grid;
      MortonHash(Grid grid): grid(grid){}
      //Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros. See [1]
      inline __host__ __device__ uint encodeMorton(const uint &i) const{
	uint x = i;
	x &= 0x3ff;
	x = (x | x << 16) & 0x30000ff;
	x = (x | x << 8)  & 0x300f00f;
	x = (x | x << 4)  & 0x30c30c3;
	x = (x | x << 2)  & 0x9249249;
	return x;
      }
      /*Fuse three 10 bit numbers in 32 bits, producing a Z order Morton hash*/
      inline __host__ __device__ uint hash(int3 cell) const{
	return encodeMorton(cell.x) | (encodeMorton(cell.y) << 1) | (encodeMorton(cell.z) << 2);
      }

      inline __host__ __device__ uint operator()(real4 pos) const{
	const int3 cell = grid.getCell(pos);
	return hash(cell);
      }

    };
    //The hash is the cell 1D index, this pattern is better than random for neighbour transverse, but worse than Morton
    struct CellIndexHash{
      Grid grid;
      CellIndexHash(Grid grid): grid(grid){}
      inline __device__ __host__ uint hash(int3 cell) const{
	return grid.getCellIndex(cell);
      }
      inline __host__ __device__ uint operator()(real4 pos) const{
	const int3 cell = grid.getCell(pos);
	return hash(cell);
      }

    };

    static int _ffs(int i)
    {
        int bit;

        if (0 == i)
            return 0;

        for (bit = 1; !(i & 1); ++bit)
            i >>= 1;
        return bit;
    }

    int clz(uint n){
      n |= (n >>  1);
      n |= (n >>  2);
      n |= (n >>  4);
      n |= (n >>  8);
      n |= (n >> 16);
      return 32-_ffs(n - (n >> 1));
    }

    template<class HashIterator>
    __global__ void assignHash(HashIterator hasher,
			       int* __restrict__ index,
			       uint* __restrict__ hash , int N){
      const int i = blockIdx.x*blockDim.x + threadIdx.x;
      if(i >= N) return;
      const uint ihash = hasher[i];
      index[i] = i;
      hash[i]  = ihash;
    }

    template<class InputIterator, class OutputIterator>
    __global__ void reorderArray(const InputIterator old,
				 OutputIterator sorted,
				 int* __restrict__ pindex, int N){
      int i = blockIdx.x*blockDim.x + threadIdx.x;
      if(i>=N) return;
      sorted[i] = old[pindex[i]];
    }

  }

  class ParticleSorter{
  public:
    ParticleSorter() = delete;
    ParticleSorter(std::shared_ptr<System> sys):sys(sys){};

    template<class hashType>
    void sortByKey(cub::DoubleBuffer<int> &index,
		   cub::DoubleBuffer<hashType> &hash,
		   int N, cudaStream_t st = 0, int end_bit = sizeof(uint)*8){
      if(N > maxRequestedElements){
        maxRequestedElements = N;
	cub_temp_storage_bytes = 0;
	CudaSafeCall(cub::DeviceRadixSort::SortPairs(nullptr, cub_temp_storage_bytes,
						     hash,
						     index,
						     N,
						     0, end_bit,
						     st));
      }
      auto alloc = sys->getTemporaryDeviceAllocator<char>();
      std::shared_ptr<char> d_temp_storage(alloc.allocate(cub_temp_storage_bytes),
					   [=](char* ptr){ alloc.deallocate(ptr);});
      void* d_temp_storage_ptr = d_temp_storage.get();
      CudaSafeCall(cub::DeviceRadixSort::SortPairs(d_temp_storage_ptr, cub_temp_storage_bytes,
						   hash,
						   index,
						   N,
						   0, end_bit,
						   st));
    }
    //Return the most significant bit of an unsigned integral type
    template<typename T> inline int msb(T n){
      static_assert(std::is_integral<T>::value && !std::is_signed<T>::value,
		    "msb<T>(): T must be an unsigned integral type.");
      for(T i = std::numeric_limits<T>::digits - 1, mask = 1 << i;
	   i >= 0;
	   --i, mask >>= 1){
	if((n & mask) != 0) return i;
	}
      return 0;
    }

    template<class HashIterator>
    void updateOrderWithCustomHash(HashIterator &hasher,
				   uint N, uint maxHash = std::numeric_limits<uint>::max(),
				   cudaStream_t st = 0){
      sys->log<System::DEBUG1>("[ParticleSorter] Updating with custom hash iterator");
      sys->log<System::DEBUG1>("[ParticleSorter] Assigning hash to %d elements", N);
      sys->log<System::DEBUG2>("[ParticleSorter] Maximum hash 0x%x ( dec: %u ) last bit: %d",
			       maxHash, maxHash, 32-Sorter::clz(maxHash) );
      try{
	tryToUpdateOrderWithCustomHash(hasher, N, maxHash, st);
      }
      catch(...){
	sys->log<System::ERROR>("ParticleSorter raised an exception in updateOrderWithCustomHash");
	throw;
      }
    }

    template<class CellHasher = Sorter::MortonHash, class InputIterator>
    void updateOrderByCellHash(InputIterator pos, uint N, Box box, int3 cellDim, cudaStream_t st = 0){
      Grid grid(box, cellDim);
      CellHasher hasher(grid);
      auto hashIterator = thrust::make_transform_iterator(pos, hasher);
      auto maxHash = hasher.hash(cellDim-1);
      this->updateOrderWithCustomHash(hashIterator, N, maxHash, st);
    }

    void updateOrderById(int *id, int N, cudaStream_t st = 0){
      try{
	tryToUpdateOrderById(id, N, st);
      }
      catch(...){
	sys->log<System::ERROR>("Exception raised in ParticleSorter::updateOrderById");
	throw;
      }
    }

    //WARNING: _unsorted and _sorted cannot be aliased!
    template<class InputIterator, class OutputIterator>
    void applyCurrentOrder(InputIterator d_property_unsorted, OutputIterator d_property_sorted,
			   int N, cudaStream_t st = 0){
      int Nthreads=128;
      int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);
      Sorter::reorderArray<<<Nblocks, Nthreads, 0, st>>>(d_property_unsorted,
							 d_property_sorted,
							 thrust::raw_pointer_cast(index.data()),
							 N);
      CudaCheckError();
    }

    int * getSortedIndexArray(int N){
      try{
	return tryToGetSortedIndexArray(N);
      }
      catch(...){
	sys->log<System::ERROR>("Exception raised in ParticleSorter::getSortedIndexArray");
	throw;
      }
    }

    uint * getSortedHashes(){
      return thrust::raw_pointer_cast(hash.data());
    }

    int * getIndexArrayById(int * id, int N, cudaStream_t st = 0){
      try{
	return tryToGetIndexArrayById(id, N, st);
      }
      catch(...){
	sys->log<System::ERROR>("Exception raised in ParticleSorter::getIndexArrayById");
	throw;
      }
    }

  private:
    bool init = false;
    bool originalOrderNeedsUpdate = true;

    int maxRequestedElements = 0;
    size_t cub_temp_storage_bytes = 0;

    thrust::device_vector<int>  original_index;
    thrust::device_vector<int>  index, index_alt;
    thrust::device_vector<uint> hash, hash_alt;

    std::shared_ptr<System> sys;

    void tryToUpdateOrderById(int *id, int N, cudaStream_t st){
      original_index.resize(N);
      index_alt.resize(N);
      hash.resize(N);
      hash_alt.resize(N);
      cub::CountingInputIterator<int> ci(0);
      thrust::copy(thrust::cuda::par, ci, ci+N, original_index.begin());
      int* d_hash = (int*)thrust::raw_pointer_cast(hash.data());
      CudaSafeCall(cudaMemcpyAsync(d_hash, id, N*sizeof(int), cudaMemcpyDeviceToDevice,st));
      auto db_index = cub::DoubleBuffer<int>(thrust::raw_pointer_cast(original_index.data()),
					     thrust::raw_pointer_cast(index_alt.data()));
      auto db_hash  = cub::DoubleBuffer<int>(d_hash, (int*)thrust::raw_pointer_cast(hash_alt.data()));
      this->sortByKey(db_index, db_hash, N, st);
      if(db_index.selector)
	original_index.swap(index_alt);
    }

    template<class HashIterator>
    void tryToUpdateOrderWithCustomHash(HashIterator &hasher, uint N, uint maxHash, cudaStream_t st){
      init = true;
      hash.resize(N);
      index.resize(N);
      int Nthreads=128;
      int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);
      Sorter::assignHash<<<Nblocks, Nthreads, 0, st>>>(hasher,
						       thrust::raw_pointer_cast(index.data()),
						       thrust::raw_pointer_cast(hash.data()),
						       N);
      CudaCheckError();
      hash_alt.resize(N);
      index_alt.resize(N);
      auto db_index = cub::DoubleBuffer<int>(thrust::raw_pointer_cast(index.data()),
					     thrust::raw_pointer_cast(index_alt.data()));
      auto db_hash  = cub::DoubleBuffer<uint>(thrust::raw_pointer_cast(hash.data()),
					      thrust::raw_pointer_cast(hash_alt.data()));
      //Cub just needs this endbit at least
      int maxbit = 32-Sorter::clz(maxHash);
      maxbit = std::min(maxbit, 32);
      this->sortByKey(db_index, db_hash, N, st, maxbit);
      CudaCheckError();
      //Sometimes CUB will not swap the references in the DoubleBuffer
      if(db_index.selector)
	index.swap(index_alt);
      if(db_hash.selector)
	hash.swap(hash_alt);
      originalOrderNeedsUpdate = true;
    }

    int * tryToGetSortedIndexArray(int N){
      int lastN = index.size();
      if(lastN != N){
	cub::CountingInputIterator<int> ci(lastN);
        index.resize(N);
	thrust::copy(thrust::cuda::par, ci, ci+(N-lastN), index.begin()+lastN);
      }
      return thrust::raw_pointer_cast(index.data());
    }

    int * tryToGetIndexArrayById(int * id, int N, cudaStream_t st = 0){
      if(!init) return id;
      if(originalOrderNeedsUpdate){
	this->updateOrderById(id, N, st);
	originalOrderNeedsUpdate = false;
      }
      int lastN = original_index.size();
      if(lastN != N){
	original_index.resize(N);
	thrust::copy(thrust::cuda::par, id, id+(N-lastN), original_index.begin()+lastN);
      }
      return thrust::raw_pointer_cast(original_index.data());
    }

  };
}
#endif
