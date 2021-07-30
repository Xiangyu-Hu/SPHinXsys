/*
References:

[1] http://developer.download.nvidia.com/assets/cuda/files/particles.pdf

 */
#ifndef CELLLIST_CUH
#define CELLLIST_CUH

#include"ParticleData.cuh"
#include"ParticleGroup.cuh"
#include"ParticleSorter.cuh"
#include"Box.cuh"
#include"Grid.cuh"
#include"System.h"
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include<cub/cub.cuh>
#include<limits>

namespace gpu{
  namespace CellList_ns{
    template<class InputIterator>
    __global__ void fillCellList(InputIterator sortPos,
				 uint *cellStart, int *cellEnd,
				 uint currentValidCell,
				 int *errorFlag,
				 int N, Grid grid){
      uint id = blockIdx.x*blockDim.x + threadIdx.x;
      if(id<N){
	uint icell, icell2;
	icell = grid.getCellIndex(grid.getCell(make_real3(sortPos[id])));
	if(id>0){ /*Shared memory target VVV*/
	  icell2 = grid.getCellIndex(grid.getCell(make_real3(sortPos[id-1])));
	}
	else{
	  icell2 = 0;
	}
	const int ncells = grid.getNumberCells();
	if(icell>=ncells or icell2>=ncells){
	  errorFlag[0] = 1;
	  return;
	}
	if(icell != icell2 or id == 0){
	  cellStart[icell] = id+currentValidCell;
	  if(id>0)
	    cellEnd[icell2] = id;
	}
	if(id == N-1) cellEnd[icell] = N;
      }
    }

    template<class NeighbourContainer>
    __global__ void fillNeighbourList(const real4* sortPos,
				      const int *groupIndex,
				      NeighbourContainer ni,
				      int *neighbourList, int*nNeighbours,
				      int maxNeighbours,
				      real cutOff2,
				      int N, Box box,
				      int* tooManyNeighboursFlag){
      int id = blockIdx.x*blockDim.x + threadIdx.x;
      if(id>=N) return;
      int nneigh = 0;
      const int offset = id*maxNeighbours;
#if CUB_PTX_ARCH < 300
      constexpr auto cubModifier = cub::LOAD_DEFAULT;
#else
      constexpr auto cubModifier = cub::LOAD_LDG;
#endif
      const real3 pi = make_real3(cub::ThreadLoad<cubModifier>(sortPos + id));
      ni.set(id);
      auto it = ni.begin();
      while(it){
	auto n = *it++;
	const int cur_j = n.getInternalIndex();
	const real3 pj = make_real3(cub::ThreadLoad<cubModifier>(sortPos + cur_j));
	const real3 rij = box.apply_pbc(pj-pi);
	if(dot(rij, rij) <= cutOff2){
	  nneigh++;
	  if(nneigh>=maxNeighbours){
	    atomicMax(tooManyNeighboursFlag, nneigh);
	    return;
	  }
	  neighbourList[offset + nneigh-1] = cur_j;
	}
      }
      //Include self interactions
      neighbourList[offset + nneigh] = id;
      nNeighbours[id] = nneigh+1;
    }
  }

  class CellList{
  protected:
    thrust::device_vector<uint> cellStart;
    uint currentValidCell;
    int currentValidCell_counter;
    thrust::device_vector<int>  cellEnd;

    Grid grid;

    thrust::device_vector<int>  errorFlags;

    thrust::device_vector<int> neighbourList, numberNeighbours;
    thrust::device_vector<real4> sortPos;
    int maxNeighboursPerParticle;
    ParticleSorter ps;

    shared_ptr<ParticleData> pd;
    shared_ptr<ParticleGroup> pg;
    shared_ptr<System> sys;

    bool force_next_update = true;
    bool rebuildNlist;

    real3 currentCutOff;

    connection numParticlesChangedConnection, posWriteConnection;
    cudaEvent_t event;


    void handleNumParticlesChanged(int Nnew){
      sys->log<System::DEBUG>("[CellList] Number particles changed signal handled.");
      int numberParticles = pg->getNumberParticles();
      if(neighbourList.size()){
	neighbourList.resize(numberParticles*maxNeighboursPerParticle);
	numberNeighbours.resize(numberParticles);
      }
      currentValidCell_counter = -1;
      force_next_update = true;
    }

    void handlePosWriteRequested(){
      sys->log<System::DEBUG1>("[CellList] Issuing a list update after positions were written to.");
      force_next_update = true;
    }

    struct NeighbourListOffsetFunctor{
      NeighbourListOffsetFunctor(int str, int* groupIndex):
	stride(str), groupIndex(groupIndex){}
      int stride;
      int *groupIndex;
      inline __host__ __device__ int operator()(const int &index) const{
	return groupIndex[index]*stride;
      }
    };

    using CountingIterator = cub::CountingInputIterator<int>;
    using StrideIterator = cub::TransformInputIterator<int, NeighbourListOffsetFunctor, CountingIterator>;

  public:

    CellList(shared_ptr<ParticleData> pd,
	     shared_ptr<System> sys):
      CellList(pd, std::make_shared<ParticleGroup>(pd, sys), sys){}

    CellList(shared_ptr<ParticleData> pd,
	     shared_ptr<ParticleGroup> pg,
	     shared_ptr<System> sys): pd(pd), pg(pg), sys(sys), ps(sys), currentCutOff(real3()){
      sys->log<System::MESSAGE>("[CellList] Created");

      maxNeighboursPerParticle = 32;

      pd->getNumParticlesChangedSignal()->connect([this](int Nnew){this->handleNumParticlesChanged(Nnew);});
      pd->getPosWriteRequestedSignal()->connect([this](){this->handlePosWriteRequested();});

      CudaSafeCall(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));

      currentValidCell_counter = -1;
      CudaCheckError();
    }

    ~CellList(){
      sys->log<System::DEBUG>("[CellList] Destroyed");
      numParticlesChangedConnection.disconnect();
      posWriteConnection.disconnect();
    }

    struct NeighbourListData{
      int * neighbourList;
      int *numberNeighbours;
      StrideIterator particleStride = StrideIterator(CountingIterator(0), NeighbourListOffsetFunctor(0, nullptr));
    };

    NeighbourListData getNeighbourList(cudaStream_t st = 0){
      if(currentCutOff.x != currentCutOff.y or
	 currentCutOff.x != currentCutOff.z or
	 currentCutOff.z != currentCutOff.y){
	sys->log<System::ERROR>("[CellList] Invalid cutoff in getNeighbourList: %f %f %f",
				currentCutOff.x, currentCutOff.y, currentCutOff.z);
	throw std::invalid_argument("Cannot use NeighbourList with a different cutOff in each direction");
      }
      const int numberParticles = pg->getNumberParticles();
      if(rebuildNlist){
	rebuildNeighbourList(st);
      }
      else{
	auto pos = pd->getPos(access::location::gpu, access::mode::read);
	auto posGroupIterator = pg->getPropertyInputIterator(pos.raw(), access::location::gpu);
	ps.applyCurrentOrder(posGroupIterator, sortPos.begin(), numberParticles, st);
      }
      NeighbourListData nl;
      nl.neighbourList = thrust::raw_pointer_cast(neighbourList.data());
      nl.numberNeighbours = thrust::raw_pointer_cast(numberNeighbours.data());
      nl.particleStride = StrideIterator(CountingIterator(0),
					 NeighbourListOffsetFunctor(maxNeighboursPerParticle,
								    ps.getSortedIndexArray(numberParticles)));
      return nl;
    }

    bool needsRebuild(Box box, real3 cutOff){
      if(force_next_update) return true;
      if(cutOff.x != currentCutOff.x) return true;
      if(cutOff.y != currentCutOff.y) return true;
      if(cutOff.z != currentCutOff.z) return true;
      if(box != grid.box) return true;
      return false;
    }

    void updateNeighbourList(Grid in_grid, real3 cutOff, cudaStream_t st = 0){
      sys->log<System::DEBUG1>("[CellList] Updating list");
      if(!isGridAndCutOffValid(in_grid, cutOff))
	throw std::runtime_error("CellList encountered an invalid grid and/or cutoff");
      if(needsRebuild(in_grid.box, cutOff) == false)
	return;
      else
	currentCutOff = cutOff;
      pd->hintSortByHash(in_grid.box, cutOff);
      force_next_update = false;
      this->grid = in_grid;
      sys->log<System::DEBUG2>("[CellList] Using %d %d %d cells", grid.cellDim.x, grid.cellDim.y, grid.cellDim.z);
      resizeCellListToCurrentGrid();
      updateCurrentValidCell();
      updateOrderAndStoreInSortPos(st);
      fillCellList(st);
      CudaCheckError();
      rebuildNlist = true;
    }

    void updateNeighbourList(Box box, real cutOff, cudaStream_t st = 0){
      updateNeighbourList(box, make_real3(cutOff), st);
    }

    void updateNeighbourList(Box box, real3 cutOff, cudaStream_t st = 0){
      Grid a_grid = Grid(box, cutOff);
      int3 cellDim = a_grid.cellDim;
      if(cellDim.x < 3) cellDim.x = 3;
      if(cellDim.y < 3) cellDim.y = 3;
      if(box.boxSize.z > real(0.0) && cellDim.z < 3) cellDim.z = 3;
      a_grid = Grid(box, cellDim);
      updateNeighbourList(a_grid, cutOff, st);
    }

    //This accesor function is part of CellList only, not part of the NeighbourList interface
    //They allow to obtain a reference to the cell list structures to use them outside
    struct CellListData{
      //[all particles in cell 0, all particles in cell 1,..., all particles in cell ncells]
      //cellStart[i] stores the index of the first particle in cell i (in internal index)
      //cellEnd[i] stores the last particle in cell i (in internal index)
      //So the number of particles in cell i is cellEnd[i]-cellStart[i]
      const uint * cellStart;
      const int  * cellEnd;
      const real4 *sortPos;   //Particle positions in internal index
      const int* groupIndex; //Transformation between internal indexes and group indexes
      Grid grid;
      uint VALID_CELL;
    };
    CellListData getCellList(){
      this->updateNeighbourList(grid, currentCutOff);
      CellListData cl;
      try{
	cl.cellStart   =  thrust::raw_pointer_cast(cellStart.data());
	cl.cellEnd     =  thrust::raw_pointer_cast(cellEnd.data());
	cl.sortPos     =  thrust::raw_pointer_cast(sortPos.data());
      }
      catch(thrust::system_error &e){
	sys->log<System::CRITICAL>("[CellList] Thrust could not access cellList arrays with error: %s", e.what());
      }
      int numberParticles = pg->getNumberParticles();
      cl.groupIndex =  ps.getSortedIndexArray(numberParticles);
      cl.grid = grid;
      cl.VALID_CELL = currentValidCell;
      return cl;
    }

    class NeighbourContainer; //forward declaration for befriending
  private:
    class NeighbourIterator; //forward declaration for befriending

    //Neighbour is a small accesor for NeighbourIterator
    //Represents a particle, you can ask for its index and position
    struct Neighbour{
      __device__ Neighbour(const Neighbour &other):
      internal_i(other.internal_i){
	groupIndex = other.groupIndex;
	sortPos = other.sortPos;
      }
      //Index in the internal sorted index of the cell list
      __device__ int getInternalIndex(){return internal_i;}
      //Index in the particle group
      __device__ int getGroupIndex(){return groupIndex[internal_i];}
      __device__ real4 getPos(){return cub::ThreadLoad<cub::LOAD_LDG>(sortPos+internal_i);}

    private:
      int internal_i;
      const int* groupIndex;
      const real4* sortPos;
      friend class NeighbourIterator;
      __device__ Neighbour(int i, const int* gi, const real4* sp):
	internal_i(i), groupIndex(gi), sortPos(sp){}
    };

    //This forward iterator must be constructed by NeighbourContainer,
    class NeighbourIterator:
      public thrust::iterator_adaptor<
 NeighbourIterator,
   int,   Neighbour,
   thrust::any_system_tag,
   thrust::forward_device_iterator_tag,
   Neighbour,   int
   >{
      friend class thrust::iterator_core_access;

      int j; //Current neighbour index
      CellListData nl;

      int ci; //Current cell

      int3 celli; //Cell of particle i
      int lastParticle; //Index of last particle in current cell

      uint VALID_CELL;
      //Take j to the start of the next cell and return true, if no more cells remain then return false
      __device__ bool nextcell(){
	const bool is2D = nl.grid.cellDim.z<=1;
	if(ci >= (is2D?9:27)) return false;
	bool empty = true;
	do{
	  int3 cellj = celli;
	  cellj.x += ci%3-1;
	  cellj.y += (ci/3)%3-1;
	  cellj.z =  is2D?0:(celli.z+ci/9-1);
	  cellj = nl.grid.pbc_cell(cellj);
	  const bool isPeriodicCellInNonPeriodicBox = (!nl.grid.box.isPeriodicX() and abs(cellj.x-celli.x)>1) or
	    (!nl.grid.box.isPeriodicY() and abs(cellj.y-celli.y)>1) or
	    (!nl.grid.box.isPeriodicZ() and abs(cellj.z-celli.z)>1);
	  if(!isPeriodicCellInNonPeriodicBox){
	    const int icellj = nl.grid.getCellIndex(cellj);
	    const uint cs = nl.cellStart[icellj];
	    empty = cs<VALID_CELL;
	    lastParticle = empty?-1:nl.cellEnd[icellj];
	    j = empty?-1:int(cs-VALID_CELL);
	  }
	  ci++;
	  if(ci >= (is2D?9:27)) return !empty;
	}while(empty);
	return true;
      }

      //Take j to the next neighbour
      __device__  void increment(){
	if(++j == lastParticle) j = nextcell()?j:-1;
      }

      __device__ Neighbour dereference() const{
        return Neighbour(j, nl.groupIndex, nl.sortPos);
      }

    public:
      //j==-1 means there are no more neighbours and the iterator is invalidated
      __device__  operator bool(){ return j!= -1;}

    private:
      friend class NeighbourContainer;
      __device__ NeighbourIterator(int i, CellListData nl, bool begin):
        j(-2), nl(nl), ci(0), lastParticle(-1), VALID_CELL(nl.VALID_CELL){
	if(begin){
	  celli = nl.grid.getCell(make_real3(cub::ThreadLoad<cub::LOAD_LDG>(nl.sortPos+i)));
	  increment();
	}
	else j = -1;
      }
      //enables searching the neighbors of a given point: probe
      __device__ NeighbourIterator(real3 probe, CellListData nl, bool begin):
        j(-2), nl(nl), ci(0), lastParticle(-1), VALID_CELL(nl.VALID_CELL){
        if(begin){
          celli = nl.grid.getCell(probe);
          increment();
        }
        else j = -1;
      }
    };

  public:
    //This is a pseudocontainer which only purpose is to provide begin() and end() NeighbourIterators for a certain particle
    struct NeighbourContainer{
      int my_i = -1;
      CellListData nl;
      NeighbourContainer(CellListData nl): nl(nl){}
      __device__ void set(int i){this->my_i = i;}
      __device__ NeighbourIterator begin(){return NeighbourIterator(my_i, nl, true);}
      //to search neighbors for a given spatial point, e.g. a porbe
      __device__ NeighbourIterator begin(real3 probe){return NeighbourIterator(probe, nl, true);}
      __device__ NeighbourIterator end(){  return NeighbourIterator(my_i, nl, false);}
    };

    NeighbourContainer getNeighbourContainer(){
      auto nl = getCellList();
      return NeighbourContainer(nl);
    }

    const real4* getPositionIterator(){
      return thrust::raw_pointer_cast(sortPos.data());
    }
    const int* getGroupIndexIterator(){
      auto nl = getCellList();
      return nl.groupIndex;
    }

  private:

    bool isGridAndCutOffValid(Grid in_grid, real3 cutOff){
      if(in_grid.cellDim.x < 3 or
       	 in_grid.cellDim.y < 3 or
       	 (in_grid.cellDim.z < 3 and in_grid.box.boxSize.z != real(0.0))){
       	sys->log<System::ERROR>("[CellList] I cannot work with less than 3 cells per dimension!");
	return false;
      }
      //In the case of 3 cells per direction all the particles are checked anyway
      if(in_grid.cellDim.x != 3 or
       	 in_grid.cellDim.y != 3 or
       	 (in_grid.cellDim.z != 3 and in_grid.box.boxSize.z != real(0.0))){

	if(in_grid.cellSize.x < cutOff.x or
	   in_grid.cellSize.y < cutOff.y or
	   (in_grid.cellSize.z < cutOff.z and in_grid.cellSize.z>1)){
	  sys->log<System::ERROR>("[CellList] The cell size cannot be smaller than the cut off.");
	  return false;
	}
      }
      return true;
    }

    void tryToResizeCellListToCurrentGrid(){
      const int ncells = grid.getNumberCells();
      if(cellStart.size()!= ncells){
	currentValidCell_counter = -1;
	cellStart.resize(ncells);
      }
      if(cellEnd.size()!= ncells) cellEnd.resize(ncells);
      CudaCheckError();
    }

    void resizeCellListToCurrentGrid(){
      try{
	tryToResizeCellListToCurrentGrid();
      }
      catch(...){
	sys->log<System::ERROR>("[CellList] Raised exception at cell list resize");
	throw;
      }
    }

    void updateCurrentValidCell(){
      const int numberParticles = pg->getNumberParticles();
      const bool isCounterUninitialized = (currentValidCell_counter < 0);
      const ullint nextStepMaximumValue = ullint(numberParticles)*(currentValidCell_counter+2);
      constexpr ullint maximumStorableValue = ullint(std::numeric_limits<uint>::max())-1ull;
      const bool nextStepOverflows  = (nextStepMaximumValue >= maximumStorableValue);
      if(isCounterUninitialized or nextStepOverflows){
	currentValidCell = numberParticles;
	currentValidCell_counter = 1;
	const int ncells = grid.getNumberCells();
	const int Nthreads = 512;
	const int Nblocks= ncells/Nthreads + ((ncells%Nthreads)?1:0);
	auto it = thrust::make_counting_iterator<int>(0);
	fillWithGPU<<<Nblocks, Nthreads>>>(thrust::raw_pointer_cast(cellStart.data()),
					   it, 0, ncells);
	CudaCheckError();
      }
      else{
	currentValidCell_counter++;
	currentValidCell = uint(numberParticles)*currentValidCell_counter;
      }
    }

    void updateOrderAndStoreInSortPos(cudaStream_t st){
      const int numberParticles = pg->getNumberParticles();
      auto pos = pd->getPos(access::location::gpu, access::mode::read);
      auto posGroupIterator = pg->getPropertyInputIterator(pos.begin(), access::location::gpu);
      ps.updateOrderByCellHash<Sorter::MortonHash>(posGroupIterator,
						   numberParticles,
						   grid.box, grid.cellDim, st);
      CudaCheckError();
      sortPos.resize(numberParticles);
      ps.applyCurrentOrder(posGroupIterator, thrust::raw_pointer_cast(sortPos.data()), numberParticles, st);
      CudaCheckError();
    }

    void fillCellList(cudaStream_t st){
      const int numberParticles = pg->getNumberParticles();
      const int Nthreads = 512;
      int Nblocks = numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);
      sys->log<System::DEBUG3>("[CellList] fill Cell List, currentValidCell: %d", currentValidCell);
      int h_errorFlag = 0;
      errorFlags.resize(1);
      int *d_errorFlag = thrust::raw_pointer_cast(errorFlags.data());
      CudaSafeCall(cudaMemcpyAsync(d_errorFlag, &h_errorFlag, sizeof(int), cudaMemcpyHostToDevice, st));
      CellList_ns::fillCellList<<<Nblocks, Nthreads, 0, st>>>(thrust::raw_pointer_cast(sortPos.data()),
							      thrust::raw_pointer_cast(cellStart.data()),
							      thrust::raw_pointer_cast(cellEnd.data()),
							      currentValidCell,
							      d_errorFlag,
							      numberParticles,
							      grid);
      CudaSafeCall(cudaMemcpyAsync(&h_errorFlag, d_errorFlag, sizeof(int), cudaMemcpyDeviceToHost, st));
      CudaSafeCall(cudaEventRecord(event, st));
      CudaSafeCall(cudaEventSynchronize(event));
      if(h_errorFlag > 0){
	sys->log<System::ERROR>("[CellList] NaN positions found during construction");
	throw std::overflow_error("CellList encountered NaN positions");
      }
      CudaCheckError();
    }

    void rebuildNeighbourList(cudaStream_t st){
      const int numberParticles = pg->getNumberParticles();
      neighbourList.resize(numberParticles*maxNeighboursPerParticle);
      numberNeighbours.resize(numberParticles);
      rebuildNlist = false;
      int Nthreads=128;
      int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);
      int flag;
      do{
	flag = 0;
	auto neighbourList_ptr = thrust::raw_pointer_cast(neighbourList.data());
	auto numberNeighbours_ptr = thrust::raw_pointer_cast(numberNeighbours.data());
	errorFlags.resize(1);
	int* d_tooManyNeighboursFlag = thrust::raw_pointer_cast(errorFlags.data());
	sys->log<System::DEBUG3>("[CellList] fill Neighbour List");
	CellList_ns::fillNeighbourList<<<Nblocks, Nthreads, 0, st>>>(thrust::raw_pointer_cast(sortPos.data()),
								     this->getGroupIndexIterator(),
								     this->getNeighbourContainer(),
								     neighbourList_ptr, numberNeighbours_ptr,
								     maxNeighboursPerParticle,
								     currentCutOff.x*currentCutOff.x,
								     numberParticles, grid.box,
								     d_tooManyNeighboursFlag);
	CudaSafeCall(cudaMemcpyAsync(&flag, d_tooManyNeighboursFlag, sizeof(int), cudaMemcpyDeviceToHost, st));
	CudaSafeCall(cudaEventRecord(event, st));
	CudaSafeCall(cudaEventSynchronize(event));
	if(flag != 0){
	  this->maxNeighboursPerParticle += 32;
	  sys->log<System::DEBUG>("[CellList] Resizing list to %d neighbours per particle",
				  maxNeighboursPerParticle);
	  int zero = 0;
	  CudaSafeCall(cudaMemcpyAsync(d_tooManyNeighboursFlag, &zero, sizeof(int), cudaMemcpyHostToDevice, st));
	  neighbourList.resize(numberParticles*maxNeighboursPerParticle);
	}
      }while(flag!=0);
    }

  };
}
#endif


