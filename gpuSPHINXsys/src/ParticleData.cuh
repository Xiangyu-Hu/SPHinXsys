#ifndef PARTICLEDATA_CUH
#define PARTICLEDATA_CUH
#include"System.h"

#include"Property.cuh"
#include"ParticleSorter.cuh"
#include"vector.cuh"

#include<include/nod/nod.hpp>
#include<boost/preprocessor.hpp>
#include<boost/preprocessor/stringize.hpp>
#include<boost/preprocessor/seq/for_each.hpp>
#include<boost/preprocessor/tuple/elem.hpp>
#include<thrust/device_vector.h>
#include <thrust/system_error.h>

//List here all the properties with this syntax:
/*       ((PropertyName, propertyName, TYPE))				\      */
//The preprocessor ensures that they are included wherever is needed
#define ALL_PROPERTIES_LIST ((Pos, pos, real4))        \
                            ((Id, id, int))            \
                            ((Mass, mass, real))       \
                            ((Vol, vol, real))         \
                            ((Force, force, real4))    \
                            ((Energy, energy, real))   \
                            ((Vel, vel, real3))        \
                            ((Vel_tv, vel_tv, real3))  \
                            ((F_Pb, f_Pb, real3))      \
                            ((Rho, rho, real))         \
                            ((Rho0, rho0, real))       \
                            ((Drho, drho, real))       \
                            ((Sigma0, sigma0, real))   \
                            ((Pressure, pressure, real))

namespace gpu{

  template<class T>
  using signal = typename nod::unsafe_signal<T>;

  using connection = nod::connection;

  //Get the Name (first letter capital) from a tuple in the property list
#define PROPNAME_CAPS(tuple) BOOST_PP_TUPLE_ELEM(3, 0 ,tuple)
  //Get the name (no capital) from a tuple in the property list
#define PROPNAME(tuple) BOOST_PP_TUPLE_ELEM(3, 1 ,tuple)
  //Get the type from a tuple in the property list
#define PROPTYPE(tuple) BOOST_PP_TUPLE_ELEM(3, 2 ,tuple)

//This macro iterates through all properties applying some macro
#define PROPERTY_LOOP(macro)  BOOST_PP_SEQ_FOR_EACH(macro, _, ALL_PROPERTIES_LIST)



  class ParticleData{
  public:
    //Hints to ParticleData about how to perform different task. Mainly how to sort the particles.
    struct Hints{
      bool orderByHash = false;
      Box hash_box = Box(make_real3(128));
      real3 hash_cutOff = make_real3(10.0);
      bool orderByType = false;

    };

  private:
    shared_ptr<System> sys;
#define DECLARE_PROPERTIES_T(type, name) Property<type> name;
#define DECLARE_PROPERTIES(r,data, tuple) DECLARE_PROPERTIES_T(PROPTYPE(tuple), PROPNAME(tuple))

    //Declare all property containers
    PROPERTY_LOOP(DECLARE_PROPERTIES)

    int numberParticles;
    shared_ptr<signal<void(void)>> reorderSignal = std::make_shared<signal<void(void)>>();
    shared_ptr<signal<void(int)>> numParticlesChangedSignal = std::make_shared<signal<void(int)>>();

//Declare write access signals for all properties
#define DECLARE_SIGNAL_PROPERTIES_T(type, name) shared_ptr<signal<void(void)>> BOOST_PP_CAT(name,WriteRequestedSignal = std::make_shared<signal<void(void)>>();)
#define DECLARE_SIGNAL_PROPERTIES(r,data, tuple) DECLARE_SIGNAL_PROPERTIES_T(PROPTYPE(tuple), PROPNAME(tuple))
    //Declare all property write signals
    PROPERTY_LOOP(DECLARE_SIGNAL_PROPERTIES)



    std::shared_ptr<ParticleSorter> particle_sorter;
    thrust::host_vector<int> originalOrderIndexCPU;
    bool originalOrderIndexCPUNeedsUpdate;
    Hints hints;

  public:
    ParticleData() = delete;
    ParticleData(int numberParticles, shared_ptr<System> sys);
    ~ParticleData(){
      sys->log<System::DEBUG>("[ParticleData] Destroyed");
    }


    //Generate getters for all properties except ID
#define GET_PROPERTY_T(Name,name)  GET_PROPERTY_R(Name,name)
#define GET_PROPERTY_R(Name, name)					\
  inline auto get ## Name(access::location dev, access::mode mode) -> decltype(name.data(dev,mode)){ \
    if(!name.isAllocated()) name.resize(numberParticles);		\
    if(!name.isAllocated() or mode==access::mode::write or mode==access::mode::readwrite){ \
    (*name ## WriteRequestedSignal)();				\
    }									\
      return name.data(dev,mode);	                            	\
    }									\

#define GET_PROPERTY(r, data, tuple) GET_PROPERTY_T(PROPNAME_CAPS(tuple), PROPNAME(tuple))

    //Define getProperty() functions for all properties in list
    PROPERTY_LOOP(GET_PROPERTY)

        //Generate getters for all properties except ID
#define GET_PROPERTY_IF_ALLOC_T(Name,name)  GET_PROPERTY_IF_ALLOC_R(Name,name)
#define GET_PROPERTY_IF_ALLOC_R(Name, name)					\
    inline auto get ## Name ## IfAllocated(access::location dev, access::mode mode) -> decltype(name.data(dev,mode)){ \
      if(!name.isAllocated()){                    \
	decltype(name.data(dev,mode)) tmp;        \
	return tmp;	                          \
      }						  \
      return this->get ## Name(dev,mode);	  \
    }						  \

#define GET_PROPERTY_IF_ALLOC(r, data, tuple) GET_PROPERTY_IF_ALLOC_T(PROPNAME_CAPS(tuple), PROPNAME(tuple))

    //Define getProperty() functions for all properties in list
    PROPERTY_LOOP(GET_PROPERTY_IF_ALLOC)

    //Generate isPropAllocated for all properties
#define IS_ALLOCATED_T(Name, name) IS_ALLOCATED_R(Name, name)
#define IS_ALLOCATED_R(Name, name)					\
    inline bool is##Name##Allocated(){return name.isAllocated();}	\

#define IS_ALLOCATED(r, data, tuple) IS_ALLOCATED_T(PROPNAME_CAPS(tuple), PROPNAME(tuple))

    PROPERTY_LOOP(IS_ALLOCATED)

    void sortParticles();

    const int * getIdOrderedIndices(access::location dev){
      sys->log<System::DEBUG5>("[ParticleData] Id order requested for %d (0=cpu, 1=gpu)", dev);
      auto id = getId(access::location::gpu, access::mode::read);
      int *sortedIndex = particle_sorter->getIndexArrayById(id.raw(), numberParticles);
      sys->log<System::DEBUG6>("[ParticleData] Id reorder completed.");
      if(dev == access::location::gpu){
	return sortedIndex;
      }
      else{
	if(originalOrderIndexCPUNeedsUpdate){
	  sys->log<System::DEBUG1>("[ParticleData] Updating CPU original order array");
	  originalOrderIndexCPU.resize(numberParticles);
	  int * sortedIndexCPU = thrust::raw_pointer_cast(originalOrderIndexCPU.data());
	  CudaSafeCall(cudaMemcpy(sortedIndexCPU,
				  sortedIndex,
				  numberParticles*sizeof(int),
				  cudaMemcpyDeviceToHost));
	  originalOrderIndexCPUNeedsUpdate = false;
	  return sortedIndexCPU;
	}
	else{
	  return thrust::raw_pointer_cast(originalOrderIndexCPU.data());
	}
      }
    }

    template<class InputIterator, class OutputIterator>
    void applyCurrentOrder(InputIterator in, OutputIterator out, int numElements){
      particle_sorter->applyCurrentOrder(in, out, numElements);
    }

    const int * getCurrentOrderIndexArray(){
      return particle_sorter->getSortedIndexArray(numberParticles);
    }

    void changeNumParticles(int Nnew);

    int getNumParticles(){ return this->numberParticles;}

    shared_ptr<signal<void(void)>> getReorderSignal(){
      sys->log<System::DEBUG>("[ParticleData] Reorder signal requested");
      return this->reorderSignal;
    }

#define GET_PROPERTY_SIGNAL_T(Name,name)  GET_PROPERTY_SIGNAL_R(Name,name)
#define GET_PROPERTY_SIGNAL_R(Name, name)				\
    inline shared_ptr<signal<void(void)>> get ## Name ## WriteRequestedSignal(){ \
      return this->name ## WriteRequestedSignal;			\
    }
#define GET_PROPERTY_SIGNAL(r, data, tuple) GET_PROPERTY_SIGNAL_T(PROPNAME_CAPS(tuple), PROPNAME(tuple))
    PROPERTY_LOOP(GET_PROPERTY_SIGNAL)

    void emitReorder(){
      sys->log<System::DEBUG>("[ParticleData] Emitting reorder signal...");
      (*this->reorderSignal)();
    }

    shared_ptr<signal<void(int)>> getNumParticlesChangedSignal(){
      return this->numParticlesChangedSignal;
    }


    void hintSortByHash(Box hash_box, real3 hash_cutOff){
      hints.orderByHash = true;
      hints.hash_box = hash_box;
      hints.hash_cutOff = hash_cutOff;

    }

  private:

    void emitNumParticlesChanged(int Nnew){
      (*numParticlesChangedSignal)(Nnew);
    }

  };


#define INIT_PROPERTIES_T(NAME, name) ,  name(BOOST_PP_STRINGIZE(NAME), sys)
#define INIT_PROPERTIES(r,data, tuple) INIT_PROPERTIES_T(PROPNAME_CAPS(tuple), PROPNAME(tuple))

  ParticleData::ParticleData(int numberParticles, shared_ptr<System> sys):
    numberParticles(numberParticles),
    originalOrderIndexCPUNeedsUpdate(true),
    sys(sys)
    PROPERTY_LOOP(INIT_PROPERTIES)
  {
    sys->log<System::MESSAGE>("[ParticleData] Created with %d particles.", numberParticles);

    id.resize(numberParticles);
    CudaCheckError();

    auto id_prop = id.data(access::location::gpu, access::mode::write);

    cub::CountingInputIterator<int> ci(0);
    thrust::copy(thrust::cuda::par,
		 ci, ci + numberParticles,
		 id_prop.begin());

    particle_sorter = std::make_shared<ParticleSorter>(sys);
  }

  //Sort the particles to improve a certain kind of access pattern.
  void ParticleData::sortParticles(){
    sys->log<System::DEBUG>("[ParticleData] Sorting particles...");

    {
      auto posPtr     = pos.data(access::gpu, access::read);
      if(hints.orderByHash || !hints.orderByType){
	int3 cellDim = make_int3(hints.hash_box.boxSize/hints.hash_cutOff);
	particle_sorter->updateOrderByCellHash(posPtr.raw(), numberParticles, hints.hash_box, cellDim);
      }

    }
  //This macro reorders to the newest order a property given its name
#define APPLY_CURRENT_ORDER(r, data, tuple) APPLY_CURRENT_ORDER_R(PROPNAME(tuple))
#define APPLY_CURRENT_ORDER_R(name) {					\
      if(name.isAllocated()){						\
	auto devicePtr     = name.data(access::gpu, access::write);	\
	auto device_altPtr = name.getAltGPUBuffer();			\
	particle_sorter->applyCurrentOrder(devicePtr.raw(), device_altPtr, numberParticles); \
	name.swapInternalBuffers();						\
      }									\
    }
    //Apply current order to all allocated properties. See APPLY_CURRENT_ORDER macro
    PROPERTY_LOOP(APPLY_CURRENT_ORDER)

    originalOrderIndexCPUNeedsUpdate = true;
    this->emitReorder();
  }



  void ParticleData::changeNumParticles(int Nnew){
    sys->log<System::CRITICAL>("[ParticleData] CHANGE PARTICLES FUNCTIONALITY NOT IMPLEMENTED YET!!!");
    sys->log<System::DEBUG>("[ParticleData] Adding/Removing particles...");
    this->numberParticles = Nnew;
    pos.resize(Nnew);
#define RESIZE_PROPERTY_R(name) {if(this->name.isAllocated()){this->name.resize(this->numberParticles);}}
#define RESIZE_PROPERTY(r, data, tuple) RESIZE_PROPERTY_R(PROPNAME(tuple))

    PROPERTY_LOOP(RESIZE_PROPERTY)

    originalOrderIndexCPUNeedsUpdate = true;
    this->emitNumParticlesChanged(Nnew);
  }
}

#undef ALL_PROPERTIES_LIST
#undef PROPNAME_CAPS
#undef PROPNAME
#undef PROPTYPE
#undef PROPERTY_LOOP
#undef DECLARE_PROPERTIES_T
#undef DECLARE_PROPERTIES
#undef DECLARE_SIGNAL_PROPERTIES_T
#undef DECLARE_SIGNAL_PROPERTIES
#undef GET_PROPERTY_T
#undef GET_PROPERTY_R
#undef GET_PROPERTY
#undef GET_PROPERTY_SIGNAL_T
#undef GET_PROPERTY_SIGNAL_R
#undef GET_PROPERTY_SIGNAL
#undef IS_ALLOCATED_T
#undef IS_ALLOCATED_R
#undef IS_ALLOCATED
#undef GET_PROPERTY_IF_ALLOC
#undef GET_PROPERTY_IF_ALLOC_T
#undef GET_PROPERTY_IF_ALLOC_R
#undef APPLY_CURRENT_ORDER
#undef APPLY_CURRENT_ORDER_R
#undef RESIZE_PROPERTY_R
#undef RESIZE_PROPERTY




#endif
