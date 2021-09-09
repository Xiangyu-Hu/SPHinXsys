#ifndef MY_MEMORY_POOL_H
#define MY_MEMORY_POOL_H

#define TBB_PREVIEW_MEMORY_POOL 1

#include "tbb/memory_pool.h"
#include "tbb/enumerable_thread_specific.h"

#include <list>

using namespace tbb;
//-------------------------------------------------------------------------------------------------
//my memory pool
//-------------------------------------------------------------------------------------------------
template<class T>
class MyMemoryPool {
	T sample;

#ifdef __EMSCRIPTEN__
	std::allocator<T> my_pool;				//memory pool
	std::list<T> data_list;					//list of all nodes allocated
#else
	tbb::memory_pool< std::allocator<T> > my_pool;				//memory pool
	typedef tbb::memory_pool_allocator<T> pool_allocator_t;		//memory allocator
	std::list<T, pool_allocator_t> data_list;					//list of all nodes allocated
#endif
	std::list<T*> free_list;									//list of all free nodes

public:

	//constructor
#ifdef __EMSCRIPTEN__
	MyMemoryPool() : data_list(my_pool) {};
#else
	MyMemoryPool() : data_list(pool_allocator_t(my_pool)) {};
#endif
	//deconstructor
	~MyMemoryPool() {
		//my_pool.recycle();
	};
	//prepare an avaliable node
	T* malloc()
	{
		if (free_list.empty()) {
			data_list.push_back(sample);
			return (&data_list.back());
		}
		else
		{
			T* result = free_list.front();
			free_list.pop_front();
			return result;
		}
	};
	//relinquish an unused node
	void free(T* ptr)
	{
		free_list.push_back(ptr);
	};
	//return the total number of nodes allocated
	int capicity()
	{
		return data_list.size();
	};
	//return the number of current available nodes
	int available_node()
	{
		return free_list.size();
	};
};

#endif //MY_MEMORY_POOL_H
