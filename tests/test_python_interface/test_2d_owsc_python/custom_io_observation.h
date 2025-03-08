#ifndef CUSTOM_IO_OBSERVATION_H
#define CUSTOM_IO_OBSERVATION_H

#include "io_observation.h"

namespace SPH
{
template <class LocalReduceMethodType>
class ExtendedReducedQuantityRecording : public ReducedQuantityRecording<LocalReduceMethodType>
{
public:
    // Inherit constructors from the base class
    using ReducedQuantityRecording<LocalReduceMethodType>::ReducedQuantityRecording;

    // Function to directly return the result of reduce_method_.exec()
    typename LocalReduceMethodType::ReturnType getReducedQuantity() 
    {
        return this->reduce_method_.exec();
    }
};

} // namespace SPH 
#endif // CUSTOM_IO_OBSERVATION_H
