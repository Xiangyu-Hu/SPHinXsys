#include "sph_system.h"

#include "all_body_relations.h"
#include "elastic_dynamics.h"

namespace SPH
{
//=================================================================================================//
SimulationContext &SPHSystem::getSimulationContext()
{
    if (simulation_ctx_ptr_keeper_.getPtr() == nullptr)
    {
        std::cout << "\n Error: Simulation context not setup yet! \n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *simulation_ctx_ptr_keeper_.getPtr();
}
//=================================================================================================//
void SPHSystem::initializeSystemCellLinkedLists()
{
    for (auto &body : real_bodies_)
    {
        DynamicCast<RealBody>(this, body)->updateCellLinkedList();
    }
}
//=================================================================================================//
void SPHSystem::initializeSystemConfigurations()
{
    for (auto &body : sph_bodies_)
    {
        StdVec<SPHRelation *> &body_relations = body->getBodyRelations();
        for (size_t i = 0; i < body_relations.size(); i++)
        {
            body_relations[i]->updateConfiguration();
        }
    }
}
//=================================================================================================//
Real SPHSystem::getSmallestTimeStepAmongSolidBodies(Real CFL)
{
    Real dt = MaxReal;
    for (size_t i = 0; i < solid_bodies_.size(); i++)
    {
        ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(*solid_bodies_[i], CFL);
        Real dt_temp = computing_time_step_size.exec();
        if (dt_temp < dt)
            dt = dt_temp;
    }
    return dt;
}
//=================================================================================================//
} // namespace SPH
