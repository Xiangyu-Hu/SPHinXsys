#pragma once

#include "udf_P_refinement.h"

namespace SPH
{

namespace fluid_dynamics
{

namespace udf
{

    template <class DataDelegationType>
    template <class BaseRelationType>
    P_refinement_GetVelocityGradient<DataDelegationType>::
        P_refinement_GetVelocityGradient(BaseRelationType& base_relation)
        : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
        velocity_gradient_only_P_(this->particles_->template registerStateVariableData<Matd>("VelocityGradientInnerOnlyP")),
        k_gradient_only_P_(this->particles_->template registerStateVariableData<Vecd>("TurbulentKineticEnergyGradientOnlyP")),
        omega_gradient_only_P_(this->particles_->template registerStateVariableData<Vecd>("TurbulentSpecificDissipationGradientOnlyP")),
        Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
        vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
        is_near_wall_P1_(this->particles_->template getVariableDataByName<int>("IsNearWallP1")),
        is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
        turbu_k_(this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
        turbu_omega_(this->particles_->template getVariableDataByName<Real>("TurbulentSpecificDissipation")) {}

}

}

}
