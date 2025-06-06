#ifndef IMPLICIT_DISSIPATION_HPP
#define IMPLICIT_DISSIPATION_HPP

#include "implicit_dissipation.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
QuantityTensorProductAverage<DataType, DynamicsIdentifier>::
    QuantityTensorProductAverage(DynamicsIdentifier &identifier,
                                 const std::string &variable_name1,
                                 const std::string &variable_name2)
    : BaseReduceDynamics(identifier),
      dv_variable1_(this->particles_->template getVariableByName<DataType>(variable_name1)),
      dv_variable2_(this->particles_->template getVariableByName<DataType>(variable_name1)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
QuantityTensorProductAverage<DataType, DynamicsIdentifier>::ReduceKernel::
    ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable1_(encloser.dv_variable1_->DelegatedData(ex_policy)),
      variable2_(encloser.dv_variable2_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
DissipativeTransform<Inner<DissipationType, Parameters...>>::
    DissipativeTransform(Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDissipationType(inner_relation, variable_name),
      dv_transformed_(this->particles_->template registerStateVariableOnly<DataType>(
          "Transformed" + variable_name)) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
DissipativeTransform<Inner<DissipationType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDissipationType::InteractKernel(ex_policy, encloser),
      transformed_(encloser.dv_transformed_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
void DissipativeTransform<Inner<DissipationType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    // compute the error and parameters
    DataType result = -this->variable_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        DataType pair_difference = this->variable_[index_i] - this->variable_[index_j];
        result += 2.0 * this->dis_coeff_(index_i, index_j) * pair_difference *
                  this->dW_ij(index_i, index_j) * this->Vol_[index_j] * dt /
                  this->vec_r_ij(index_i, index_j).norm();
    }
    transformed_[index_i] = result;
}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
DissipationRHS<DataType, DynamicsIdentifier>::
    DissipationRHS(DynamicsIdentifier &identifier, const std::string &variable_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_variable_(this->particles_->template registerStateVariableOnly<DataType>(variable_name)),
      dv_old_state_(this->particles_->template registerStateVariableOnlyFrom<DataType>(
          "Old" + variable_name, variable_name)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
DissipationRHS<DataType, DynamicsIdentifier>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
      old_state_(encloser.dv_old_state_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
void DissipationRHS<DataType, DynamicsIdentifier>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    old_state_[index_i] = variable_[index_i];
}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
FullDissipationResidue<DataType, DynamicsIdentifier>::
    FullDissipationResidue(DynamicsIdentifier &identifier, const std::string &variable_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_residue_(this->particles_->template registerStateVariableOnly<DataType>("Residue" + variable_name)),
      dv_old_state_(this->particles_->template getVariableByName<DataType>("Old" + variable_name)),
      dv_transformed_(this->particles_->template getVariableByName<DataType>("Transformed" + variable_name)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
FullDissipationResidue<DataType, DynamicsIdentifier>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : residue_(encloser.dv_residue_->DelegatedData(ex_policy)),
      old_state_(encloser.dv_old_state_->DelegatedData(ex_policy)),
      transformed_(encloser.dv_transformed_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
void FullDissipationResidue<DataType, DynamicsIdentifier>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    residue_[index_i] = old_state_[index_i] - transformed_[index_i];
}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
UpdateDissipationResidue<DataType, DynamicsIdentifier>::
    UpdateDissipationResidue(DynamicsIdentifier &identifier,
                             const std::string &variable_name,
                             SingularVariable<TensorProductType<DataType>> *sv_search_depth)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_residue_(this->particles_->template registerStateVariableOnly<DataType>("Residue" + variable_name)),
      dv_transformed_(this->particles_->template getVariableByName<DataType>(
          "Transformed" + dv_residue_->Name())),
      sv_search_depth_(sv_search_depth) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
UpdateDissipationResidue<DataType, DynamicsIdentifier>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : residue_(encloser.dv_residue_->DelegatedData(ex_policy)),
      transformed_(encloser.dv_transformed_->DelegatedData(ex_policy)),
      search_depth_(encloser.sv_search_depth_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
void UpdateDissipationResidue<DataType, DynamicsIdentifier>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    residue_[index_i] -= (*search_depth_) * transformed_[index_i];
}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
UpdateDissipationSolution<DataType, DynamicsIdentifier>::
    UpdateDissipationSolution(DynamicsIdentifier &identifier,
                              const std::string &variable_name,
                              SingularVariable<TensorProductType<DataType>> *sv_search_depth)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_residue_(this->particles_->template registerStateVariableOnly<DataType>("Residue" + variable_name)),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)),
      sv_search_depth_(sv_search_depth) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
UpdateDissipationSolution<DataType, DynamicsIdentifier>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : residue_(encloser.dv_residue_->DelegatedData(ex_policy)),
      variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
      search_depth_(encloser.sv_search_depth_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
void UpdateDissipationSolution<DataType, DynamicsIdentifier>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    variable_[index_i] += (*search_depth_) * residue_[index_i];
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename... Parameters>
ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, Parameters...>>::
    ImplicitDissipation(RelationType<Parameters...> &first_relation,
                        const std::string &variable_name, Real convergence_criteria)
    : BaseDynamics<void>(),
      convergence_criteria_(convergence_criteria),
      dynamics_identifier_(first_relation.getDynamicsIdentifier()),
      sv_search_depth_("SearchDepth" + variable_name),
      dissipation_rhs_(dynamics_identifier_, variable_name),
      transformed_variable_(first_relation, variable_name),
      full_residue_(dynamics_identifier_, variable_name),
      transformed_residue_(first_relation, "Residue" + variable_name),
      update_residue_(dynamics_identifier_, variable_name, &sv_search_depth_),
      residue_average_(dynamics_identifier_, variable_name),
      transformed_residue_average_(dynamics_identifier_, variable_name),
      update_dissipation_solution_(dynamics_identifier_, variable_name, &sv_search_depth_) {}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename... Parameters>
template <typename... ControlParameters, typename... RelationParameters, typename... Args>
auto &ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, Parameters...>>::
    addContactInteraction(Contact<RelationParameters...> &contact_relation, Args &&...args)
{
    transformed_variable_.addContactInteraction(contact_relation, std::forward<Args>(args)...);
    transformed_residue_.addContactInteraction(contact_relation, std::forward<Args>(args)...);
    return *this;
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename... Parameters>
void ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, Parameters...>>::
    exec(Real dt)
{
    Real residue_norm = MaxReal;
    UnsignedInt iteration_count = 0;
    dissipation_rhs_.exec();
    while (residue_norm > convergence_criteria_)
    {
        if (iteration_count % 10 == 0)
        {
            transformed_variable_.exec(dt);
            full_residue_.exec();
        }
        else
        {
            update_residue_.exec();
        }

        TensorProductType<DataType> tensor_residue_average = residue_average_.exec();
        residue_norm = math::sqrt(getSquaredNorm(tensor_residue_average));
        transformed_residue_.exec(dt);
        sv_search_depth_.setValue(tensor_residue_average * getInverse(transformed_residue_average_.exec()));
        update_dissipation_solution_.exec();
        ++iteration_count;
        std::cout << "Iteration: " << iteration_count << ", Residue Norm: " << residue_norm << std::endl;
        if (iteration_count > 100)
        {
            std::cout << "Implicit dissipation did not converge within 100 iterations." << std::endl;
            break;
        }
    }
}
//=================================================================================================//
} // namespace SPH
#endif // IMPLICIT_DISSIPATION_HPP
