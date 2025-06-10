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
      dv_variable2_(this->particles_->template getVariableByName<DataType>(variable_name2)) {}
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
    // the transform matrix should be symmetric and positive definite
    DataType result = this->variable_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        DataType pair_difference = this->variable_[index_i] - this->variable_[index_j];
        result -= 2.0 * this->dis_coeff_(index_i, index_j) * pair_difference *
                  this->dW_ij(index_i, index_j) * this->Vol_[index_j] * dt /
                  this->vec_r_ij(index_i, index_j).norm();
    }
    transformed_[index_i] = result;
}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
DissipativeTransform<Contact<Dirichlet<DissipationType>, Parameters...>>::
    DissipativeTransform(Contact<Parameters...> &contact_relation, const std::string &variable_name)
    : BaseDissipationType(contact_relation, variable_name),
      dv_transformed_(this->particles_->template registerStateVariableOnly<DataType>(
          "Transformed" + variable_name))
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        contact_dv_variable_.push_back(
            this->contact_particles_[k]->template getVariableByName<DataType>(variable_name));
    }
}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
DissipativeTransform<Contact<Dirichlet<DissipationType>, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseDissipationType::InteractKernel(ex_policy, encloser, contact_index),
      transformed_(encloser.dv_transformed_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_variable_(encloser.contact_dv_variable_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DissipationType, typename... Parameters>
void DissipativeTransform<Contact<Dirichlet<DissipationType>, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    DataType result = this->zero_value_;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        DataType pair_difference = this->variable_[index_i] - contact_variable_[index_j];
        result -= 2.0 * this->dis_coeff_(index_i, index_j) * pair_difference *
                  this->dW_ij(index_i, index_j) * contact_Vol_[index_j] * dt /
                  this->vec_r_ij(index_i, index_j).norm();
    }
    transformed_[index_i] += result;
}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
DissipationRHS<DataType, DynamicsIdentifier>::
    DissipationRHS(DynamicsIdentifier &identifier, const std::string &variable_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)),
      dv_rhs_(this->particles_->template registerStateVariableOnly<DataType>("RHS" + variable_name)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
DissipationRHS<DataType, DynamicsIdentifier>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
      rhs_(encloser.dv_rhs_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
FullDissipationResidue<DataType, DynamicsIdentifier>::
    FullDissipationResidue(DynamicsIdentifier &identifier, const std::string &variable_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_residue_(this->particles_->template registerStateVariableOnly<DataType>("Residue" + variable_name)),
      dv_rhs_(this->particles_->template getVariableByName<DataType>("RHS" + variable_name)),
      dv_transformed_(this->particles_->template getVariableByName<DataType>("Transformed" + variable_name)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
FullDissipationResidue<DataType, DynamicsIdentifier>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : residue_(encloser.dv_residue_->DelegatedData(ex_policy)),
      rhs_(encloser.dv_rhs_->DelegatedData(ex_policy)),
      transformed_(encloser.dv_transformed_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
UpdateDissipationSearchDirection<DataType, DynamicsIdentifier>::
    UpdateDissipationSearchDirection(DynamicsIdentifier &identifier,
                                     const std::string &variable_name,
                                     SingularVariable<TensorProductType<DataType>> *sv_residue_ratio)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_residue_(this->particles_->template getVariableByName<DataType>("Residue" + variable_name)),
      dv_search_direction_(
          this->particles_->template getVariableByName<DataType>("SearchDirection" + variable_name)),
      sv_residue_ratio_(sv_residue_ratio) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
UpdateDissipationSearchDirection<DataType, DynamicsIdentifier>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : residue_(encloser.dv_residue_->DelegatedData(ex_policy)),
      search_direction_(encloser.dv_search_direction_->DelegatedData(ex_policy)),
      residue_ratio_(encloser.sv_residue_ratio_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
UpdateDissipationSolution<DataType, DynamicsIdentifier>::
    UpdateDissipationSolution(DynamicsIdentifier &identifier,
                              const std::string &variable_name,
                              SingularVariable<TensorProductType<DataType>> *sv_search_depth)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_search_direction_(this->particles_->template registerStateVariableOnly<DataType>(
          "SearchDirection" + variable_name)),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)),
      sv_search_depth_(sv_search_depth) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
UpdateDissipationSolution<DataType, DynamicsIdentifier>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : search_direction_(encloser.dv_search_direction_->DelegatedData(ex_policy)),
      variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
      search_depth_(encloser.sv_search_depth_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename... Parameters>
ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, Parameters...>>::
    ImplicitDissipation(RelationType<Parameters...> &first_relation,
                        const std::string &variable_name, Real convergence_criteria)
    : BaseDynamics<void>(), variable_name_(variable_name),
      convergence_criteria_(convergence_criteria),
      identifier_(first_relation.getDynamicsIdentifier()),
      sv_search_depth_("SearchDepth" + variable_name),
      sv_residue_ratio_("ResidueRatio" + variable_name),
      dissipation_rhs_(identifier_, variable_name),
      transformed_variable_(first_relation, variable_name),
      full_residue_(identifier_, variable_name),
      initial_search_direction_(
          identifier_, "SearchDirection" + variable_name, "Residue" + variable_name),
      transformed_search_direction_(first_relation, "SearchDirection" + variable_name),
      update_search_direction_(identifier_, variable_name, &sv_residue_ratio_),
      residue_average_(identifier_, variable_name),
      transformed_search_direction_average_(identifier_, variable_name),
      update_dissipation_solution_(identifier_, variable_name, &sv_search_depth_) {}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename... Parameters>
template <typename... ControlParameters, typename... RelationParameters, typename... Args>
auto &ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, Parameters...>>::
    addPostContactInteraction(Contact<RelationParameters...> &contact_relation, Args &&...args)
{
    for (auto &identifier : contact_relation.getContactIdentifiers())
    {
        initialization_methods_.push_back(
            initialization_ptrs_keeper_.createPtr<
                StateDynamics<ExecutionPolicy, VariableAssignment<SPHBody, ConstantValue<DataType>>>>(
                *identifier, "SearchDirection" + variable_name_, ZeroData<DataType>::value));
    }
    transformed_variable_.template addPostContactInteraction<ControlParameters...>(
        contact_relation, std::forward<Args>(args)...);
    transformed_search_direction_.template addPostContactInteraction<ControlParameters...>(
        contact_relation, std::forward<Args>(args)...);
    return *this;
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename... Parameters>
void ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, Parameters...>>::
    initializeImplicitDissipation()
{
    for (auto &initialization : initialization_methods_)
    {
        initialization->exec();
    }
}
//=================================================================================================//
template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename... Parameters>
void ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, Parameters...>>::
    exec(Real dt)
{
    // initial residue and search direction
    dissipation_rhs_.exec();
    transformed_variable_.exec(dt);
    full_residue_.exec();
    TensorProductType<DataType> present_residue_average = residue_average_.exec();
    Real residue_norm = math::sqrt(getSquaredNorm(present_residue_average));
    initial_search_direction_.exec();

    UnsignedInt iteration_count = 0;
    while (residue_norm > convergence_criteria_)
    {
        transformed_search_direction_.exec(dt);
        TensorProductType<DataType> transformed_average = transformed_search_direction_average_.exec();
        sv_search_depth_.setValue(present_residue_average * getInverse(transformed_average));
        update_dissipation_solution_.exec();

        // update residue and search direction
        TensorProductType<DataType> previous_residue_average = present_residue_average;
        transformed_variable_.exec(dt);
        full_residue_.exec();
        present_residue_average = residue_average_.exec();
        residue_norm = math::sqrt(getSquaredNorm(present_residue_average));
        sv_residue_ratio_.setValue(present_residue_average * getInverse(previous_residue_average));
        update_search_direction_.exec();

        ++iteration_count;
        if (iteration_count > 1000)
        {
            std::cout << "Implicit dissipation did not converge within 1000 iterations." << std::endl;
            break;
        }
    }
}
//=================================================================================================//
} // namespace SPH
#endif // IMPLICIT_DISSIPATION_HPP
