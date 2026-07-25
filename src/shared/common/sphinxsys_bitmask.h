/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	sphinxsys_bitmask.h
 * @brief 	Data container for large vector, e.g. particle data.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef SPHINXSYS_BITMASK_H
#define SPHINXSYS_BITMASK_H

#include "sphinxsys_variable.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace SPH
{
class GroupManager
{
  public:
    GroupManager(DiscreteVariable<UnsignedInt> *group_variable)
        : group_variable_(group_variable), next_bit_(0) {}
    GroupManager(const GroupManager &) = delete;
    GroupManager &operator=(const GroupManager &) = delete;

    std::string GroupVariableName() const { return group_variable_->Name(); }
    DataView<UnsignedInt> GroupVariableDataView() const { return group_variable_->getDataView(); }
    bool hasGroup(const std::string &name) const { return name_to_bit_.count(name) > 0; }
    std::size_t numGroups() const { return name_to_bit_.size(); }

    UnsignedInt registerGroup(const std::string &name)
    {
        if (name_to_bit_.count(name))
        {
            return getGroupMask(name);
        }
        if (next_bit_ >= UnsignedIntBits)
        {
            std::runtime_error(
                "Maximum groups reached (<" + std::to_string(UnsignedIntBits) + ").");
        }
        UnsignedInt bit = (UnsignedInt(1) << next_bit_);
        name_to_bit_[name] = bit;
        ++next_bit_;

        return bit;
    }

    // Helper: get mask for a single group (throws if not found)
    UnsignedInt getGroupMask(const std::string &name) const
    {
        auto it = name_to_bit_.find(name);
        if (it == name_to_bit_.end())
            throw std::runtime_error("Group not found: " + name);
        return it->second;
    }

    // Return the OR of masks for all given group names (union).
    UnsignedInt unionMask(const std::vector<std::string> &names) const
    {
        UnsignedInt combined = 0;
        for (auto &n : names)
            combined |= getGroupMask(n);
        return combined;
    }

    class AddMaskKernel
    {
      public:
        template <class ExecutionPolicy>
        AddMaskKernel(const ExecutionPolicy &ex_policy, GroupManager &manager, UnsignedInt group_mask)
            : group_variable_(manager.group_variable_->DelegatedDataView(ex_policy)),
              group_mask_(group_mask) {}

        void operator()(UnsignedInt index)
        {
            group_variable_[index] |= group_mask_;
        }

      private:
        DataView<UnsignedInt> group_variable_;
        UnsignedInt group_mask_;
    };

    class RemoveMaskKernel
    {
      public:
        template <class ExecutionPolicy>
        RemoveMaskKernel(const ExecutionPolicy &ex_policy, GroupManager &manager, UnsignedInt group_mask)
            : group_variable_(manager.group_variable_->DelegatedDataView(ex_policy)),
              group_mask_(group_mask) {}

        void operator()(UnsignedInt index)
        {
            group_variable_[index] &= ~group_mask_;
        }

      private:
        DataView<UnsignedInt> group_variable_;
        UnsignedInt group_mask_;
    };

    class CheckMaskKernel
    {
      public:
        template <class ExecutionPolicy>
        CheckMaskKernel(const ExecutionPolicy &ex_policy, GroupManager &manager, UnsignedInt group_mask)
            : group_variable_(manager.group_variable_->DelegatedDataView(ex_policy)),
              group_mask_(group_mask) {}

        bool operator()(UnsignedInt index)
        {
            return (group_variable_[index] & group_mask_) != 0;
        }

      private:
        DataView<UnsignedInt> group_variable_;
        UnsignedInt group_mask_;
    };

    template <class ExecutionPolicy>
    AddMaskKernel createAddMaskKernel(const ExecutionPolicy &ex_policy, const std::string &name)
    {
        UnsignedInt group_mask = registerGroup(name);
        return AddMaskKernel(ex_policy, *this, group_mask);
    }

  private:
    DiscreteVariable<UnsignedInt> *group_variable_;
    std::unordered_map<std::string, UnsignedInt> name_to_bit_;
    int next_bit_;
};
} // namespace SPH
#endif // SPHINXSYS_BITMASK_H