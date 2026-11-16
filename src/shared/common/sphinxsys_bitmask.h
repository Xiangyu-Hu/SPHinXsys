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
 * @brief 	tbd.
 * @author	Xiangyu Hu
 */
#ifndef SPHINXSYS_BITMASK_H
#define SPHINXSYS_BITMASK_H

#include "sphinxsys_variable.h"

#include <string>
#include <unordered_map>

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
    UnsignedInt numGroups() const { return name_to_bit_.size(); }
    bool hasDerivedGroup(const std::string &derived_name) const { return derived_groups_.count(derived_name) > 0; }

    UnsignedInt registerGroup(const std::string &name)
    {
        if (name_to_bit_.count(name))
        {
            return getGroupMask(name);
        }
        if (next_bit_ >= UnsignedIntBits)
        {
            throw std::runtime_error(
                "Maximum groups reached (<" + std::to_string(UnsignedIntBits) + ").");
        }
        UnsignedInt bit = (UnsignedInt(1) << next_bit_);
        name_to_bit_[name] = bit;
        ++next_bit_;

        return bit;
    }

    UnsignedInt getGroupMask(const std::string &name) const
    {
        auto it = name_to_bit_.find(name);
        if (it == name_to_bit_.end())
            throw std::runtime_error("Group not found: " + name);
        return it->second;
    }

    void registerDerivedGroup(const std::string &derived_name)
    {
        if (!hasDerivedGroup(derived_name))
        {
            derived_groups_[derived_name] = 0;
        }
    }

    void addToDerivedGroup(const std::string &derived_name, const std::string &name)
    {
        validateDerivedGroup(derived_name);
        derived_groups_[derived_name] |= getGroupMask(name);
    }

    void removeFromDerivedGroup(const std::string &derived_name, const std::string &name)
    {
        validateDerivedGroup(derived_name);
        derived_groups_[derived_name] &= ~getGroupMask(name);
    }

    void setDerivedGroup(const std::string &derived_name, UnsignedInt mask)
    {
        validateDerivedGroup(derived_name);
        derived_groups_[derived_name] = mask;
    }

    UnsignedInt getDerivedGroup(const std::string &derived_name) const
    {
        validateDerivedGroup(derived_name);
        return derived_groups_.at(derived_name);
    }

    class MaskKernel
    {
      public:
        template <class ExecutionPolicy>
        MaskKernel(const ExecutionPolicy &ex_policy, GroupManager &manager, UnsignedInt mask)
            : group_variable_(manager.group_variable_->DelegatedDataView(ex_policy)),
              mask_(mask) {}

        void add(UnsignedInt index)
        {
            group_variable_[index] |= mask_;
        }

        void remove(UnsignedInt index)
        {
            group_variable_[index] &= ~mask_;
        }

        void set(UnsignedInt index)
        {
            group_variable_[index] = mask_;
        }

        bool check(UnsignedInt index) const
        {
            return (group_variable_[index] & mask_) != 0;
        }

      private:
        DataView<UnsignedInt> group_variable_;
        UnsignedInt mask_;
    };

    MaskKernel createHostMaskKernel(const std::string &name)
    {
        UnsignedInt mask = registerGroup(name);
        return MaskKernel(seq, *this, mask);
    }

  private:
    DiscreteVariable<UnsignedInt> *group_variable_;
    std::unordered_map<std::string, UnsignedInt> name_to_bit_;
    std::unordered_map<std::string, UnsignedInt> derived_groups_;
    int next_bit_;

    void validateDerivedGroup(const std::string &derived_name) const
    {
        if (!hasDerivedGroup(derived_name))
            throw std::runtime_error("Derived group not found: " + derived_name);
    }
};
} // namespace SPH
#endif // SPHINXSYS_BITMASK_H