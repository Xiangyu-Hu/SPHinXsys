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
 * @file sphinxsys_entity.h
 * @brief Here gives the definitions of the entity and component system.
 * @details These variables are those discretized in spaces and time.
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_ENTITY_H
#define SPHINXSYS_ENTITY_H

#include "base_data_type.h"
#include "ownership.h"

#include <memory>
#include <typeindex>
#include <unordered_map>

namespace SPH
{

class Entity
{
  public:
    explicit Entity(const std::string &name) : name_(name) {};
    ~Entity() {};
    std::string Name() const { return name_; };

  protected:
    const std::string name_;
};

class EntityManager
{
    std::unordered_map<std::type_index, std::unordered_map<size_t, Entity *>> all_entities_;
    size_t entity_count_ = 0;

  public:
    EntityManager() = default;
    ~EntityManager() {};

    template <typename T>
    T *addEntity(T *entity)
    {
        static_assert(std::is_base_of<Entity, T>::value, "T must derive from Entity");
        T *existing_entity = findEntityByName<T>(entity->Name());
        if (existing_entity == nullptr)
        {
            entity_count_++;
            all_entities_[typeid(T)][entity_count_] = entity;
        }
        return entity;
    }

    template <typename T>
    T &getEntityByName(const std::string &name)
    {
        T *entity = findEntityByName<T>(name);
        if (entity != nullptr)
        {
            return *entity;
        }
        throw std::runtime_error(type_name<T>() + ": " + name + " not found");
    }

    template <typename T>
    std::vector<T *> entitiesWith()
    {
        std::vector<T *> result;
        for (const auto &[entity_id, entity] : all_entities_[typeid(T)])
        {
            result.push_back(static_cast<T *>(entity));
        }
        return result;
    }

  protected:
    template <typename T>
    T *findEntityByName(const std::string &name)
    {
        for (const auto &[entity_id, entity] : all_entities_[typeid(T)])
        {
            if (entity->Name() == name)
            {
                return static_cast<T *>(entity);
            }
        }
        return nullptr;
    }
};
} // namespace SPH
#endif // SPHINXSYS_ENTITY_H
