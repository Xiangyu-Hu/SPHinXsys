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
#include <algorithm>
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
    std::unordered_map<std::type_index, std::unordered_map<std::string, void *>> all_entities_;

  public:
    EntityManager() = default;
    ~EntityManager() {};

    /** Remove all registered entities (non-owning registry reset). */
    void clear() { all_entities_.clear(); }

    template <typename T>
    T *addEntity(const std::string &name, T *entity)
    {
        T *existing_entity = findEntityByName<T>(name);
        if (existing_entity == nullptr)
        {
            all_entities_[typeid(T)][name] = entity;
            return entity;
        }
        return existing_entity;
    }

    template <typename T>
    T *addEntityOrThrow(const std::string &name, T *entity)
    {
        T *existing_entity = findEntityByName<T>(name);
        if (existing_entity != nullptr)
        {
            throw std::runtime_error(std::string(type_name<T>()) + ": duplicated entity name '" + name + "'");
        }
        all_entities_[typeid(T)][name] = entity;
        return entity;
    }

    template <typename T>
    bool removeEntity(const std::string &name)
    {
        auto type_it = all_entities_.find(typeid(T));
        if (type_it == all_entities_.end())
            return false;

        return type_it->second.erase(name) > 0;
    }

    template <typename T>
    bool hasEntity(const std::string &name) const
    {
        return findEntityByName<T>(name) != nullptr;
    }

    template <typename T>
    T *tryGetEntityByName(const std::string &name)
    {
        return findEntityByName<T>(name);
    }

    template <typename T>
    const T *tryGetEntityByName(const std::string &name) const
    {
        return findEntityByName<T>(name);
    }

    template <typename T>
    T &getEntityByName(const std::string &name)
    {
        T *entity = tryGetEntityByName<T>(name);
        if (entity != nullptr)
        {
            return *entity;
        }
        throw std::runtime_error(std::string(type_name<T>()) + ": " + name + " not found");
    }

    template <typename T>
    std::vector<T *> entitiesWith()
    {
        std::vector<T *> result;
        auto type_it = all_entities_.find(typeid(T));
        if (type_it == all_entities_.end())
            return result;

        for (const auto &[name, entity] : type_it->second)
        {
            result.push_back(static_cast<T *>(entity));
        }
        return result;
    }

    template <typename T>
    std::vector<T *> entitiesWithSorted()
    {
        std::vector<std::pair<std::string, T *>> named_entities;
        auto type_it = all_entities_.find(typeid(T));
        if (type_it == all_entities_.end())
            return {};

        named_entities.reserve(type_it->second.size());
        for (const auto &[name, entity] : type_it->second)
        {
            named_entities.emplace_back(name, static_cast<T *>(entity));
        }

        std::sort(named_entities.begin(), named_entities.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });

        std::vector<T *> result;
        result.reserve(named_entities.size());
        for (const auto &entry : named_entities)
        {
            result.push_back(entry.second);
        }
        return result;
    }

  protected:
    template <typename T>
    T *findEntityByName(const std::string &name)
    {
        auto type_it = all_entities_.find(typeid(T));
        if (type_it == all_entities_.end())
            return nullptr;

        auto entity_it = type_it->second.find(name);
        if (entity_it == type_it->second.end())
            return nullptr;

        return static_cast<T *>(entity_it->second);
    }

    template <typename T>
    const T *findEntityByName(const std::string &name) const
    {
        auto type_it = all_entities_.find(typeid(T));
        if (type_it == all_entities_.end())
            return nullptr;

        auto entity_it = type_it->second.find(name);
        if (entity_it == type_it->second.end())
            return nullptr;

        return static_cast<const T *>(entity_it->second);
    }
};
} // namespace SPH
#endif // SPHINXSYS_ENTITY_H
