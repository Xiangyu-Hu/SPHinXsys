
#ifndef COMPOSITE_MATERIAL_H
#define COMPOSITE_MATERIAL_H

#include "elastic_solid.h"

#include <general_dynamics.h>

namespace SPH
{
/**@class CompositeMaterial*/
class CompositeMaterial : public ElasticSolid
{
  protected:
    UniquePtrsKeeper<ElasticSolid> composite_ptrs_keeper_;
    StdVec<ElasticSolid *> composite_materials_;

  public:
    StdLargeVec<int> material_id_;

    explicit CompositeMaterial(Real rho0) : ElasticSolid(rho0)
    {
        material_type_name_ = "CompositeMaterial";
    };
    virtual ~CompositeMaterial(){};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override
    {
        return composite_materials_[material_id_[particle_index_i]]->StressPK2(deformation, particle_index_i);
    };

    Real CompositeDensity(size_t particle_index_i)
    {
        return composite_materials_[material_id_[particle_index_i]]->ReferenceDensity();
    };

    virtual Matd StressCauchy(Matd &almansi_strain, Matd &F, size_t particle_index_i) override
    {
        return Matd::Identity();
    };

    virtual Real VolumetricKirchhoff(Real J) override
    {
        return 0.0;
    };

    virtual std::string getRelevantStressMeasureName() override { return "PK2"; };

    template <class ElasticSolidType, typename... Args>
    void add(Args &&...args)
    {
        ElasticSolid *added_material =
            composite_ptrs_keeper_.createPtr<ElasticSolidType>(std::forward<Args>(args)...);
        composite_materials_.push_back(added_material);
        c0_ = SMAX(c0_, added_material->ReferenceSoundSpeed());
        setContactStiffness(c0_);
    };
};

/**
 * @class MaterialIdInitialization
 */
class MaterialIdInitialization
    : public LocalDynamics,
      public GeneralDataDelegateSimple
{
  public:
    explicit MaterialIdInitialization(SPHBody &sph_body);

  protected:
    StdLargeVec<int> &material_id_;
    StdLargeVec<Vecd> &pos0_;
};

/**
 * @class ActiveModelSolid
 */
class ActiveModelSolid : public SaintVenantKirchhoffSolid
{
  public:
    explicit ActiveModelSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
        : SaintVenantKirchhoffSolid(rho0, youngs_modulus, poisson_ratio)
    {
        material_type_name_ = "ActiveModelSolid";
    };
    virtual ~ActiveModelSolid(){};

    /** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i);
};
} // namespace SPH
#endif // COMPOSITE_MATERIAL_H
