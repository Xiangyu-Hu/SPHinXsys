
#ifndef COMPOSITE_MATERIAL_H
#define COMPOSITE_MATERIAL_H

#include "elastic_solid.h"

#include <general_dynamics.h>

namespace SPH
{
/**@class CompositeSolid*/
class CompositeSolid : public ElasticSolid
{
  protected:
    UniquePtrsKeeper<ElasticSolid> composite_ptrs_keeper_;
    StdVec<ElasticSolid *> composite_materials_;

  public:
    StdLargeVec<int> material_id_;

    explicit CompositeSolid(Real rho0) : ElasticSolid(rho0)
    {
        material_type_name_ = "CompositeSolid";
    };
    virtual ~CompositeSolid(){};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Matd StressPK2(Matd &deformation, size_t index_i) override
    {
        return composite_materials_[material_id_[index_i]]->StressPK2(deformation, index_i);
    };

    virtual Matd StressPK1(Matd &deformation, size_t index_i) override
    {
        return composite_materials_[material_id_[index_i]]->StressPK1(deformation, index_i);
    };

    Real CompositeDensity(size_t index_i)
    {
        return composite_materials_[material_id_[index_i]]->ReferenceDensity();
    };

    virtual Matd StressCauchy(Matd &almansi_strain, Matd &F, size_t index_i) override
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
} // namespace SPH
#endif // COMPOSITE_MATERIAL_H
