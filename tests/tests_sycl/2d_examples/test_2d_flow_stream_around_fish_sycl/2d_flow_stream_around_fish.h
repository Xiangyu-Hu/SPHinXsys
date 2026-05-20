
#include "2d_fish_and_bones.h"
#include "active_model.h"
#include "sphinxsys.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.8;                                /**< Channel length. */
Real DH = 0.4;                                /**< Channel height. */
Real particle_spacing_ref = 0.0025;           /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0;         /**< Boundary width, determined by specific layer of boundary particles. */

Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;

Real cx = 0.3 * DL;           /**< Center of fish in x direction. */
Real cy = DH / 2;             /**< Center of fish in y direction. */
Real fish_length = 0.2;       /**< Length of fish. */
Real fish_thickness = 0.03;   /**< The maximum fish thickness. */
Real muscle_thickness = 0.02; /**< The maximum fish thickness. */
Real head_length = 0.03;      /**< Length of fish bone. */
Real bone_thickness = 0.003;  /**< Length of fish bone. */
Real fish_shape_resolution = particle_spacing_ref * 0.5;
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                /**< Density. */
Real U_f = 1.0;                      /**< freestream velocity. */
Real c_f = 10.0 * U_f;               /**< Speed of sound. */
Real Re = 30000.0;                   /**< Reynolds number. */
Real mu_f = rho0_f * U_f * 0.3 / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 1050.0;
Real Youngs_modulus1 = 0.8e6;
Real Youngs_modulus2 = 0.5e6;
Real Youngs_modulus3 = 1.1e6;
Real poisson = 0.49;
//----------------------------------------------------------------------
//	SPH bodies with cases dependent geometries.
//----------------------------------------------------------------------
class FishBody : public MultiPolygonShape
{

  public:
    explicit FishBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> fish_shape = CreatFishShape(cx, cy, fish_length, fish_shape_resolution);
        multi_polygon_.addPolygon(fish_shape, GeometricOps::add);
    }
};
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon fish(CreatFishShape(cx, cy, fish_length, fish_shape_resolution));
        subtract<MultiPolygonShape>(fish);
    }
};
//----------------------------------------------------------------------
//	Define case dependent bodies material, constraint and boundary conditions.
//----------------------------------------------------------------------
struct FreeStreamVelocity
{
    Real u_ref_, t_ref_;

    // Default constructor for EmitterInflowConditionCK (passes no extra args).
    FreeStreamVelocity() : u_ref_(U_f), t_ref_(2.0) {}

    // Constructor required by old FreeStreamVelocityCorrection (passes *this as BoundaryConditionType).
    template <class BoundaryConditionType>
    explicit FreeStreamVelocity(BoundaryConditionType &) : u_ref_(U_f), t_ref_(2.0) {}

    // Interface for EmitterInflowConditionCK: returns scalar axis velocity.
    Real getAxisVelocity(const Vecd &position, const Real &current_axis_velocity, Real time)
    {
        Real time_factor = time / t_ref_;
        return time_factor < 1.0 ? 0.5 * u_ref_ * (1.0 - std::cos(Pi * time_factor)) : u_ref_;
    }

    // Interface for FreeStreamVelocityCorrection: returns full velocity vector.
    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();
        Real time_factor = current_time / t_ref_;
        target_velocity[0] = time_factor < 1.0 ? 0.5 * u_ref_ * (1.0 - std::cos(Pi * time_factor)) : u_ref_;
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Case dependent composite material
//----------------------------------------------------------------------
class FishBodyComposite : public CompositeSolid
{
  public:
    FishBodyComposite() : CompositeSolid(rho0_s)
    {
        add<ActiveModelSolid>(rho0_s, Youngs_modulus1, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus2, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus3, poisson);
    };

    /**
     * @class ConstituteKernel
     * @brief GPU-safe stress dispatch for FishBodyComposite.
     *
     * Replaces the CPU runtime dispatch (virtual call through pointer vector)
     * with compile-time if/else on material_id_[index_i].
     * Each sub-kernel is already GPU-safe via its own ConstituteKernel.
     *
     * material_id == 0 -> ActiveModelSolid   (muscle)
     * material_id == 1 -> SVK material 1     (bone)
     * material_id == 2 -> SVK material 2     (tissue)
     */
    class ConstituteKernel
    {
      public:
        template <typename ExecutionPolicy>
        ConstituteKernel(const ExecutionPolicy &ex_policy, FishBodyComposite &encloser)
            : material_id_(encloser.dv_material_id_->DelegatedData(ex_policy)),
              active_kernel_(ex_policy,
                  static_cast<ActiveModelSolid &>(*encloser.composite_materials_[0])),
              svk1_kernel_(ex_policy,
                  static_cast<SaintVenantKirchhoffSolid &>(*encloser.composite_materials_[1])),
              svk2_kernel_(ex_policy,
                  static_cast<SaintVenantKirchhoffSolid &>(*encloser.composite_materials_[2])) {}

        inline Matd StressPK1(const Matd &F, size_t index_i)
        {
            if (material_id_[index_i] == 0)
                return active_kernel_.StressPK1(F, index_i);
            else if (material_id_[index_i] == 1)
                return svk1_kernel_.StressPK1(F, index_i);
            else
                return svk2_kernel_.StressPK1(F, index_i);
        }

        inline Real VolumetricKirchhoff(Real J) { return 0.0; }

      protected:
        int *material_id_;
        ActiveModelSolid::ConstituteKernel active_kernel_;
        SaintVenantKirchhoffSolid::ConstituteKernel svk1_kernel_;
        SaintVenantKirchhoffSolid::ConstituteKernel svk2_kernel_;
    };
};
//----------------------------------------------------------------------
//	Case dependent initialization material ids
//----------------------------------------------------------------------
Real a1 = 1.22 * fish_thickness / fish_length;
Real a2 = 3.19 * fish_thickness / fish_length / fish_length;
Real a3 = -15.73 * fish_thickness / pow(fish_length, 3);
Real a4 = 21.87 * fish_thickness / pow(fish_length, 4);
Real a5 = -10.55 * fish_thickness / pow(fish_length, 5);

Real b1 = 1.22 * muscle_thickness / fish_length;
Real b2 = 3.19 * muscle_thickness / fish_length / fish_length;
Real b3 = -15.73 * muscle_thickness / pow(fish_length, 3);
Real b4 = 21.87 * muscle_thickness / pow(fish_length, 4);
Real b5 = -10.55 * muscle_thickness / pow(fish_length, 5);

class InitializeDisplacementCK : public LocalDynamics
{
  public:
    explicit InitializeDisplacementCK(SolidBody &solid_body)
        : LocalDynamics(solid_body),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_pos_temp_(particles_->registerStateVariable<Vecd>("TemporaryPosition")) {}

    struct UpdateKernel
    {
        template <typename ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, InitializeDisplacementCK &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              pos_temp_(encloser.dv_pos_temp_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0) { pos_temp_[index_i] = pos_[index_i]; }

      protected:
        Vecd *pos_, *pos_temp_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos_temp_;
};

class UpdateAverageVelocityAndAccelerationCK : public LocalDynamics
{
  public:
    explicit UpdateAverageVelocityAndAccelerationCK(SolidBody &solid_body)
        : LocalDynamics(solid_body),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_pos_temp_(particles_->getVariableByName<Vecd>("TemporaryPosition")),
          dv_vel_ave_(particles_->registerStateVariable<Vecd>("AverageVelocity")),
          dv_acc_ave_(particles_->registerStateVariable<Vecd>("AverageAcceleration")) {}

    struct UpdateKernel
    {
        template <typename ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, UpdateAverageVelocityAndAccelerationCK &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              pos_temp_(encloser.dv_pos_temp_->DelegatedData(ex_policy)),
              vel_ave_(encloser.dv_vel_ave_->DelegatedData(ex_policy)),
              acc_ave_(encloser.dv_acc_ave_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            Vecd updated_vel_ave = (pos_[index_i] - pos_temp_[index_i]) / (dt + Eps);
            acc_ave_[index_i] = (updated_vel_ave - vel_ave_[index_i]) / (dt + Eps);
            vel_ave_[index_i] = updated_vel_ave;
        }

      protected:
        Vecd *pos_, *pos_temp_, *vel_ave_, *acc_ave_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos_temp_, *dv_vel_ave_, *dv_acc_ave_;
};

class ZeroForceCK : public LocalDynamics
{
  public:
    explicit ZeroForceCK(SolidBody &solid_body)
        : LocalDynamics(solid_body),
          dv_force_(particles_->getVariableByName<Vecd>("Force")),
          dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")) {}

    struct UpdateKernel
    {
        template <typename ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, ZeroForceCK &encloser)
            : force_(encloser.dv_force_->DelegatedData(ex_policy)),
              force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            force_[index_i] = Vecd::Zero();
            force_prior_[index_i] = Vecd::Zero();
        }

      protected:
        Vecd *force_, *force_prior_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_force_, *dv_force_prior_;
};

class FishMaterialInitialization : public MaterialIdInitialization
{
  public:
    explicit FishMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body) {};

    void update(size_t index_i, Real dt = 0.0)
    {
        Real x = pos_[index_i][0] - cx;
        Real y = pos_[index_i][1];

        Real y1 = a1 * pow(x, 0 + 1) + a2 * pow(x, 1 + 1) + a3 * pow(x, 2 + 1) + a4 * pow(x, 3 + 1) + a5 * pow(x, 4 + 1);
        if (x <= (fish_length - head_length) && y > (y1 - 0.004 + cy) && y > (cy + bone_thickness / 2))
        {
            material_id_[index_i] = 0; // region for muscle
        }
        else if (x <= (fish_length - head_length) && y < (-y1 + 0.004 + cy) && y < (cy - bone_thickness / 2))
        {
            material_id_[index_i] = 0; // region for muscle
        }
        else if ((x > (fish_length - head_length)) || ((y < (cy + bone_thickness / 2)) && (y > (cy - bone_thickness / 2))))
        {
            material_id_[index_i] = 2;
        }
        else
        {
            material_id_[index_i] = 1;
        }
    };

    class UpdateKernel
    {
      public:
        template <typename ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, FishMaterialInitialization &encloser)
            : material_id_(encloser.dv_material_id_->DelegatedData(ex_policy)),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              cx_(cx), cy_(cy), fish_length_(fish_length), head_length_(head_length),
              bone_thickness_(bone_thickness),
              a1_(a1), a2_(a2), a3_(a3), a4_(a4), a5_(a5) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            Real x = pos_[index_i][0] - cx_;
            Real y = pos_[index_i][1];

            Real y1 = a1_ * math::pow(x, Real(1)) + a2_ * math::pow(x, Real(2)) +
                      a3_ * math::pow(x, Real(3)) + a4_ * math::pow(x, Real(4)) +
                      a5_ * math::pow(x, Real(5));
            if (x <= (fish_length_ - head_length_) && y > (y1 - 0.004 + cy_) && y > (cy_ + bone_thickness_ / 2))
            {
                material_id_[index_i] = 0;
            }
            else if (x <= (fish_length_ - head_length_) && y < (-y1 + 0.004 + cy_) && y < (cy_ - bone_thickness_ / 2))
            {
                material_id_[index_i] = 0;
            }
            else if ((x > (fish_length_ - head_length_)) || ((y < (cy_ + bone_thickness_ / 2)) && (y > (cy_ - bone_thickness_ / 2))))
            {
                material_id_[index_i] = 2;
            }
            else
            {
                material_id_[index_i] = 1;
            }
        }

      protected:
        int *material_id_;
        Vecd *pos_;
        Real cx_, cy_, fish_length_, head_length_, bone_thickness_;
        Real a1_, a2_, a3_, a4_, a5_;
    };
};
//----------------------------------------------------------------------
//	imposing active strain to fish muscle
//----------------------------------------------------------------------
class ImposingActiveStrain : public LocalDynamics
{
  public:
    explicit ImposingActiveStrain(SolidBody &solid_body)
        : LocalDynamics(solid_body),
          sv_physical_time_(&sph_system_->svPhysicalTime()),
          dv_material_id_(particles_->getVariableByName<int>("MaterialID")),
          dv_pos0_(particles_->getVariableByName<Vecd>("InitialPosition")),
          dv_active_strain_(particles_->getVariableByName<Matd>("ActiveStrain")) {}

    struct UpdateKernel
    {
        template <typename ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, ImposingActiveStrain &encloser)
            : physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
              material_id_(encloser.dv_material_id_->DelegatedData(ex_policy)),
              pos0_(encloser.dv_pos0_->DelegatedData(ex_policy)),
              active_strain_(encloser.dv_active_strain_->DelegatedData(ex_policy)),
              cx_(cx), cy_(cy), fish_length_(fish_length), bone_thickness_(bone_thickness) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            if (material_id_[index_i] == 0)
            {
                Real x = pos0_[index_i][0] - cx_;
                Real y = pos0_[index_i][1];

                Real Am = 0.12;
                Real frequency = 4.0;
                Real w = 2 * Pi * frequency;
                Real lambda = 3.0 * fish_length_;
                Real wave_number = 2 * Pi / lambda;
                Real hx = -(math::pow(x, Real(2)) - math::pow(fish_length_, Real(2))) / math::pow(fish_length_, Real(2));
                Real start_time = 0.2;
                Real current_time = *physical_time_;
                Real strength = 1 - math::exp(-current_time / start_time);

                Real phase_shift = y > (cy_ + bone_thickness_ / 2) ? Real(0) : Pi / 2;
                active_strain_[index_i](0, 0) =
                    -Am * hx * strength * math::pow(math::sin(w * current_time / 2 + wave_number * x / 2 + phase_shift), Real(2));
            }
        }

      protected:
        Real *physical_time_;
        int *material_id_;
        Vecd *pos0_;
        Matd *active_strain_;
        Real cx_, cy_, fish_length_, bone_thickness_;
    };

  protected:
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<int> *dv_material_id_;
    DiscreteVariable<Vecd> *dv_pos0_;
    DiscreteVariable<Matd> *dv_active_strain_;
};