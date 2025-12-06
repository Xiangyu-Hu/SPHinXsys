#ifndef WINDKESSEL_BC_H
#define WINDKESSEL_BC_H

#include "sphinxsys.h"
#include "pressure_boundary.h"

namespace SPH
{
namespace fluid_dynamics
{
class TargetOutletPressureWindkessel : public BaseLocalDynamics<BodyPartByCell>
{
  public:
    explicit TargetOutletPressureWindkessel(AlignedBoxByCell& aligned_box_part)
        : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
          part_id_(aligned_box_part.getPartID()),
          Rp_(0.0), C_(0.0), Rd_(0.0), delta_t_(0.0), 
          Q_n_(0.0), Q_0_(0.0), p_n_(80*133.32), p_0_(80*133.32),
          flow_rate_(*(this->particles_->registerSingularVariable<Real>("FlowRate" + std::to_string(part_id_ - 1))->Data())),
          current_flow_rate_(0.0), previous_flow_rate_(0.0),
          M_n_(0.0), current_mass_flow_rate_(0.0), previous_mass_flow_rate_(0.0),
          acc_mass_flow_rate_(*(this->particles_->registerSingularVariable<Real>("AccMassFlowRate" + std::to_string(part_id_ - 1))->Data())),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime"))
    {};
    virtual ~TargetOutletPressureWindkessel(){};

    void setWindkesselParams(Real Rp, Real C, Real Rd, Real dt)
    {
        Rp_ = Rp;
        C_ = C;
        Rd_ = Rd;
        delta_t_ = dt;
    }

    void updateNextPressure()
    {
        getFlowRate();

        Q_n_ = current_flow_rate_ / delta_t_;
        M_n_ = current_mass_flow_rate_ / delta_t_;

        Real dp_dt = - p_0_ / (C_ * Rd_) + (Rp_ + Rd_) * Q_n_ / (C_ * Rd_) + Rp_ * (Q_n_ - Q_0_) / (delta_t_ + TinyReal);
        Real p_star = p_0_ + dp_dt * delta_t_;
        Real dp_dt_star = - p_star / (C_ * Rd_) + (Rp_ + Rd_) * Q_n_ / (C_ * Rd_) + Rp_ * (Q_n_ - Q_0_) / (delta_t_ + TinyReal);
        p_n_ = p_0_ + 0.5 * delta_t_ * (dp_dt + dp_dt_star);

        std::cout << "p_n_ = " << p_n_ / 133.32 << " mmHg" << std::endl;

        writeOutletPressureData();
        writeOutletFlowRateData();
    }

    Real operator()(Real p, Real current_time)
    {
        return p_n_ - 80*133.32;
    }

  protected:
    int part_id_;
    Real Rp_, C_, Rd_, delta_t_;
    Real Q_n_, Q_0_;
    Real p_n_, p_0_;
    Real &flow_rate_, current_flow_rate_, previous_flow_rate_;
    Real M_n_;
    Real current_mass_flow_rate_, previous_mass_flow_rate_, &acc_mass_flow_rate_;
    Real *physical_time_;

    void getFlowRate()
    {
        Q_0_ = Q_n_;
        p_0_ = p_n_;
        current_flow_rate_ = flow_rate_ - previous_flow_rate_;
        previous_flow_rate_ = flow_rate_;

        current_mass_flow_rate_ = acc_mass_flow_rate_ - previous_mass_flow_rate_;
        previous_mass_flow_rate_ = acc_mass_flow_rate_;
    }

    void writeOutletPressureData()
    {
        std::string output_folder = "./output";
        std::string filefullpath = output_folder + "/" + std::to_string(part_id_ - 1) + "_windkessel_outlet_pressure.dat";
        std::ofstream out_file(filefullpath.c_str(), std::ios::app);
        out_file << *physical_time_ << "   " << p_n_ <<  "\n";
        out_file.close();
    }

    void writeOutletFlowRateData()
    {
        std::string output_folder = "./output";
        std::string filefullpath = output_folder + "/" + std::to_string(part_id_ - 1) + "_volume_flow_rate.dat";
        std::ofstream out_file(filefullpath.c_str(), std::ios::app);
        out_file << *physical_time_ << "   " << Q_n_ <<  "\n";
        out_file.close();

        std::string filefullpath_mass = output_folder + "/" + std::to_string(part_id_ - 1) + "_mass_flow_rate.dat";
        std::ofstream out_file_mass(filefullpath_mass.c_str(), std::ios::app);
        out_file_mass << *physical_time_ << "   " << M_n_ <<  "\n";
        out_file_mass.close();
    }
};

using WindkesselBoundaryCondition = PressureCondition<TargetOutletPressureWindkessel>;

//----------------------------------------------------------------------
//	Windkessel buffer
//----------------------------------------------------------------------
template <typename TargetPressure, class ExecutionPolicy = ParallelPolicy>
class BidirectionalBufferWindkessel
{
  protected:
    TargetPressure target_pressure_;

    class TagBufferParticles : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        TagBufferParticles(AlignedBoxByCell &aligned_box_part)
            : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
              part_id_(aligned_box_part.getPartID()),
              pos_(particles_->getVariableDataByName<Vecd>("Position")),
              aligned_box_(aligned_box_part.getAlignedBox()),
              buffer_indicator_(particles_->registerStateVariableData<int>("BufferIndicator"))
        {
            particles_->addEvolvingVariable<int>("BufferIndicator");
        };
        virtual ~TagBufferParticles() {};

        virtual void update(size_t index_i, Real dt = 0.0)
        {
            if (aligned_box_.checkInBounds(pos_[index_i]))
            {
                buffer_indicator_[index_i] = part_id_;
            }
        };

      protected:
        int part_id_;
        Vecd *pos_;
        AlignedBox &aligned_box_;
        int *buffer_indicator_;
    };

    class Injection : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        Injection(AlignedBoxByCell &aligned_box_part, ParticleBuffer<Base> &particle_buffer,
                  TargetPressure &target_pressure)
            : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
              part_id_(aligned_box_part.getPartID()),
              particle_buffer_(particle_buffer),
              aligned_box_(aligned_box_part.getAlignedBox()),
              fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
              pos_(particles_->getVariableDataByName<Vecd>("Position")),
              rho_(particles_->getVariableDataByName<Real>("Density")),
              p_(particles_->getVariableDataByName<Real>("Pressure")),
              Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
              previous_surface_indicator_(particles_->getVariableDataByName<int>("PreviousSurfaceIndicator")),
              buffer_indicator_(particles_->getVariableDataByName<int>("BufferIndicator")),
              upper_bound_fringe_(0.5 * sph_body_->getSPHBodyResolutionRef()),
              physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")),
              flow_rate_(*(this->particles_->template getSingularVariableByName<Real>("FlowRate" + std::to_string(part_id_ - 1))->Data())),
              acc_mass_flow_rate_(*(this->particles_->template getSingularVariableByName<Real>("AccMassFlowRate" + std::to_string(part_id_ - 1))->Data())),
              target_pressure_(target_pressure)
              {
            particle_buffer_.checkParticlesReserved();
        };
        virtual ~Injection() {};

        void update(size_t index_i, Real dt = 0.0)
        {
            if (!aligned_box_.checkInBounds(pos_[index_i]))
            {
                if (aligned_box_.checkUpperBound(pos_[index_i], upper_bound_fringe_) &&
                    buffer_indicator_[index_i] == part_id_ &&
                    index_i < particles_->TotalRealParticles())
                {
                    mutex_switch.lock();
                    particle_buffer_.checkEnoughBuffer(*particles_);
                    size_t new_particle_index = particles_->createRealParticleFrom(index_i);
                    buffer_indicator_[new_particle_index] = 0;

                    /** Periodic bounding. */
                    pos_[index_i] = aligned_box_.getUpperPeriodic(pos_[index_i]);
                    Real sound_speed = fluid_.getSoundSpeed(rho_[index_i]);
                    p_[index_i] = target_pressure_(p_[index_i], *physical_time_);
                    rho_[index_i] = p_[index_i] / pow(sound_speed, 2.0) + fluid_.ReferenceDensity();
                    previous_surface_indicator_[index_i] = 1;
                    mutex_switch.unlock();

                    flow_rate_ -= Vol_[index_i];
                    acc_mass_flow_rate_ -= Vol_[index_i] * rho_[index_i];
                }
            }
        }

      protected:
        int part_id_;
        std::mutex mutex_switch;
        ParticleBuffer<Base> &particle_buffer_;
        AlignedBox &aligned_box_;
        Fluid &fluid_;
        Vecd *pos_;
        Real *rho_, *p_, *Vol_;
        int *previous_surface_indicator_, *buffer_indicator_;
        Real upper_bound_fringe_;
        Real *physical_time_;
        Real &flow_rate_, &acc_mass_flow_rate_;

      private:
        TargetPressure &target_pressure_;
    };

    class Deletion : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        Deletion(AlignedBoxByCell &aligned_box_part)
            : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
              part_id_(aligned_box_part.getPartID()),
              aligned_box_(aligned_box_part.getAlignedBox()),
              pos_(particles_->getVariableDataByName<Vecd>("Position")),
              rho_(particles_->getVariableDataByName<Real>("Density")),
              Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
              buffer_indicator_(particles_->getVariableDataByName<int>("BufferIndicator")),
              flow_rate_(*(this->particles_->template getSingularVariableByName<Real>("FlowRate" + std::to_string(part_id_ - 1))->Data())),
              acc_mass_flow_rate_(*(this->particles_->template getSingularVariableByName<Real>("AccMassFlowRate" + std::to_string(part_id_ - 1))->Data())) {};
        virtual ~Deletion() {};

        void update(size_t index_i, Real dt = 0.0)
        {
            if (!aligned_box_.checkInBounds(pos_[index_i]))
            {
                mutex_switch.lock();
                while (aligned_box_.checkLowerBound(pos_[index_i]) &&
                       buffer_indicator_[index_i] == part_id_ &&
                       index_i < particles_->TotalRealParticles())
                {
                    particles_->switchToBufferParticle(index_i);
                    flow_rate_ += Vol_[index_i];
                    acc_mass_flow_rate_ += Vol_[index_i] * rho_[index_i];
                }
                mutex_switch.unlock();
            }
        }

      protected:
        int part_id_;
        std::mutex mutex_switch;
        AlignedBox &aligned_box_;
        Vecd *pos_;
        Real *rho_, *Vol_;
        int *buffer_indicator_;
        Real &flow_rate_, &acc_mass_flow_rate_;
    };

  public:
    BidirectionalBufferWindkessel(AlignedBoxByCell &aligned_box_part, ParticleBuffer<Base> &particle_buffer)
        : target_pressure_(TargetPressure(aligned_box_part)),
          tag_buffer_particles(aligned_box_part),
          injection(aligned_box_part, particle_buffer, target_pressure_),
          deletion(aligned_box_part) {};
    virtual ~BidirectionalBufferWindkessel() {};

    SimpleDynamics<TagBufferParticles, ExecutionPolicy> tag_buffer_particles;
    SimpleDynamics<Injection, ExecutionPolicy> injection;
    SimpleDynamics<Deletion, ExecutionPolicy> deletion;
};

using WindkesselOutletBidirectionalBuffer = BidirectionalBufferWindkessel<TargetOutletPressureWindkessel>;

class TotalVelocityNormVal
    : public BaseLocalDynamicsReduce<ReduceSum<Real>, BodyPartByCell>
{
  protected:
    Vecd *vel_;
    AlignedBox &aligned_box_;
    const int alignment_axis_;
    Transform &transform_;

  public:
      explicit TotalVelocityNormVal(AlignedBoxByCell& aligned_box_part)
          : BaseLocalDynamicsReduce<ReduceSum<Real>, BodyPartByCell>(aligned_box_part),
          vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
          aligned_box_(aligned_box_part.getAlignedBox()),
          alignment_axis_(aligned_box_.AlignmentAxis()),
          transform_(aligned_box_.getTransform()) {};

    virtual ~TotalVelocityNormVal(){};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd frame_velocity = Vecd::Zero();
        frame_velocity[alignment_axis_] = transform_.xformBaseVecToFrame(vel_[index_i])[alignment_axis_];
        return frame_velocity[alignment_axis_];
    }
};

template <class ReduceSumType>
class AreaAverageFlowRate : public ReduceSumType
{
  public:
    template <class DynamicsIdentifier, typename... Args>
    explicit AreaAverageFlowRate(DynamicsIdentifier &identifier, Real outlet_area,  Args &&...args)
          : ReduceSumType(identifier, std::forward<Args>(args)...),
            part_id_(identifier.getPartID()),
            tansient_flow_rate_(*(this->particles_->template registerSingularVariable<Real>("TransientVolumeFlowRate" + std::to_string(part_id_ - 1))->Data())),
            outlet_area_(outlet_area), physical_time_(this->sph_system_->template getSystemVariableDataByName<Real>("PhysicalTime")) {};
    virtual ~AreaAverageFlowRate(){};

    virtual Real outputResult(Real reduced_value) override
    {
        Real average_velocity_norm = ReduceSumType::outputResult(reduced_value) / Real(this->identifier_->SizeOfLoopRange());
        tansient_flow_rate_ = average_velocity_norm * outlet_area_;

        std::string output_folder = "./output";
        std::string filefullpath = output_folder + "/" + std::to_string(part_id_ - 1) + "_transient_VolumeFlowRate.dat";
        std::ofstream out_file(filefullpath.c_str(), std::ios::app);
        out_file << *physical_time_ << "   " << tansient_flow_rate_ <<  "\n";
        out_file.close();

        return tansient_flow_rate_;
    }

  private:
    int part_id_;
    Real &tansient_flow_rate_;
    Real outlet_area_;
    Real *physical_time_;
};

using SectionTransientFlowRate = AreaAverageFlowRate<TotalVelocityNormVal>;

} // namespace fluid_dynamics
} // namespace SPH
#endif // WINDKESSEL_BC_H