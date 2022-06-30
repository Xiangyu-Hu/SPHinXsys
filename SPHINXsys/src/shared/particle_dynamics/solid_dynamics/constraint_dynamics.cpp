/**
 * @file 	constraint_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "constraint_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		ConstrainSolidBodyRegion::
			ConstrainSolidBodyRegion(SPHBody &sph_body, BodyPartByParticle &body_part)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidDataSimple(sph_body),
			  pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			  n_(particles_->n_), n_0_(particles_->n_0_),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_)
		{
		}
		//=================================================================================================//
		void ConstrainSolidBodyRegion::Update(size_t index_i, Real dt)
		{
			Vecd pos_0 = pos_0_[index_i];
			Vecd pos_n = pos_n_[index_i];
			Vecd vel_n = vel_n_[index_i];
			Vecd dvel_dt = dvel_dt_[index_i];

			pos_n_[index_i] = getDisplacement(pos_0, pos_n);
			vel_n_[index_i] = getVelocity(pos_0, pos_n, vel_n);
			dvel_dt_[index_i] = getAcceleration(pos_0, pos_n, dvel_dt);
			/** the average values are prescirbed also. */
			vel_ave_[index_i] = vel_n_[index_i];
			dvel_dt_ave_[index_i] = dvel_dt_[index_i];
		}
		//=================================================================================================//
		ConstrainSolidBodySurfaceRegion::
			ConstrainSolidBodySurfaceRegion(SPHBody &body, BodyPartByParticle &body_part)
			: PartSimpleDynamicsByParticle(body, body_part), SolidDataSimple(body),
			  pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  apply_constrain_to_particle_(StdLargeVec<bool>(pos_0_.size(), false))
		{
			// get the surface layer of particles
			BodySurface surface_layer(body);
			// select which particles the spring is applied to
			// if the particle is in the surface layer, the force is applied
			for (size_t particle_i : surface_layer.body_part_particles_)
				apply_constrain_to_particle_[particle_i] = true;
		}
		//=================================================================================================//
		void ConstrainSolidBodySurfaceRegion::Update(size_t index_i, Real dt)
		{
			if (apply_constrain_to_particle_[index_i])
			{
				pos_n_[index_i] = pos_0_[index_i];
				vel_n_[index_i] = Vecd(0);
				dvel_dt_[index_i] = Vecd(0);
			}
		}
		//=================================================================================================//
		PositionSolidBody::
			PositionSolidBody(SPHBody &sph_body, BodyPartByParticle &body_part,
							  Real start_time, Real end_time, Vecd pos_end_center)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidDataSimple(sph_body),
			  pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_),
			  start_time_(start_time), end_time_(end_time), pos_end_center_(pos_end_center)
		{
			BoundingBox bounds = sph_body.getBodyDomainBounds();
			pos_0_center_ = (bounds.first + bounds.second) * 0.5;
			translation_ = pos_end_center_ - pos_0_center_;
		}
		//=================================================================================================//
		Vecd PositionSolidBody::getDisplacement(size_t index_i, Real dt)
		{
			Vecd displacement;
			try
			{
				// displacement from the initial position
				Vecd pos_final = pos_0_[index_i] + translation_;
				displacement = (pos_final - pos_n_[index_i]) * dt /
							   (end_time_ - GlobalStaticVariables::physical_time_);
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("PositionSolidBody::getDisplacement: particle index out of bounds") + std::to_string(index_i));
			}
			return displacement;
		}
		//=================================================================================================//
		void PositionSolidBody::Update(size_t index_i, Real dt)
		{
			try
			{
				// only apply in the defined time period
				if (GlobalStaticVariables::physical_time_ >= start_time_ &&
					GlobalStaticVariables::physical_time_ <= end_time_)
				{
					pos_n_[index_i] = pos_n_[index_i] + getDisplacement(index_i, dt); // displacement from the initial position
					vel_n_[index_i] = getVelocity();
				}
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("PositionSolidBody::Update: particle index out of bounds") + std::to_string(index_i));
			}
		}
		//=================================================================================================//
		PositionScaleSolidBody::
			PositionScaleSolidBody(SPHBody &sph_body, BodyPartByParticle &body_part,
								   Real start_time, Real end_time, Real end_scale)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidDataSimple(sph_body),
			  pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_),
			  start_time_(start_time), end_time_(end_time), end_scale_(end_scale)
		{
			BoundingBox bounds = sph_body.getBodyDomainBounds();
			pos_0_center_ = (bounds.first + bounds.second) * 0.5;
		}
		//=================================================================================================//
		Vecd PositionScaleSolidBody::getDisplacement(size_t index_i, Real dt)
		{
			Vecd displacement(0);
			try
			{
				// displacement from the initial position
				Vecd pos_final = pos_0_center_ + end_scale_ * (pos_0_[index_i] - pos_0_center_);
				displacement = (pos_final - pos_n_[index_i]) * dt /
							   (end_time_ - GlobalStaticVariables::physical_time_);
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("PositionScaleSolidBody::getDisplacement: particle index out of bounds") + std::to_string(index_i));
			}
			return displacement;
		}
		//=================================================================================================//
		void PositionScaleSolidBody::Update(size_t index_i, Real dt)
		{
			try
			{
				// only apply in the defined time period
				if (GlobalStaticVariables::physical_time_ >= start_time_ &&
					GlobalStaticVariables::physical_time_ <= end_time_)
				{
					pos_n_[index_i] = pos_n_[index_i] + getDisplacement(index_i, dt); // displacement from the initial position
					vel_n_[index_i] = getVelocity();
				}
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("PositionScaleSolidBody::Update: particle index out of bounds") + std::to_string(index_i));
			}
		}
		//=================================================================================================//
		bool checkIfPointInBoundingBox(Vecd point, BoundingBox &bbox)
		{
			if (point.size() >= 3 && bbox.first.size() >= 3 && bbox.second.size() >= 3)
			{
				return point[0] >= bbox.first[0] && point[0] <= bbox.second[0] &&
					   point[1] >= bbox.first[1] && point[1] <= bbox.second[1] &&
					   point[2] >= bbox.first[2] && point[2] <= bbox.second[2];
			}
			if (point.size() >= 2 && bbox.first.size() >= 2 && bbox.second.size() >= 2)
			{
				return point[0] >= bbox.first[0] && point[0] <= bbox.second[0] &&
					   point[1] >= bbox.first[1] && point[1] <= bbox.second[1];
			}
			if (point.size() >= 1 && bbox.first.size() >= 1 && bbox.second.size() >= 1)
			{
				return point[0] >= bbox.first[0] && point[0] <= bbox.second[0];
			}
			throw std::runtime_error(std::string("checkIfPointInBoundingBox: Vecd point or BoundingBox& bbox has a dimension of <1"));
			return false;
		}
		//=================================================================================================//
		TranslateSolidBody::
			TranslateSolidBody(SPHBody &sph_body, BodyPartByParticle &body_part, Real start_time, Real end_time, Vecd translation)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidDataSimple(sph_body),
			  pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_), pos_end_({}),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  start_time_(start_time), end_time_(end_time), translation_(translation)
		{
			// record the particle positions that should be reached at end time
			for (size_t index_i = 0; index_i < pos_n_.size(); index_i++)
			{
				pos_end_.push_back(pos_n_[index_i] + translation_);
			}
		}
		//=================================================================================================//
		Vecd TranslateSolidBody::getDisplacement(size_t index_i, Real dt)
		{
			Vecd displacement(0);
			try
			{
				displacement = (pos_end_[index_i] - pos_n_[index_i]) * dt / (end_time_ - GlobalStaticVariables::physical_time_);
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("TranslateSolidBody::getDisplacement: particle index out of bounds") + std::to_string(index_i));
			}
			return displacement;
		}
		//=================================================================================================//
		void TranslateSolidBody::Update(size_t index_i, Real dt)
		{
			// only apply in the defined time period
			if (GlobalStaticVariables::physical_time_ >= start_time_ && GlobalStaticVariables::physical_time_ <= end_time_)
			{
				try
				{
					pos_n_[index_i] = pos_n_[index_i] + 0.5 * getDisplacement(index_i, dt); // displacement from the initial position, 0.5x because it's executed twice
					vel_n_[index_i] = getVelocity();
				}
				catch (std::out_of_range &e)
				{
					throw std::runtime_error(std::string("TranslateSolidBody::Update: particle index out of bounds") + std::to_string(index_i));
				}
			}
		}
		//=================================================================================================//
		TranslateSolidBodyPart::
			TranslateSolidBodyPart(SPHBody &sph_body, BodyPartByParticle &body_part, Real start_time, Real end_time, Vecd translation, BoundingBox bbox)
			: TranslateSolidBody(sph_body, body_part, start_time, end_time, translation), bbox_(bbox)
		{
		}
		//=================================================================================================//
		void TranslateSolidBodyPart::Update(size_t index_i, Real dt)
		{
			try
			{
				Vecd point = pos_0_[index_i];
				if (checkIfPointInBoundingBox(point, bbox_))
				{
					if (GlobalStaticVariables::physical_time_ >= start_time_ && GlobalStaticVariables::physical_time_ <= end_time_)
					{
						vel_n_[index_i] = getDisplacement(index_i, dt) / dt;
					}
					else
					{
						vel_n_[index_i] = 0;
						dvel_dt_[index_i] = 0;
					}
				}
			}
			catch (std::out_of_range &e)
			{
				throw std::runtime_error(std::string("TranslateSolidBodyPart::Update: particle index out of bounds") + std::to_string(index_i));
			}
		}
		//=================================================================================================//
		SoftConstrainSolidBodyRegion::
			SoftConstrainSolidBodyRegion(BaseBodyRelationInner &inner_relation, BodyPartByParticle &body_part)
			: PartInteractionDynamicsByParticleWithUpdate(*inner_relation.sph_body_, body_part),
			  SolidDataInner(inner_relation),
			  Vol_(particles_->Vol_),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_)
		{
			particles_->registerAVariable(vel_temp_, "TemporaryVelocity");
			particles_->registerAVariable(dvel_dt_temp_, "TemporaryAcceleration");
		}
		//=================================================================================================//
		void SoftConstrainSolidBodyRegion::Interaction(size_t index_i, Real dt)
		{
			Real ttl_weight(Eps);
			Vecd vel_i = vel_n_[index_i];
			Vecd dvel_dt_i = dvel_dt_[index_i];

			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real weight_j = inner_neighborhood.W_ij_[n] * Vol_[index_j];

				ttl_weight += weight_j;
				vel_i += vel_n_[index_j] * weight_j;
				dvel_dt_i += dvel_dt_[index_j] * weight_j;
			}

			vel_temp_[index_i] = vel_i / ttl_weight;
			dvel_dt_temp_[index_i] = dvel_dt_i / ttl_weight;
		}
		//=================================================================================================//
		void SoftConstrainSolidBodyRegion::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] = vel_temp_[index_i];
			dvel_dt_[index_i] = dvel_dt_temp_[index_i];
			/** the average values are prescirbed also. */
			vel_ave_[index_i] = vel_n_[index_i];
			dvel_dt_ave_[index_i] = dvel_dt_[index_i];
		}
		//=================================================================================================//
		ClampConstrainSolidBodyRegion::
			ClampConstrainSolidBodyRegion(BaseBodyRelationInner &inner_relation, BodyPartByParticle &body_part)
			: ParticleDynamics<void>(*inner_relation.sph_body_),
			  constrianing_(ConstrainSolidBodyRegion(*inner_relation.sph_body_, body_part)),
			  softing_(SoftConstrainSolidBodyRegion(inner_relation, body_part)) {}
		//=================================================================================================//
		void ClampConstrainSolidBodyRegion::exec(Real dt)
		{
			constrianing_.exec();
			softing_.exec();
		}
		//=================================================================================================//
		void ClampConstrainSolidBodyRegion::parallel_exec(Real dt)
		{
			constrianing_.parallel_exec();
			softing_.parallel_exec();
		}
		//=================================================================================================//
		ConstrainSolidBodyMassCenter::
			ConstrainSolidBodyMassCenter(SPHBody &sph_body, Vecd constrain_direction)
			: ParticleDynamicsSimple(sph_body), SolidDataSimple(sph_body),
			  correction_matrix_(Matd(1.0)), vel_n_(particles_->vel_n_),
			  compute_total_momentum_(sph_body, "Velocity")
		{
			for (int i = 0; i != Dimensions; ++i)
				correction_matrix_[i][i] = constrain_direction[i];
			BodySummation<Real> compute_total_mass_(sph_body, "Mass");
			total_mass_ = compute_total_mass_.parallel_exec();
		}
		//=================================================================================================//
		void ConstrainSolidBodyMassCenter::setupDynamics(Real dt)
		{
			velocity_correction_ =
				correction_matrix_ * compute_total_momentum_.parallel_exec(dt) / total_mass_;
		}
		//=================================================================================================//
		void ConstrainSolidBodyMassCenter::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] -= velocity_correction_;
		}
		//=================================================================================================//
		ConstrainSolidBodyPartBySimBody::
			ConstrainSolidBodyPartBySimBody(SolidBody &solid_body,
											SolidBodyPartForSimbody &body_part,
											SimTK::MultibodySystem &MBsystem,
											SimTK::MobilizedBody &mobod,
											SimTK::Force::DiscreteForces &force_on_bodies,
											SimTK::RungeKuttaMersonIntegrator &integ)
			: ConstrainSolidBodyRegion(solid_body, body_part),
			  MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			initial_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
		void ConstrainSolidBodyPartBySimBody::setupDynamics(Real dt)
		{
			body_->setNewlyUpdated();
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
		}
		//=================================================================================================//
		TotalForceOnSolidBodyPartForSimBody::
			TotalForceOnSolidBodyPartForSimBody(SolidBody &solid_body,
												SolidBodyPartForSimbody &body_part,
												SimTK::MultibodySystem &MBsystem,
												SimTK::MobilizedBody &mobod,
												SimTK::Force::DiscreteForces &force_on_bodies,
												SimTK::RungeKuttaMersonIntegrator &integ)
			: PartDynamicsByParticleReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>(solid_body, body_part),
			  SolidDataSimple(solid_body), mass_(particles_->mass_),
			  force_from_fluid_(particles_->force_from_fluid_), dvel_dt_prior_(particles_->dvel_dt_prior_),
			  pos_n_(particles_->pos_n_),
			  MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			initial_reference_ = SpatialVec(Vec3(0), Vec3(0));
		}
		//=================================================================================================//
		void TotalForceOnSolidBodyPartForSimBody::SetupReduce()
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
	}
}
