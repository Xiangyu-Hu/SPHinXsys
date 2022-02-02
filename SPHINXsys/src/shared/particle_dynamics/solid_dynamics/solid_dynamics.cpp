/**
 * @file 	solid_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "solid_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SolidBody &solid_body, Real CFL)
			: ParticleDynamicsReduce<Real, ReduceMin>(solid_body),
			  ElasticSolidDataSimple(solid_body), CFL_(CFL),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength())
		{
			initial_reference_ = DBL_MAX;
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			Real sound_speed = material_->ReferenceSoundSpeed();
			return CFL_ * SMIN(sqrt(smoothing_length_ / (dvel_dt_[index_i].norm() + TinyReal)),
							   smoothing_length_ / (sound_speed + vel_n_[index_i].norm()));
		}
		//=================================================================================================//
		CorrectConfiguration::
			CorrectConfiguration(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  SolidDataInner(inner_relation),
			  Vol_(particles_->Vol_), B_(particles_->B_) {}
		//=================================================================================================//
		void CorrectConfiguration::Interaction(size_t index_i, Real dt)
		{
			Matd local_configuration(Eps); // a small number added to diagonal to avoid divide zero
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				local_configuration -= Vol_[index_j] * SimTK::outer(r_ji, gradw_ij);
			}
			B_[index_i] = SimTK::inverse(local_configuration);
		}
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
			for (size_t particle_i: surface_layer.body_part_particles_) apply_constrain_to_particle_[particle_i] = true;
		}
		//=================================================================================================//
		void ConstrainSolidBodySurfaceRegion::Update(size_t index_i, Real dt)
		{
			if(apply_constrain_to_particle_[index_i])
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
			catch (out_of_range &e)
			{
				throw runtime_error(string("PositionSolidBody::getDisplacement: particle index out of bounds") + to_string(index_i));
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
			catch (out_of_range &e)
			{
				throw runtime_error(string("PositionSolidBody::Update: particle index out of bounds") + to_string(index_i));
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
			catch (out_of_range &e)
			{
				throw runtime_error(string("PositionScaleSolidBody::getDisplacement: particle index out of bounds") + to_string(index_i));
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
			catch (out_of_range &e)
			{
				throw runtime_error(string("PositionScaleSolidBody::Update: particle index out of bounds") + to_string(index_i));
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
			throw runtime_error(string("checkIfPointInBoundingBox: Vecd point or BoundingBox& bbox has a dimension of <1"));
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
			catch (out_of_range &e)
			{
				throw runtime_error(string("TranslateSolidBody::getDisplacement: particle index out of bounds") + to_string(index_i));
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
				catch (out_of_range &e)
				{
					throw runtime_error(string("TranslateSolidBody::Update: particle index out of bounds") + to_string(index_i));
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
			catch (out_of_range &e)
			{
				throw runtime_error(string("TranslateSolidBodyPart::Update: particle index out of bounds") + to_string(index_i));
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
			particles_->registerAVariable<indexVector, Vecd>(vel_temp_, "TemporaryVelocity");
			particles_->registerAVariable<indexVector, Vecd>(dvel_dt_temp_, "TemporaryAcceleration");
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
			BodySummation<indexScalar, Real> compute_total_mass_(sph_body, "Mass");
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
		ImposeExternalForce::
			ImposeExternalForce(SolidBody &solid_body, SolidBodyPartForSimbody &body_part)
			: PartSimpleDynamicsByParticle(solid_body, body_part), SolidDataSimple(solid_body),
			  pos_0_(particles_->pos_0_), vel_n_(particles_->vel_n_),
			  vel_ave_(particles_->vel_ave_) {}
		//=================================================================================================//
		void ImposeExternalForce::Update(size_t index_i, Real dt)
		{
			Vecd induced_acceleration = getAcceleration(pos_0_[index_i]);
			vel_n_[index_i] += induced_acceleration * dt;
			vel_ave_[index_i] = vel_n_[index_i];
		}
		//=================================================================================================//
		SpringDamperConstraintParticleWise::
			SpringDamperConstraintParticleWise(SolidBody &solid_body, Vecd stiffness, Real damping_ratio)
			: ParticleDynamicsSimple(solid_body), SolidDataSimple(solid_body),
			  pos_n_(particles_->pos_n_),
			  pos_0_(particles_->pos_0_),
			  vel_n_(particles_->vel_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_)
		{
			// scale stiffness and damping by mass here, so it's not necessary in each iteration
			stiffness_ = stiffness / std::accumulate(particles_->mass_.begin(), particles_->mass_.end(), 0.0);
			damping_coeff_ = stiffness_ * damping_ratio;

			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		Vecd SpringDamperConstraintParticleWise::getSpringForce(size_t index_i, Vecd &disp)
		{
			Vecd spring_force(0);
			for (int i = 0; i < disp.size(); i++)
			{
				spring_force[i] = -stiffness_[i] * disp[i];
			}
			return spring_force;
		}
		//=================================================================================================//
		Vecd SpringDamperConstraintParticleWise::getDampingForce(size_t index_i)
		{
			Vecd damping_force(0);
			for (int i = 0; i < vel_n_[index_i].size(); i++)
			{
				damping_force[i] = -damping_coeff_[i] * vel_n_[index_i][i];
			}
			return damping_force;
		}
		//=================================================================================================//
		void SpringDamperConstraintParticleWise::Update(size_t index_i, Real dt)
		{
			Vecd delta_x = pos_n_[index_i] - pos_0_[index_i];
			dvel_dt_prior_[index_i] += getSpringForce(index_i, delta_x);
			dvel_dt_prior_[index_i] += getDampingForce(index_i);
		}
		//=================================================================================================//
		SpringNormalOnSurfaceParticles::
			SpringNormalOnSurfaceParticles(SolidBody &solid_body, BodyPartByParticle &body_part, bool outer_surface, Vecd source_point, Real stiffness, Real damping_ratio)
			: PartSimpleDynamicsByParticle(solid_body, body_part), SolidDataSimple(solid_body),
			  pos_n_(particles_->pos_n_),
			  pos_0_(particles_->pos_0_),
			  n_(particles_->n_),
			  n_0_(particles_->n_0_),
			  vel_n_(particles_->vel_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  mass_(particles_->mass_),
			  apply_spring_force_to_particle_(StdLargeVec<bool>(pos_0_.size(), false))
		{
			// get the surface layer of particles
			BodySurface surface_layer(solid_body);
			// select which particles the spring is applied to
			for (size_t particle_i : surface_layer.body_part_particles_)
			{
				// vector to the source point from the particle
				Vecd vector_to_particle = source_point - pos_0_[particle_i];
				// normal of the particle
				Vecd normal = n_0_[particle_i];
				// get the cos of the angle between the vector and the normal
				Real cos_theta = getCosineOfAngleBetweenTwoVectors(vector_to_particle, normal);
				// if outer surface, the normals close an angle greater than 90°
				// if the angle is greater than 90°, we apply the spring force to the surface particle
				Real epsilon = 1e-6; // to ignore exactly perpendicular surfaces
				if (outer_surface && cos_theta < -epsilon)
				{
					apply_spring_force_to_particle_[particle_i] = true;
				}
				// if not outer surface, it's inner surface, meaning the normals close an angle smaller than 90°
				// if the angle is less than 90°, we apply the spring force to the surface particle
				if (!outer_surface && cos_theta > epsilon)
				{
					apply_spring_force_to_particle_[particle_i] = true;
				}
			}
			// scale stiffness and damping by area here, so it's not necessary in each iteration
			// we take the area of the first particle, assuming they are uniform
			Real area = std::pow(particles_->Vol_[0], 2.0 / 3.0);
			stiffness_ = stiffness * area;
			damping_coeff_ = stiffness_ * damping_ratio;

			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		Vecd SpringNormalOnSurfaceParticles::getSpringForce(size_t index_i, Vecd disp)
		{
			// normal of the particle
			Vecd normal = particles_->n_0_[index_i];
			// get the normal portion of the displacement, which is parallel to the normal of particles, meaning it is the normal vector * scalar
			Vecd normal_disp = getVectorProjectionOfVector(disp, normal);

			Vecd spring_force_vector = -stiffness_ * normal_disp;

			return spring_force_vector;
		}
		//=================================================================================================//
		Vecd SpringNormalOnSurfaceParticles::getDampingForce(size_t index_i)
		{
			// normal of the particle
			Vecd normal = particles_->n_0_[index_i];
			//velocity of the particle
			Vecd velocity_n = vel_n_[index_i];
			// get the normal portion of the velocity, which is parallel to the normal of particles, meaning it is the normal vector * scalar
			Vecd normal_vel = getVectorProjectionOfVector(velocity_n, normal);

			Vecd damping_force_vector = -damping_coeff_ * normal_vel;

			return damping_force_vector;
		}
		//=================================================================================================//
		void SpringNormalOnSurfaceParticles::Update(size_t index_i, Real dt)
		{
			try
			{
				if (apply_spring_force_to_particle_[index_i])
				{
					Vecd delta_x = pos_n_[index_i] - pos_0_[index_i];
					dvel_dt_prior_[index_i] += getSpringForce(index_i, delta_x) / mass_[index_i];
					dvel_dt_prior_[index_i] += getDampingForce(index_i) / mass_[index_i];
				}
			}
			catch (out_of_range &e)
			{
				throw runtime_error(string("SpringNormalOnSurfaceParticles::Update: particle index out of bounds") + to_string(index_i));
			}
		}
		//=================================================================================================//
		SpringOnSurfaceParticles::
			SpringOnSurfaceParticles(SolidBody &body, Real stiffness, Real damping_ratio)
			: ParticleDynamicsSimple(body), SolidDataSimple(body),
			  pos_n_(particles_->pos_n_),
			  pos_0_(particles_->pos_0_),
			  vel_n_(particles_->vel_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  mass_(particles_->mass_),
			  apply_spring_force_to_particle_(StdLargeVec<bool>(pos_0_.size(), false))
		{
			// get the surface layer of particles
			BodySurface surface_layer(body);
			// select which particles the spring is applied to
			// if the particle is in the surface layer, the force is applied
			for (size_t particle_i: surface_layer.body_part_particles_) apply_spring_force_to_particle_[particle_i] = true;

			// scale stiffness and damping by area here, so it's not necessary in each iteration
			// we take the area of the first particle, assuming they are uniform
			Real area = std::pow(particles_->Vol_[0], 2.0 / 3.0);
			stiffness_ = stiffness * area;
			damping_coeff_ = stiffness_ * damping_ratio;

			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		void SpringOnSurfaceParticles::Update(size_t index_i, Real dt)
		{
			try{
				if (apply_spring_force_to_particle_[index_i])
				{
					dvel_dt_prior_[index_i] += -stiffness_ * (pos_n_[index_i] - pos_0_[index_i]) / mass_[index_i];
					dvel_dt_prior_[index_i] += -damping_coeff_ * vel_n_[index_i]  / mass_[index_i];
				}
			}
				catch(out_of_range& e){
				throw runtime_error(string("SpringOnSurfaceParticles::Update: particle index out of bounds") + to_string(index_i));
			}
		}
		//=================================================================================================//
		AccelerationForBodyPartInBoundingBox::
			AccelerationForBodyPartInBoundingBox(SolidBody &solid_body, BoundingBox &bounding_box, Vecd acceleration)
			: ParticleDynamicsSimple(solid_body), SolidDataSimple(solid_body),
			  pos_n_(particles_->pos_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  bounding_box_(bounding_box),
			  acceleration_(acceleration)
		{
			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		void AccelerationForBodyPartInBoundingBox::Update(size_t index_i, Real dt)
		{
			try
			{
				Vecd point = pos_n_[index_i];
				if (checkIfPointInBoundingBox(point, bounding_box_))
				{
					dvel_dt_prior_[index_i] += acceleration_;
				}
			}
			catch (out_of_range &e)
			{
				throw runtime_error(string("AccelerationForBodyPartInBoundingBox::Update: particle index out of bounds") + to_string(index_i));
			}
		}
		//=================================================================================================//
		ForceInBodyRegion::
			ForceInBodyRegion(SPHBody &sph_body, BodyPartByParticle &body_part, Vecd force, Real end_time)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidDataSimple(sph_body),
			  pos_0_(particles_->pos_0_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  acceleration_(0),
			  end_time_(end_time),
			  apply_force_to_particle_(StdLargeVec<bool>(pos_0_.size(), false))
		{
			Real mass_in_region = 0.0;
			for(size_t i = 0; i < body_part_particles_.size(); i++)
			{
					int particle_ID = body_part_particles_[i];
					mass_in_region += particles_->mass_[particle_ID];
					apply_force_to_particle_[particle_ID] = true;
			};
			// calculate acceleration: force / mass of particles in region
			acceleration_ = force / mass_in_region;
			// set ghost particles to zero
			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		void ForceInBodyRegion::Update(size_t index_i, Real dt)
		{
			try
			{
				Real time_factor = SMIN(GlobalStaticVariables::physical_time_ / end_time_, 1.0);
				dvel_dt_prior_[index_i] = acceleration_ * time_factor;
			}
			catch (out_of_range &e)
			{
				throw runtime_error(string("ForceInBodyRegion::Update: particle index out of bounds") + to_string(index_i));
			}
		}
		//=================================================================================================//
		SurfacePressureFromSource::
			SurfacePressureFromSource(SPHBody &sph_body, BodyPartByParticle &body_part, Vecd source_point, StdVec<array<Real, 2>> pressure_over_time)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidDataSimple(sph_body),
			  pos_0_(particles_->pos_0_),
			  n_(particles_->n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  mass_(particles_->mass_),
			  pressure_over_time_(pressure_over_time),
			  // set apply_pressure_to_particle_ to false for each particle
			  apply_pressure_to_particle_(StdLargeVec<bool>(pos_0_.size(), false))
		{
			// get the surface layer of particles
			BodySurface surface_layer(sph_body);
			// select which particles the pressure is applied to
			for (size_t particle_i : surface_layer.body_part_particles_)
			{
				// vector to the source point from the particle
				Vecd vector_to_particle = source_point - particles_->pos_0_[particle_i];
				// normal of the particle
				Vecd normal = particles_->n_0_[particle_i];
				// get the cos of the angle between the vector and the normal
				Real cos_theta = getCosineOfAngleBetweenTwoVectors(vector_to_particle, normal);
				// if the angle is less than 90°, we apply the pressure to the surface particle
				// ignore exactly perpendicular surfaces
				if (cos_theta > 1e-6)
					apply_pressure_to_particle_[particle_i] = true;
			}

			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		Real SurfacePressureFromSource::getPressure()
		{
			// check if we have reached the max time, if yes, return the last pressure
			bool max_time_reached = GlobalStaticVariables::physical_time_ > pressure_over_time_[pressure_over_time_.size() - 1][0];
			if (max_time_reached)
				return pressure_over_time_[pressure_over_time_.size() - 1][1];

			int interval = 0;
			// find out the interval
			for (size_t i = 0; i < pressure_over_time_.size(); i++)
			{
				if (GlobalStaticVariables::physical_time_ < pressure_over_time_[i][0])
				{
					interval = i;
					break;
				}
			}
			// interval has to be at least 1
			if (interval < 1)
				throw runtime_error(string("SurfacePressureFromSource::getPressure(): pressure_over_time input not correct, should start with {0.0, 0.0}"));
			// scale the pressure to the current time
			Real t_0 = pressure_over_time_[interval - 1][0];
			Real t_1 = pressure_over_time_[interval][0];
			Real p_0 = pressure_over_time_[interval - 1][1];
			Real p_1 = pressure_over_time_[interval][1];

			return p_0 + (p_1 - p_0) * (GlobalStaticVariables::physical_time_ - t_0) / (t_1 - t_0);
		}
		//=================================================================================================//
		void SurfacePressureFromSource::Update(size_t index_i, Real dt)
		{
			try
			{
				if (apply_pressure_to_particle_[index_i])
				{
					// get the surface area of the particle, assuming it has a cubic volume
					// acceleration is particle force / particle mass
					Real area = std::pow(particles_->Vol_[index_i], 2.0 / 3.0);
					Real acc_from_pressure = getPressure() * area / mass_[index_i];
					// vector is made by multiplying it with the surface normal
					// add the acceleration to the particle
					dvel_dt_prior_[index_i] += (-1.0) * n_[index_i] * acc_from_pressure;
				}
			}
			catch (out_of_range &e)
			{
				throw runtime_error(string("SurfacePressureFromSource::Update: particle index out of bounds") + to_string(index_i));
			}
		}
		//=================================================================================================//
		ElasticDynamicsInitialCondition::
			ElasticDynamicsInitialCondition(SolidBody &solid_body)
			: ParticleDynamicsSimple(solid_body),
			  ElasticSolidDataSimple(solid_body),
			  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_)
		{
		}
		//=================================================================================================//
		UpdateElasticNormalDirection::
			UpdateElasticNormalDirection(SolidBody &solid_body)
			: ParticleDynamicsSimple(solid_body),
			  ElasticSolidDataSimple(solid_body),
			  n_(particles_->n_), n_0_(particles_->n_0_), F_(particles_->F_)
		{
		}
		//=================================================================================================//
		DeformationGradientTensorBySummation::
			DeformationGradientTensorBySummation(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  ElasticSolidDataInner(inner_relation),
			  Vol_(particles_->Vol_), pos_n_(particles_->pos_n_),
			  B_(particles_->B_), F_(particles_->F_)
		{
		}
		//=================================================================================================//
		void DeformationGradientTensorBySummation::Interaction(size_t index_i, Real dt)
		{
			Vecd &pos_n_i = pos_n_[index_i];

			Matd deformation(0.0);
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation -= Vol_[index_j] * SimTK::outer((pos_n_i - pos_n_[index_j]), gradw_ij);
			}

			F_[index_i] = B_[index_i] * deformation;
		}
		//=================================================================================================//
		BaseElasticRelaxation::
			BaseElasticRelaxation(BaseBodyRelationInner &inner_relation)
			: ParticleDynamics1Level(*inner_relation.sph_body_),
			  ElasticSolidDataInner(inner_relation), Vol_(particles_->Vol_),
			  rho_n_(particles_->rho_n_), mass_(particles_->mass_),
			  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_) {}
		//=================================================================================================//
		StressRelaxationFirstHalf::
			StressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation)
			: BaseElasticRelaxation(inner_relation),
			  dvel_dt_prior_(particles_->dvel_dt_prior_), force_from_fluid_(particles_->force_from_fluid_),
			  stress_PK1_(particles_->stress_PK1_)
		{
			rho0_ = material_->ReferenceDensity();
			inv_rho0_ = 1.0 / rho0_;
			smoothing_length_ = sph_adaptation_->ReferenceSmoothingLength();
			numerical_dissipation_factor_ = 0.25;
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] = rho0_ / det(F_[index_i]);
			//obtain the first Piola-Kirchhoff stress from the second Piola-Kirchhoff stress
			//it seems using reproducing correction here increases convergence rate near the free surface
			stress_PK1_[index_i] = F_[index_i] * material_->ConstitutiveRelation(F_[index_i], index_i) * B_[index_i];
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Interaction(size_t index_i, Real dt)
		{
			//including gravity and force from fluid
			Vecd acceleration = dvel_dt_prior_[index_i] + force_from_fluid_[index_i] / mass_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];
				Real dim_r_ij_1 = Dimensions / r_ij;
				Vecd pos_jump = pos_n_[index_i] - pos_n_[index_j];
				Vecd vel_jump = vel_n_[index_i] - vel_n_[index_j];
				Real strain_rate = SimTK::dot(pos_jump, vel_jump) * dim_r_ij_1 * dim_r_ij_1;
				Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
				Matd numerical_stress_ij =
					0.5 * (F_[index_i] + F_[index_j]) * material_->PairNumericalDamping(strain_rate, smoothing_length_);
				acceleration += (stress_PK1_[index_i] + stress_PK1_[index_j] + numerical_dissipation_factor_ * weight * numerical_stress_ij) *
								inner_neighborhood.dW_ij_[n] * e_ij * Vol_[index_j] * inv_rho0_;
			}

			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		void StressRelaxationFirstHalf::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += dvel_dt_[index_i] * dt;
		}
		//=================================================================================================//
		KirchhoffParticleStressRelaxationFirstHalf::
			KirchhoffParticleStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation)
			: StressRelaxationFirstHalf(inner_relation){};
		//=================================================================================================//
		void KirchhoffParticleStressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] = rho0_ / det(F_[index_i]);
			Real J = det(F_[index_i]);
			Real one_over_J = 1.0 / J;
			rho_n_[index_i] = rho0_ * one_over_J;
			Real J_to_minus_2_over_dimension = pow(one_over_J, 2.0 * one_over_dimensions_);
			Matd normalized_b = (F_[index_i] * ~F_[index_i]) * J_to_minus_2_over_dimension;
			Matd deviatoric_b = normalized_b - Matd(1.0) * normalized_b.trace() * one_over_dimensions_;
			Matd inverse_F_T = ~SimTK::inverse(F_[index_i]);
			//obtain the first Piola-Kirchhoff stress from the Kirchhoff stress
			//it seems using reproducing correction here increases convergence rate
			//near the free surface however, this correction is not used for the numerical disspation
			stress_PK1_[index_i] = (Matd(1.0) * material_->VolumetricKirchhoff(J) +
									material_->DeviatoricKirchhoff(deviatoric_b)) *
								   inverse_F_T * B_[index_i];
		}
		//=================================================================================================//
		KirchhoffStressRelaxationFirstHalf::
			KirchhoffStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation)
			: StressRelaxationFirstHalf(inner_relation)
		{
			particles_->registerAVariable<indexScalar, Real>(J_to_minus_2_over_dimension_, "DeterminantTerm");
			particles_->registerAVariable<indexMatrix, Matd>(stress_on_particle_, "StressOnParticle");
			particles_->registerAVariable<indexMatrix, Matd>(inverse_F_T_, "InverseTransposedDeformation");
		};
		//=================================================================================================//
		void KirchhoffStressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			Real J = det(F_[index_i]);
			//throw an exception if the determinant becomes negative
			if (J <= 0) throw std::runtime_error(std::string("Determinant of F_ became negative! SPHBody: ") + this->body_->getBodyName()
				+ " particle ID: " + std::to_string(index_i));
			Real one_over_J = 1.0 / J;
			rho_n_[index_i] = rho0_ * one_over_J;
			J_to_minus_2_over_dimension_[index_i] = pow(one_over_J, 2.0 * one_over_dimensions_);
			inverse_F_T_[index_i] = ~SimTK::inverse(F_[index_i]);
			stress_on_particle_[index_i] =
				inverse_F_T_[index_i] * (material_->VolumetricKirchhoff(J) - 
					correction_factor_ * material_->ShearModulus() * J_to_minus_2_over_dimension_[index_i] * (F_[index_i] * ~F_[index_i]).trace() * one_over_dimensions_) +
				material_->NumericalDampingLeftCauchy(F_[index_i], dF_dt_[index_i], smoothing_length_, index_i) * inverse_F_T_[index_i];
			stress_PK1_[index_i] = F_[index_i] * material_->ConstitutiveRelation(F_[index_i], index_i);
			
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					if (std::isnan(stress_on_particle_[index_i][i][j])){
						throw std::runtime_error(std::string("stress_on_particle_ is Not A Number"));
					}
				}
			}
		}
		//=================================================================================================//
		void KirchhoffStressRelaxationFirstHalf::Interaction(size_t index_i, Real dt)
		{
			//including gravity and force from fluid
			Vecd acceleration = dvel_dt_prior_[index_i] + force_from_fluid_[index_i] / mass_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd shear_force_ij = correction_factor_ * material_->ShearModulus() *
									  (J_to_minus_2_over_dimension_[index_i] + J_to_minus_2_over_dimension_[index_j]) *
									  (pos_n_[index_i] - pos_n_[index_j]) / inner_neighborhood.r_ij_[n];
				acceleration += ((stress_on_particle_[index_i] + stress_on_particle_[index_j]) * inner_neighborhood.e_ij_[n] + shear_force_ij) *
								inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inv_rho0_;
			}
			dvel_dt_[index_i] = acceleration;

			for (int i = 0; i < 3; i++){
				if (std::isnan(acceleration[i])){
					throw std::runtime_error(std::string("acceleration is Not A Number! SPHBody: ") + this->body_->getBodyName() 
						+ " particle ID: " + std::to_string(index_i));
				}
			}
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Interaction(size_t index_i, Real dt)
		{
			const Vecd &vel_n_i = vel_n_[index_i];

			Matd deformation_gradient_change_rate(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation_gradient_change_rate -=
					Vol_[index_j] * SimTK::outer((vel_n_i - vel_n_[index_j]), gradw_ij);
			}

			dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
		}
		//=================================================================================================//
		void StressRelaxationSecondHalf::Update(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
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
			  SolidDataSimple(solid_body),
			  force_from_fluid_(particles_->force_from_fluid_), contact_force_(particles_->contact_force_),
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
