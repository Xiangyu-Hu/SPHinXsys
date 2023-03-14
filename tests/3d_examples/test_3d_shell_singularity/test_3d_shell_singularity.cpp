/**
 * @file 	test_3d_sphere_compression.cpp
 * @brief 	Shell verification  incl. refinement study
 * @details Circular plastic shell verification case with relaxed shell particles
 * @author 	Bence Rochlitz
 * @ref 	ANSYS Workbench Verification Manual, Release 15.0, November 2013, VMMECH051: Bending of a Circular Plate Using Axisymmetric Elements
 */

#include "sphinxsys.h"
#include <numeric>
#include <gtest/gtest.h>

using namespace SPH;

class ShellPlateParticleGenerator : public SurfaceParticleGenerator
{
	const StdVec<Vec3d>& pos_0_;
	const StdVec<Vec3d>& n_0_;
	const double particle_area_;
	const double thickness_;
public:
	ShellPlateParticleGenerator(SPHBody &sph_body, const StdVec<Vec3d>& pos_0, const StdVec<Vec3d>& n_0, const double particle_area, const double thickness)
		: SurfaceParticleGenerator(sph_body),
		pos_0_(pos_0),
		n_0_(n_0),
		particle_area_(particle_area),
		thickness_(thickness)
	{};
	virtual void initializeGeometricVariables() override
	{
		for (size_t i=0; i<pos_0_.size(); ++i)
		{
			initializePositionAndVolumetricMeasure(pos_0_[i], particle_area_);
			initializeSurfaceProperties(n_0_[i], thickness_);
		}
	}
};

template<typename VectorType>
BoundingBox get_particles_bounding_box(const VectorType& pos_0)
{
    Vec3d lower(pos_0[0]);
    Vec3d upper(pos_0[0]);
    for (const auto& pos: pos_0){
        for (int i=0; i<3; i++){
            if (lower[i] > pos[i]) lower[i] = pos[i];
            if (upper[i] < pos[i]) upper[i] = pos[i];
        }
    }
    return BoundingBox(lower, upper);
}

void read_obj_vertices(StdVec<Vec3d>& pos_0, StdVec<Vec3d>& n_0, const std::string& file_name)
{
	std::cout << "read_obj_vertices started" << std::endl;

    std::ifstream myfile(file_name, std::ios_base::in);
	if(!myfile.is_open()) throw std::runtime_error("read_obj_vertices: file doesn't exist");

	Vec3d particle(0);
	unsigned int count = 0;
	double value = 0;

	std::vector<double> vec;
	vec.reserve(6);

    while (myfile >> value)
    {
        vec.push_back(value);
		++count;
		if (vec.size() % 6 == 0)
		{
			pos_0.emplace_back(vec[0], vec[1], vec[2]);
			n_0.emplace_back(vec[3], vec[4], vec[5]);
			vec.clear();
		}
    }

	std::cout << "read_obj_vertices finished" << std::endl;
}

const double get_physical_viscosity_general(const double rho, const double youngs_modulus, const double length_scale, const double shape_constant = 0.4)
{
	// the physical viscosity is defined in the paper of prof. Hu
	// https://arxiv.org/pdf/2103.08932.pdf
	// physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
	// beta: shape constant (0.4 for beam)
	return shape_constant/ 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

namespace unit_system
{// Units: g, mm, ms, N, MPa

	// MASS	LENGTH	TIME	FORCE	STRESS	ENERGY	DENSITY		YOUNG's		35MPH 56.33KMPH		GRAVITY
	// g	mm		ms		N		MPa		N-mm	7.83e-03	2.07e+05	15.65				9.806e-03
	// https://www.dynasupport.com/howtos/general/consistent-units
	
	constexpr double mass_unit = 1e3; // kg to g
	constexpr double length_unit = 1e3; // m to mm
	constexpr double time_unit = 1e3; // s to ms
	constexpr double pressure_unit = 1e-6; // Pa to MPa
	constexpr double density_unit = 1e-6; // kg/m^3 to g/mm^3
	constexpr double acceleration_unit = 1e-3; // m/s^2 to mm/ms^2
}

void shell_bending()
{
	// main geometric parameters
	const double scale = 1;
	const double thickness = 0.2*scale; // 1 mm

	// boundary condition parameters
	const double acc_max = 20 * unit_system::acceleration_unit;
	const Vec3d center = Vec3d(16, 20, 27.5) * scale;
	const double height_fix = 8 * scale;
	const double height_pressure_min = 20 * scale;
	const double end_time = 0.1 * unit_system::time_unit;

	// resolution
	const double dp = thickness;
	const double total_area = 858.7; // measured in Meshmixer
	std::cout << "total_area: " << total_area << std::endl;
	// material
	const double rho = 1e3 * unit_system::density_unit;
	const double E = 4e5 * unit_system::pressure_unit;
	const double mu = 0.3;
	auto material = makeShared<LinearElasticSolid>(rho, E, mu);
	double physical_viscosity = 7e3;
	std::cout << "physical_viscosity: " << physical_viscosity << std::endl;
	physical_viscosity = get_physical_viscosity_general(rho, E, thickness);
	std::cout << "physical_viscosity: " << physical_viscosity << std::endl;

	// system bounding box
	BoundingBox bb_system;

	// generating particles from predefined positions from obj file
	StdVec<Vec3d> obj_vertices, obj_normals;
	read_obj_vertices(obj_vertices, obj_normals, "input/rose_pos_normals_0_4_mm.txt");

	std::for_each(obj_vertices.begin(), obj_vertices.end(), [&](Vec3d &vec){ vec *= scale; });
	const double particle_area = total_area / obj_vertices.size();
	// find out BoundingBox
	bb_system = get_particles_bounding_box(obj_vertices);
	std::cout << "bb_system.first_: " << bb_system.first_ << std::endl;
	std::cout << "bb_system.second_: " << bb_system.second_ << std::endl;

	// shell
	auto shell_shape = makeShared<ComplexShape>("shell_shape"); // keep all data for parameter study

	// starting the actual simulation
	SPHSystem system(bb_system, dp);
	SolidBody shell_body(system, shell_shape);
	shell_body.defineParticlesWithMaterial<ShellParticles>(material.get());
	shell_body.generateParticles<ShellPlateParticleGenerator>(obj_vertices, obj_normals, particle_area, thickness);
	auto shell_particles = dynamic_cast<ShellParticles*>(&shell_body.getBaseParticles());
	// output_opening
	IOEnvironment io_env(system, false);
	shell_body.addBodyStateForRecording<Vec3d>("NormalDirection");
	BodyStatesRecordingToVtp vtp_output(io_env, {shell_body});
	vtp_output.writeToFile(0);

	// methods
	InnerRelation shell_body_inner(shell_body);
	SimpleDynamics<TimeStepInitialization> initialize_external_force(shell_body);
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(shell_body_inner);
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(shell_body);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(shell_body_inner, 3, true);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(shell_body_inner);

	// pressure boundary condition
	auto apply_pressure = [&](ElasticSolidParticles& particles)
	{
		const double time_ratio = GlobalStaticVariables::physical_time_/end_time;
		const double acc_temp = acc_max * time_ratio;

		// force application
		for(size_t i=0; i<particles.total_real_particles_; ++i)
		{
			if (particles.pos0_[i][1] < height_pressure_min)
				continue;
			
			particles.acc_prior_[i] += acc_temp * particles.n_[i];
		}
	};

	BodyPartByParticle constrained_edges(shell_body, "constrained_edges");
	auto constrained_edge_ids = [&]()
	{// brute force finding the edges
		IndexVector ids;
		for (size_t i = 0; i < shell_body.getBaseParticles().pos_.size(); ++i)
			if (shell_body.getBaseParticles().pos_[i][1] < height_fix)
				ids.push_back(i);
		return ids;
	}();
	constrained_edges.body_part_particles_ = constrained_edge_ids;

	SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constrain_holder(constrained_edges);

	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
		shell_velocity_damping(0.2, shell_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
		shell_rotation_damping(0.2, shell_body_inner, "AngularVelocity", physical_viscosity);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.parallel_exec();

	{// tests on initialization
		{// checking particle distances - avoid bugs of reading file
			double min_rij = Infinity;
			double max_rij = 0;
			for (size_t i = 0; i < shell_particles->pos0_.size(); ++i)
			{
				Neighborhood &inner_neighborhood = shell_body_inner.inner_configuration_[i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					const double r_ij = inner_neighborhood.r_ij_[n];
					min_rij = std::min(min_rij, r_ij);
					max_rij = std::max(max_rij, r_ij);
				}
			}
			std::cout << "min_rij: " << min_rij << std::endl;
			std::cout << "max_rij: " << max_rij << std::endl;
			EXPECT_GT(min_rij, dp/2);
		}

		// test volume
		const double total_volume = std::accumulate(shell_particles->Vol_.begin(), shell_particles->Vol_.end(), 0.0);
		std::cout << "total_volume: " << total_volume << std::endl;
		const double total_mass = std::accumulate(shell_particles->mass_.begin(), shell_particles->mass_.end(), 0.0);
		std::cout << "total_mass: " << total_mass << std::endl;
		EXPECT_FLOAT_EQ(total_volume, total_area);
		EXPECT_FLOAT_EQ(total_mass, total_area*rho);
	}

	/**
	 * From here the time stepping begins.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	int ite = 0;
	const double output_period = end_time / 25.0;
	double dt = 0.0;
	tick_count t1 = tick_count::now();
	/**
	 * Main loop
	 */
	double max_dt = 0.0;
	// recording - not pushed to GitHub due to lack of matplotlib there
	StdVec<double> time, max_displacement, center_deflection;
	try
	{
		while (GlobalStaticVariables::physical_time_ < end_time)
		{
			double integral_time = 0.0;
			while (integral_time < output_period)
			{
				if (ite % 1000 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							<< GlobalStaticVariables::physical_time_ << "	dt: "
							<< dt << "\n";
				}

				initialize_external_force.parallel_exec(dt);

				apply_pressure(*shell_particles);

				dt = computing_time_step_size.parallel_exec();
				{// checking for excessive time step reduction
					if (dt > max_dt) max_dt = dt;
					if (dt < max_dt/1e3) throw std::runtime_error("time step decreased too much, iteration: " + std::to_string(ite));
				}
				
				stress_relaxation_first_half.parallel_exec(dt);
				constrain_holder.parallel_exec();
				shell_velocity_damping.parallel_exec(dt);
				shell_rotation_damping.parallel_exec(dt);
				constrain_holder.parallel_exec();
				stress_relaxation_second_half.parallel_exec(dt);

				++ite;
				integral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				shell_body.updateCellLinkedList();

				{// checking if any position has become nan
					for (const auto& pos: shell_body.getBaseParticles().pos_)
						if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
							throw std::runtime_error("position has become nan, iteration: " + std::to_string(ite));
				}
			}
			{// output data
				std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
				vtp_output.writeToFile(ite);
			}
			{// recording - not pushed to GitHub due to lack of matplotlib there
				time.push_back(GlobalStaticVariables::physical_time_);
				max_displacement.push_back(shell_particles->getMaxDisplacement());
			}
		}
		tick_count::interval_t tt = tick_count::now()-t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
		std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
	}
	catch(const std::exception& e)
	{
		std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
		vtp_output.writeToFile(ite);
		throw std::runtime_error(e.what());
	}
}

TEST(shell_bending, half_sphere)
{
	fs::remove_all("output");
	fs::create_directory("output");

	EXPECT_NO_THROW(shell_bending());
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
