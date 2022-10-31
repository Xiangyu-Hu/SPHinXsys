/**
 * @file 	3d_roof_analytical.cpp
 * @brief 	Shell verificaiton  incl. refinement study
 * @details Roof shell verification case with relaxed shell particles
 * @author 	Bence Rochlitz
 * @ref 	ANSYS Workbench Verification Manual, Release 15.0, November 2013, VMMECH069: Barrel Vault Roof Under Self Weight
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>

using namespace SPH;

Real to_rad(Real angle){return angle*Pi/180;}
static const int simtk_res = 1;

void relax_shell(RealBody& plate_body, Real thickness, Real level_set_refinement_ratio)
{
	// BUG: apparently only works if dp > thickness, otherwise ShellNormalDirectionPrediction::correctNormalDirection() throws error

	InnerRelation imported_model_inner(plate_body);
	SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(plate_body);
	relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(imported_model_inner, thickness, level_set_refinement_ratio);
	relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(imported_model_inner, thickness);	
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.mid_surface_bounding_.parallel_exec();
	plate_body.updateCellLinkedList();
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
		}
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
	}
	shell_normal_prediction.exec();
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;
}

class TimeDependentExternalForce : public Gravity
{
	Real time_to_full_external_force_;
public:
	explicit TimeDependentExternalForce(Vec3d external_force, Real time_to_full_external_force)
		: Gravity(external_force),
		time_to_full_external_force_(time_to_full_external_force)
		{};
	virtual Vec3d InducedAcceleration(Vec3d &position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_external_force_
				   ? current_time * global_acceleration_ / time_to_full_external_force_
				   : global_acceleration_;
	}
};

class ShellRoofParticleGenerator : public SurfaceParticleGenerator
{
	const StdVec<Vec3d>& pos_0_;
	const Vec3d center_;
	const Real dp_square_;
	const Real thickness_;
public:
	explicit ShellRoofParticleGenerator(SPHBody &sph_body, const StdVec<Vec3d>& pos_0, const Vec3d& center, Real dp, Real thickness)
		: SurfaceParticleGenerator(sph_body),
		pos_0_(pos_0),
		center_(center),
		dp_square_(dp*dp),
		thickness_(thickness)
		{};
	virtual void initializeGeometricVariables() override
	{
		for (const auto& pos: pos_0_)
		{
			// creating the normal direction - z coordinate is always zero
			Vec3d center_to_pos = pos-center_;
			center_to_pos[2] = 0;
			initializePositionAndVolumetricMeasure(pos, dp_square_);
			initializeSurfaceProperties(center_to_pos.normalize(), thickness_);
		}
	}
};

template<typename VectorType>
SPH::BoundingBox get_particles_bounding_box(const VectorType& pos_0)
{
    SPH::Vec3d lower(pos_0[0]);
    SPH::Vec3d upper(pos_0[0]);
    for (const auto& pos: pos_0){
        for (int i=0; i<3; i++){
            if (lower[i] > pos[i]) lower[i] = pos[i];
            if (upper[i] < pos[i]) upper[i] = pos[i];
        }
    }
    return SPH::BoundingBox(lower, upper);
}

StdVec<Vec3d> read_obj_vertices(const std::string& file_name)
{
	std::cout << "read_obj_vertices started" << std::endl;
	StdVec<Vec3d> pos_0;
    std::ifstream myfile(file_name, std::ios_base::in);
	if(!myfile.is_open()) throw std::runtime_error("read_obj_vertices: file doesn't exist");
    while (!myfile.eof())
    {
		Real x, y, z;
		myfile >> x;
		myfile >> y;
		myfile >> z;
		pos_0.push_back({x,y,z});
    }
	std::cout << "read_obj_vertices finished" << std::endl;
	return pos_0;
}

int main(int ac, char *av[])
{
	// main geometric parameters
	Vec3d tangential_vec(1,0,0);
	Vec3d radial_vec(0,1,0);
	Vec3d length_vec(0,0,1);
	unsigned int tangential_axis = 0;
	unsigned int radial_axis = 1;
	unsigned int length_axis = 2;
	Real radius = 25;
	Real length = 50;
	Real thickness = 0.25;
	Real teta = 40;
	Real arc = radius*to_rad(teta);
	Vec3d center(0,-radius,0);
	// resolution
	Real dp = 1;
	// material
	Real rho = 36.7347;
	Real E = 4.32e8;
	Real mu = 0.3;
	auto material = makeShared<SaintVenantKirchhoffSolid>(rho, E, mu);
	Real physical_viscosity = 7e3; // where is this value coming from?
	std::cout << "physical_viscosity: " << physical_viscosity << std::endl;
	// gravity
	Vec3d gravity = -9.8066*radial_vec;
	Real time_to_full_external_force = 0.1;
	// system bounding box
	BoundingBox bb_system;

	// Option A: generating particles from plate and warping
	// BUG: edge particle density is too high, ShellNormalDirectionPrediction randomly throws error or not
	auto shell_particles_warped_plate = [&]()
	{// generate particle positions
		// 1. to create the roof geometry we create the initial shell particles based on a plate
		// 2. then we warp the plate according to the given radius and angle

		// 1. Plate
		// fake thickness for fast levelset generation - we just need particle positions
		Real thickness_temp = dp*2;
		Real level_set_refinement_ratio = dp / (0.1 * thickness_temp);
		// plate dimensions
		Vec3d plate_dim(2*arc+dp, thickness_temp, length+dp);
		// shape
		auto plate_shape = makeShared<ComplexShape>("plate_shape");
		plate_shape->add<TriangleMeshShapeBrick>(plate_dim/2, simtk_res, Vec3d(0));
		SPHSystem system(plate_shape->getBounds(), dp);
		// body and particles
		RealBody plate_body(system, plate_shape);
		plate_body.defineBodyLevelSetShape(level_set_refinement_ratio)->correctLevelSetSign();
		plate_body.defineParticlesWithMaterial<ShellParticles>(material.get());
		plate_body.generateParticles<ThickSurfaceParticleGeneratorLattice>(thickness_temp);
		relax_shell(plate_body, thickness_temp, level_set_refinement_ratio);
		// output
		IOEnvironment io_env(system);
		BodyStatesRecordingToVtp vtp_output(io_env, {plate_body});
		vtp_output.writeToFile(0);

		// 2. Warping
		// creating a curved shell from the plate based on the arc length at each particle
		auto warping_transform = [&](const Vec3d& pos)
		{
			// take the x coordinate as arc length
			Real angle = pos[tangential_axis]/radius;
			// transfer x, y coordinates, keep z as is
			Real x_new = radius*std::sin(angle);
			Real y_new = radius*std::cos(angle)-radius;
			return Vec3d(x_new, y_new, pos[length_axis]);
		};
		StdVec<Vec3d> pos_0;
		pos_0.reserve(plate_body.getBaseParticles().pos_.size());
		for (const auto& pos: plate_body.getBaseParticles().pos_) pos_0.push_back(warping_transform(pos));
		// update bb_system and return
		bb_system = get_particles_bounding_box(pos_0);
		return pos_0;
	};

	// Option B: generating particles from stl
	// BUG: edge particle density is too high
	auto shell_particles_stl = [&]()
	{
		Real thickness_temp = 4; // 2 or 4
		std::string stl_path = "input/shell_50mm_80d_" + std::to_string(int(thickness_temp)) + "mm.stl";
		Real level_set_refinement_ratio = dp / (0.1 * thickness_temp);
		// shape
		auto plate_shape = makeShared<ComplexShape>("plate_shape_stl");
		plate_shape->add<TriangleMeshShapeSTL>(stl_path, Vec3d(0), 1);
		bb_system = plate_shape->getBounds();
		SPHSystem system(bb_system, dp);
		// body and particles
		RealBody plate_body(system, plate_shape);
		plate_body.defineBodyLevelSetShape(level_set_refinement_ratio)->correctLevelSetSign();
		plate_body.defineParticlesWithMaterial<ShellParticles>(material.get());
		plate_body.generateParticles<ThickSurfaceParticleGeneratorLattice>(thickness_temp);
		relax_shell(plate_body, thickness_temp, level_set_refinement_ratio);
		// output
		IOEnvironment io_env(system, false);
		BodyStatesRecordingToVtp vtp_output(io_env, {plate_body});
		vtp_output.writeToFile(0);

		return plate_body.getBaseParticles().pos_;
	};

	// Option C: generating particles from predefined positions from obj file
	StdVec<Vec3d> obj_vertices = read_obj_vertices("input/shell_50mm_80d_1mm.txt"); // dp = 1
	// transform to flat plate
	auto flatten_transform = [&](StdVec<Vec3d>& pos_0)
	{// editing input vector - non-const input!
		for (auto& pos: pos_0)
		{
			// get the arc length at the x coordinate and set it as new x coordinate
			Real angle = std::asin(pos[tangential_axis]/radius);
			pos[tangential_axis] = angle*radius;
			// y will be 0
			pos[radial_axis] = 0;
			// z stays as is
		}
	};
	// flatten_transform(obj_vertices);
	// find out BoundingBox
	BoundingBox obj_vertices_bb = get_particles_bounding_box(obj_vertices); // store this
	bb_system = obj_vertices_bb;
	// just making sure nothing leaves the bounding box
	bb_system.first += Vec3d(-5);
	bb_system.second += Vec3d(5);
	std::cout << "bb_system.first: " << bb_system.first << std::endl;
	std::cout << "bb_system.second: " << bb_system.second << std::endl;

	// shell
	auto shell_shape = makeShared<ComplexShape>("shell_shape");

	// starting the actual simulation
	SPHSystem system(bb_system, dp);
	SolidBody shell_body(system, shell_shape);
	shell_body.defineParticlesWithMaterial<ShellParticles>(material.get());
	shell_body.generateParticles<ShellRoofParticleGenerator>(obj_vertices, center, dp, thickness);
	auto shell_particles = dynamic_cast<ShellParticles*>(&shell_body.getBaseParticles());
	// output
	IOEnvironment io_env(system, true);
	shell_body.addBodyStateForRecording<Vec3d>("NormalDirection");
	BodyStatesRecordingToVtp vtp_output(io_env, {shell_body});
	vtp_output.writeToFile(0);

	// methods
	InnerRelation shell_body_inner(shell_body);
	SimpleDynamics<TimeStepInitialization> initialize_external_force(shell_body, makeShared<Gravity>(gravity));
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(shell_body_inner);
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(shell_body);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(shell_body_inner, 3, false);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(shell_body_inner);

	BodyPartByParticle constrained_edges(shell_body, "constrained_edges");
	auto constrained_edge_ids = [&]()
	{// brute force finding the edges
		IndexVector ids;
		for (size_t i = 0; i < shell_body.getBaseParticles().pos_.size(); ++i)
			if (shell_body.getBaseParticles().pos_[i][length_axis] < obj_vertices_bb.first[length_axis]+dp/2 ||
				shell_body.getBaseParticles().pos_[i][length_axis] > obj_vertices_bb.second[length_axis]-dp/2)
					ids.push_back(i);
		return ids;
	}();
	constrained_edges.body_part_particles_ = constrained_edge_ids;

	SimpleDynamics<solid_dynamics::FixedInAxisDirection, BodyPartByParticle> constrain_holder(constrained_edges, length_vec);
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
		shell_position_damping(0.2, shell_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
		shell_rotation_damping(0.2, shell_body_inner, "AngularVelocity", physical_viscosity);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.parallel_exec();

	/**
	 * From here the time stepping begins.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	int ite = 0;
	Real end_time = 2.0;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	tick_count t1 = tick_count::now();
	/**
	 * Main loop
	 */
	Real max_dt = 0.0;
	try
	{
		while (GlobalStaticVariables::physical_time_ < end_time)
		{
			Real integral_time = 0.0;
			while (integral_time < output_period)
			{
				if (ite % 1000 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							<< GlobalStaticVariables::physical_time_ << "	dt: "
							<< dt << "\n";
				}

				initialize_external_force.parallel_exec(dt);

				dt = 0.2 * computing_time_step_size.parallel_exec();
				{// checking for excessive time step reduction
					if (dt > max_dt) max_dt = dt;
					if (dt < max_dt/1e3) throw std::runtime_error("time step decreased too much");
				}
				
				stress_relaxation_first_half.parallel_exec(dt);
				constrain_holder.parallel_exec();
				// shell_velocity_damping.parallel_exec(dt);
				// shell_rotation_damping.parallel_exec(dt);
				constrain_holder.parallel_exec();
				stress_relaxation_second_half.parallel_exec(dt);

				++ite;
				integral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				// shell_body.updateCellLinkedList();

				{// checking if any position has become nan
					// BUG: damping throws nan error it seems - regardless of particle generation method
					for (const auto& pos: shell_body.getBaseParticles().pos_)
						if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
							throw std::runtime_error("position has become nan");
				}
			}
			{// output displacement
				std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
			}
			vtp_output.writeToFile(ite);
		}
		tick_count::interval_t tt = tick_count::now()-t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		{// output displacement
			auto shell_particles = dynamic_cast<ElasticSolidParticles*>(&shell_body.getBaseParticles());
			std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
		}
		vtp_output.writeToFile(ite);
	}

	return 0;
}
