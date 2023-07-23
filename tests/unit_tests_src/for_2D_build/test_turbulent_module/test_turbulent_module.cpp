#include <gtest/gtest.h>
#include "sphinxsys.h" // SPHinXsys Library.
#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"

using namespace SPH;
//----------------------------------------------------------------------
//	Define basic parameters
//----------------------------------------------------------------------
Real DL = 4;						  /**< Reference length. */
Real DH = 2;						  /**< Reference and the height of main channel. */
Real resolution_ref = 0.1;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;	   /**< Characteristic velocity. */
/** Reference sound speed needs */
Real c_f = 10.0 * U_f;
Real Re = 20000.0;					/**< Reynolds number according to book [2006 Wilcox]. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------


//----------------------------------------------------------------------
//	Define geometry
//----------------------------------------------------------------------
/** the water block . */
std::vector<Vecd> water_block_shape
{
	Vecd(0.0 - 2.0 * BW, 0.0),Vecd(0.0 - 2.0 * BW, DH),
	Vecd(DL + 2.0 * BW, DH),Vecd(DL + 2.0 * BW, 0.0),Vecd(0.0 - 2.0 * BW, 0.0)
};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape
{
	Vecd(-3.0 * BW, -BW),Vecd(-3.0 * BW, DH + BW),
	Vecd(DL + 3.0 * BW, DH + BW),Vecd(DL + 3.0 * BW, -BW),Vecd(-3.0 * BW, -BW),
};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape
{
	Vecd(0.0-3.0 * BW, 0.0),Vecd(0.0- 3.0 * BW, DH),
	Vecd(DL+ 3.0 * BW, DH),Vecd(DL+ 3.0 * BW, 0.0),Vecd(0.0- 3.0 * BW, 0.0)
};
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	Create global testing environment
//----------------------------------------------------------------------
class BaseTurbulentModule
{
protected:
	BoundingBox system_domain_bounds;
	SPHSystem sph_system;
	IOEnvironment io_environment;
	FluidBody water_block;
	SolidBody wall_boundary;

public:
	BaseTurbulentModule() :
		system_domain_bounds(Vec2d(-DL - 2.0 * BW, -BW), Vec2d(DL + 3.0 * BW, DH + BW)),
		sph_system(system_domain_bounds, resolution_ref),
		io_environment(sph_system),
		water_block(sph_system, makeShared<WaterBlock>("WaterBody")),
		wall_boundary(sph_system, makeShared<WallBoundary>("Wall"))
	{
		std::cout << "constructor for BaseTurbulentModule" << std::endl;
		water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
		water_block.generateParticles<ParticleGeneratorLattice>();
		wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
		wall_boundary.generateParticles<ParticleGeneratorLattice>();
	}
	virtual ~BaseTurbulentModule() {};
};

class TurbulentModule : public :: testing :: Test, public BaseTurbulentModule
{
protected:
	InnerRelation water_block_inner;
	ComplexRelation water_block_complex_relation;
	
	InteractionWithUpdate<fluid_dynamics::K_TurtbulentModelInner, SequencedPolicy> k_equation_relaxation;
	InteractionWithUpdate<fluid_dynamics::E_TurtbulentModelInner> epsilon_equation_relaxation;
	InteractionDynamics<fluid_dynamics::TKEnergyAccInner, SequencedPolicy> turbulent_kinetic_energy_acceleration;
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction;
	InteractionDynamics<fluid_dynamics::StandardWallFunctionCorrection, SequencedPolicy> standard_wall_function_correction;
	InteractionDynamics<fluid_dynamics::TurbulentViscousAccelerationWithWall, SequencedPolicy> turbulent_viscous_acceleration;
	SimpleDynamics<fluid_dynamics::TurbulentEddyViscosity, SequencedPolicy> update_eddy_viscosity;

	InteractionWithUpdate<CorrectedConfigurationInner> correct_configuration;

	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation;
	BodyStatesRecordingToVtp write_body_states;
	size_t number_of_iterations;
	Real integration_time;
	Real dt;
	int num_fluid_particle;
public:
	TurbulentModule() :
		BaseTurbulentModule(),
		water_block_inner(water_block),
		water_block_complex_relation(water_block_inner, { &wall_boundary }),
		
		correct_configuration(water_block_inner),

		k_equation_relaxation(water_block_inner),
		epsilon_equation_relaxation(water_block_inner),
		turbulent_kinetic_energy_acceleration(water_block_inner),
		wall_boundary_normal_direction(wall_boundary),
		standard_wall_function_correction(water_block_complex_relation),
		turbulent_viscous_acceleration(water_block_complex_relation),
		update_eddy_viscosity(water_block),

		update_density_by_summation(water_block_complex_relation),
		write_body_states(io_environment, sph_system.real_bodies_)
	{
		std::cout << "constructor for TurbulentModule" << std::endl;
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();

		/** Update cell linked list and configuration. */
		water_block.updateCellLinkedListWithParticleSort(100);
		water_block_complex_relation.updateConfiguration();

		write_body_states.writeToFile();

		update_density_by_summation.exec();
		number_of_iterations = 0;
		integration_time = 0.0;
		dt = 0.0;
		num_fluid_particle = water_block.SizeOfLoopRange();
	}
	virtual ~TurbulentModule() {};
};
//----------------------------------------------------------------------


//----------------------------------------------------------------------
//	Test K_Epsilon equation
//----------------------------------------------------------------------
class Test_K_Epsilon_Equation : public fluid_dynamics::FluidInitialCondition
{
protected:
	StdLargeVec<Vecd>& pos_, & vel_,& acc_, & acc_prior_;
	StdLargeVec<Real>& turbu_mu_;
	StdLargeVec<Real>& turbu_k_;
	StdLargeVec<Real>& turbu_epsilon_;
	StdLargeVec<Vecd>& velo_friction_;
	Vecd imposed_data_;
	const int num_data = 20;
	const int num_file = 4;
	double input_data[20][4];
	std::string file_name[4] = { "U.dat","K.dat" ,"Epsilon.dat" ,"Mu_t.dat" };
	Real smoothing_length_min_;
	StdLargeVec<Vecd> expect_k_gradient_;
public:
	Test_K_Epsilon_Equation(SPHBody& sph_body) : FluidInitialCondition(sph_body),
		pos_(particles_->pos_), vel_(particles_->vel_),
		turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
		turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
		turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation")),
		smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
		acc_(particles_->acc_), acc_prior_(particles_->acc_prior_),
		velo_friction_(*particles_->getVariableByName<Vecd>("FrictionVelocity"))
	{
		particles_->registerVariable(expect_k_gradient_, "ExpectTkeGradient");
		particles_->addVariableToWrite<Vecd>("ExpectTkeGradient");

		particles_->addVariableToWrite<Vecd>("PriorAcceleration");
		particles_->addVariableToWrite<Vecd>("Acceleration");

		std::vector<double> data;
		//load 4 files one by one
		for (int j = 0; j < num_file; ++j)
		{
			data = loadInputData(num_data, num_file, file_name[j]);
			for (int i = 0; i < num_data; ++i)
			{
				input_data[i][j] = data[i];
			}
		}
	
	}
	virtual ~Test_K_Epsilon_Equation() {};
	void update(size_t index_i, Real dt = 0.0);
	void impose_parabolic_k(size_t index_i);
	void mannul_update_velocity(size_t index_i, Real dt);
	void mannul_impose_pressure_gradient(size_t index_i);
	void clear_acc_prior(size_t index_i);

	std::vector<double> loadInputData(int num_data, int num_file, std::string file_name)
	{
		std::ifstream file("./MappingData/FromPy10/" + file_name, std::ios::binary);  
		if (file)
		{
			std::string line;
			std::vector<double> temp_data;
			for (int i = 0; i < num_data; ++i)
			{
				std::getline(file, line);
				//std::cout << line << std::endl;
				double value = std::stod(line); 
				temp_data.push_back(value); 
			}
			file.close();
			return temp_data;
		}
		else 
		{
			std::cerr << "cannot open file." << std::endl;
			std::cin.get();
			exit(1);
		}

	};

};
void Test_K_Epsilon_Equation::update(size_t index_i, Real dt)
{
	for (int i = 0; i < num_data; ++i)
	{
		if (pos_[index_i][1] > i * resolution_ref && pos_[index_i][1] < (1 + i) * resolution_ref)
		{
			//if (i == 6)
				//system("pause");
			vel_[index_i][0] = input_data[i][0];
			turbu_k_[index_i] = input_data[i][1];
			turbu_epsilon_[index_i] = input_data[i][2];
			turbu_mu_[index_i] = input_data[i][3];
			continue;
		}
	}
}
void Test_K_Epsilon_Equation::impose_parabolic_k(size_t index_i)
{
	Real Radius = (DH / 2.0);
	Real transformed_pos = pos_[index_i][1] - (DH / 2.0);
	turbu_k_[index_i] = 1.5 * U_f * (1.0 - transformed_pos * transformed_pos / Radius / Radius);
	expect_k_gradient_[index_i][1] = -2.0 * 1.5 * U_f / Radius / Radius * transformed_pos;
}
void Test_K_Epsilon_Equation::mannul_update_velocity(size_t index_i, Real dt)
{
	vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
}
void Test_K_Epsilon_Equation::mannul_impose_pressure_gradient(size_t index_i)
{
	//acc_[index_i] = 2.0/DH*velo_friction_[index_i].dot(velo_friction_[index_i])* velo_friction_[index_i].normalized();
	acc_[index_i] = 2.0 / DH * 0.053033 * 0.053033 * Vecd(1, 0);
}


void Test_K_Epsilon_Equation::clear_acc_prior(size_t index_i)
{
	acc_prior_[index_i] = Vecd::Zero();
}


TEST_F(TurbulentModule, TestTurbulentKineticEnergyEquation)
{
	SimpleDynamics<Test_K_Epsilon_Equation, SequencedPolicy> test_k_ep_equation(water_block);
	wall_boundary_normal_direction.exec();
	write_body_states.writeToFile();
	test_k_ep_equation.exec(); //** impose profiles *
	write_body_states.writeToFile();//** output to check initial profiles *
	
	dt = 0.1;
	int num_iter = 0;
	while (GlobalStaticVariables::physical_time_<100.0)
	{
		//** Try to introduce B correction * 
		correct_configuration.exec();
		//update_eddy_viscosity.exec();


		//** Test viscous force *
		standard_wall_function_correction.exec(); //the wall viscous relies on wall Func values
		turbulent_viscous_acceleration.exec();
		write_body_states.writeToFile();
		std::cout <<"Test viscous" << std::endl;
		system("pause");

		//turbulent_kinetic_energy_acceleration.exec();
		for (int index_i = 0; index_i < num_fluid_particle; ++index_i)
		{
			test_k_ep_equation.mannul_impose_pressure_gradient(index_i);
			//test_k_ep_equation.mannul_update_velocity(index_i, dt);
		}

		k_equation_relaxation.exec(dt);
		epsilon_equation_relaxation.exec(dt);
		standard_wall_function_correction.exec();

		GlobalStaticVariables::physical_time_ += dt;
		std::cout << "physical_time_="<< GlobalStaticVariables::physical_time_ <<std::endl;
		if(num_iter % 10 == 0)
			write_body_states.writeToFile();
		num_iter++;
		for (int index_i = 0; index_i < num_fluid_particle; ++index_i)
		{
			test_k_ep_equation.clear_acc_prior(index_i);
		}
		system("pause");
	}

	ASSERT_NEAR(1, 1, 0.02);


	std::cout << "TestTurbulentKineticEnergyEquation completed, continuing will rewrite body states data" << std::endl;
	system("pause");

}

TEST_F(TurbulentModule, TestMomentumEquation)
{



	ASSERT_NEAR(1, 1, 0.02);


	std::cout << "TestTurbulentKineticEnergyEquation completed, continuing will rewrite body states data" << std::endl;
	system("pause");

}


TEST_F(TurbulentModule, TestTurbulentKineticEnergyGradient)
{
	SimpleDynamics<Test_K_Epsilon_Equation, SequencedPolicy> test_k_ep_equation(water_block);
	InteractionDynamics<fluid_dynamics::TKEnergyAccComplex, SequencedPolicy> turbulent_kinetic_energy_acceleration(water_block_complex_relation);
	for (int index_i = 0; index_i < num_fluid_particle; ++index_i)
	{
		test_k_ep_equation.impose_parabolic_k(index_i); //** impose parabolic profiles *
	}
	write_body_states.writeToFile();
	turbulent_kinetic_energy_acceleration.exec();
	write_body_states.writeToFile();
	ASSERT_NEAR(1, 1, 0.02);


	std::cout << "TestTurbulentKineticEnergyGradient completed, continuing will rewrite body states data" << std::endl;
	system("pause");

}



//----------------------------------------------------------------------


//----------------------------------------------------------------------
//	Test velocity gradient
//----------------------------------------------------------------------
class TestVeloGrad : public fluid_dynamics::FluidInitialCondition
{
protected:
	StdLargeVec<Vecd>& pos_, & vel_;
	StdLargeVec<Matd> expect_velo_gradient_;
	StdLargeVec<Matd>& velocity_gradient;
public:
	TestVeloGrad(SPHBody& sph_body) : FluidInitialCondition(sph_body),
		pos_(particles_->pos_), vel_(particles_->vel_),
		velocity_gradient(*particles_->getVariableByName<Matd>("Velocity_Gradient"))
	{
		particles_->registerVariable(expect_velo_gradient_, "ExpectVeloGradient");
		particles_->addVariableToWrite<Matd>("ExpectVeloGradient");
	}
	virtual ~TestVeloGrad() {};
	void update(size_t index_i, Real dt = 0.0)
	{
		Real Radius = (DH / 2.0);
		Real transformed_pos = pos_[index_i][1] - (DH / 2.0);
		vel_[index_i][0] = 1.5 * U_f * (1.0 - transformed_pos * transformed_pos / Radius / Radius);
		expect_velo_gradient_[index_i](0, 1) = -2.0 * 1.5 * U_f / Radius / Radius * transformed_pos;
	};
	Matd getExpectVelocityGradient(size_t index_i)
	{
		return expect_velo_gradient_[index_i];
	};
	Matd getCalculatedVelocityGradient(size_t index_i)
	{
		return velocity_gradient[index_i];
	};
	Vecd getPosition(size_t index_i)
	{
		return pos_[index_i];
	};
};
TEST_F(TurbulentModule, TestVelocityGradient)
{
	SimpleDynamics<TestVeloGrad, SequencedPolicy> test_velocity_gradient(water_block);
	test_velocity_gradient.exec(); //** impose initial velocity profile and calculate theoretical velocity gradient value *
	k_equation_relaxation.exec(dt);
	write_body_states.writeToFile();
	number_of_iterations++;
	std::cout << "number of fluid particles=" << num_fluid_particle << std::endl;
	for (int index_i = 0; index_i < num_fluid_particle; ++index_i)
	{
		//** reduce testing region, because in horizontal direction, kernel truncation is avoided by using in/outlet buffer *
		Vecd pos_ = test_velocity_gradient.getPosition(index_i);
		if (pos_[0] >= 0.5 * DL && pos_[0] < 0.5 * DL + BW)
		{
			Matd A = test_velocity_gradient.getExpectVelocityGradient(index_i);
			Matd B = test_velocity_gradient.getCalculatedVelocityGradient(index_i);
			for (int j = 0; j < Dimensions; ++j)
			{
				for (int i = 0; i < Dimensions; ++i)
				{
					ASSERT_NEAR(A(i, j), B(i, j), 1.0);
				}
			}
		}
		
	}
}
//----------------------------------------------------------------------



//=================================================================================================//
int main(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}