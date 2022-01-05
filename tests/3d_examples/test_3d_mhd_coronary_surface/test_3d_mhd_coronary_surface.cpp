
#include "implicit_surface.h"
#include "complex_implicit_shape.h"
#include "sphinxsys.h"

#include <memory>

using namespace SPH;


class CoronaryBody : public SolidBody
{
public:
	CoronaryBody(SPHSystem& system, std::string body_name)
		: SolidBody(system, body_name)
	{
		ImplicitSurface  coronary;
		coronary.loadFromMHDImge("./input/coronary_out_resample.mhd");

		body_shape_.add<ComplexImplicitShape>(this, coronary);

	}
};


int main(int argc, char **argv)
{
    //----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
    BoundingBox boundingBox(std::make_pair(Vec3d{-1,-1,-1}, Vec3d(511,511,351)));
	SPHSystem system(boundingBox, 1);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = true;
	/** output environment. */
	In_Output 	in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	CoronaryBody inputbody(system, "coronary");
	SolidParticles inputbody_particles(inputbody);
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp inputbody_recording_to_vtu(in_output, { inputbody });
	MeshRecordingToPlt 	cell_linked_list_recording(in_output, inputbody, inputbody.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner inputbody_inner(inputbody);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	relax_dynamics::RelaxationStepInner relaxation_step_inner(inputbody_inner, true);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	relaxation_step_inner.surface_bounding_.parallel_exec();
	inputbody.updateCellLinkedList();
	//----------------------------------------------------------------------
	//	First output before the simulation.
	//----------------------------------------------------------------------
	inputbody_recording_to_vtu.writeToFile(0.0);
	cell_linked_list_recording.writeToFile(0.0);
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 100)
	{
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the airfoil N = " << ite_p << "\n";
			inputbody_recording_to_vtu.writeToFile(ite_p);
		}
	}
	std::cout << "The physics relaxation process of airfoil finish !" << std::endl;

	return 0;
}