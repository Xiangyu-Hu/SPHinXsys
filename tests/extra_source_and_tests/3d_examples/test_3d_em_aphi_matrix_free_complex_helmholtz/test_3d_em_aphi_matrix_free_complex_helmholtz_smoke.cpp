#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_solver.h"

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::matrix_free;

int main()
{
    StdVec<Complex> field{Complex(1.0, -0.5), Complex(-0.25, 0.75)};
    StdVec<Complex> rhs(2, Complex(0.0, 0.0));
    StdVec<Complex> reaction_coefficient(2, Complex(1.0, 0.0));
    StdVec<ScalarLaplaceEdge> edges{ScalarLaplaceEdge(0, 1, 0.5)};

    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(field.size());

    ScalarComplexHelmholtzSolverParameters parameters;
    parameters.max_iterations_ = 12;
    parameters.relaxation_factor_ = 0.5;
    parameters.absolute_tolerance_ = 1.0e-8;

    const ScalarComplexHelmholtzSolverState state = solveScalarComplexHelmholtz(
        field, rhs, reaction_coefficient, residuals, parameters,
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            accumulateScalarLaplaceResidualsFromEdges(current_field, edges, current_residuals);
        });

    if (state.iterations_ == 0)
    {
        return 1;
    }

    if (state.current_residual_l2_ > state.initial_residual_l2_)
    {
        return 2;
    }

    return 0;
}
