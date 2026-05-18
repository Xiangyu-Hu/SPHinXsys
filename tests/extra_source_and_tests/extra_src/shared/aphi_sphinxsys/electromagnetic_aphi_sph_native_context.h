#ifndef ELECTROMAGNETIC_APHI_SPH_NATIVE_CONTEXT_H
#define ELECTROMAGNETIC_APHI_SPH_NATIVE_CONTEXT_H

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_assembly.h"
#include "inner_body_relation.h"
#include "relation_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

struct MatrixFreeAPhiSphNativeContext;

/** When true, coupled residual/solver paths use SPH-native CK operators (requires a live context). */
bool matrixFreeAPhiUseSphNativeOperators();

/** When true, Helmholtz/Jacobi inner solves use SPH-native Laplace residuals (falls back to RESIDUALS env). */
bool matrixFreeAPhiUseSphNativeHelmholtz();

bool useMatrixFreeAPhiSphNativeHelmholtz(MatrixFreeAPhiSphNativeContext *sph_native_context);

/**
 * Holds RealBody and one shared_ck Inner relation for lazily-built SPH-native operators.
 * Graph-based Jacobi / Helmholtz solves remain on the validated graph branch.
 */
struct MatrixFreeAPhiSphNativeContext
{
    MatrixFreeAPhiSphNativeContext(RealBody &body, Inner<> &ck_inner_relation);

    sph::SPHAPhiOperatorAssembly &operatorAssembly();
    const sph::SPHAPhiOperatorAssembly &operatorAssembly() const;
    const MatrixFreeAPhiSphNativeDiagnostics &diagnostics() const;

    RealBody &body_;
    Inner<> &ck_inner_relation_;
    std::unique_ptr<sph::SPHAPhiOperatorAssembly> assembly_;
};

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_sph_native_context.hpp"

#endif // ELECTROMAGNETIC_APHI_SPH_NATIVE_CONTEXT_H
