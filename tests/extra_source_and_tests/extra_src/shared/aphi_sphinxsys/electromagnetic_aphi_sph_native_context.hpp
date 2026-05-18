#ifndef ELECTROMAGNETIC_APHI_SPH_NATIVE_CONTEXT_HPP
#define ELECTROMAGNETIC_APHI_SPH_NATIVE_CONTEXT_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_sph_native_context.h"

#include <cstdlib>

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

inline bool matrixFreeAPhiUseSphNativeOperators()
{
    const char *value = std::getenv("EM_APHI_USE_SPH_NATIVE_RESIDUALS");
    return value != nullptr && std::atoi(value) != 0;
}

inline bool matrixFreeAPhiUseSphNativeHelmholtz()
{
    const char *value = std::getenv("EM_APHI_USE_SPH_NATIVE_HELMHOLTZ");
    if (value != nullptr)
    {
        return std::atoi(value) != 0;
    }
    return matrixFreeAPhiUseSphNativeOperators();
}

inline bool useMatrixFreeAPhiSphNativePath(MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    return sph_native_context != nullptr && matrixFreeAPhiUseSphNativeOperators();
}

inline bool useMatrixFreeAPhiSphNativeHelmholtz(MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    return sph_native_context != nullptr && matrixFreeAPhiUseSphNativeHelmholtz();
}

inline MatrixFreeAPhiSphNativeContext::MatrixFreeAPhiSphNativeContext(RealBody &body, Inner<> &ck_inner_relation)
    : body_(body), ck_inner_relation_(ck_inner_relation)
{
}

inline sph::SPHAPhiOperatorAssembly &MatrixFreeAPhiSphNativeContext::operatorAssembly()
{
    if (!assembly_)
    {
        assembly_ = std::make_unique<sph::SPHAPhiOperatorAssembly>(body_, ck_inner_relation_);
    }
    return *assembly_;
}

inline const sph::SPHAPhiOperatorAssembly &MatrixFreeAPhiSphNativeContext::operatorAssembly() const
{
    return const_cast<MatrixFreeAPhiSphNativeContext *>(this)->operatorAssembly();
}

inline const MatrixFreeAPhiSphNativeDiagnostics &MatrixFreeAPhiSphNativeContext::diagnostics() const
{
    return operatorAssembly().diagnostics();
}

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_SPH_NATIVE_CONTEXT_HPP
