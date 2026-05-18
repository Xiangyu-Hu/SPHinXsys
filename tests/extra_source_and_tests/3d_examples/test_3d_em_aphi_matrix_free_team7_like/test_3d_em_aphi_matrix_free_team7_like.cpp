/**
 * @file test_3d_em_aphi_matrix_free_team7_like.cpp
 * @brief Matrix-free TEAM7-like three-region baseline for the matrix-free A-phi to Joule-heating chain.
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_solver.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_native_context.h"

#include <cstdlib>
#include <memory>
#include <iomanip>
#include <string>
#include <Eigen/Dense>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::matrix_free;

namespace
{
using ThermalRelaxationInner =
    DiffusionRelaxationRK2<DiffusionRelaxation<Inner<KernelGradientInner>, IsotropicDiffusion>>;

class Team7LikeMatrixFreeBoxShape : public ComplexShape
{
  public:
    Team7LikeMatrixFreeBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};
class TemperatureInitialization : public LocalDynamics
{
  public:
    explicit TemperatureInitialization(SPHBody &sph_body, Real initial_temperature)
        : LocalDynamics(sph_body),
          initial_temperature_(initial_temperature),
          temperature_(particles_->registerStateVariableData<Real>("Temperature"))
    {
        particles_->addVariableToWrite<Real>("Temperature");
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        temperature_[index_i] = initial_temperature_;
    }

  protected:
    Real initial_temperature_;
    Real *temperature_;
};

class JouleHeatingVariableInitialization : public LocalDynamics
{
  public:
    explicit JouleHeatingVariableInitialization(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          joule_heat_source_(particles_->registerStateVariableData<Real>("JouleHeatSource")),
          temperature_change_rate_by_joule_(particles_->registerStateVariableData<Real>("TemperatureChangeRateByJoule"))
    {
        particles_->addVariableToWrite<Real>("JouleHeatSource");
        particles_->addVariableToWrite<Real>("TemperatureChangeRateByJoule");
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        joule_heat_source_[index_i] = 0.0;
        temperature_change_rate_by_joule_[index_i] = 0.0;
    }

  protected:
    Real *joule_heat_source_;
    Real *temperature_change_rate_by_joule_;
};

class AddJouleHeatToTemperature : public LocalDynamics
{
  public:
    explicit AddJouleHeatToTemperature(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          temperature_(particles_->getVariableDataByName<Real>("Temperature")),
          temperature_change_rate_by_joule_(particles_->getVariableDataByName<Real>("TemperatureChangeRateByJoule"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        temperature_[index_i] += dt * temperature_change_rate_by_joule_[index_i];
    }

  protected:
    Real *temperature_;
    Real *temperature_change_rate_by_joule_;
};

bool inside_box_region(const Vecd &position,
                       Real xmin, Real xmax,
                       Real ymin, Real ymax,
                       Real zmin, Real zmax)
{
    return position[0] >= xmin && position[0] <= xmax &&
           position[1] >= ymin && position[1] <= ymax &&
           position[2] >= zmin && position[2] <= zmax;
}

Real coordinate_from_fraction(Real fraction, Real extent)
{
    return fraction * extent;
}

Real smooth_box_profile(const Vecd &position,
                        Real xmin, Real xmax,
                        Real ymin, Real ymax,
                        Real zmin, Real zmax,
                        Real amplitude)
{
    if (!inside_box_region(position, xmin, xmax, ymin, ymax, zmin, zmax))
    {
        return 0.0;
    }

    const Real xi = (position[0] - xmin) / (xmax - xmin + TinyReal);
    const Real eta = (position[1] - ymin) / (ymax - ymin + TinyReal);
    const Real zeta = (position[2] - zmin) / (zmax - zmin + TinyReal);
    return amplitude * std::sin(Pi * xi) * std::sin(Pi * eta) * std::sin(Pi * zeta);
}

struct PhysicalRegionSummary
{
    size_t particles = 0;
    Real mean_abs_a = 0.0;
    Real max_abs_a = 0.0;
    Real mean_abs_phi = 0.0;
    Real max_abs_phi = 0.0;
    Real mean_abs_e = 0.0;
    Real max_abs_e = 0.0;
    Real mean_abs_j = 0.0;
    Real max_abs_j = 0.0;
    Real mean_joule = 0.0;
    Real max_joule = 0.0;
    Real mean_temperature_delta = 0.0;
    Real max_temperature_delta = 0.0;
};

void finalize_region_summary(PhysicalRegionSummary &summary)
{
    if (summary.particles == 0)
    {
        return;
    }
    const Real inv = 1.0 / static_cast<Real>(summary.particles);
    summary.mean_abs_a *= inv;
    summary.mean_abs_phi *= inv;
    summary.mean_abs_e *= inv;
    summary.mean_abs_j *= inv;
    summary.mean_joule *= inv;
    summary.mean_temperature_delta *= inv;
}

struct Team7LikeVisualizationParticleFields
{
    Real *region_id = nullptr;
    Real *is_conductor = nullptr;
    Real *is_coil = nullptr;
    Real *is_source = nullptr;
    Real *phi_real = nullptr;
    Real *phi_imag = nullptr;
    Real *abs_phi = nullptr;
    Real *ax_real = nullptr;
    Real *ax_imag = nullptr;
    Real *ay_real = nullptr;
    Real *ay_imag = nullptr;
    Real *az_real = nullptr;
    Real *az_imag = nullptr;
    Real *abs_a = nullptr;
    Real *ex_real = nullptr;
    Real *ex_imag = nullptr;
    Real *ey_real = nullptr;
    Real *ey_imag = nullptr;
    Real *ez_real = nullptr;
    Real *ez_imag = nullptr;
    Real *abs_e = nullptr;
    Real *jx_real = nullptr;
    Real *jx_imag = nullptr;
    Real *jy_real = nullptr;
    Real *jy_imag = nullptr;
    Real *jz_real = nullptr;
    Real *jz_imag = nullptr;
    Real *abs_j = nullptr;
    Real *joule_density = nullptr;
    Real *temperature_delta = nullptr;
};

Team7LikeVisualizationParticleFields registerTeam7LikeVisualizationVariables(BaseParticles &particles)
{
    Team7LikeVisualizationParticleFields fields;
    auto register_real = [&](const std::string &name) -> Real *
    {
        Real *data = particles.registerStateVariableData<Real>(name, Real(0));
        particles.addVariableToWrite<Real>(name);
        return data;
    };
    fields.region_id = register_real("RegionId");
    fields.is_conductor = register_real("IsConductor");
    fields.is_coil = register_real("IsCoil");
    fields.is_source = register_real("IsSource");
    fields.phi_real = register_real("PhiReal");
    fields.phi_imag = register_real("PhiImag");
    fields.abs_phi = register_real("AbsPhi");
    fields.ax_real = register_real("AxReal");
    fields.ax_imag = register_real("AxImag");
    fields.ay_real = register_real("AyReal");
    fields.ay_imag = register_real("AyImag");
    fields.az_real = register_real("AzReal");
    fields.az_imag = register_real("AzImag");
    fields.abs_a = register_real("AbsA");
    fields.ex_real = register_real("ExReal");
    fields.ex_imag = register_real("ExImag");
    fields.ey_real = register_real("EyReal");
    fields.ey_imag = register_real("EyImag");
    fields.ez_real = register_real("EzReal");
    fields.ez_imag = register_real("EzImag");
    fields.abs_e = register_real("AbsE");
    fields.jx_real = register_real("JxReal");
    fields.jx_imag = register_real("JxImag");
    fields.jy_real = register_real("JyReal");
    fields.jy_imag = register_real("JyImag");
    fields.jz_real = register_real("JzReal");
    fields.jz_imag = register_real("JzImag");
    fields.abs_j = register_real("AbsJ");
    fields.joule_density = register_real("JouleDensity");
    fields.temperature_delta = register_real("TemperatureDelta");
    return fields;
}

void registerTeam7LikeVisualizationWithBodyStatesRecording(SolidBody &body,
                                                             BodyStatesRecordingToVtp &write_states,
                                                             const Team7LikeVisualizationParticleFields &fields)
{
    write_states.addToWrite<Real>(body, "RegionId");
    write_states.addToWrite<Real>(body, "IsConductor");
    write_states.addToWrite<Real>(body, "IsCoil");
    write_states.addToWrite<Real>(body, "IsSource");
    write_states.addToWrite<Real>(body, "PhiReal");
    write_states.addToWrite<Real>(body, "PhiImag");
    write_states.addToWrite<Real>(body, "AbsPhi");
    write_states.addToWrite<Real>(body, "AxReal");
    write_states.addToWrite<Real>(body, "AxImag");
    write_states.addToWrite<Real>(body, "AyReal");
    write_states.addToWrite<Real>(body, "AyImag");
    write_states.addToWrite<Real>(body, "AzReal");
    write_states.addToWrite<Real>(body, "AzImag");
    write_states.addToWrite<Real>(body, "AbsA");
    write_states.addToWrite<Real>(body, "ExReal");
    write_states.addToWrite<Real>(body, "ExImag");
    write_states.addToWrite<Real>(body, "EyReal");
    write_states.addToWrite<Real>(body, "EyImag");
    write_states.addToWrite<Real>(body, "EzReal");
    write_states.addToWrite<Real>(body, "EzImag");
    write_states.addToWrite<Real>(body, "AbsE");
    write_states.addToWrite<Real>(body, "JxReal");
    write_states.addToWrite<Real>(body, "JxImag");
    write_states.addToWrite<Real>(body, "JyReal");
    write_states.addToWrite<Real>(body, "JyImag");
    write_states.addToWrite<Real>(body, "JzReal");
    write_states.addToWrite<Real>(body, "JzImag");
    write_states.addToWrite<Real>(body, "AbsJ");
    write_states.addToWrite<Real>(body, "JouleDensity");
    write_states.addToWrite<Real>(body, "TemperatureDelta");
    write_states.addToWrite<Real>(body, "Temperature");
    write_states.addToWrite<Real>(body, "JouleHeatSource");
    (void)fields;
}

void populateTeam7LikeVisualizationVariables(const Team7LikeVisualizationParticleFields &viz,
                                             const MatrixFreeAPhiFields &fields,
                                             const StdVec<Vec3c> &electric_field,
                                             const StdVec<Vec3c> &current_density,
                                             const StdVec<Real> &joule_density,
                                             const StdVec<Real> &temperature_delta,
                                             const StdVec<bool> &particle_is_conductor,
                                             const StdVec<bool> &particle_is_coil,
                                             const StdVec<bool> &particle_is_source)
{
    const size_t n = fields.ax.size();
    for (size_t i = 0; i != n; ++i)
    {
        Real region_id = Real(0);
        if (particle_is_conductor[i])
        {
            region_id = Real(1);
        }
        else if (particle_is_coil[i])
        {
            region_id = Real(2);
        }
        viz.region_id[i] = region_id;
        viz.is_conductor[i] = particle_is_conductor[i] ? Real(1) : Real(0);
        viz.is_coil[i] = particle_is_coil[i] ? Real(1) : Real(0);
        viz.is_source[i] = particle_is_source[i] ? Real(1) : Real(0);
        viz.phi_real[i] = fields.phi[i].real();
        viz.phi_imag[i] = fields.phi[i].imag();
        viz.abs_phi[i] = std::abs(fields.phi[i]);
        viz.ax_real[i] = fields.ax[i].real();
        viz.ax_imag[i] = fields.ax[i].imag();
        viz.ay_real[i] = fields.ay[i].real();
        viz.ay_imag[i] = fields.ay[i].imag();
        viz.az_real[i] = fields.az[i].real();
        viz.az_imag[i] = fields.az[i].imag();
        viz.abs_a[i] = std::sqrt(std::norm(fields.ax[i]) + std::norm(fields.ay[i]) + std::norm(fields.az[i]));
        viz.ex_real[i] = electric_field[i][0].real();
        viz.ex_imag[i] = electric_field[i][0].imag();
        viz.ey_real[i] = electric_field[i][1].real();
        viz.ey_imag[i] = electric_field[i][1].imag();
        viz.ez_real[i] = electric_field[i][2].real();
        viz.ez_imag[i] = electric_field[i][2].imag();
        viz.abs_e[i] = static_cast<Real>(electric_field[i].norm());
        viz.jx_real[i] = current_density[i][0].real();
        viz.jx_imag[i] = current_density[i][0].imag();
        viz.jy_real[i] = current_density[i][1].real();
        viz.jy_imag[i] = current_density[i][1].imag();
        viz.jz_real[i] = current_density[i][2].real();
        viz.jz_imag[i] = current_density[i][2].imag();
        viz.abs_j[i] = static_cast<Real>(current_density[i].norm());
        viz.joule_density[i] = joule_density[i];
        viz.temperature_delta[i] = temperature_delta[i];
    }
}

Real computeRealFieldAverage(const StdVec<Real> &field)
{
    if (field.empty())
    {
        return 0.0;
    }
    Real sum = 0.0;
    for (Real value : field)
    {
        sum += value;
    }
    return sum / static_cast<Real>(field.size());
}

Real computeL2Error(const StdVec<Complex> &solution, const StdVec<Complex> &reference)
{
    if (solution.empty())
    {
        return 0.0;
    }

    Real squared_sum = 0.0;
    for (size_t i = 0; i != solution.size(); ++i)
    {
        squared_sum += std::norm(solution[i] - reference[i]);
    }
    return std::sqrt(squared_sum / static_cast<Real>(solution.size()));
}

Real computeComplexFieldL2NormLocal(const StdVec<Complex> &field)
{
    if (field.empty())
    {
        return 0.0;
    }

    Real squared_sum = 0.0;
    for (const Complex &value : field)
    {
        squared_sum += std::norm(value);
    }
    return std::sqrt(squared_sum / static_cast<Real>(field.size()));
}

StdVec<Complex> computeDivergenceOfVectorFieldLocal(const MatrixFreePairwiseGraph &graph,
                                                    const MatrixFreeAPhiFields &fields)
{
    const size_t number_of_particles = fields.ax.size();
    const StdVec<Vec3c> grad_ax = applyMatrixFreeGradient(graph, fields.ax);
    const StdVec<Vec3c> grad_ay = applyMatrixFreeGradient(graph, fields.ay);
    const StdVec<Vec3c> grad_az = applyMatrixFreeGradient(graph, fields.az);
    StdVec<Complex> divergence(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        divergence[i] = grad_ax[i][0] + grad_ay[i][1] + grad_az[i][2];
    }
    return divergence;
}

Real computeMaxError(const StdVec<Complex> &solution, const StdVec<Complex> &reference)
{
    Real max_error = 0.0;
    for (size_t i = 0; i != solution.size(); ++i)
    {
        max_error = SMAX(max_error, std::abs(solution[i] - reference[i]));
    }
    return max_error;
}


Real computeVectorFieldL2Error(const StdVec<Vec3c> &solution, const StdVec<Vec3c> &reference)
{
    if (solution.size() != reference.size() || solution.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (size_t i = 0; i != solution.size(); ++i)
    {
        for (int axis = 0; axis != Dimensions; ++axis)
        {
            squared_sum += std::norm(solution[i][axis] - reference[i][axis]);
        }
    }
    return std::sqrt(squared_sum / static_cast<Real>(solution.size()));
}

Real computeRealFieldL2Error(const StdVec<Real> &solution, const StdVec<Real> &reference)
{
    if (solution.size() != reference.size() || solution.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (size_t i = 0; i != solution.size(); ++i)
    {
        const Real difference = solution[i] - reference[i];
        squared_sum += difference * difference;
    }
    return std::sqrt(squared_sum / static_cast<Real>(solution.size()));
}

StdVec<Real> computeJouleHeatingDensity(const StdVec<Vec3c> &electric_field,
                                        const StdVec<Vec3c> &current_density)
{
    const size_t number_of_particles = electric_field.size();
    StdVec<Real> joule_density(number_of_particles, 0.0);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        Complex power_density(0.0, 0.0);
        for (int axis = 0; axis != Dimensions; ++axis)
        {
            power_density += electric_field[i][axis] * std::conj(current_density[i][axis]);
        }
        joule_density[i] = 0.5 * power_density.real();
    }
    return joule_density;
}


Real computeRealFieldL2Norm(const StdVec<Real> &field)
{
    if (field.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (const Real value : field)
    {
        squared_sum += value * value;
    }
    return std::sqrt(squared_sum / static_cast<Real>(field.size()));
}

Real computeWeightedCentroidX(const StdVec<Real> &weights, const Vecd *positions)
{
    if (weights.empty())
    {
        return 0.0;
    }
    Real weighted_sum = 0.0;
    Real total_weight = 0.0;
    for (size_t i = 0; i != weights.size(); ++i)
    {
        const Real w = SMAX(weights[i], Real(0.0));
        weighted_sum += w * positions[i][0];
        total_weight += w;
    }
    return total_weight <= TinyReal ? 0.0 : weighted_sum / total_weight;
}

int findMirrorIndex(const Vecd *positions, size_t number_of_particles, size_t target_index,
                    Real center_y, Real tolerance)
{
    const Real mirror_y = 2.0 * center_y - positions[target_index][1];
    for (size_t candidate = 0; candidate != number_of_particles; ++candidate)
    {
        if (std::abs(positions[candidate][0] - positions[target_index][0]) <= tolerance &&
            std::abs(positions[candidate][1] - mirror_y) <= tolerance &&
            std::abs(positions[candidate][2] - positions[target_index][2]) <= tolerance)
        {
            return static_cast<int>(candidate);
        }
    }
    return -1;
}

Real computeComplexFieldMirrorSymmetryError(const StdVec<Complex> &field, const Vecd *positions,
                                            size_t number_of_particles, Real center_y, Real tolerance)
{
    Real squared_sum = 0.0;
    size_t pair_count = 0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        if (positions[i][1] > center_y + tolerance)
        {
            continue;
        }
        const int mirror_index = findMirrorIndex(positions, number_of_particles, i, center_y, tolerance);
        if (mirror_index < 0)
        {
            continue;
        }
        squared_sum += std::norm(field[i] - field[static_cast<size_t>(mirror_index)]);
        ++pair_count;
    }
    return pair_count == 0 ? 0.0 : std::sqrt(squared_sum / static_cast<Real>(pair_count));
}

Real computeComplexFieldMirrorAntisymmetryError(const StdVec<Complex> &field, const Vecd *positions,
                                                size_t number_of_particles, Real center_y, Real tolerance)
{
    Real squared_sum = 0.0;
    size_t pair_count = 0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        if (positions[i][1] > center_y + tolerance)
        {
            continue;
        }
        const int mirror_index = findMirrorIndex(positions, number_of_particles, i, center_y, tolerance);
        if (mirror_index < 0)
        {
            continue;
        }
        squared_sum += std::norm(field[i] + field[static_cast<size_t>(mirror_index)]);
        ++pair_count;
    }
    return pair_count == 0 ? 0.0 : std::sqrt(squared_sum / static_cast<Real>(pair_count));
}

StdVec<Complex> makeMeanCenteredCopy(const StdVec<Complex> &field)
{
    StdVec<Complex> centered = field;
    removeMeanOffset(centered);
    return centered;
}

Real computeRealFieldMirrorSymmetryError(const StdVec<Real> &field, const Vecd *positions,
                                         size_t number_of_particles, Real center_y, Real tolerance)
{
    Real squared_sum = 0.0;
    size_t pair_count = 0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        if (positions[i][1] > center_y + tolerance)
        {
            continue;
        }
        const int mirror_index = findMirrorIndex(positions, number_of_particles, i, center_y, tolerance);
        if (mirror_index < 0)
        {
            continue;
        }
        const Real difference = field[i] - field[static_cast<size_t>(mirror_index)];
        squared_sum += difference * difference;
        ++pair_count;
    }
    return pair_count == 0 ? 0.0 : std::sqrt(squared_sum / static_cast<Real>(pair_count));
}

size_t readEnvSizeT(const char *name, size_t default_value)
{
    const char *value = std::getenv(name);
    return value == nullptr ? default_value : static_cast<size_t>(std::strtoull(value, nullptr, 10));
}

Real readEnvReal(const char *name, Real default_value)
{
    const char *value = std::getenv(name);
    return value == nullptr ? default_value : static_cast<Real>(std::strtod(value, nullptr));
}

bool readEnvBool(const char *name, bool default_value)
{
    const char *value = std::getenv(name);
    if (value == nullptr)
    {
        return default_value;
    }
    return std::atoi(value) != 0;
}

int readEnvInt(const char *name, int default_value)
{
    const char *value = std::getenv(name);
    if (value == nullptr || value[0] == '\0')
    {
        return default_value;
    }
    return std::atoi(value);
}

std::string readEnvString(const char *name, const std::string &default_value)
{
    const char *value = std::getenv(name);
    return value == nullptr ? default_value : std::string(value);
}

void enforceReference(StdVec<Complex> &field, size_t reference_index, Complex reference_value)
{
    if (field.empty() || reference_index >= field.size())
    {
        return;
    }
    const Complex offset = field[reference_index] - reference_value;
    for (Complex &value : field)
    {
        value -= offset;
    }
}

ScalarComplexHelmholtzSolverState solveDiscretePhiForSourceFreeCase(const MatrixFreePairwiseGraph &graph,
                                                                    const StdVec<Real> &sigma,
                                                                    const StdVec<Complex> &rhs,
                                                                    StdVec<Complex> &phi)
{
    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(phi.size());
    ScalarComplexHelmholtzSolverParameters parameters;
    parameters.max_iterations_ = 8000;
    parameters.relaxation_factor_ = 1.0;
    parameters.absolute_tolerance_ = 1.0e-6;
    parameters.diagonal_regularization_ = 1.0e-12;

    const StdVec<Complex> zero_reaction(phi.size(), Complex(0.0, 0.0));
    return solveScalarComplexHelmholtz(
        phi, rhs, zero_reaction, residuals, parameters,
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            accumulateScalarLaplaceResidualsFromClearedGraph(graph, current_field, sigma, current_residuals);
        });
}

Real computeResidualL2FromDifference(const StdVec<Complex> &lhs, const StdVec<Complex> &rhs)
{
    if (lhs.empty())
    {
        return 0.0;
    }

    Real squared_sum = 0.0;
    for (size_t i = 0; i != lhs.size(); ++i)
    {
        squared_sum += std::norm(lhs[i] - rhs[i]);
    }
    return std::sqrt(squared_sum / static_cast<Real>(lhs.size()));
}

StdVec<Complex> solveReferenceGaugeChi(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &rhs,
                                       Real &reference_residual_l2)
{
    using DenseComplexMatrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
    using DenseComplexVector = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;

    const size_t number_of_particles = rhs.size();
    DenseComplexMatrix matrix(number_of_particles, number_of_particles);
    DenseComplexVector rhs_vector(number_of_particles);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        rhs_vector(static_cast<Eigen::Index>(i)) = rhs[i];
    }

    for (size_t j = 0; j != number_of_particles; ++j)
    {
        StdVec<Complex> basis(number_of_particles, Complex(0.0, 0.0));
        basis[j] = Complex(1.0, 0.0);
        const StdVec<Complex> column = applyScalarDivergenceOfGradientFromGraph(graph, basis);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            matrix(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = column[i];
        }
    }

    DenseComplexVector solution = matrix.completeOrthogonalDecomposition().solve(rhs_vector);
    StdVec<Complex> chi(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        chi[i] = solution(static_cast<Eigen::Index>(i));
    }
    removeMeanOffset(chi);

    const StdVec<Complex> operator_response = applyScalarDivergenceOfGradientFromGraph(graph, chi);
    reference_residual_l2 = computeResidualL2FromDifference(operator_response, rhs);
    return chi;
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = readEnvReal("EM_APHI_STAGGERED_SMOKE_DP", 0.05);
    const bool write_team7_vtp = readEnvBool("EM_APHI_TEAM7LIKE_WRITE_VTP", true);
    const std::string case_mode = readEnvString("EM_APHI_STAGGERED_SMOKE_CASE", "coulomb_variable_sigma_forced_response");
    const bool use_sycl_jacobi = readEnvBool("EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI", false);
    const bool use_sycl_gradient = readEnvBool("EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT", false);
    const bool use_sycl_harmonic_gradient =
        readEnvBool("EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT", use_sycl_gradient);
    const bool use_sycl_laplace_residual = readEnvBool("EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL", false);
    setMatrixFreeJacobiUseSycl(use_sycl_jacobi);
    setMatrixFreeGradientUseSycl(use_sycl_gradient);
    setMatrixFreeHarmonicGradientUseSycl(use_sycl_harmonic_gradient);
    setMatrixFreeLaplaceResidualUseSycl(use_sycl_laplace_residual);
#if SPHINXSYS_USE_SYCL
    {
        const int la_backend = readEnvInt("EM_APHI_MATRIX_FREE_LA_BACKEND", 0);
        setMatrixFreeAPhiLinearAlgebraBackend(la_backend == 1 ? MatrixFreeAPhiLinearAlgebraBackendKind::ReservedDeviceKrylov
                                                             : MatrixFreeAPhiLinearAlgebraBackendKind::HostControlledWithSyclPrimitives);
    }
    if (readEnvBool("EM_APHI_MATRIX_FREE_SYCL_POLICY_SELFTEST", false))
    {
        const MatrixFreeSyclPolicySnapshot snap = captureMatrixFreeSyclPolicySnapshot();
        setMatrixFreeGradientUseSycl(!snap.gradient);
        setMatrixFreeHarmonicGradientUseSycl(!snap.harmonic_gradient);
        setMatrixFreeLaplaceResidualUseSycl(!snap.laplace_residual);
        setMatrixFreeJacobiUseSycl(!snap.jacobi);
        applyMatrixFreeSyclPolicySnapshot(snap);
        if (matrixFreeGradientUseSycl() != snap.gradient || matrixFreeHarmonicGradientUseSycl() != snap.harmonic_gradient ||
            matrixFreeLaplaceResidualUseSycl() != snap.laplace_residual || matrixFreeJacobiUseSycl() != snap.jacobi)
        {
            std::cout << "matrix_free_sycl_policy_selftest_fail=1\n";
            return 1;
        }
        std::cout << "matrix_free_sycl_policy_selftest_pass=1\n";
    }
#endif
    const Real body_length = readEnvReal("EM_APHI_TEAM7LIKE_LENGTH", 1.20);
    const Real body_height = readEnvReal("EM_APHI_TEAM7LIKE_HEIGHT", 1.00);
    const Real body_width = readEnvReal("EM_APHI_TEAM7LIKE_WIDTH", 0.30);
    const Real boundary_width = readEnvReal("EM_APHI_TEAM7LIKE_BOUNDARY_WIDTH", 3.0 * dp_0);
    const Vecd body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    SolidBody body(sph_system,
                   makeShared<Team7LikeMatrixFreeBoxShape>("MatrixFreeTeam7LikeBody", body_center, body_halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    const Real thermal_diffusivity = readEnvReal("EM_APHI_JOULE_THERMAL_DIFFUSIVITY", 0.02);
    const Real rho_cp_thermal = readEnvReal("EM_APHI_JOULE_RHO_CP", 1.0);
    const Real initial_temperature = readEnvReal("EM_APHI_JOULE_INITIAL_TEMPERATURE", 300.0);
    const Real thermal_end_time = readEnvReal("EM_APHI_JOULE_END_TIME", 0.1);
    body.defineClosure<Solid, IsotropicDiffusion>(Solid(), ConstructArgs(std::string("Temperature"), thermal_diffusivity));
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(body);
    Inner<> inner_ck(body);
    BodyStatesRecordingToVtp write_states(sph_system);
    SimpleDynamics<TemperatureInitialization> temperature_initialization(body, initial_temperature);
    SimpleDynamics<JouleHeatingVariableInitialization> joule_heating_variable_initialization(body);
    GetDiffusionTimeStepSize get_thermal_time_step(body);
    ThermalRelaxationInner thermal_relaxation(inner_relation);
    SimpleDynamics<AddJouleHeatToTemperature> add_joule_heat_to_temperature(body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = body.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    const Team7LikeVisualizationParticleFields visualization_fields = registerTeam7LikeVisualizationVariables(particles);
    registerTeam7LikeVisualizationWithBodyStatesRecording(body, write_states, visualization_fields);
    temperature_initialization.exec();
    joule_heating_variable_initialization.exec();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real *volumetric_measure = particles.getVariableDataByName<Real>("VolumetricMeasure");

    MatrixFreeAPhiParameters parameters;
    parameters.angular_frequency = readEnvReal("EM_APHI_STAGGERED_SMOKE_OMEGA", 5.0);
    parameters.interface_contrast_threshold = readEnvReal("EM_APHI_INTERFACE_SIGMA_RATIO_THRESHOLD", 10.0);
    parameters.phi_reference_index = 0;
    parameters.phi_reference_value = Complex(0.0, 0.0);

    MatrixFreeAPhiDiscreteView discrete_view{
        number_of_particles,
        body.getSPHAdaptation().ReferenceSmoothingLength(),
        positions,
        volumetric_measure,
        nullptr,
        &inner_relation.inner_configuration_,
        nullptr};

    const MatrixFreePairwiseGraph graph = buildMatrixFreePairwiseGraph(discrete_view, parameters);
#if SPHINXSYS_USE_SYCL
    matrixFreeAPhiSyclBindTimestepResources(graph, number_of_particles);
#endif

    const bool use_sph_native_residuals = readEnvBool("EM_APHI_USE_SPH_NATIVE_RESIDUALS", false);
    std::unique_ptr<MatrixFreeAPhiSphNativeContext> sph_native_context_storage;
    if (use_sph_native_residuals)
    {
        sph_native_context_storage = std::make_unique<MatrixFreeAPhiSphNativeContext>(body, inner_ck);
    }
    MatrixFreeAPhiSphNativeContext *sph_native_context = sph_native_context_storage.get();

    const Real sigma_air = readEnvReal("EM_APHI_TEAM7LIKE_SIGMA_AIR", 1.0e-4);
    const Real sigma_conductor = readEnvReal("EM_APHI_TEAM7LIKE_SIGMA_CONDUCTOR", 1.0);
    const Real sigma_coil = readEnvReal("EM_APHI_TEAM7LIKE_SIGMA_COIL", 1.0e-4);
    const Real nu_air = readEnvReal("EM_APHI_TEAM7LIKE_NU_AIR", 1.0);
    const Real nu_conductor = readEnvReal("EM_APHI_TEAM7LIKE_NU_CONDUCTOR", 1.0);
    const Real nu_coil = readEnvReal("EM_APHI_TEAM7LIKE_NU_COIL", 1.0);
    const Real coil_source_amplitude = readEnvReal("EM_APHI_TEAM7LIKE_COIL_J0", 0.8);

    const Real conductor_xmin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_XMIN_FRACTION", 0.52);
    const Real conductor_xmax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_XMAX_FRACTION", 0.68);
    const Real conductor_ymin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_YMIN_FRACTION", 0.30);
    const Real conductor_ymax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_YMAX_FRACTION", 0.70);
    const Real conductor_zmin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_ZMIN_FRACTION", 0.20);
    const Real conductor_zmax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_ZMAX_FRACTION", 0.80);

    const Real coil_xmin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_XMIN_FRACTION", 0.24);
    const Real coil_xmax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_XMAX_FRACTION", 0.38);
    const Real coil_ymin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_YMIN_FRACTION", 0.15);
    const Real coil_ymax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_YMAX_FRACTION", 0.85);
    const Real coil_zmin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_ZMIN_FRACTION", 0.20);
    const Real coil_zmax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_ZMAX_FRACTION", 0.80);

    const Real conductor_xmin = coordinate_from_fraction(conductor_xmin_fraction, body_length);
    const Real conductor_xmax = coordinate_from_fraction(conductor_xmax_fraction, body_length);
    const Real conductor_ymin = coordinate_from_fraction(conductor_ymin_fraction, body_height);
    const Real conductor_ymax = coordinate_from_fraction(conductor_ymax_fraction, body_height);
    const Real conductor_zmin = coordinate_from_fraction(conductor_zmin_fraction, body_width);
    const Real conductor_zmax = coordinate_from_fraction(conductor_zmax_fraction, body_width);

    const Real coil_xmin = coordinate_from_fraction(coil_xmin_fraction, body_length);
    const Real coil_xmax = coordinate_from_fraction(coil_xmax_fraction, body_length);
    const Real coil_ymin = coordinate_from_fraction(coil_ymin_fraction, body_height);
    const Real coil_ymax = coordinate_from_fraction(coil_ymax_fraction, body_height);
    const Real coil_zmin = coordinate_from_fraction(coil_zmin_fraction, body_width);
    const Real coil_zmax = coordinate_from_fraction(coil_zmax_fraction, body_width);

    StdVec<Real> sigma(number_of_particles, sigma_air);
    StdVec<Real> nu(number_of_particles, nu_air);
    StdVec<bool> is_conductor(number_of_particles, false);
    StdVec<bool> is_coil(number_of_particles, false);
    StdVec<bool> is_source(number_of_particles, false);
    size_t conductor_particles = 0;
    size_t coil_particles = 0;
    size_t air_particles = 0;
    size_t source_particles = 0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const bool conductor_region = inside_box_region(
            positions[i], conductor_xmin, conductor_xmax, conductor_ymin, conductor_ymax, conductor_zmin, conductor_zmax);
        const bool coil_region = inside_box_region(
            positions[i], coil_xmin, coil_xmax, coil_ymin, coil_ymax, coil_zmin, coil_zmax);
        is_conductor[i] = conductor_region;
        is_coil[i] = coil_region;
        if (conductor_region)
        {
            sigma[i] = sigma_conductor;
            nu[i] = nu_conductor;
            ++conductor_particles;
        }
        else if (coil_region)
        {
            sigma[i] = sigma_coil;
            nu[i] = nu_coil;
            ++coil_particles;
        }
        else
        {
            ++air_particles;
        }
    }

    MatrixFreeAPhiFields exact_fields;
    exact_fields.ax.resize(number_of_particles, Complex(0.0, 0.0));
    exact_fields.ay.resize(number_of_particles, Complex(0.0, 0.0));
    exact_fields.az.resize(number_of_particles, Complex(0.0, 0.0));
    exact_fields.phi.resize(number_of_particles, Complex(0.0, 0.0));

    const Complex transverse_amplitude(1.0, 0.1);
    const Complex coupled_ax_amplitude(0.25, -0.08);
    const Complex driven_phi_amplitude(0.12, -0.03);
    const bool has_discrete_manufactured_reference = case_mode != "coulomb_variable_sigma_forced_response";
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const Real x = positions[i][0];
        const Real phase_x = Pi * x / body_length;
        if (case_mode == "coupled_source_free")
        {
            exact_fields.ax[i] = coupled_ax_amplitude * std::sin(phase_x);
        }
        else if (case_mode != "coulomb_variable_sigma_forced_response")
        {
            exact_fields.ay[i] = transverse_amplitude * std::sin(phase_x);
            if (case_mode == "coulomb_variable_sigma_driven")
            {
                const Real y = positions[i][1];
                const Real phase_y = Pi * y / body_height;
                exact_fields.phi[i] = driven_phi_amplitude * std::cos(phase_x) * std::sin(phase_y);
            }
        }
    }
    enforceReference(exact_fields.phi, parameters.phi_reference_index, parameters.phi_reference_value);

    Real reference_phi_build_residual = 0.0;
    if (case_mode == "coupled_source_free" || case_mode == "coulomb_variable_sigma_source_free")
    {
        // StdVec<Complex> sigma_ax(number_of_particles, Complex(0.0, 0.0));
        // StdVec<Complex> sigma_ay(number_of_particles, Complex(0.0, 0.0));
        // StdVec<Complex> sigma_az(number_of_particles, Complex(0.0, 0.0));
        // for (size_t i = 0; i != number_of_particles; ++i)
        // {
        //     sigma_ax[i] = sigma[i] * exact_fields.ax[i];
        //     sigma_ay[i] = sigma[i] * exact_fields.ay[i];
        //     sigma_az[i] = sigma[i] * exact_fields.az[i];
        // }
        // const StdVec<Vec3c> grad_sigma_ax = applyMatrixFreeGradient(graph, sigma_ax);
        // const StdVec<Vec3c> grad_sigma_ay = applyMatrixFreeGradient(graph, sigma_ay);
        // const StdVec<Vec3c> grad_sigma_az = applyMatrixFreeGradient(graph, sigma_az);

        const StdVec<Vec3c> sigma_grad_ax =
            applyMatrixFreeHarmonicWeightedGradient(graph, exact_fields.ax, sigma);
        const StdVec<Vec3c> sigma_grad_ay =
            applyMatrixFreeHarmonicWeightedGradient(graph, exact_fields.ay, sigma);
        const StdVec<Vec3c> sigma_grad_az =
            applyMatrixFreeHarmonicWeightedGradient(graph, exact_fields.az, sigma);


        StdVec<Complex> rhs_phi(number_of_particles, Complex(0.0, 0.0));
        const Complex imaginary_unit(0.0, 1.0);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            //const Complex divergence_sigma_a = grad_sigma_ax[i][0] + grad_sigma_ay[i][1] + grad_sigma_az[i][2];
            const Complex divergence_sigma_a = sigma_grad_ax[i][0] + sigma_grad_ay[i][1] + sigma_grad_az[i][2];

            rhs_phi[i] = imaginary_unit * parameters.angular_frequency * divergence_sigma_a;
        }
        ScalarComplexHelmholtzSolverState exact_phi_state =
            solveDiscretePhiForSourceFreeCase(graph, sigma, rhs_phi, exact_fields.phi);
        enforceReference(exact_fields.phi, parameters.phi_reference_index, parameters.phi_reference_value);
        reference_phi_build_residual = exact_phi_state.current_residual_l2_;
        const Real acceptable_reference_residual = 1.0e-5;
        if (!exact_phi_state.converged_ && exact_phi_state.current_residual_l2_ > acceptable_reference_residual)
        {
            std::cerr << "failed_to_build_coupled_source_free_phi_reference"
                      << " iterations=" << exact_phi_state.iterations_
                      << " residual_l2=" << exact_phi_state.current_residual_l2_ << std::endl;
            return 1;
        }
    }

    const StdVec<Complex> laplace_ax = applyScalarNegativeLaplaceFromGraph(graph, exact_fields.ax, nu);
    const StdVec<Complex> laplace_ay = applyScalarNegativeLaplaceFromGraph(graph, exact_fields.ay, nu);
    const StdVec<Complex> laplace_az = applyScalarNegativeLaplaceFromGraph(graph, exact_fields.az, nu);
    const StdVec<Complex> laplace_phi_exact = applyScalarNegativeLaplaceFromGraph(graph, exact_fields.phi, sigma);
    const StdVec<Vec3c> grad_phi = applyMatrixFreeGradient(graph, exact_fields.phi);
    const StdVec<Vec3c> sigma_grad_phi =
        applyMatrixFreeHarmonicWeightedGradient(graph, exact_fields.phi, sigma);

    MatrixFreeAPhiSources sources;
    sources.source_ax.resize(number_of_particles, Complex(0.0, 0.0));
    sources.source_ay.resize(number_of_particles, Complex(0.0, 0.0));
    sources.source_az.resize(number_of_particles, Complex(0.0, 0.0));
    sources.source_phi.resize(number_of_particles, Complex(0.0, 0.0));

    const Complex imaginary_unit(0.0, 1.0);
    // StdVec<Complex> sigma_ax(number_of_particles, Complex(0.0, 0.0));
    // StdVec<Complex> sigma_ay(number_of_particles, Complex(0.0, 0.0));
    // StdVec<Complex> sigma_az(number_of_particles, Complex(0.0, 0.0));
    // for (size_t i = 0; i != number_of_particles; ++i)
    // {
    //     sigma_ax[i] = sigma[i] * exact_fields.ax[i];
    //     sigma_ay[i] = sigma[i] * exact_fields.ay[i];
    //     sigma_az[i] = sigma[i] * exact_fields.az[i];
    // }
    // const StdVec<Vec3c> grad_sigma_ax_for_source = applyMatrixFreeGradient(graph, sigma_ax);
    // const StdVec<Vec3c> grad_sigma_ay_for_source = applyMatrixFreeGradient(graph, sigma_ay);
    // const StdVec<Vec3c> grad_sigma_az_for_source = applyMatrixFreeGradient(graph, sigma_az);

    const StdVec<Vec3c> sigma_grad_ax_for_source =
        applyMatrixFreeHarmonicWeightedGradient(graph, exact_fields.ax, sigma);
    const StdVec<Vec3c> sigma_grad_ay_for_source =
        applyMatrixFreeHarmonicWeightedGradient(graph, exact_fields.ay, sigma);
    const StdVec<Vec3c> sigma_grad_az_for_source =
        applyMatrixFreeHarmonicWeightedGradient(graph, exact_fields.az, sigma);

    const Complex forced_source_amplitude(coil_source_amplitude, 0.0);
    const std::string forced_profile = readEnvString("EM_APHI_FORCED_PROFILE", "sin");
    const Real forced_gaussian_center = readEnvReal("EM_APHI_FORCED_CENTER_X", 0.5 * (coil_xmin + coil_xmax));
    const Real forced_gaussian_width = readEnvReal("EM_APHI_FORCED_WIDTH", 0.5 * (coil_xmax - coil_xmin));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const Complex divergence_sigma_a =
            sigma_grad_ax_for_source[i][0] + sigma_grad_ay_for_source[i][1] + sigma_grad_az_for_source[i][2];


        sources.source_ax[i] = laplace_ax[i] + imaginary_unit * parameters.angular_frequency * sigma[i] * exact_fields.ax[i] +
                               sigma_grad_phi[i][0];
        sources.source_ay[i] = laplace_ay[i] + imaginary_unit * parameters.angular_frequency * sigma[i] * exact_fields.ay[i] +
                               sigma_grad_phi[i][1];
        sources.source_az[i] = laplace_az[i] + imaginary_unit * parameters.angular_frequency * sigma[i] * exact_fields.az[i] +
                               sigma_grad_phi[i][2];
        //sources.source_phi[i] = Complex(0.0, 0.0);
        sources.source_phi[i] =
            laplace_phi_exact[i] - imaginary_unit * parameters.angular_frequency * divergence_sigma_a;

    }

    if (case_mode == "coulomb_variable_sigma_forced_response")
    {
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            Real profile_weight = 0.0;
            bool use_eigen_matched_z_source = false;
            if (is_coil[i])
            {
                if (forced_profile == "eigen_z_sin")
                {
                    const Real xi = (positions[i][0] - coil_xmin) / (coil_xmax - coil_xmin + TinyReal);
                    const Real eta = (positions[i][1] - coil_ymin) / (coil_ymax - coil_ymin + TinyReal);
                    const Real zeta = (positions[i][2] - coil_zmin) / (coil_zmax - coil_zmin + TinyReal);
                    profile_weight = std::sin(Pi * xi) * std::sin(Pi * eta) * std::sin(Pi * zeta);
                    use_eigen_matched_z_source = true;
                }
                else if (forced_profile == "gaussian")
                {
                    const Real width = SMAX(forced_gaussian_width, static_cast<Real>(0.1) * dp_0);
                    const Real normalized_distance = (positions[i][0] - forced_gaussian_center) / width;
                    const Real yz_taper = smooth_box_profile(positions[i], coil_xmin, coil_xmax, coil_ymin, coil_ymax,
                                                             coil_zmin, coil_zmax, 1.0);
                    profile_weight = std::exp(-normalized_distance * normalized_distance) * yz_taper;
                }
                else
                {
                    profile_weight = smooth_box_profile(positions[i], coil_xmin, coil_xmax, coil_ymin, coil_ymax,
                                                        coil_zmin, coil_zmax, 1.0);
                }
            }
            sources.source_ax[i] = Complex(0.0, 0.0);
            sources.source_ay[i] = use_eigen_matched_z_source ? Complex(0.0, 0.0)
                                                              : forced_source_amplitude * profile_weight;
            sources.source_az[i] = use_eigen_matched_z_source ? forced_source_amplitude * profile_weight
                                                              : Complex(0.0, 0.0);
            sources.source_phi[i] = Complex(0.0, 0.0);
            is_source[i] = std::abs(profile_weight) > TinyReal;
            if (is_source[i])
            {
                ++source_particles;
            }
        }
    }

    MatrixFreeAPhiSolverParameters solver_parameters;
    const size_t default_outer_iterations =
        case_mode == "coupled_source_free" || case_mode == "coulomb_variable_sigma_source_free" ||
                case_mode == "coulomb_variable_sigma_driven" ||
                case_mode == "coulomb_variable_sigma_forced_response"
            ? 160
            : 80;
    solver_parameters.max_outer_iterations =
        readEnvSizeT("EM_APHI_STAGGERED_SMOKE_OUTER_ITERS", default_outer_iterations);
    solver_parameters.a_component_solver.max_iterations_ = readEnvSizeT("EM_APHI_STAGGERED_SMOKE_A_ITERS", 800);
    solver_parameters.a_component_solver.absolute_tolerance_ =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_A_TOL", 1.0e-6);
    solver_parameters.phi_solver.max_iterations_ = readEnvSizeT("EM_APHI_STAGGERED_SMOKE_PHI_ITERS", 800);
    solver_parameters.phi_solver.absolute_tolerance_ = readEnvReal("EM_APHI_STAGGERED_SMOKE_PHI_TOL", 1.0e-6);
    solver_parameters.enable_gauge_projection = readEnvBool("EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE", false);
    const std::string gauge_operator_mode = readEnvString("EM_APHI_STAGGERED_SMOKE_GAUGE_OPERATOR", "consistent");
    solver_parameters.use_operator_consistent_gauge_projection = gauge_operator_mode != "laplace";
    if (case_mode == "coupled_source_free" && solver_parameters.enable_gauge_projection)
    {
        std::cerr << "coupled_source_free_does_not_match_coulomb_gauge_in_quasi_1d" << std::endl;
        return 2;
    }
    solver_parameters.residual_tolerance = readEnvReal("EM_APHI_STAGGERED_SMOKE_RESIDUAL_TOL", 1.0e-5);
    solver_parameters.divergence_tolerance = readEnvReal("EM_APHI_STAGGERED_SMOKE_DIVERGENCE_TOL", 1.0e-5);
    solver_parameters.enable_gauge_penalty = readEnvBool("EM_APHI_ENABLE_GAUGE_PENALTY", false);
    solver_parameters.gauge_penalty_coefficient =
        readEnvReal("EM_APHI_GAUGE_PENALTY_COEFF", solver_parameters.enable_gauge_penalty ? 1.0 : 0.0);
    solver_parameters.gauge_penalty_ramp_iterations =
        readEnvSizeT("EM_APHI_GAUGE_PENALTY_RAMP_ITERS", 0);
    solver_parameters.gauge_penalty_initial_ratio =
        readEnvReal("EM_APHI_GAUGE_PENALTY_INITIAL_RATIO", 1.0);
    solver_parameters.outer_relaxation_factor =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_OUTER_RELAX", solver_parameters.enable_gauge_penalty ? 0.8 : 1.0);
    solver_parameters.update_tolerance =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_UPDATE_TOL", 1.0e-7);
    solver_parameters.residual_growth_guard_start_iteration =
        readEnvSizeT("EM_APHI_STAGGERED_SMOKE_GROWTH_GUARD_DELAY", 4);
    solver_parameters.residual_growth_limit =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_GROWTH_LIMIT", 20.0);
    solver_parameters.source_ramp_iterations =
        readEnvSizeT("EM_APHI_SOURCE_RAMP_ITERS", 20);
    solver_parameters.source_initial_ratio =
        readEnvReal("EM_APHI_SOURCE_INITIAL_RATIO", 0.1);
    const size_t load_step_count = readEnvSizeT("EM_APHI_TEAM7LIKE_LOAD_STEPS", 1);
    const Real load_step_initial_ratio =
        readEnvReal("EM_APHI_TEAM7LIKE_LOAD_INITIAL_RATIO", 1.0);

    if (solver_parameters.enable_gauge_penalty)
    {
        MatrixFreeAPhiFields exact_penalty_fields = exact_fields;
        const StdVec<Complex> exact_divergence_a = computeDivergenceOfVectorFieldLocal(graph, exact_penalty_fields);
        const StdVec<Vec3c> exact_gauge_penalty_gradient = applyMatrixFreeGradient(graph, exact_divergence_a);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            sources.source_ax[i] += solver_parameters.gauge_penalty_coefficient * exact_gauge_penalty_gradient[i][0];
            sources.source_ay[i] += solver_parameters.gauge_penalty_coefficient * exact_gauge_penalty_gradient[i][1];
            sources.source_az[i] += solver_parameters.gauge_penalty_coefficient * exact_gauge_penalty_gradient[i][2];
        }
    }

    if (solver_parameters.enable_gauge_projection)
    {
        solver_parameters.gauge_solver.max_iterations_ =
            readEnvSizeT("EM_APHI_STAGGERED_SMOKE_GAUGE_ITERS",
                         solver_parameters.use_operator_consistent_gauge_projection ? 12000 : 500);
        solver_parameters.gauge_solver.relaxation_factor_ =
            readEnvReal("EM_APHI_STAGGERED_SMOKE_GAUGE_RELAX",
                        solver_parameters.use_operator_consistent_gauge_projection ? 5.0e-4 : 1.0);
        solver_parameters.gauge_solver.diagonal_regularization_ =
            readEnvReal("EM_APHI_STAGGERED_SMOKE_GAUGE_DIAG_REG",
                        solver_parameters.use_operator_consistent_gauge_projection ? 1.0e-6 : 1.0e-12);
        solver_parameters.gauge_solver.absolute_tolerance_ =
            readEnvReal("EM_APHI_STAGGERED_SMOKE_GAUGE_TOL", 1.0e-6);
    }

    const Real source_ax_l2 = computeComplexFieldL2NormLocal(sources.source_ax);
    const Real source_ay_l2 = computeComplexFieldL2NormLocal(sources.source_ay);
    const Real source_phi_l2 = computeComplexFieldL2NormLocal(sources.source_phi);

    const bool gauge_only_analysis = readEnvBool("EM_APHI_STAGGERED_SMOKE_GAUGE_ONLY", false);

    if (gauge_only_analysis && !solver_parameters.enable_gauge_projection)
    {
        std::cerr << "gauge_only_requires_enable_gauge" << std::endl;
        return 4;
    }

    const StdVec<Complex> exact_divergence_a_field = computeDivergenceOfVectorFieldLocal(graph, exact_fields);
    const Real exact_divergence_a_l2 =
        has_discrete_manufactured_reference ? computeComplexFieldL2NormLocal(exact_divergence_a_field) : Real(-1.0);

    Real gauge_seed_div_delta_l2 = 0.0;
    Real gauge_seed_operator_response_l2 = 0.0;
    Real gauge_seed_match_l2 = 0.0;
    Real gauge_seed_match_neg_l2 = 0.0;

    MatrixFreeAPhiFields fields;
    fields.ax.resize(number_of_particles, Complex(0.0, 0.0));
    fields.ay.resize(number_of_particles, Complex(0.0, 0.0));
    fields.az.resize(number_of_particles, Complex(0.0, 0.0));
    fields.phi.resize(number_of_particles, Complex(0.0, 0.0));

    if (solver_parameters.enable_gauge_projection)
    {
        const Real seed_amplitude = readEnvReal("EM_APHI_STAGGERED_SMOKE_GAUGE_SEED", 0.2);
        StdVec<Complex> chi_seed(number_of_particles, Complex(0.0, 0.0));
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            const Real x = positions[i][0];
            chi_seed[i] = seed_amplitude * std::cos(Pi * x / body_length);
        }
        const StdVec<Vec3c> grad_chi_seed = applyMatrixFreeGradient(graph, chi_seed);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            fields.ax[i] = exact_fields.ax[i] + grad_chi_seed[i][0];
            fields.ay[i] = exact_fields.ay[i] + grad_chi_seed[i][1];
            fields.az[i] = exact_fields.az[i] + grad_chi_seed[i][2];
            fields.phi[i] = exact_fields.phi[i] - imaginary_unit * parameters.angular_frequency * chi_seed[i];
        }

        const StdVec<Complex> seeded_divergence_a = computeDivergenceOfVectorFieldLocal(graph, fields);
        const StdVec<Complex> operator_response = applyScalarDivergenceOfGradientFromGraph(graph, chi_seed);
        StdVec<Complex> seeded_divergence_delta(number_of_particles, Complex(0.0, 0.0));
        StdVec<Complex> operator_response_negative(number_of_particles, Complex(0.0, 0.0));
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            seeded_divergence_delta[i] = seeded_divergence_a[i] - exact_divergence_a_field[i];
            operator_response_negative[i] = -operator_response[i];
        }
        gauge_seed_div_delta_l2 = computeComplexFieldL2NormLocal(seeded_divergence_delta);
        gauge_seed_operator_response_l2 = computeComplexFieldL2NormLocal(operator_response);
        gauge_seed_match_l2 = computeL2Error(seeded_divergence_delta, operator_response);
        gauge_seed_match_neg_l2 = computeL2Error(seeded_divergence_delta, operator_response_negative);
    }

    if (gauge_only_analysis)
    {
        MatrixFreeAPhiFields pregauge_fields = fields;
        MatrixFreeAPhiSolverParameters pregauge_solver = solver_parameters;
        pregauge_solver.enable_gauge_projection = false;
        pregauge_solver.max_outer_iterations =
            readEnvSizeT("EM_APHI_STAGGERED_SMOKE_PREGAUGE_OUTER_ITERS", default_outer_iterations);

        const MatrixFreeAPhiSolverState pregauge_state = solveMatrixFreeAPhiStaggered(
            graph, pregauge_fields, sigma, nu, parameters, sources, pregauge_solver, sph_native_context);
        const StdVec<Complex> pregauge_divergence_a_field = computeDivergenceOfVectorFieldLocal(graph, pregauge_fields);
        const Real pregauge_divergence_a_error_l2 = computeL2Error(pregauge_divergence_a_field, exact_divergence_a_field);
        const MatrixFreeAPhiGaugeProjectionResult isolated_gauge =
            applyMatrixFreeAPhiGaugeProjectionStep(graph, pregauge_fields, sigma, parameters,
                                                  solver_parameters.gauge_solver,
                                                  solver_parameters.remove_gauge_mean_offset, true,
                                                  solver_parameters.use_operator_consistent_gauge_projection,
                                                  sph_native_context);
        const MatrixFreeAPhiResiduals post_gauge_residuals = evaluateMatrixFreeAPhiResiduals(
            graph, pregauge_fields, sigma, nu, parameters, sources, solver_parameters.enable_gauge_penalty,
            solver_parameters.gauge_penalty_coefficient, sph_native_context);
        const StdVec<Complex> post_gauge_divergence_a_field = computeDivergenceOfVectorFieldLocal(graph, pregauge_fields);
        const Real post_gauge_divergence_a_error_l2 = computeL2Error(post_gauge_divergence_a_field, exact_divergence_a_field);

        StdVec<Complex> reference_rhs = computeDivergenceOfVectorFieldLocal(graph, pregauge_fields);
        removeMeanOffset(reference_rhs);
        Real reference_gauge_solver_residual_l2 = 0.0;
        const StdVec<Complex> reference_chi = solveReferenceGaugeChi(graph, reference_rhs, reference_gauge_solver_residual_l2);
        MatrixFreeAPhiFields reference_gauge_fields = pregauge_fields;
        const StdVec<Vec3c> reference_grad_chi = applyMatrixFreeGradient(graph, reference_chi);
        for (size_t i = 0; i != reference_gauge_fields.ax.size(); ++i)
        {
            reference_gauge_fields.ax[i] -= reference_grad_chi[i][0];
            reference_gauge_fields.ay[i] -= reference_grad_chi[i][1];
            reference_gauge_fields.az[i] -= reference_grad_chi[i][2];
            reference_gauge_fields.phi[i] += imaginary_unit * parameters.angular_frequency * reference_chi[i];
        }
        enforceReference(reference_gauge_fields.phi, parameters.phi_reference_index, parameters.phi_reference_value);
        const MatrixFreeAPhiResiduals reference_post_gauge_residuals = evaluateMatrixFreeAPhiResiduals(
            graph, reference_gauge_fields, sigma, nu, parameters, sources, solver_parameters.enable_gauge_penalty,
            solver_parameters.gauge_penalty_coefficient, sph_native_context);
        const StdVec<Complex> reference_divergence_a = computeDivergenceOfVectorFieldLocal(graph, reference_gauge_fields);
        const Real reference_gauge_div_a_after_l2 = computeComplexFieldL2NormLocal(reference_divergence_a);
        const Real reference_gauge_div_a_error_l2 = computeL2Error(reference_divergence_a, exact_divergence_a_field);
        const Real reference_gauge_chi_l2 = computeComplexFieldL2NormLocal(reference_chi);

        const Real ax_l2_error = computeL2Error(pregauge_fields.ax, exact_fields.ax);
        const Real ay_l2_error = computeL2Error(pregauge_fields.ay, exact_fields.ay);
        const Real phi_l2_error = computeL2Error(pregauge_fields.phi, exact_fields.phi);

        std::cout << std::setprecision(12)
                  << "matrix_free_aphi_gauge_isolation_smoke"
                  << " case_mode=" << case_mode
                  << " sycl_jacobi=" << (use_sycl_jacobi ? 1 : 0)
                  << " sycl_gradient=" << (use_sycl_gradient ? 1 : 0)
                  << " sycl_harmonic_gradient=" << (use_sycl_harmonic_gradient ? 1 : 0)
                  << " sycl_laplace_residual=" << (use_sycl_laplace_residual ? 1 : 0)
                  << " dp=" << dp_0
                  << " points=" << number_of_particles
                  << " gauge_operator_mode=" << gauge_operator_mode
                  << " gauge_penalty_enabled=" << (solver_parameters.enable_gauge_penalty ? 1 : 0)
                  << " gauge_penalty_coeff=" << solver_parameters.gauge_penalty_coefficient
                  << " effective_gauge_penalty_coeff=" << pregauge_state.effective_gauge_penalty_coefficient
                  << " outer_relaxation_factor=" << solver_parameters.outer_relaxation_factor
                  << " update_tolerance=" << solver_parameters.update_tolerance
              << " growth_guard_delay=" << solver_parameters.residual_growth_guard_start_iteration
              << " growth_limit=" << solver_parameters.residual_growth_limit
              << " source_ramp_iters=" << solver_parameters.source_ramp_iterations
              << " source_initial_ratio=" << solver_parameters.source_initial_ratio
              << " growth_guard_delay=" << solver_parameters.residual_growth_guard_start_iteration
              << " growth_limit=" << solver_parameters.residual_growth_limit
                  << " gauge_seed_div_delta_l2=" << gauge_seed_div_delta_l2
                  << " gauge_seed_operator_response_l2=" << gauge_seed_operator_response_l2
                  << " gauge_seed_match_l2=" << gauge_seed_match_l2
                  << " gauge_seed_match_neg_l2=" << gauge_seed_match_neg_l2
                  << " pregauge_outer_iterations=" << pregauge_state.outer_iterations
                  << " pregauge_converged=" << (pregauge_state.converged ? 1 : 0)
                  << " pregauge_residual_ax_l2=" << pregauge_state.residuals.residual_ax_l2
                  << " pregauge_residual_ay_l2=" << pregauge_state.residuals.residual_ay_l2
                  << " pregauge_residual_phi_l2=" << pregauge_state.residuals.residual_phi_l2
                  << " pregauge_divergence_a_l2=" << pregauge_state.residuals.divergence_a_l2
                  << " pregauge_divergence_j_l2=" << pregauge_state.residuals.divergence_j_l2
                  << " exact_divergence_a_l2=" << exact_divergence_a_l2
                  << " pregauge_divergence_a_error_l2=" << pregauge_divergence_a_error_l2
                  << " gauge_solver_iterations=" << isolated_gauge.solver_state.iterations_
                  << " gauge_solver_converged=" << (isolated_gauge.solver_state.converged_ ? 1 : 0)
                  << " gauge_solver_residual_l2=" << isolated_gauge.solver_state.current_residual_l2_
                  << " gauge_solver_min_diagonal_abs=" << isolated_gauge.solver_state.current_min_diagonal_abs_
                  << " gauge_solver_max_diagonal_abs=" << isolated_gauge.solver_state.current_max_diagonal_abs_
                  << " gauge_solver_nonfinite_diagonal_count=" << isolated_gauge.solver_state.current_nonfinite_diagonal_count_
                  << " gauge_rhs_mean_abs_before_projection=" << isolated_gauge.diagnostics.gauge_rhs_mean_abs_before_projection
                  << " gauge_rhs_mean_abs_after_projection=" << isolated_gauge.diagnostics.gauge_rhs_mean_abs_after_projection
                  << " gauge_chi_mean_abs_after_solve=" << isolated_gauge.diagnostics.chi_mean_abs_after_solve
                  << " gauge_chi_l2=" << isolated_gauge.diagnostics.chi_l2
                  << " gauge_grad_chi_l2=" << isolated_gauge.diagnostics.grad_chi_l2
                  << " gauge_chi_max_abs=" << isolated_gauge.diagnostics.chi_max_abs
                  << " gauge_phi_ref_after_phi_abs=" << pregauge_state.gauge_diagnostics.phi_reference_offset_after_phi_solve_abs
                  << " gauge_phi_ref_after_update_abs=" << isolated_gauge.diagnostics.phi_reference_offset_after_gauge_update_abs
                  << " gauge_phi_ref_after_final_abs=" << isolated_gauge.diagnostics.phi_reference_offset_after_final_reference_abs
                  << " gauge_div_a_before_l2=" << isolated_gauge.diagnostics.divergence_a_before_l2
                  << " gauge_div_a_after_raw_l2=" << isolated_gauge.diagnostics.divergence_a_after_raw_l2
                  << " gauge_div_a_after_final_l2=" << isolated_gauge.diagnostics.divergence_a_after_final_l2
                  << " gauge_div_j_before_l2=" << isolated_gauge.diagnostics.divergence_j_before_l2
                  << " gauge_div_j_after_raw_l2=" << isolated_gauge.diagnostics.divergence_j_after_raw_l2
                  << " gauge_div_j_after_final_l2=" << isolated_gauge.diagnostics.divergence_j_after_final_l2
                  << " gauge_e_before_l2=" << isolated_gauge.diagnostics.electric_field_before_l2
                  << " gauge_e_after_raw_l2=" << isolated_gauge.diagnostics.electric_field_after_raw_l2
                  << " gauge_e_after_final_l2=" << isolated_gauge.diagnostics.electric_field_after_final_l2
                  << " gauge_e_change_raw_l2=" << isolated_gauge.diagnostics.electric_field_change_raw_l2
                  << " gauge_e_change_final_l2=" << isolated_gauge.diagnostics.electric_field_change_final_l2
                  << " gauge_j_before_l2=" << isolated_gauge.diagnostics.current_density_before_l2
                  << " gauge_j_after_raw_l2=" << isolated_gauge.diagnostics.current_density_after_raw_l2
                  << " gauge_j_after_final_l2=" << isolated_gauge.diagnostics.current_density_after_final_l2
                  << " gauge_j_change_raw_l2=" << isolated_gauge.diagnostics.current_density_change_raw_l2
                  << " gauge_j_change_final_l2=" << isolated_gauge.diagnostics.current_density_change_final_l2
                  << " post_gauge_residual_ax_l2=" << post_gauge_residuals.residual_ax_l2
                  << " post_gauge_residual_ay_l2=" << post_gauge_residuals.residual_ay_l2
                  << " post_gauge_residual_phi_l2=" << post_gauge_residuals.residual_phi_l2
                  << " post_gauge_divergence_a_l2=" << post_gauge_residuals.divergence_a_l2
                  << " post_gauge_divergence_j_l2=" << post_gauge_residuals.divergence_j_l2
                  << " post_gauge_divergence_a_error_l2=" << post_gauge_divergence_a_error_l2
                  << " reference_gauge_solver_residual_l2=" << reference_gauge_solver_residual_l2
                  << " reference_gauge_chi_l2=" << reference_gauge_chi_l2
                  << " reference_gauge_div_a_after_l2=" << reference_gauge_div_a_after_l2
                  << " reference_gauge_div_a_error_l2=" << reference_gauge_div_a_error_l2
                  << " reference_post_gauge_residual_ax_l2=" << reference_post_gauge_residuals.residual_ax_l2
                  << " reference_post_gauge_residual_phi_l2=" << reference_post_gauge_residuals.residual_phi_l2
                  << " ax_l2_error=" << ax_l2_error
                  << " ay_l2_error=" << ay_l2_error
                  << " phi_l2_error=" << phi_l2_error
                  << std::endl;

        return pregauge_state.outer_iterations == 0 ? 1 : 0;
    }

    MatrixFreeAPhiSolverState state;
    size_t completed_load_steps = 0;
    Real effective_load_scale = 1.0;
    for (size_t load_step = 0; load_step != load_step_count; ++load_step)
    {
        if (load_step_count <= 1)
        {
            effective_load_scale = Real(1.0);
        }
        else
        {
            const Real load_fraction = static_cast<Real>(load_step) /
                                       static_cast<Real>(load_step_count - 1);
            effective_load_scale = load_step_initial_ratio +
                                   (1.0 - load_step_initial_ratio) * load_fraction;
        }

        MatrixFreeAPhiSources staged_sources;
        staged_sources.source_ax.resize(number_of_particles);
        staged_sources.source_ay.resize(number_of_particles);
        staged_sources.source_az.resize(number_of_particles);
        staged_sources.source_phi.resize(number_of_particles);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            staged_sources.source_ax[i] = effective_load_scale * sources.source_ax[i];
            staged_sources.source_ay[i] = effective_load_scale * sources.source_ay[i];
            staged_sources.source_az[i] = effective_load_scale * sources.source_az[i];
            staged_sources.source_phi[i] = effective_load_scale * sources.source_phi[i];
        }

        state = solveMatrixFreeAPhiStaggered(graph, fields, sigma, nu, parameters, staged_sources, solver_parameters,
                                             sph_native_context);
        if (!state.converged)
        {
            break;
        }
        completed_load_steps = load_step + 1;
    }

    const StdVec<Complex> final_divergence_a_field = computeDivergenceOfVectorFieldLocal(graph, fields);
    const Real final_divergence_a_error_l2 =
        has_discrete_manufactured_reference ? computeL2Error(final_divergence_a_field, exact_divergence_a_field) : Real(-1.0);
    const Real ax_l2_error =
        has_discrete_manufactured_reference ? computeL2Error(fields.ax, exact_fields.ax) : Real(-1.0);
    const Real ay_l2_error =
        has_discrete_manufactured_reference ? computeL2Error(fields.ay, exact_fields.ay) : Real(-1.0);
    const Real phi_l2_error =
        has_discrete_manufactured_reference ? computeL2Error(fields.phi, exact_fields.phi) : Real(-1.0);
    const Real ax_max_error =
        has_discrete_manufactured_reference ? computeMaxError(fields.ax, exact_fields.ax) : Real(-1.0);
    const Real ay_max_error =
        has_discrete_manufactured_reference ? computeMaxError(fields.ay, exact_fields.ay) : Real(-1.0);
    const Real phi_max_error =
        has_discrete_manufactured_reference ? computeMaxError(fields.phi, exact_fields.phi) : Real(-1.0);

    const StdVec<Vec3c> exact_electric_field =
        has_discrete_manufactured_reference ? computeElectricField(graph, exact_fields, parameters) : StdVec<Vec3c>();
    const StdVec<Vec3c> exact_current_density =
        has_discrete_manufactured_reference ? computeCurrentDensity(graph, exact_fields, sigma, parameters) : StdVec<Vec3c>();
    const StdVec<Real> exact_joule_density =
        has_discrete_manufactured_reference ? computeJouleHeatingDensity(exact_electric_field, exact_current_density) : StdVec<Real>();
    const StdVec<Vec3c> final_electric_field = computeElectricField(graph, fields, parameters);
    const StdVec<Vec3c> final_current_density = computeCurrentDensity(graph, fields, sigma, parameters);
    const StdVec<Real> final_joule_density = computeJouleHeatingDensity(final_electric_field, final_current_density);
    const Real electric_l2_error =
        has_discrete_manufactured_reference ? computeVectorFieldL2Error(final_electric_field, exact_electric_field) : Real(-1.0);
    const Real current_l2_error =
        has_discrete_manufactured_reference ? computeVectorFieldL2Error(final_current_density, exact_current_density) : Real(-1.0);
    const Real joule_l2_error =
        has_discrete_manufactured_reference ? computeRealFieldL2Error(final_joule_density, exact_joule_density) : Real(-1.0);
    const Real electric_field_l2 = computeVectorComplexFieldL2Norm(final_electric_field);
    const Real current_density_l2 = computeVectorComplexFieldL2Norm(final_current_density);
    const Real joule_density_l2 = computeRealFieldL2Norm(final_joule_density);
    const Real joule_density_max = final_joule_density.empty() ? 0.0 : *std::max_element(final_joule_density.begin(), final_joule_density.end());
    const Real symmetry_tolerance = 0.25 * dp_0;
    const Real ay_mirror_symmetry_l2 =
        computeComplexFieldMirrorSymmetryError(fields.ay, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real phi_mirror_symmetry_l2 =
        computeComplexFieldMirrorSymmetryError(fields.phi, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real phi_mirror_antisymmetry_l2 =
        computeComplexFieldMirrorAntisymmetryError(fields.phi, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const StdVec<Complex> centered_phi = makeMeanCenteredCopy(fields.phi);
    const Real phi_centered_mirror_symmetry_l2 =
        computeComplexFieldMirrorSymmetryError(centered_phi, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real phi_centered_mirror_antisymmetry_l2 =
        computeComplexFieldMirrorAntisymmetryError(centered_phi, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real joule_mirror_symmetry_l2 =
        computeRealFieldMirrorSymmetryError(final_joule_density, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);

    Real *temperature_data = particles.getVariableDataByName<Real>("Temperature");
    Real *joule_heat_source_data = particles.getVariableDataByName<Real>("JouleHeatSource");
    Real *temperature_change_rate_by_joule_data = particles.getVariableDataByName<Real>("TemperatureChangeRateByJoule");
    StdVec<Real> initial_temperature_field(number_of_particles, initial_temperature);
    Real average_joule_density = computeRealFieldAverage(final_joule_density);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        joule_heat_source_data[i] = final_joule_density[i];
        temperature_change_rate_by_joule_data[i] = final_joule_density[i] / rho_cp_thermal;
    }

    Real thermal_time = 0.0;
    size_t thermal_steps = 0;
    while (thermal_time < thermal_end_time - TinyReal)
    {
        Real dt = SMIN(get_thermal_time_step.exec(), thermal_end_time - thermal_time);
        thermal_relaxation.exec(dt);
        add_joule_heat_to_temperature.exec(dt);
        thermal_time += dt;
        ++thermal_steps;
    }

    StdVec<Real> final_temperature_field(number_of_particles, initial_temperature);
    StdVec<Real> temperature_delta_field(number_of_particles, 0.0);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        final_temperature_field[i] = temperature_data[i];
        temperature_delta_field[i] = temperature_data[i] - initial_temperature;
    }
    const Real temperature_average = computeRealFieldAverage(final_temperature_field);
    const Real temperature_delta_average = computeRealFieldAverage(temperature_delta_field);
    const Real temperature_min = final_temperature_field.empty() ? 0.0 : *std::min_element(final_temperature_field.begin(), final_temperature_field.end());
    const Real temperature_max = final_temperature_field.empty() ? 0.0 : *std::max_element(final_temperature_field.begin(), final_temperature_field.end());
    const Real temperature_mirror_symmetry_l2 =
        computeRealFieldMirrorSymmetryError(final_temperature_field, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real temperature_delta_mirror_symmetry_l2 =
        computeRealFieldMirrorSymmetryError(temperature_delta_field, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real expected_temperature_delta_average = average_joule_density * thermal_time / rho_cp_thermal;
    const Real temperature_delta_average_error = std::abs(temperature_delta_average - expected_temperature_delta_average);

    PhysicalRegionSummary conductor_summary;
    PhysicalRegionSummary coil_summary;
    PhysicalRegionSummary air_summary;
    PhysicalRegionSummary source_summary;
    StdVec<Real> source_weights(number_of_particles, 0.0);
    StdVec<Real> conductor_joule_weights(number_of_particles, 0.0);
    StdVec<Real> conductor_temperature_delta_weights(number_of_particles, 0.0);

    auto accumulate_summary = [&](PhysicalRegionSummary &summary, size_t i)
    {
        const Real abs_a = std::sqrt(std::norm(fields.ax[i]) + std::norm(fields.ay[i]) + std::norm(fields.az[i]));
        const Real abs_phi = std::abs(fields.phi[i]);
        const Real abs_e = static_cast<Real>(final_electric_field[i].norm());
        const Real abs_j = static_cast<Real>(final_current_density[i].norm());
        const Real joule = final_joule_density[i];
        const Real temperature_delta = temperature_delta_field[i];
        summary.particles++;
        summary.mean_abs_a += abs_a;
        summary.max_abs_a = SMAX(summary.max_abs_a, abs_a);
        summary.mean_abs_phi += abs_phi;
        summary.max_abs_phi = SMAX(summary.max_abs_phi, abs_phi);
        summary.mean_abs_e += abs_e;
        summary.max_abs_e = SMAX(summary.max_abs_e, abs_e);
        summary.mean_abs_j += abs_j;
        summary.max_abs_j = SMAX(summary.max_abs_j, abs_j);
        summary.mean_joule += joule;
        summary.max_joule = SMAX(summary.max_joule, joule);
        summary.mean_temperature_delta += temperature_delta;
        summary.max_temperature_delta = SMAX(summary.max_temperature_delta, temperature_delta);
    };

    for (size_t i = 0; i != number_of_particles; ++i)
    {
        source_weights[i] = std::abs(sources.source_ax[i]) + std::abs(sources.source_ay[i]) +
                           std::abs(sources.source_az[i]) + std::abs(sources.source_phi[i]);
        if (is_conductor[i])
        {
            accumulate_summary(conductor_summary, i);
            conductor_joule_weights[i] = final_joule_density[i];
            conductor_temperature_delta_weights[i] = temperature_delta_field[i];
        }
        else if (is_coil[i])
        {
            accumulate_summary(coil_summary, i);
        }
        else
        {
            accumulate_summary(air_summary, i);
        }
        if (is_source[i])
        {
            accumulate_summary(source_summary, i);
        }
    }

    finalize_region_summary(conductor_summary);
    finalize_region_summary(coil_summary);
    finalize_region_summary(air_summary);
    finalize_region_summary(source_summary);

    if (write_team7_vtp)
    {
        populateTeam7LikeVisualizationVariables(visualization_fields, fields, final_electric_field, final_current_density,
                                                final_joule_density, temperature_delta_field, is_conductor, is_coil,
                                                is_source);
        body.setNewlyUpdated();
        write_states.writeToFile(0.0);
    }

    const Real source_centroid_x = computeWeightedCentroidX(source_weights, positions);
    const Real conductor_joule_centroid_x = computeWeightedCentroidX(conductor_joule_weights, positions);
    const Real conductor_temperature_delta_centroid_x = computeWeightedCentroidX(conductor_temperature_delta_weights, positions);
    const Real source_centroid_shift = source_centroid_x - 0.5 * body_length;
    const Real conductor_joule_centroid_shift = conductor_joule_centroid_x - 0.5 * body_length;
    const Real conductor_temperature_delta_centroid_shift = conductor_temperature_delta_centroid_x - 0.5 * body_length;
    const Real conductor_joule_relative_to_source_shift = conductor_joule_centroid_x - source_centroid_x;
    const Real conductor_temperature_delta_relative_to_source_shift = conductor_temperature_delta_centroid_x - source_centroid_x;
    const Real conductor_to_air_joule_ratio =
        air_summary.mean_joule <= TinyReal ? 0.0 : conductor_summary.mean_joule / air_summary.mean_joule;
    const Real conductor_to_air_temperature_delta_ratio =
        air_summary.mean_temperature_delta <= TinyReal ? 0.0 : conductor_summary.mean_temperature_delta / air_summary.mean_temperature_delta;
    const bool require_converged = readEnvBool("EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED", false);
    const Real validation_max_residual_ay = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY", -1.0);
    const Real validation_max_divergence_a_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR", -1.0);
    const Real validation_max_ay_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR", -1.0);
    const Real validation_max_phi_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR", -1.0);
    const Real validation_max_e_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_E_ERROR", -1.0);
    const Real validation_max_j_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_J_ERROR", -1.0);
    const Real validation_max_joule_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_JOULE_ERROR", -1.0);
    const Real validation_max_ay_mirror = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_AY_MIRROR", -1.0);
    const Real validation_max_phi_mirror = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_MIRROR", -1.0);
    const Real validation_max_joule_mirror = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_JOULE_MIRROR", -1.0);
    const Real validation_max_phi_antimirror = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_ANTIMIRROR", -1.0);
    const Real validation_max_phi_centered_mirror =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_MIRROR", -1.0);
    const Real validation_max_phi_centered_antimirror =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_ANTIMIRROR", -1.0);
    const Real validation_max_temperature_mirror = readEnvReal("EM_APHI_JOULE_MAX_TEMPERATURE_MIRROR", -1.0);
    const Real validation_max_temperature_delta_mirror = readEnvReal("EM_APHI_JOULE_MAX_TEMPERATURE_DELTA_MIRROR", -1.0);
    const Real validation_max_temperature_delta_average_error = readEnvReal("EM_APHI_JOULE_MAX_TEMPERATURE_DELTA_AVG_ERROR", -1.0);
    const Real validation_min_joule_ratio = readEnvReal("EM_APHI_TEAM7LIKE_MIN_CONDUCTOR_TO_AIR_JOULE_RATIO", -1.0);
    const Real validation_min_temperature_ratio = readEnvReal("EM_APHI_TEAM7LIKE_MIN_CONDUCTOR_TO_AIR_TEMPERATURE_RATIO", -1.0);
    const Real validation_min_joule_shift = readEnvReal("EM_APHI_TEAM7LIKE_MIN_CONDUCTOR_JOULE_SHIFT", -1.0);
    const Real validation_min_temperature_shift = readEnvReal("EM_APHI_TEAM7LIKE_MIN_CONDUCTOR_TEMPERATURE_SHIFT", -1.0);

    bool validation_pass = true;
    if (require_converged && !state.converged)
    {
        validation_pass = false;
    }
    if (validation_max_residual_ay >= Real(0.0) && state.residuals.residual_ay_l2 > validation_max_residual_ay)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_divergence_a_error >= Real(0.0) &&
        final_divergence_a_error_l2 > validation_max_divergence_a_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_ay_error >= Real(0.0) &&
        ay_l2_error > validation_max_ay_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_phi_error >= Real(0.0) &&
        phi_l2_error > validation_max_phi_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_e_error >= Real(0.0) &&
        electric_l2_error > validation_max_e_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_j_error >= Real(0.0) &&
        current_l2_error > validation_max_j_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_joule_error >= Real(0.0) &&
        joule_l2_error > validation_max_joule_error)
    {
        validation_pass = false;
    }
    if (validation_max_ay_mirror >= Real(0.0) && ay_mirror_symmetry_l2 > validation_max_ay_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_phi_mirror >= Real(0.0) && phi_mirror_symmetry_l2 > validation_max_phi_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_joule_mirror >= Real(0.0) && joule_mirror_symmetry_l2 > validation_max_joule_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_phi_antimirror >= Real(0.0) &&
        phi_mirror_antisymmetry_l2 > validation_max_phi_antimirror)
    {
        validation_pass = false;
    }
    if (validation_max_phi_centered_mirror >= Real(0.0) &&
        phi_centered_mirror_symmetry_l2 > validation_max_phi_centered_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_phi_centered_antimirror >= Real(0.0) &&
        phi_centered_mirror_antisymmetry_l2 > validation_max_phi_centered_antimirror)
    {
        validation_pass = false;
    }
    if (validation_max_temperature_mirror >= Real(0.0) &&
        temperature_mirror_symmetry_l2 > validation_max_temperature_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_temperature_delta_mirror >= Real(0.0) &&
        temperature_delta_mirror_symmetry_l2 > validation_max_temperature_delta_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_temperature_delta_average_error >= Real(0.0) &&
        temperature_delta_average_error > validation_max_temperature_delta_average_error)
    {
        validation_pass = false;
    }
    if (validation_min_joule_ratio >= Real(0.0) &&
        conductor_to_air_joule_ratio < validation_min_joule_ratio)
    {
        validation_pass = false;
    }
    if (validation_min_temperature_ratio >= Real(0.0) &&
        conductor_to_air_temperature_delta_ratio < validation_min_temperature_ratio)
    {
        validation_pass = false;
    }
    if (validation_min_joule_shift >= Real(0.0) &&
        conductor_joule_relative_to_source_shift < validation_min_joule_shift)
    {
        validation_pass = false;
    }
    if (validation_min_temperature_shift >= Real(0.0) &&
        conductor_temperature_delta_relative_to_source_shift < validation_min_temperature_shift)
    {
        validation_pass = false;
    }

    int sph_native_device_policy = 0;
    size_t sph_native_standard_gradient_calls = 0;
    size_t sph_native_laplace_helmholtz_calls = 0;
    size_t sph_native_divergence_of_gradient_calls = 0;
    size_t sph_native_sigma_grad_phi_calls = 0;
    size_t sph_native_harmonic_weighted_gradient_calls = 0;
    size_t sph_native_standard_divergence_calls = 0;
    size_t sph_native_harmonic_divergence_calls = 0;
    double sph_native_standard_gradient_seconds = 0.0;
    double sph_native_laplace_helmholtz_seconds = 0.0;
    double sph_native_divergence_of_gradient_seconds = 0.0;
    double sph_native_sigma_grad_phi_seconds = 0.0;
    double sph_native_harmonic_weighted_gradient_seconds = 0.0;
    double sph_native_standard_divergence_seconds = 0.0;
    double sph_native_harmonic_divergence_seconds = 0.0;
    if (sph_native_context != nullptr)
    {
        const auto &sph_native_diagnostics = sph_native_context->diagnostics();
        sph_native_device_policy = sph_native_diagnostics.uses_device_policy_ ? 1 : 0;
        sph_native_standard_gradient_calls = sph_native_diagnostics.standard_gradient_calls_;
        sph_native_laplace_helmholtz_calls = sph_native_diagnostics.laplace_helmholtz_calls_;
        sph_native_divergence_of_gradient_calls = sph_native_diagnostics.divergence_of_gradient_calls_;
        sph_native_sigma_grad_phi_calls = sph_native_diagnostics.sigma_grad_phi_calls_;
        sph_native_harmonic_weighted_gradient_calls = sph_native_diagnostics.harmonic_weighted_gradient_calls_;
        sph_native_standard_divergence_calls = sph_native_diagnostics.standard_divergence_calls_;
        sph_native_harmonic_divergence_calls = sph_native_diagnostics.harmonic_divergence_calls_;
        sph_native_standard_gradient_seconds = sph_native_diagnostics.standard_gradient_seconds_;
        sph_native_laplace_helmholtz_seconds = sph_native_diagnostics.laplace_helmholtz_seconds_;
        sph_native_divergence_of_gradient_seconds = sph_native_diagnostics.divergence_of_gradient_seconds_;
        sph_native_sigma_grad_phi_seconds = sph_native_diagnostics.sigma_grad_phi_seconds_;
        sph_native_harmonic_weighted_gradient_seconds = sph_native_diagnostics.harmonic_weighted_gradient_seconds_;
        sph_native_standard_divergence_seconds = sph_native_diagnostics.standard_divergence_seconds_;
        sph_native_harmonic_divergence_seconds = sph_native_diagnostics.harmonic_divergence_seconds_;
    }

    std::cout << std::setprecision(12)
              << "matrix_free_aphi_team7_like"
              << " case_mode=" << case_mode
              << " sycl_jacobi=" << (use_sycl_jacobi ? 1 : 0)
              << " sycl_gradient=" << (use_sycl_gradient ? 1 : 0)
              << " sycl_harmonic_gradient=" << (use_sycl_harmonic_gradient ? 1 : 0)
              << " sycl_laplace_residual=" << (use_sycl_laplace_residual ? 1 : 0)
              << " sycl_graph_uploads=" << matrixFreeAPhiSyclGraphWorkspaceUploadCount()
              << " sycl_value_field_uploads=" << matrixFreeAPhiSyclValueFieldUploadCount()
              << " sycl_value_field_upload_skips=" << matrixFreeAPhiSyclValueFieldUploadSkipCount()
              << " sycl_value_coefficient_uploads=" << matrixFreeAPhiSyclValueCoefficientUploadCount()
              << " sycl_value_laplace_uploads=" << matrixFreeAPhiSyclValueLaplaceUploadCount()
              << " sycl_value_gradient_downloads=" << matrixFreeAPhiSyclValueGradientDownloadCount()
              << " sycl_value_laplace_downloads=" << matrixFreeAPhiSyclValueLaplaceDownloadCount()
              << " sycl_jacobi_reaction_uploads=" << matrixFreeJacobiReactionUploadCount()
              << " sycl_fused_graph_helmholtz_solves=" << matrixFreeFusedGraphHelmholtzSolveCount()
              << " sycl_fused_graph_helmholtz_metric_downloads=" << matrixFreeFusedGraphHelmholtzMetricDownloadCount()
              << " sycl_fused_graph_helmholtz_metric_scalar_downloads="
              << matrixFreeFusedGraphHelmholtzMetricScalarDownloadCount()
              << " sycl_fused_graph_helmholtz_rhs_uploads=" << matrixFreeFusedGraphHelmholtzRhsUploadCount()
              << " sycl_fused_graph_helmholtz_reaction_uploads=" << matrixFreeFusedGraphHelmholtzReactionUploadCount()
              << " dp=" << dp_0
              << " points=" << number_of_particles
              << " sph_native_residuals=" << (use_sph_native_residuals ? 1 : 0)
              << " sph_native_helmholtz=" << (matrixFreeAPhiUseSphNativeHelmholtz() ? 1 : 0)
              << " sph_native_device_policy=" << sph_native_device_policy
              << " sph_native_standard_gradient_calls=" << sph_native_standard_gradient_calls
              << " sph_native_laplace_helmholtz_calls=" << sph_native_laplace_helmholtz_calls
              << " sph_native_divergence_of_gradient_calls=" << sph_native_divergence_of_gradient_calls
              << " sph_native_sigma_grad_phi_calls=" << sph_native_sigma_grad_phi_calls
              << " sph_native_harmonic_weighted_gradient_calls=" << sph_native_harmonic_weighted_gradient_calls
              << " sph_native_standard_divergence_calls=" << sph_native_standard_divergence_calls
              << " sph_native_harmonic_divergence_calls=" << sph_native_harmonic_divergence_calls
              << " sph_native_standard_gradient_seconds=" << sph_native_standard_gradient_seconds
              << " sph_native_laplace_helmholtz_seconds=" << sph_native_laplace_helmholtz_seconds
              << " sph_native_divergence_of_gradient_seconds=" << sph_native_divergence_of_gradient_seconds
              << " sph_native_sigma_grad_phi_seconds=" << sph_native_sigma_grad_phi_seconds
              << " sph_native_harmonic_weighted_gradient_seconds=" << sph_native_harmonic_weighted_gradient_seconds
              << " sph_native_standard_divergence_seconds=" << sph_native_standard_divergence_seconds
              << " sph_native_harmonic_divergence_seconds=" << sph_native_harmonic_divergence_seconds
              << " operator_backend=" << (use_sph_native_residuals ? "sph_native_ck" : "graph_cpu_sycl_jacobi")
              << " has_discrete_reference=" << (has_discrete_manufactured_reference ? 1 : 0)
              << " gauge_operator_mode=" << gauge_operator_mode
              << " gauge_enabled=" << (solver_parameters.enable_gauge_projection ? 1 : 0)
              << " gauge_penalty_enabled=" << (solver_parameters.enable_gauge_penalty ? 1 : 0)
              << " gauge_penalty_coeff=" << solver_parameters.gauge_penalty_coefficient
              << " effective_gauge_penalty_coeff=" << state.effective_gauge_penalty_coefficient
              << " effective_source_scale=" << state.effective_source_scale
              << " outer_relaxation_factor=" << solver_parameters.outer_relaxation_factor
              << " update_tolerance=" << solver_parameters.update_tolerance
              << " growth_guard_delay=" << solver_parameters.residual_growth_guard_start_iteration
              << " growth_limit=" << solver_parameters.residual_growth_limit
              << " source_ramp_iters=" << solver_parameters.source_ramp_iterations
              << " source_initial_ratio=" << solver_parameters.source_initial_ratio
              << " load_steps=" << load_step_count
              << " load_initial_ratio=" << load_step_initial_ratio
              << " completed_load_steps=" << completed_load_steps
              << " effective_load_scale=" << effective_load_scale
              << " field_update_l2=" << state.field_update_l2
              << " relative_field_update_l2=" << state.relative_field_update_l2
              << " outer_iterations=" << state.outer_iterations
              << " converged=" << (state.converged ? 1 : 0)
              << " residual_ax_l2=" << state.residuals.residual_ax_l2
              << " residual_ay_l2=" << state.residuals.residual_ay_l2
              << " residual_phi_l2=" << state.residuals.residual_phi_l2
              << " divergence_a_l2=" << state.residuals.divergence_a_l2
              << " divergence_j_l2=" << state.residuals.divergence_j_l2
              << " interface_sigma_ratio_threshold=" << parameters.interface_contrast_threshold
              << " interface_edge_count=" << state.residuals.interface_edge_count
              << " high_contrast_edge_count=" << state.residuals.high_contrast_edge_count
              << " high_contrast_particle_count=" << state.residuals.high_contrast_particle_count
              << " max_interface_sigma_ratio=" << state.residuals.max_interface_sigma_ratio
              << " high_contrast_sigma_grad_phi_l2=" << state.residuals.high_contrast_sigma_grad_phi_l2
              << " high_contrast_div_sigma_a_l2=" << state.residuals.high_contrast_div_sigma_a_l2
              << " high_contrast_residual_a_l2=" << state.residuals.high_contrast_residual_a_l2
              << " high_contrast_residual_phi_l2=" << state.residuals.high_contrast_residual_phi_l2
              << " exact_divergence_a_l2=" << exact_divergence_a_l2
              << " divergence_a_error_l2=" << final_divergence_a_error_l2
              << " ax_l2_error=" << ax_l2_error
              << " ay_l2_error=" << ay_l2_error
              << " phi_l2_error=" << phi_l2_error
              << " electric_l2_error=" << electric_l2_error
              << " current_l2_error=" << current_l2_error
              << " joule_l2_error=" << joule_l2_error
              << " electric_field_l2=" << electric_field_l2
              << " current_density_l2=" << current_density_l2
              << " joule_density_l2=" << joule_density_l2
              << " joule_density_max=" << joule_density_max
              << " thermal_time=" << thermal_time
              << " thermal_steps=" << thermal_steps
              << " temperature_average=" << temperature_average
              << " temperature_min=" << temperature_min
              << " temperature_max=" << temperature_max
              << " temperature_delta_average=" << temperature_delta_average
              << " expected_temperature_delta_average=" << expected_temperature_delta_average
              << " temperature_delta_average_error=" << temperature_delta_average_error
              << " conductor_particles=" << conductor_particles
              << " coil_particles=" << coil_particles
              << " air_particles=" << air_particles
              << " source_particles=" << source_particles
              << " source_centroid_x=" << source_centroid_x
              << " conductor_joule_centroid_x=" << conductor_joule_centroid_x
              << " conductor_temperature_delta_centroid_x=" << conductor_temperature_delta_centroid_x
              << " source_centroid_shift=" << source_centroid_shift
              << " conductor_joule_centroid_shift=" << conductor_joule_centroid_shift
              << " conductor_temperature_delta_centroid_shift=" << conductor_temperature_delta_centroid_shift
              << " conductor_joule_relative_to_source_shift=" << conductor_joule_relative_to_source_shift
              << " conductor_temperature_delta_relative_to_source_shift=" << conductor_temperature_delta_relative_to_source_shift
              << " conductor_to_air_joule_ratio=" << conductor_to_air_joule_ratio
              << " conductor_to_air_temperature_delta_ratio=" << conductor_to_air_temperature_delta_ratio
              << " conductor_mean_abs_a=" << conductor_summary.mean_abs_a
              << " conductor_mean_abs_phi=" << conductor_summary.mean_abs_phi
              << " conductor_mean_abs_e=" << conductor_summary.mean_abs_e
              << " conductor_mean_abs_j=" << conductor_summary.mean_abs_j
              << " conductor_mean_joule=" << conductor_summary.mean_joule
              << " conductor_max_joule=" << conductor_summary.max_joule
              << " conductor_mean_temperature_delta=" << conductor_summary.mean_temperature_delta
              << " conductor_max_temperature_delta=" << conductor_summary.max_temperature_delta
              << " coil_mean_joule=" << coil_summary.mean_joule
              << " coil_max_joule=" << coil_summary.max_joule
              << " air_mean_joule=" << air_summary.mean_joule
              << " air_max_joule=" << air_summary.max_joule
              << " source_mean_joule=" << source_summary.mean_joule
              << " source_max_joule=" << source_summary.max_joule
              << " temperature_mirror_symmetry_l2=" << temperature_mirror_symmetry_l2
              << " temperature_delta_mirror_symmetry_l2=" << temperature_delta_mirror_symmetry_l2
              << " ay_mirror_symmetry_l2=" << ay_mirror_symmetry_l2
              << " phi_mirror_symmetry_l2=" << phi_mirror_symmetry_l2
              << " phi_mirror_antisymmetry_l2=" << phi_mirror_antisymmetry_l2
              << " phi_centered_mirror_symmetry_l2=" << phi_centered_mirror_symmetry_l2
              << " phi_centered_mirror_antisymmetry_l2=" << phi_centered_mirror_antisymmetry_l2
              << " joule_mirror_symmetry_l2=" << joule_mirror_symmetry_l2
              << " ax_max_error=" << ax_max_error
              << " ay_max_error=" << ay_max_error
              << " phi_max_error=" << phi_max_error
              << " source_ax_l2=" << source_ax_l2
              << " source_ay_l2=" << source_ay_l2
              << " source_phi_l2=" << source_phi_l2
              << " forced_profile=" << forced_profile
              << " forced_center_x=" << forced_gaussian_center
              << " forced_width=" << forced_gaussian_width
              << " reference_phi_build_residual=" << reference_phi_build_residual
              << " gauge_applied=" << (state.gauge_diagnostics.applied ? 1 : 0)
              << " gauge_solver_min_diagonal_abs=" << state.gauge_state.current_min_diagonal_abs_
              << " gauge_solver_max_diagonal_abs=" << state.gauge_state.current_max_diagonal_abs_
              << " gauge_solver_nonfinite_diagonal_count=" << state.gauge_state.current_nonfinite_diagonal_count_
              << " gauge_rhs_mean_abs_before_projection=" << state.gauge_diagnostics.gauge_rhs_mean_abs_before_projection
              << " gauge_rhs_mean_abs_after_projection=" << state.gauge_diagnostics.gauge_rhs_mean_abs_after_projection
              << " gauge_chi_mean_abs_after_solve=" << state.gauge_diagnostics.chi_mean_abs_after_solve
              << " gauge_chi_l2=" << state.gauge_diagnostics.chi_l2
              << " gauge_grad_chi_l2=" << state.gauge_diagnostics.grad_chi_l2
              << " gauge_chi_max_abs=" << state.gauge_diagnostics.chi_max_abs
              << " gauge_phi_ref_after_phi_abs=" << state.gauge_diagnostics.phi_reference_offset_after_phi_solve_abs
              << " gauge_phi_ref_after_update_abs=" << state.gauge_diagnostics.phi_reference_offset_after_gauge_update_abs
              << " gauge_phi_ref_after_final_abs=" << state.gauge_diagnostics.phi_reference_offset_after_final_reference_abs
              << " gauge_div_a_before_l2=" << state.gauge_diagnostics.divergence_a_before_l2
              << " gauge_div_a_after_raw_l2=" << state.gauge_diagnostics.divergence_a_after_raw_l2
              << " gauge_div_a_after_final_l2=" << state.gauge_diagnostics.divergence_a_after_final_l2
              << " gauge_div_j_before_l2=" << state.gauge_diagnostics.divergence_j_before_l2
              << " gauge_div_j_after_raw_l2=" << state.gauge_diagnostics.divergence_j_after_raw_l2
              << " gauge_div_j_after_final_l2=" << state.gauge_diagnostics.divergence_j_after_final_l2
              << " gauge_e_before_l2=" << state.gauge_diagnostics.electric_field_before_l2
              << " gauge_e_after_raw_l2=" << state.gauge_diagnostics.electric_field_after_raw_l2
              << " gauge_e_after_final_l2=" << state.gauge_diagnostics.electric_field_after_final_l2
              << " gauge_e_change_raw_l2=" << state.gauge_diagnostics.electric_field_change_raw_l2
              << " gauge_e_change_final_l2=" << state.gauge_diagnostics.electric_field_change_final_l2
              << " gauge_j_before_l2=" << state.gauge_diagnostics.current_density_before_l2
              << " gauge_j_after_raw_l2=" << state.gauge_diagnostics.current_density_after_raw_l2
              << " gauge_j_after_final_l2=" << state.gauge_diagnostics.current_density_after_final_l2
              << " gauge_j_change_raw_l2=" << state.gauge_diagnostics.current_density_change_raw_l2
              << " gauge_j_change_final_l2=" << state.gauge_diagnostics.current_density_change_final_l2
              << " validation_require_converged=" << (require_converged ? 1 : 0)
              << " validation_max_residual_ay=" << validation_max_residual_ay
              << " validation_max_diva_error=" << validation_max_divergence_a_error
              << " validation_max_ay_error=" << validation_max_ay_error
              << " validation_max_phi_error=" << validation_max_phi_error
              << " validation_max_e_error=" << validation_max_e_error
              << " validation_max_j_error=" << validation_max_j_error
              << " validation_max_joule_error=" << validation_max_joule_error
              << " validation_max_ay_mirror=" << validation_max_ay_mirror
              << " validation_max_phi_mirror=" << validation_max_phi_mirror
              << " validation_max_joule_mirror=" << validation_max_joule_mirror
              << " validation_max_phi_antimirror=" << validation_max_phi_antimirror
              << " validation_max_phi_centered_mirror=" << validation_max_phi_centered_mirror
              << " validation_max_phi_centered_antimirror=" << validation_max_phi_centered_antimirror
              << " validation_max_temperature_mirror=" << validation_max_temperature_mirror
              << " validation_max_temperature_delta_mirror=" << validation_max_temperature_delta_mirror
              << " validation_max_temperature_delta_avg_error=" << validation_max_temperature_delta_average_error
              << " validation_min_joule_ratio=" << validation_min_joule_ratio
              << " validation_min_temperature_ratio=" << validation_min_temperature_ratio
              << " validation_min_joule_shift=" << validation_min_joule_shift
              << " validation_min_temperature_shift=" << validation_min_temperature_shift
              << " validation_pass=" << (validation_pass ? 1 : 0)
              << " write_team7_vtp=" << (write_team7_vtp ? 1 : 0)
              << std::endl;

    if (state.outer_iterations == 0)
    {
        return 1;
    }
    return validation_pass ? 0 : 4;
}
