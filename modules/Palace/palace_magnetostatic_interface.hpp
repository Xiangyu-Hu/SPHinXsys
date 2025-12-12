// palace_magnetostatic_interface.hpp
#pragma once

#include <memory>
#include <string>
#include <vector>

#include <mpi.h>
#include <mfem.hpp>

// Palace headers required in the interface
#include "fem/mesh.hpp"                    // palace::Mesh
#include "drivers/magnetostaticsolver.hpp" // palace::MagnetostaticSolver
#include "linalg/vector.hpp"               // palace::Vector
#include "utils/communication.hpp"         // palace::Mpi
#include "utils/iodata.hpp"                // palace::IoData

namespace palace
{
class CurlCurlOperator;  // forward declaration
}

namespace sphinxsys_palace
{

// Wrapper for a single Palace magnetostatic simulation.
//
// Responsibilities:
//   - Load JSON configuration (IoData).
//   - Initialize device, Hypre, SLEPc/PETSc.
//   - Build mesh hierarchy (ParMesh -> palace::Mesh).
//   - Construct and run palace::MagnetostaticSolver.
//   - On the final adapted mesh, rebuild CurlCurlOperator and compute
//     B = curl(A) for each current source boundary.
//   - Store B both as true-dof vectors (palace::Vector) and as
//     mfem::ParGridFunction for convenient access and postprocessing.
//
// IMPORTANT:
//   - We assume palace::Mpi::Init(argc, argv) has already been called in main().
class MagnetostaticCase
{
public:
    // json_file: Path to Palace JSON configuration file.
    // comm     : MPI communicator to use (typically MPI_COMM_WORLD).
    // verbose  : If true, prints some SphinxSys-Palace info messages.
    explicit MagnetostaticCase(const std::string &json_file,
                               MPI_Comm comm = MPI_COMM_WORLD,
                               bool verbose = false);

    MagnetostaticCase(const MagnetostaticCase &) = delete;
    MagnetostaticCase &operator=(const MagnetostaticCase &) = delete;

    MagnetostaticCase(MagnetostaticCase &&) noexcept = default;
    MagnetostaticCase &operator=(MagnetostaticCase &&) noexcept = default;

    ~MagnetostaticCase();

    // Initialize Palace objects:
    //   - Read IoData from JSON.
    //   - Create output directory.
    //   - Configure MFEM device + libCEED backend.
    //   - Initialize Hypre and (optionally) SLEPc/PETSc.
    //   - Build mesh hierarchy (ParMesh -> palace::Mesh).
    void Initialize();

    // Solve the magnetostatic problem and compute B fields:
    //   - Run MagnetostaticSolver::SolveEstimateMarkRefine(mesh_).
    //   - On the final mesh, build CurlCurlOperator and, for each source,
    //       * assemble excitation RHS,
    //       * solve K A = RHS,
    //       * compute B = Curl * A,
    //       * store B as true-dof vector and as ParGridFunction.
    //   - Print Palace timing summary and save solver metadata.
    //   - Finalize libCEED and SLEPc if enabled.
    void Solve();

    // Convenience method: equivalent to { Initialize(); Solve(); }.
    void Run();

    // --------- Accessors for coupling with SPHinXsys ---------

    // Access Palace mesh (usually there is only one mesh).
    palace::Mesh &GetMesh(std::size_t i = 0);
    const palace::Mesh &GetMesh(std::size_t i = 0) const;

    // Direct access to underlying Palace solver, e.g., for metadata queries.
    palace::MagnetostaticSolver &GetSolver();
    const palace::MagnetostaticSolver &GetSolver() const;

    // Check whether Initialize() has completed successfully.
    bool IsInitialized() const { return initialized_; }

    // --------- Accessors for magnetic field results ---------

    // Number of surface-current source boundaries (i.e., number of independent B fields).
    int GetNumSources() const;

    // Get B field for a given source as a true-dof vector in the RT space.
    // You must call Solve() or Run() before calling this.
    const palace::Vector &GetBTrueDofs(int source_index = 0) const;

    // Get B field for a given source as an mfem::ParGridFunction defined on the RT space
    // of the final adapted mesh. This is convenient for pointwise evaluation and VTK output.
    const mfem::ParGridFunction &GetBField(int source_index = 0) const;

    // Save B field of a given source in ParaView (VTK) format (.pvtu + .vtu).
    // prefix     : base file name prefix (e.g. "rings_B").
    // output_dir : directory where files are written (e.g. "postpro/sphinxsys").
    void SaveBFieldParaView(const std::string &prefix,
                            const std::string &output_dir,
                            int source_index = 0) const;

private:
    // Internal helper: after Palace driver has finished its adaptive solve on mesh_,
    // rebuild CurlCurlOperator on the final mesh and compute B = Curl(A) for each
    // source boundary. Store the results in b_fields_true_ and b_fields_grid_.
    void ComputeMagneticFields();

private:
    // Input configuration
    std::string json_file_;
    MPI_Comm    comm_;
    bool        verbose_;

    // Initialization status
    bool initialized_ = false;

    // MPI information
    bool root_ = false;
    int  size_ = 0;

    // Thread / device information
    int omp_threads_ = 0;
    int ngpu_        = 0;

    // Core Palace objects
    std::unique_ptr<palace::IoData>              iodata_;
    std::vector<std::unique_ptr<palace::Mesh>>   mesh_;
    std::unique_ptr<palace::MagnetostaticSolver> solver_;

    // Curl-curl operator on the final adapted mesh, used only in the interface layer
    // to reconstruct A and B fields, without modifying Palace core code.
    std::unique_ptr<palace::CurlCurlOperator> curl_op_;

    // Magnetic field storage:
    //   - One B(true-dof vector) per source in the RT space.
    //   - One ParGridFunction per source for convenient access/postprocessing.
    std::vector<palace::Vector> b_fields_true_;
    std::vector<std::unique_ptr<mfem::ParGridFunction>> b_fields_grid_;
};

// C-style helper wrapper for test cases (e.g., test_3d_palace_rings_M):
//   - Assumes palace::Mpi::Init(argc, argv) has already been called.
//   - Returns 0 on success, non-zero on failure.
int RunMagnetostaticCase(const std::string &json_file,
                         MPI_Comm comm = MPI_COMM_WORLD,
                         bool verbose = false);

}  // namespace sphinxsys_palace
