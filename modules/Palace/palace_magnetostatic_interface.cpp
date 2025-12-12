// palace_magnetostatic_interface.cpp
#include "palace_magnetostatic_interface.hpp"

#include <sstream>
#include <cstring>
#include <stdexcept>

#include <mfem.hpp>

#include "drivers/magnetostaticsolver.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "linalg/hypre.hpp"
#include "linalg/ksp.hpp"
#include "linalg/slepc.hpp"
#include "models/curlcurloperator.hpp"
#include "utils/device.hpp"
#include "utils/geodata.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

#if defined(MFEM_USE_STRUMPACK)
#include <StrumpackConfig.hpp>
#endif

using namespace palace;

namespace sphinxsys_palace
{

// Optional: expose Palace git tag if available.
static const char *GetPalaceGitTag()
{
#if defined(PALACE_GIT_COMMIT)
  static const char *commit = PALACE_GIT_COMMIT_ID;
#else
  static const char *commit = "UNKNOWN";
#endif
  return commit;
}

// Configure MFEM device string from Palace Device enum.
static std::string ConfigureDevice(Device device)
{
  std::string device_str;
  switch (device)
  {
    case Device::CPU:
      device_str = "cpu";
      break;
    case Device::GPU:
#if defined(MFEM_USE_CUDA)
      device_str = "cuda";
#elif defined(MFEM_USE_HIP)
      device_str = "hip";
#else
      Mpi::Warning("Palace must be built with either CUDA or HIP support for GPU device "
                   "usage, reverting to CPU!\n");
      device_str = "cpu";
#endif
      break;
    case Device::DEBUG:
      device_str = "cpu,debug";
      break;
  }
#if defined(MFEM_USE_OPENMP)
  device_str += ",omp";
#endif
  return device_str;
}

// Configure libCEED backend, similar to Palace main.cpp (JIT dir left empty here).
static void ConfigureCeedBackend(const std::string &ceed_backend)
{
  std::string default_ceed_backend;
  if (mfem::Device::Allows(mfem::Backend::CUDA_MASK))
  {
    default_ceed_backend = "/gpu/cuda/magma";
  }
  else if (mfem::Device::Allows(mfem::Backend::HIP_MASK))
  {
    default_ceed_backend = "/gpu/hip/magma";
  }
  else if (mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE))
  {
    default_ceed_backend = "/cpu/self/ref/serial";
  }
  else
  {
    default_ceed_backend = "/cpu/self";
  }

  const std::string &backend =
      !ceed_backend.empty() ? ceed_backend : default_ceed_backend;

  ceed::Initialize(backend.c_str(), /* JIT source dir */ "");

  // Optional safety check: verify the backend actually used by libCEED.
  std::string ceed_resource = ceed::Print();
  if (backend.compare(0, backend.length(), ceed_resource, 0, backend.length()))
  {
    Mpi::Warning(
        "libCEED is not using the requested backend!\nRequested \"{}\", got \"{}\"!\n",
        backend, ceed_resource);
  }
}

// ---------------- MagnetostaticCase implementation ----------------

MagnetostaticCase::MagnetostaticCase(const std::string &json_file,
                                     MPI_Comm comm,
                                     bool verbose)
    : json_file_(json_file),
      comm_(comm),
      verbose_(verbose)
{
}

MagnetostaticCase::~MagnetostaticCase() = default;

void MagnetostaticCase::Initialize()
{
  if (initialized_)
  {
    return;
  }

  // We assume palace::Mpi::Init(argc, argv) has already been called in the outer code.
  MPI_Comm world_comm = Mpi::World();
  root_ = Mpi::Root(world_comm);
  size_ = Mpi::Size(world_comm);

  BlockTimer bt(Timer::INIT);

  if (verbose_ && root_)
  {
    Mpi::Print(world_comm,
               "[SphinxSys-Palace] MagnetostaticCase::Initialize using config: {}\n",
               json_file_);
  }

  // 1) Read IoData from JSON configuration.
  iodata_ = std::make_unique<IoData>(json_file_.c_str(), false);

  // 2) Create output directory.
  MakeOutputFolder(*iodata_, world_comm);

  // 3) Configure MFEM device, OpenMP threads and available GPUs.
  omp_threads_ = utils::ConfigureOmp();
  ngpu_        = utils::GetDeviceCount();

  mfem::Device device(ConfigureDevice(iodata_->solver.device),
                      utils::GetDeviceId(world_comm, ngpu_));
  ConfigureCeedBackend(iodata_->solver.ceed_backend);
#if defined(PALACE_WITH_GPU_AWARE_MPI)
  device.SetGPUAwareMPI(true);
#endif

  // 4) Initialize Hypre and optionally SLEPc/PETSc.
  hypre::Initialize();
#if defined(PALACE_WITH_SLEPC)
  {
    int argc = 0;
    char **argv = nullptr;
    slepc::Initialize(argc, argv, nullptr, nullptr);
    if (PETSC_COMM_WORLD != world_comm)
    {
      Mpi::Print(world_comm, "Error: Problem during MPI initialization!\n\n");
      throw std::runtime_error("SLEPc/PETSc MPI communicator mismatch!");
    }
  }
#endif

  // 5) Build mesh: read ParMesh, nondimensionalize, refine, then wrap in palace::Mesh.
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    mfem_mesh.push_back(mesh::ReadMesh(*iodata_, world_comm));

    iodata_->NondimensionalizeInputs(*mfem_mesh[0]);
    mesh::RefineMesh(*iodata_, mfem_mesh);

    for (auto &m : mfem_mesh)
    {
      mesh_.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }

  // 6) Construct the Palace MagnetostaticSolver driver.
  solver_ = std::make_unique<MagnetostaticSolver>(
      *iodata_, root_, size_, omp_threads_, GetPalaceGitTag());

  initialized_ = true;
}

void MagnetostaticCase::Solve()
{
  if (!initialized_)
  {
    Initialize();
  }

  MPI_Comm world_comm = Mpi::World();

  if (verbose_ && root_)
  {
    Mpi::Print(world_comm, "[SphinxSys-Palace] MagnetostaticCase::Solve\n");
  }

  // 1) Call Palace driver: adaptive solve / estimate / refine.
  solver_->SolveEstimateMarkRefine(mesh_);

  // 2) On the final mesh, compute B fields and build ParGridFunction representations.
  ComputeMagneticFields();

  // 3) Print timing summary and save metadata.
  BlockTimer::Print(world_comm);
  solver_->SaveMetadata(BlockTimer::GlobalTimer());

  // 4) Finalize libCEED resources.
  ceed::Finalize();

#if defined(PALACE_WITH_SLEPC)
  slepc::Finalize();
#endif
}

void MagnetostaticCase::Run()
{
  Initialize();
  Solve();
}

// --------- Accessors ---------

palace::Mesh &MagnetostaticCase::GetMesh(std::size_t i)
{
  return *mesh_.at(i);
}

const palace::Mesh &MagnetostaticCase::GetMesh(std::size_t i) const
{
  return *mesh_.at(i);
}

palace::MagnetostaticSolver &MagnetostaticCase::GetSolver()
{
  return *solver_;
}

const palace::MagnetostaticSolver &MagnetostaticCase::GetSolver() const
{
  return *solver_;
}

// --------- Magnetic field computation ---------

void MagnetostaticCase::ComputeMagneticFields()
{
  MPI_Comm world_comm = Mpi::World();

  if (mesh_.empty())
  {
    throw std::runtime_error(
        "MagnetostaticCase::ComputeMagneticFields: no mesh!");
  }

  if (verbose_ && root_)
  {
    Mpi::Print(world_comm,
               "[SphinxSys-Palace] MagnetostaticCase::ComputeMagneticFields\n");
  }

  // Build curl-curl operator on the final adapted mesh and keep it alive as a member,
  // so that its finite element spaces remain valid for ParGridFunction objects.
  curl_op_ = std::make_unique<CurlCurlOperator>(*iodata_, mesh_);
  CurlCurlOperator &curlcurl_op = *curl_op_;

  auto K = curlcurl_op.GetStiffnessMatrix();
  const auto &Curl = curlcurl_op.GetCurlMatrix();

  // Linear solver for the magnetostatic system.
  KspSolver ksp(*iodata_, curlcurl_op.GetNDSpaces(), &curlcurl_op.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  // Number of current source boundaries.
  int n_step = static_cast<int>(curlcurl_op.GetSurfaceCurrentOp().Size());
  MFEM_VERIFY(n_step > 0,
              "MagnetostaticCase::ComputeMagneticFields: no surface current sources!");

  Vector RHS(Curl.Width()), B(Curl.Height());
  std::vector<Vector> A(n_step);

  b_fields_true_.clear();
  b_fields_true_.resize(n_step);
  b_fields_grid_.clear();
  b_fields_grid_.resize(n_step);

  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : curlcurl_op.GetSurfaceCurrentOp())
  {
    if (verbose_ && root_)
    {
      Mpi::Print(world_comm,
                 "  Computing field for source %d (step %d/%d, elapsed = %.2e s)\n",
                 idx, step + 1, n_step,
                 Timer::Duration(Timer::Now() - t0).count());
    }

    A[step].SetSize(RHS.Size());
    A[step].UseDevice(true);
    A[step] = 0.0;

    // RHS for this source and solve for A.
    curlcurl_op.GetExcitationVector(idx, RHS);
    ksp.Mult(RHS, A[step]);

    // B = Curl * A.
    Curl.Mult(A[step], B);

    // Store B as true-dof vector.
    b_fields_true_[step].SetSize(B.Size());
    b_fields_true_[step].UseDevice(false);
    b_fields_true_[step] = B;

    // Build ParGridFunction for B in the RT space.
    auto &rt_space = curlcurl_op.GetRTSpace();          // palace::FiniteElementSpace
    mfem::ParFiniteElementSpace &rt_fes = rt_space.Get();

    if (!b_fields_grid_[step])
    {
      b_fields_grid_[step] = std::make_unique<mfem::ParGridFunction>(&rt_fes);
    }
    b_fields_grid_[step]->SetFromTrueDofs(B);

    step++;
  }

  if (verbose_ && root_)
  {
    Mpi::Print(world_comm,
               "[SphinxSys-Palace] MagnetostaticCase::ComputeMagneticFields: "
               "stored %d B-fields (true dofs + ParGridFunction).\n",
               n_step);
  }
}

int MagnetostaticCase::GetNumSources() const
{
  return static_cast<int>(b_fields_true_.size());
}

const palace::Vector &MagnetostaticCase::GetBTrueDofs(int source_index) const
{
  if (b_fields_true_.empty())
  {
    throw std::runtime_error(
        "MagnetostaticCase::GetBTrueDofs: no fields computed; call Solve() or Run() first.");
  }
  if (source_index < 0 ||
      source_index >= static_cast<int>(b_fields_true_.size()))
  {
    throw std::out_of_range(
        "MagnetostaticCase::GetBTrueDofs: invalid source_index.");
  }
  return b_fields_true_.at(static_cast<std::size_t>(source_index));
}

const mfem::ParGridFunction &
MagnetostaticCase::GetBField(int source_index) const
{
  if (b_fields_grid_.empty())
  {
    throw std::runtime_error(
        "MagnetostaticCase::GetBField: no fields computed; call Solve() or Run() first.");
  }
  if (source_index < 0 ||
      source_index >= static_cast<int>(b_fields_grid_.size()))
  {
    throw std::out_of_range(
        "MagnetostaticCase::GetBField: invalid source_index.");
  }
  const auto &ptr = b_fields_grid_.at(static_cast<std::size_t>(source_index));
  if (!ptr)
  {
    throw std::runtime_error(
        "MagnetostaticCase::GetBField: internal ParGridFunction is null.");
  }
  return *ptr;
}

void MagnetostaticCase::SaveBFieldParaView(const std::string &prefix,
                                           const std::string &output_dir,
                                           int source_index) const
{
  const mfem::ParGridFunction &B = GetBField(source_index);

  if (mesh_.empty())
  {
    throw std::runtime_error(
        "MagnetostaticCase::SaveBFieldParaView: no mesh available.");
  }

  palace::Mesh &pm = *mesh_.front();
  mfem::ParMesh &pmesh = pm;  // use Mesh -> ParMesh conversion operator

  mfem::ParaViewDataCollection dc(prefix, &pmesh);
  dc.RegisterField("B", const_cast<mfem::ParGridFunction *>(&B));

  dc.SetPrefixPath(output_dir);
  dc.SetDataFormat(mfem::VTKFormat::BINARY);
  dc.SetHighOrderOutput(true);
  dc.SetLevelsOfDetail(1);
  dc.SetTime(0.0);
  dc.SetCycle(0);

  dc.Save();
}

// --------- C-style helper wrapper ---------

int RunMagnetostaticCase(const std::string &json_file,
                         MPI_Comm comm,
                         bool verbose)
{
  try
  {
    MagnetostaticCase problem(json_file, comm, verbose);
    problem.Run();
  }
  catch (const std::exception &ex)
  {
    MPI_Comm world_comm = palace::Mpi::World();
    palace::Mpi::Print(world_comm,
                       "[SphinxSys-Palace] Magnetostatic case FAILED: {}\n",
                       ex.what());
    return 1;
  }
  return 0;
}

}  // namespace sphinxsys_palace
