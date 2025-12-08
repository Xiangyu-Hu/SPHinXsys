// modules/Palace/palace_interface.cpp

#include "palace_interface.hpp"

#include <cstring>
#include <memory>
#include <sstream>
#include <vector>
#include <iostream>

#include <mfem.hpp>

// Palace headers
#include "palace/drivers/magnetostaticsolver.hpp"
#include "palace/fem/errorindicator.hpp"
#include "palace/fem/libceed/ceed.hpp"
#include "palace/fem/mesh.hpp"
#include "palace/linalg/hypre.hpp"
#include "palace/linalg/slepc.hpp"
#include "palace/utils/communication.hpp"
#include "palace/utils/device.hpp"
#include "palace/utils/geodata.hpp"
#include "palace/utils/iodata.hpp"
#include "palace/utils/omp.hpp"
#include "palace/utils/outputdir.hpp"
#include "palace/utils/timer.hpp"

using namespace palace;

namespace sphinxsys_palace
{

// -----------------------------------------------------------------------------
// Helper functions derived from Palace main.cpp
// -----------------------------------------------------------------------------

// Return a string describing the embedded Palace version/tag.
// You can adjust this to match your Palace fork or git tag.
static const char *GetPalaceGitTag()
{
  return "sphinxsys-embedded-palace";
}


// Configure MFEM device string based on the requested Palace Device enum.
// This decides whether to use CPU, GPU (CUDA/HIP), debug device, and OpenMP.
static std::string ConfigureDevice(Device device)
{
  std::string device_str;
  switch (device)
  {
    case Device::CPU:
      device_str = "cpu";
      break;

    case Device::GPU:
    // GPU execution: choose between CUDA or HIP depending on how MFEM was built.
#if defined(MFEM_USE_CUDA)
      device_str = "cuda";
#elif defined(MFEM_USE_HIP)
      device_str = "hip";
#else
      // If MFEM was not built with GPU support, fall back to CPU.
      Mpi::Warning("Palace was built without CUDA/HIP support. Falling back to CPU.\n");
      device_str = "cpu";
#endif
      break;

    case Device::DEBUG:
      // Use MFEM's debug device (runs kernels on CPU but checks extra things).
      device_str = "cpu,debug";
      break;
  }

#if defined(MFEM_USE_OPENMP)
// Append OpenMP if MFEM was compiled with OpenMP support.
  device_str += ",omp";
#endif

  return device_str;
}

// Configure libCEED backend (e.g. /gpu/cuda/magma, /cpu/self, etc.)
// based on MFEM device capabilities and optional user override.
static void ConfigureCeedBackend(const std::string &ceed_backend)
{
  std::string default_ceed_backend;

  // Choose a sensible default backend depending on which MFEM backend is available.
  if (mfem::Device::Allows(mfem::Backend::CUDA_MASK))
    default_ceed_backend = "/gpu/cuda/magma";
  else if (mfem::Device::Allows(mfem::Backend::HIP_MASK))
    default_ceed_backend = "/gpu/hip/magma";
  else if (mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE))
    default_ceed_backend = "/cpu/self/ref/serial";
  else
    default_ceed_backend = "/cpu/self";

  // If user specified a backend in the JSON, prefer that; otherwise use default.
  const std::string &backend =
      (!ceed_backend.empty() ? ceed_backend : default_ceed_backend);


  // Initialize libCEED with the chosen backend.
  ceed::Initialize(backend.c_str(), "");


  // Check the actual backend that libCEED ended up using.
  std::string actual = ceed::Print();
  if (backend.compare(0, backend.length(), actual, 0, backend.length()))
  {
    Mpi::Warning("libCEED backend mismatch: requested \"{}\", but got \"{}\"\n",
                 backend, actual);
  }
}

// Print Palace ASCII banner on root process.
static void PrintPalaceBanner(MPI_Comm comm)
{
  Mpi::Print(comm,
             "_____________     _______\n"
             "_____   __   \\____ __   /____ ____________\n"
             "____   /_/  /  __ ` /  /  __ ` /  ___/  _ \\\n"
             "___   _____/  /_/  /  /  /_/  /  /__/  ___/\n"
             "  /__/     \\___,__/__/\\___,__/\\_____\\_____/\n\n");
}

// Print basic run-time information: git tag, MPI size, OpenMP threads,
// GPU devices, MFEM device info, and libCEED backend.
static void PrintPalaceInfo(MPI_Comm comm, int np, int nt,
                            int ngpu, mfem::Device &device)
{
  Mpi::Print(comm, "Palace Git tag: {}\n", GetPalaceGitTag());
  Mpi::Print(comm, "Running with {} MPI process{}\n",
             np, (np > 1 ? "es" : ""));

  if (nt > 0)
    Mpi::Print(comm, "OpenMP threads: {}\n", nt);

#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
  const char *gpu_name =
#if defined(MFEM_USE_CUDA)
      "CUDA";
#else
      "HIP";
#endif

  Mpi::Print(comm, "Detected {} {} device{}\n",
             ngpu, gpu_name, (ngpu != 1 ? "s" : ""));
#endif

  // Print MFEM device and libCEED backend in one line.
  std::ostringstream os;
  device.Print(os);
  os << "libCEED backend: " << ceed::Print();
  Mpi::Print(comm, "{}\n", os.str());
}


// -----------------------------------------------------------------------------
// Public interface used by SPHinXsys
//   This is the only function that SPHinXsys needs to call in order to run
//   the Palace magnetostatic case (e.g. rings.json).
// -----------------------------------------------------------------------------
int RunMagnetostaticCase(const std::string &config_file, bool verbose)
{
 
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0, size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  bool root = (rank == 0);

  // Optional user-facing message to confirm that Palace is being run.
  if (verbose && root)
    std::cout << "\nRunning Palace magnetostatic solver...\n";

  // Global initialization timer (Palace-style).
  BlockTimer bt(Timer::INIT);

  // Print basic info and Palace banner on root.
  if (verbose && root)
  {
    std::cout << "Using configuration file: " << config_file << "\n\n";
    PrintPalaceBanner(comm);
    std::cout << "Palace Git tag: " << GetPalaceGitTag() << "\n";
  }


  // ---------------- Configuration / output directory ----------------

  // Read all Palace input data (JSON config) without printing its own header.
  IoData iodata(config_file.c_str(), /*print=*/false);

  // Create output directory (e.g. results/, mesh/, etc.) based on IoData.
  MakeOutputFolder(iodata, comm);


  // ---------------- Parallelism & device configuration ----------------

  // Configure OpenMP threads (returns actual number of threads used).
  int nt   = utils::ConfigureOmp();
  int ngpu = utils::GetDeviceCount();

  // Configure MFEM device (CPU/GPU/DEBUG, with device id chosen from MPI rank).
  mfem::Device device(ConfigureDevice(iodata.solver.device),
                      utils::GetDeviceId(comm, ngpu));

  // Configure libCEED backend (GPU or CPU) based on solver settings.
  ConfigureCeedBackend(iodata.solver.ceed_backend);
#if defined(PALACE_WITH_GPU_AWARE_MPI)
 // If Palace was compiled with GPU-aware MPI, inform MFEM device.
  device.SetGPUAwareMPI(true);
#endif


  // Initialize Hypre (parallel linear algebra backend).
  hypre::Initialize();


  if (verbose && root)
    PrintPalaceInfo(comm, size, nt, ngpu, device);

  // ---------------- Solver and mesh setup ----------------

  // Create a magnetostatic solver instance (derived from BaseSolver).
  // Arguments: input data, root flag, MPI size, OpenMP threads, git tag string.
  std::unique_ptr<BaseSolver> solver =
      std::make_unique<MagnetostaticSolver>(
          iodata, root, size, nt, GetPalaceGitTag());

  // Build Palace Mesh object from MFEM ParMesh.
  std::vector<std::unique_ptr<Mesh>> mesh;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> pm;
    pm.push_back(mesh::ReadMesh(iodata, comm));

    iodata.NondimensionalizeInputs(*pm[0]);

    // Perform optional mesh refinement steps (uniform / adaptive)
    // as configured in the JSON.
    mesh::RefineMesh(iodata, pm);

    // Convert MFEM ParMesh objects to Palace Mesh objects.
    for (auto &m : pm)
      mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }

  // Main solve loop with (optional) estimate → mark → refine (AMR).
  solver->SolveEstimateMarkRefine(mesh);

  BlockTimer::Print(comm);
  solver->SaveMetadata(BlockTimer::GlobalTimer());


  // Note:
  //   Finalization of CEED, Hypre, SLEPc is intentionally commented out here.
  //   The outer application (SPHinXsys) may still need these libraries after
  //   Palace is done, so we avoid shutting them down from inside this function.
  // ceed::Finalize();
  // hypre::Finalize();
  // slepc::Finalize();

  return 0;
}

}  // namespace sphinxsys_palace
