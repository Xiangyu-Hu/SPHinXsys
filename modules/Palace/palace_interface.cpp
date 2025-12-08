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


static const char *GetPalaceGitTag()
{
  return "sphinxsys-embedded-palace";
}


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
      Mpi::Warning("Palace was built without CUDA/HIP support. Falling back to CPU.\n");
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


static void ConfigureCeedBackend(const std::string &ceed_backend)
{
  std::string default_ceed_backend;

  if (mfem::Device::Allows(mfem::Backend::CUDA_MASK))
    default_ceed_backend = "/gpu/cuda/magma";
  else if (mfem::Device::Allows(mfem::Backend::HIP_MASK))
    default_ceed_backend = "/gpu/hip/magma";
  else if (mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE))
    default_ceed_backend = "/cpu/self/ref/serial";
  else
    default_ceed_backend = "/cpu/self";

  const std::string &backend =
      (!ceed_backend.empty() ? ceed_backend : default_ceed_backend);


  ceed::Initialize(backend.c_str(), "");


  std::string actual = ceed::Print();
  if (backend.compare(0, backend.length(), actual, 0, backend.length()))
  {
    Mpi::Warning("libCEED backend mismatch: requested \"{}\", but got \"{}\"\n",
                 backend, actual);
  }
}


static void PrintPalaceBanner(MPI_Comm comm)
{
  Mpi::Print(comm,
             "_____________     _______\n"
             "_____   __   \\____ __   /____ ____________\n"
             "____   /_/  /  __ ` /  /  __ ` /  ___/  _ \\\n"
             "___   _____/  /_/  /  /  /_/  /  /__/  ___/\n"
             "  /__/     \\___,__/__/\\___,__/\\_____\\_____/\n\n");
}


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

  std::ostringstream os;
  device.Print(os);
  os << "libCEED backend: " << ceed::Print();
  Mpi::Print(comm, "{}\n", os.str());
}


// -----------------------------------------------------------------------------
// Public interface used by SphinxSys
// -----------------------------------------------------------------------------
int RunMagnetostaticCase(const std::string &config_file, bool verbose)
{
 
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0, size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  bool root = (rank == 0);

  if (verbose && root)
    std::cout << "\nRunning Palace magnetostatic solver...\n";

  BlockTimer bt(Timer::INIT);

  if (verbose && root)
  {
    std::cout << "Using configuration file: " << config_file << "\n\n";
    PrintPalaceBanner(comm);
    std::cout << "Palace Git tag: " << GetPalaceGitTag() << "\n";
  }


  IoData iodata(config_file.c_str(), /*print=*/false);
  MakeOutputFolder(iodata, comm);


  int nt   = utils::ConfigureOmp();
  int ngpu = utils::GetDeviceCount();

  mfem::Device device(ConfigureDevice(iodata.solver.device),
                      utils::GetDeviceId(comm, ngpu));

  ConfigureCeedBackend(iodata.solver.ceed_backend);
#if defined(PALACE_WITH_GPU_AWARE_MPI)
  device.SetGPUAwareMPI(true);
#endif


  hypre::Initialize();


  if (verbose && root)
    PrintPalaceInfo(comm, size, nt, ngpu, device);

  // ---- 4. Solver / Mesh ----
  std::unique_ptr<BaseSolver> solver =
      std::make_unique<MagnetostaticSolver>(
          iodata, root, size, nt, GetPalaceGitTag());

  std::vector<std::unique_ptr<Mesh>> mesh;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> pm;
    pm.push_back(mesh::ReadMesh(iodata, comm));

    iodata.NondimensionalizeInputs(*pm[0]);
    mesh::RefineMesh(iodata, pm);

    for (auto &m : pm)
      mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }

  solver->SolveEstimateMarkRefine(mesh);

  BlockTimer::Print(comm);
  solver->SaveMetadata(BlockTimer::GlobalTimer());


  // ceed::Finalize();
  // hypre::Finalize();
  // slepc::Finalize();

  return 0;
}

}  // namespace sphinxsys_palace
