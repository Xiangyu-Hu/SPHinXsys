/**
 * Sprint 4.5 Task 4: minimal Contact neighbor variable readback / device sync.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_readback_sync_helpers.h"

#include <cmath>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real probe_value = 42.0;
    const Real host_synced_value = 100.0;
    const Real host_unsynced_value = 200.0;
    const Real read_tol = 1.0e-6;

    AphiTwoBodyInterfaceCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    RegisterContactProbeVariablesCK register_left(case_setup.left_body);
    RegisterContactProbeVariablesCK register_right(case_setup.right_body);
    (void)register_left;
    (void)register_right;
    case_setup.updateRelations();

    const Real device_read = readNeighborProbeAfterDeviceWrite(case_setup, probe_value);
    const Real host_synced_read = readNeighborProbeAfterHostWrite(case_setup, host_synced_value, true);
    const Real host_unsynced_read = readNeighborProbeAfterHostWrite(case_setup, host_unsynced_value, false);

    const bool device_write_ok = std::abs(device_read - probe_value) < read_tol;
    const bool host_sync_ok = std::abs(host_synced_read - host_synced_value) < read_tol;
#if SPHINXSYS_USE_SYCL
    const bool host_unsynced_stale = std::abs(host_unsynced_read - host_synced_value) < read_tol &&
                                     std::abs(host_unsynced_read - host_unsynced_value) > read_tol;
    const bool passed = device_write_ok && host_sync_ok && host_unsynced_stale;
#else
    const bool host_unsynced_stale = true;
    const bool passed = device_write_ok && host_sync_ok;
#endif

    std::cout << "test_3d_contact_variable_readback_sync"
              << " device_read=" << device_read << " host_synced_read=" << host_synced_read
              << " host_unsynced_read=" << host_unsynced_read << " device_write_ok=" << (device_write_ok ? 1 : 0)
              << " host_sync_ok=" << (host_sync_ok ? 1 : 0)
              << " host_unsynced_stale=" << (host_unsynced_stale ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;

    return passed ? 0 : 1;
}
