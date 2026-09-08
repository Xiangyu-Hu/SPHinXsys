# Standalone test for DomainPartition

`DomainPartition` is pure geometry and topology — no MPI, no particles, no Eigen-specific
behaviour — so it can be tested without building SPHinXsys or installing Simbody/TBB/Boost.
This lives outside `src/` on purpose: `src/CMakeLists.txt` globs `shared/*.cpp` and
`shared/*.h` recursively, so a shim `data_type.h` under `src/` would shadow the real one
and the test `main()` would be linked into the library.

`shim/data_type.h` supplies minimal stand-ins for `Real`, `Vecd`, `Arrayi`, `BoundingBoxd`
and `SMAX`/`SMIN`, and `domain_partition.cpp` is compiled unmodified against it.

Run:

    c++ -std=c++17 -Wall -Ishim -I../src/shared/shared_ck/domain_decomposition \
        test_domain_partition.cpp ../src/shared/shared_ck/domain_decomposition/domain_partition.cpp \
        -o test_dp && ./test_dp

Covers: subdomain bounds, halo expansion, rank<->grid-coordinate round trip, ownership
lookup (including the upper-bound clamp and out-of-domain -1), neighbour topology counts
for corner / interior / slab decompositions, and halo target selection for face, edge and
corner particles.

The shim fixes `Dimensions = 3` and `Real = double`. A proper gtest under
`tests/unit_tests_src/shared/` should replace this once a full build environment is
available — this exists so the topology logic is verifiable today.
