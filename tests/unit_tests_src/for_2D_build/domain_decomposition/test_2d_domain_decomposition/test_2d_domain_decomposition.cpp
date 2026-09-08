/**
 * @file 	test_2d_domain_decomposition.cpp
 * @brief 	Unit tests of the decomposition geometry and of the subdomain fan-out.
 * @details These cover the parts of the domain decomposition that are independent of
 *          the backend, and that are therefore identical between a multi-GPU run and
 *          the host side debugging path: which subdomain owns a position, which
 *          positions belong to a halo band, and the semantics of the fan-out.
 *
 *          A failure here is a decomposition bug. It is worth keeping these cheap and
 *          exhaustive, because the same logic is far harder to observe once it runs on
 *          eight GPUs.
 * @author	Niki Loppi
 */
#include "domain_decomposition.h"
#include "subdomain_fan_out.h"

#include <gtest/gtest.h>

#include <functional>
#include <mutex>
#include <set>
#include <stdexcept>

using namespace SPH;
using namespace SPH::execution;

namespace
{
SlabDecomposition makeDecomposition()
{ // four subdomains over x in [0, 8), cut planes at 0, 2, 4, 6, 8
    return SlabDecomposition(BoundingBoxd(Vecd::Zero(), Vecd(8.0, 2.0)), 0.5, 4);
}
} // namespace

TEST(DomainDecomposition, SplitsAlongTheLongestAxis)
{
    SlabDecomposition decomposition = makeDecomposition();
    EXPECT_EQ(decomposition.SplitAxis(), 0);
    EXPECT_EQ(decomposition.NumberOfSubdomains(), 4);
    EXPECT_DOUBLE_EQ(decomposition.CutPlane(0), 0.0);
    EXPECT_DOUBLE_EQ(decomposition.CutPlane(4), 8.0);
}

TEST(DomainDecomposition, OwnershipIsAPartition)
{
    SlabDecomposition decomposition = makeDecomposition();
    const SubdomainMap &map = decomposition.getSubdomainMap();

    // Every position in the domain is owned by exactly one subdomain. A gap would
    // silently drop particles; an overlap would duplicate them.
    for (int i = 0; i < 8000; ++i)
    {
        const Real x = 8.0 * Real(i) / 8000.0;
        int containing = 0;
        for (int subdomain = 0; subdomain < 4; ++subdomain)
        {
            if (x >= map.cut_plane_[subdomain] && x < map.cut_plane_[subdomain + 1])
            {
                ++containing;
            }
        }
        EXPECT_EQ(containing, 1) << "at x = " << x;
        const int owner = map.subdomainOf(Vecd(x, 0.0));
        EXPECT_GE(owner, 0);
        EXPECT_LT(owner, 4);
    }
}

TEST(DomainDecomposition, PositionsOutsideTheDomainAreStillOwned)
{ // a particle that leaves the domain must never become unowned, or migration loses it
    SlabDecomposition decomposition = makeDecomposition();
    const SubdomainMap &map = decomposition.getSubdomainMap();
    EXPECT_EQ(map.subdomainOf(Vecd(-5.0, 0.0)), 0);
    EXPECT_EQ(map.subdomainOf(Vecd(99.0, 0.0)), 3);
}

TEST(DomainDecomposition, HaloBandIsTheSlabWidenedByTheCutOff)
{
    SlabDecomposition decomposition = makeDecomposition();
    const SubdomainMap &map = decomposition.getSubdomainMap();
    // subdomain 1 spans [2, 4), so its halo band is [1.5, 4.5)
    EXPECT_TRUE(map.inHaloBandOf(Vecd(1.6, 0.0), 1));
    EXPECT_FALSE(map.inHaloBandOf(Vecd(1.4, 0.0), 1));
    EXPECT_TRUE(map.inHaloBandOf(Vecd(4.4, 0.0), 1));
    EXPECT_FALSE(map.inHaloBandOf(Vecd(4.6, 0.0), 1));

    // A particle in a neighbor's halo band is still owned by us. That is what makes a
    // halo particle a copy rather than a transfer, and is the invariant that keeps
    // migration and halo exchange from fighting over the same particle.
    EXPECT_EQ(map.subdomainOf(Vecd(1.7, 0.0)), 0);
    EXPECT_TRUE(map.inHaloBandOf(Vecd(1.7, 0.0), 1));
}

TEST(DomainDecomposition, EndSubdomainsHaveOneNeighbor)
{
    SlabDecomposition decomposition = makeDecomposition();
    const SubdomainMap &map = decomposition.getSubdomainMap();
    EXPECT_EQ(map.neighborOf(0, 0), -1);
    EXPECT_EQ(map.neighborOf(0, 1), 1);
    EXPECT_EQ(map.neighborOf(3, 0), 2);
    EXPECT_EQ(map.neighborOf(3, 1), -1);
    EXPECT_EQ(SubdomainMap::oppositeSide(0), 1);
    EXPECT_EQ(SubdomainMap::oppositeSide(1), 0);
}

TEST(DomainDecomposition, RebalanceConvergesAndStaysMonotone)
{
    SlabDecomposition decomposition = makeDecomposition();
    const SubdomainMap &map = decomposition.getSubdomainMap();

    StdVec<UnsignedInt> load = {8000, 100, 100, 100}; // heavily skewed to the left
    for (int iteration = 0; iteration < 40; ++iteration)
    {
        decomposition.rebalance(load, 0.5);
        for (int subdomain = 0; subdomain < 4; ++subdomain)
        {
            EXPECT_LT(map.cut_plane_[subdomain], map.cut_plane_[subdomain + 1]);
        }
        EXPECT_DOUBLE_EQ(map.cut_plane_[0], 0.0); // outer planes stay pinned
        EXPECT_DOUBLE_EQ(map.cut_plane_[4], 8.0);

        // Feed back the load a uniform particle density would give for the new planes.
        for (int subdomain = 0; subdomain < 4; ++subdomain)
        {
            load[subdomain] =
                UnsignedInt(1000.0 * (map.cut_plane_[subdomain + 1] - map.cut_plane_[subdomain]));
        }
    }

    UnsignedInt lightest = load[0];
    UnsignedInt heaviest = load[0];
    for (int subdomain = 0; subdomain < 4; ++subdomain)
    {
        lightest = std::min(lightest, load[subdomain]);
        heaviest = std::max(heaviest, load[subdomain]);
    }
    EXPECT_LT(Real(heaviest) / Real(lightest), 1.05);
}

TEST(DomainDecomposition, RebalanceRespectsTheMinimumSlabThickness)
{ // below twice the halo width a halo would reach past the adjacent subdomain
    SlabDecomposition decomposition = makeDecomposition();
    const SubdomainMap &map = decomposition.getSubdomainMap();
    StdVec<UnsignedInt> degenerate_load = {100000, 1, 1, 1};
    for (int iteration = 0; iteration < 100; ++iteration)
    {
        decomposition.rebalance(degenerate_load, 0.9);
    }
    for (int subdomain = 0; subdomain < 4; ++subdomain)
    {
        const Real thickness = map.cut_plane_[subdomain + 1] - map.cut_plane_[subdomain];
        EXPECT_GE(thickness, decomposition.MinimumSlabThickness() - 1.0e-9);
    }
}

//----------------------------------------------------------------------
// Fan-out semantics. These are what the CK algorithms rely on when their
// exec() body is wrapped, so they are worth pinning down independently.
//----------------------------------------------------------------------
class SubdomainFanOut : public ::testing::TestWithParam<SubdomainRunner::Mode>
{
  protected:
    void SetUp() override { subdomain_runner.initialize(4, GetParam()); }
    void TearDown() override { subdomain_runner.initialize(1, SubdomainRunner::Mode::Sequential); }
};

TEST_P(SubdomainFanOut, RunsTheBodyOncePerSubdomain)
{
    std::mutex mutex;
    std::set<int> visited;
    fanOutOverSubdomains(par_multi_host,
                         [&]()
                         {
                             std::lock_guard<std::mutex> lock(mutex);
                             visited.insert(currentSubdomainID());
                         });
    EXPECT_EQ(visited.size(), 4u);
    EXPECT_EQ(currentSubdomainID(), 0); // binding restored on the calling thread
}

TEST_P(SubdomainFanOut, NestedFanOutStaysOnTheBoundSubdomain)
{ // dynamics compose: an interaction fans out, then calls exec() of its post-processes
    std::mutex mutex;
    std::set<int> visited;
    fanOutOverSubdomains(par_multi_host,
                         [&]()
                         {
                             fanOutOverSubdomains(par_multi_host,
                                                  [&]()
                                                  {
                                                      std::lock_guard<std::mutex> lock(mutex);
                                                      visited.insert(currentSubdomainID());
                                                  });
                         });
    EXPECT_EQ(visited.size(), 4u); // and not 16 visits, nor a deadlock
}

TEST_P(SubdomainFanOut, ReductionCombinesThePartialResults)
{
    const int total = reduceOverSubdomains<std::plus<int>>(
        par_multi_host, 0, [&]() { return currentSubdomainID() + 1; });
    EXPECT_EQ(total, 1 + 2 + 3 + 4);
}

TEST_P(SubdomainFanOut, ExceptionsPropagateToTheCaller)
{
    EXPECT_THROW(
        fanOutOverSubdomains(par_multi_host,
                             [&]()
                             {
                                 if (currentSubdomainID() == 2)
                                 {
                                     throw std::runtime_error("failure inside a subdomain");
                                 }
                             }),
        std::runtime_error);
}

TEST(SubdomainFanOutSingle, NonDecomposedPolicyRunsTheBodyExactlyOnce)
{
    int calls = 0;
    fanOutOverSubdomains(par_host, [&]() { ++calls; });
    EXPECT_EQ(calls, 1);
}

INSTANTIATE_TEST_SUITE_P(BothRunnerModes, SubdomainFanOut,
                         ::testing::Values(SubdomainRunner::Mode::Sequential,
                                           SubdomainRunner::Mode::Threaded));

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
