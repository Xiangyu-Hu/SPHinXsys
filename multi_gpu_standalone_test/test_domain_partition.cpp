#include "domain_partition.h"
#include <iostream>
using namespace SPH;

static int failures = 0;
static void check(bool ok, const std::string &what)
{
    if (!ok) { std::cout << "  FAIL: " << what << "\n"; ++failures; }
    else     { std::cout << "  ok:   " << what << "\n"; }
}
static Vecd V(Real x, Real y, Real z) { Vecd v; v[0]=x; v[1]=y; v[2]=z; return v; }
static Arrayi A(int x, int y, int z) { Arrayi a; a[0]=x; a[1]=y; a[2]=z; return a; }

int main()
{
    // Unit cube split 2x2x2 across 8 ranks, halo 0.1
    BoundingBoxd bounds(V(0,0,0), V(1,1,1));
    Real halo = 0.1;

    std::cout << "-- rank 0 of 2x2x2 --\n";
    DomainPartition p0(bounds, A(2,2,2), halo, 0, 8);
    check(p0.SubdomainBounds().lower_[0] == 0.0 && p0.SubdomainBounds().upper_[0] == 0.5,
          "rank 0 owns lower half in x");
    check(p0.HaloedBounds().lower_[0] == -0.1 && p0.HaloedBounds().upper_[0] == 0.6,
          "haloed bounds grown by halo width");
    check(p0.NeighbourRanks().size() == 7, "corner rank has 7 neighbours (of 26 possible)");

    // Ownership: row-major strides (1, nx, nx*ny) => rank = cx + 2*cy + 4*cz
    check(p0.ownerRank(V(0.25,0.25,0.25)) == 0, "centre of subdomain 0 -> rank 0");
    check(p0.ownerRank(V(0.75,0.25,0.25)) == 1, "+x -> rank 1");
    check(p0.ownerRank(V(0.25,0.75,0.25)) == 2, "+y -> rank 2");
    check(p0.ownerRank(V(0.25,0.25,0.75)) == 4, "+z -> rank 4");
    check(p0.ownerRank(V(0.75,0.75,0.75)) == 7, "+x+y+z -> rank 7");
    check(p0.ownerRank(V(1.0,1.0,1.0)) == 7, "upper corner clamps into last subdomain");
    check(p0.ownerRank(V(-0.5,0.5,0.5)) == -1, "outside system bounds -> -1");
    check(p0.isOwnedByThisRank(V(0.25,0.25,0.25)), "rank 0 owns its own centre");

    // Halo targets: a particle deep inside rank 0 is needed by nobody;
    // one within `halo` of the x-face is needed by the +x neighbour (rank 1).
    StdVec<int> targets;
    p0.haloTargetRanks(V(0.25,0.25,0.25), targets);
    check(targets.empty(), "interior particle sent to nobody");
    check(!p0.isHaloCandidate(V(0.25,0.25,0.25)), "interior particle is not a halo candidate");

    targets.clear();
    p0.haloTargetRanks(V(0.45,0.25,0.25), targets);
    check(targets.size() == 1 && targets[0] == 1, "near +x face -> sent to rank 1 only");
    check(p0.isHaloCandidate(V(0.45,0.25,0.25)), "near-face particle is a halo candidate");

    // Near the +x+y edge: needed by rank 1 (+x), rank 2 (+y) and rank 3 (+x+y)
    targets.clear();
    p0.haloTargetRanks(V(0.45,0.45,0.25), targets);
    check(targets.size() == 3, "near +x+y edge -> 3 target ranks");

    // Corner of the domain interior: all 7 neighbours need it
    targets.clear();
    p0.haloTargetRanks(V(0.45,0.45,0.45), targets);
    check(targets.size() == 7, "near +x+y+z corner -> all 7 neighbours");

    std::cout << "-- rank 13 of 3x3x3 (fully interior) --\n";
    DomainPartition p13(bounds, A(3,3,3), 0.05, 13, 27);
    check(p13.NeighbourRanks().size() == 26, "interior rank has all 26 neighbours");
    check(p13.ownerRank(V(0.5,0.5,0.5)) == 13, "centre belongs to rank 13");
    check(p13.isOwnedByThisRank(V(0.5,0.5,0.5)), "rank 13 owns the centre");

    std::cout << "-- 1D split 4x1x1 --\n";
    DomainPartition p1(bounds, A(4,1,1), 0.1, 1, 4);
    check(p1.NeighbourRanks().size() == 2, "middle slab has 2 neighbours");
    check(p1.ownerRank(V(0.3,0.5,0.5)) == 1, "x=0.3 -> rank 1");
    check(p1.ownerRank(V(0.9,0.5,0.5)) == 3, "x=0.9 -> rank 3");

    std::cout << (failures == 0 ? "\nALL PASSED\n" : "\nFAILURES\n");
    return failures;
}
