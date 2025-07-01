#include "write_mesh_field_to_plt.h"

namespace SPH
{
//=============================================================================================//
std::vector<Block3D> clusterActiveCells3D(const std::vector<Coord3D> &activeCells)
{
    std::unordered_set<Coord3D, Coord3DHash> activeSet(activeCells.begin(), activeCells.end());
    std::vector<Block3D> blocks;

    while (!activeSet.empty())
    {
        Coord3D start = *activeSet.begin();
        int x0 = start[0], y0 = start[1], z0 = start[2];
        int x1 = x0;

        // Extend in x
        while (activeSet.count({x1 + 1, y0, z0}))
            ++x1;

        // Extend in y
        int y1 = y0;
        bool extendY = true;
        while (extendY)
        {
            for (int x = x0; x <= x1; ++x)
            {
                if (!activeSet.count({x, y1 + 1, z0}))
                {
                    extendY = false;
                    break;
                }
            }
            if (extendY)
                ++y1;
        }

        // Extend in z
        int z1 = z0;
        bool extendZ = true;
        while (extendZ)
        {
            for (int y = y0; y <= y1; ++y)
            {
                for (int x = x0; x <= x1; ++x)
                {
                    if (!activeSet.count({x, y, z1 + 1}))
                    {
                        extendZ = false;
                        break;
                    }
                }
                if (!extendZ)
                    break;
            }
            if (extendZ)
                ++z1;
        }

        // Remove covered cells
        for (int z = z0; z <= z1; ++z)
        {
            for (int y = y0; y <= y1; ++y)
            {
                for (int x = x0; x <= x1; ++x)
                {
                    activeSet.erase({x, y, z});
                }
            }
        }

        blocks.push_back({{x0, y0, z0}, {x1, y1, z1}});
    }

    return blocks;
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//