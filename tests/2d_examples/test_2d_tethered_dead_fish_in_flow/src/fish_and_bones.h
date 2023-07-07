#ifndef FISH_AND_BONES_H
#define FISH_AND_BONES_H

#include "sphinxsys.h"

using namespace SPH;
// the fish outline is a 5th order polynomial function from paper
// DOI:https://doi.org/10.1016/j.jtbi.2016.08.025
Real outline(Real x, Real h, Real L)
{
    Real a[5];
    a[0] = 1.22 * h / L;
    a[1] = 3.19 * h / L / L;
    a[2] = -15.73 * h / pow(L, 3);
    a[3] = 21.87 * h / pow(L, 4);
    a[4] = -10.55 * h / pow(L, 5);
    Real y = 0.0;
    for (int n = 0; n < 5; n++)
    {
        y += a[n] * pow(x, n + 1);
    }
    return y;
}

// create bone shape
std::vector<Vecd> CreatFishShape(Real center_x, Real center_y, Real length, Real resolution)
{
    Real headtip = center_x;          // head position(cx, cy)
    Real tailtip = center_x + length; // tail position(cx + length, cy)

    Real h = 0.4;               /**< The maximum fish thickness. */
    Real L = tailtip - headtip; /**< The fish length. */

    int Nh = 100;
    Real Lstep = L / Nh;

    int start_index = 0;
    int end_index = Nh;

    std::vector<Vecd> pnts;
    for (int n = 0; n <= Nh; n++)
    {
        Real t = L - n * Lstep;
        Real x = tailtip - t;
        Real y = outline(t, h, L);

        if (y >= resolution)
        {
            pnts.push_back(Vecd(x, y));
        }
        else
        {
            if (x < 0.0)
            {
                start_index++;
            }
            else
            {
                end_index--;
            }
        }
    }

    // upper camber line
    std::vector<Vecd> pnts1;
    for (int n = 0; n <= end_index - start_index; n++)
    {
        pnts1.push_back(Vecd(pnts[n][0], pnts[n][1] + center_y));
    }
    // lower camber line
    for (int n = end_index - start_index; n >= 0; n--)
    {
        pnts1.push_back(Vecd(pnts1[n][0], -pnts1[n][1] + 2.0 * center_y));
    }
    pnts1.push_back(Vecd(pnts[0][0], pnts[0][1] + center_y));

    return pnts1;
}

// outline for
Real outline_zebra_fish(Real x, Real h, Real L)
{
    Real w_h = 0.04 * L;
    Real s_b = w_h;
    Real s_t = 0.95 * L;
    Real w_t = 0.01 * L;
    Real y = 0.0;
    if (x >= 0.0 && x < s_b)
    {
        y = sqrt(2.0 * w_h * x - x * x);
    }
    else
    {
        if (x >= s_t && x <= L)
        {
            y = w_t * (L - x) / (L - s_t);
        }
        else
        {
            y = w_h - (w_h - w_t) * (x - s_b) / (s_t - s_b);
        }
    }
    return y;
}

// create bone shape
std::vector<Vecd> CreatZebraFishShape(Real center_x, Real center_y, Real length)
{
    Real headtip = center_x;          // head position(cx, cy)
    Real tailtip = center_x + length; // tail position(cx + length, cy)

    Real h = 0.5;               // the maximum fish thickness
    Real L = tailtip - headtip; // the fish length

    int Nh = 100;
    Real Lstep = L / Nh;

    int start_index = 0;
    int end_index = Nh;

    std::vector<Vecd> pnts;
    for (int n = 0; n <= Nh; n++)
    {
        Real t = n * Lstep;
        Real x = headtip + t;
        Real y = outline_zebra_fish(t, h, L);

        if (y >= 1.0e-3)
            pnts.push_back(Vecd(x, y));
        else
        {
            if (x < 0.0)
            {
                start_index++;
            }
            else
            {
                end_index--;
            }
        }
    }

    // upper camber line
    std::vector<Vecd> pnts1;
    for (int n = 0; n <= end_index - start_index; n++)
    {
        pnts1.push_back(Vecd(pnts[n][0], pnts[n][1] + center_y));
    }
    // lower camber line
    for (int n = end_index - start_index; n >= 0; n--)
    {
        pnts1.push_back(Vecd(pnts1[n][0], -pnts1[n][1] + 2.0 * center_y));
    }
    pnts1.push_back(Vecd(pnts[0][0], pnts[0][1] + center_y));

    return pnts1;
}

// create bone shape
std::vector<Vecd> CreatBoneShape(int Nh, Real ra, Real rab, Real center_x, Real center_y)
{
    Real rb = rab * ra;
    Real headtip = center_x - ra;
    Real hstep = 2.0 * ra / Nh;

    std::vector<Vecd> pnts;

    pnts.push_back(Vecd(center_x - ra, center_y));
    for (int n = 1; n < Nh; n++)
    {
        Real x = headtip + n * hstep;
        pnts.push_back(Vecd(x, rb * sqrt(1 - pow((x - center_x) / ra, 2)) + center_y));
    }
    pnts.push_back(Vecd(center_x + ra, center_y));

    std::vector<Vecd> pnts2;
    for (int n = 0; n <= Nh; n++)
    {
        pnts2.push_back(pnts[n]);
    }
    for (int n = Nh - 1; n > 0; n--)
    {
        pnts2.push_back(Vecd(pnts[n][0], -pnts[n][1] + 2.0 * center_y));
    }
    pnts2.push_back(pnts[0]);

    return pnts2;
}
#endif // FISH_AND_BONES_H