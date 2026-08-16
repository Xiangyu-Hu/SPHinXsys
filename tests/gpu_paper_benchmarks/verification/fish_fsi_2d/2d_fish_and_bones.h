#ifndef FISH_AND_BONES_H
#define FISH_AND_BONES_H

#include "sphinxsys.h"

using namespace SPH;
//the fish outline is a 5th order polynomial function from paper 
//DOI:https://doi.org/10.1016/j.jtbi.2016.08.025
Real outline(Real x, Real h, Real L)
{
	Real a[5];
	a[0] = 1.22 * h / L;
	a[1] = 3.19 * h / L / L;
	a[2] = -15.73 * h / pow(L, 3);
	a[3] = 21.87 * h / pow(L, 4);
	a[4] = -10.55 * h / pow(L, 5);
	Real y = 0.0;
	for (int n = 0; n<5; n++)
	{
		y += a[n] * pow(x, n + 1);
	}
	return y;
}

std::vector<Vecd> CreatFishShape(Real center_x, Real center_y, Real length, Real resolution)
{
	Real headtip = center_x;		//head position(cx, cy)
	Real tailtip = center_x + length;	//tail position(cx + length, cy)

	Real h = 0.03;				/**< The maximum fish thickness. */
	Real L = tailtip - headtip;	/**< The fish length. */

	int Nh = 100;
	Real Lstep = L / Nh;

	int start_index = 0;
	int end_index = Nh;

	std::vector<Vecd> pnts;
	for (int n = 0; n <= Nh; n++)
	{
		Real t = L - n * Lstep;
		//Real t =n * Lstep;
		//Real x = tailtip - t;
		Real x = headtip + t;
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

	int N = end_index - start_index;
	std::vector<Vecd> pnts1;
	for (int n = 0; n <= N; n++)
	{
		pnts1.push_back(Vecd(pnts[n][0], pnts[n][1] + center_y));		
	}
	//lower camber line 
	for (int n = N; n >= 0; n--)
	{
		pnts1.push_back(Vecd(pnts1[n][0], -pnts1[n][1] + 2.0*center_y));
	}
	pnts1.push_back(Vecd(pnts[0][0], pnts[0][1] + center_y));

	return pnts1;
}


#endif //FISH_AND_BONES_H