#ifndef FISH_AND_BONES_H
#define FISH_AND_BONES_H

#include "sphinxsys.h"

using namespace SPH;
//the fish outline is a 5th order polynomial function from paper 
//DOI:https://doi.org/10.1016/j.jtbi.2016.08.025
class SphFish
{
  public:
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

	Real RevertOutline(Real x, Real h, Real L)
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
			y += a[n] * pow(0.2 - x, n + 1);
		}
		return y;
	}

	std::vector<Vecd> CreatFishShape(Real center_x, Real center_y, Real length, Real thickness, Real resolution)
	{
		Real headtip = center_x;		//head position(cx, cy)
		Real tailtip = center_x + length;	//tail position(cx + length, cy)

		Real h = thickness;				/**< The maximum fish thickness. */
		Real L = tailtip - headtip;	/**< The fish length. */

		int Nh = 100;
		Real Lstep = L / Nh;

		int start_index = 0;
		int end_index = Nh;

		std::vector<Vecd> pnts;
		for (int n = 0; n <= Nh; n++)
		{
			// Real t = L - n * Lstep;
			Real t =n * Lstep;
			Real x = tailtip - t;
			// Real x = headtip + t;
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
};

class SphFlap
{
  public:
    Real flapoutline(Real x, Real cord, Real thickness)
	{
		Real ratio = thickness / cord;

		Real a[5];
		a[0] = 0.2969 * pow(x / cord, 0.5);
		a[1] = -0.126 * pow(x / cord, 1.0);
		a[2] = -0.3516 * pow(x / cord, 2.0);
		a[3] = 0.2843 * pow(x / cord, 3.0);
		a[4] = -0.1015 * pow(x / cord, 4.0);
		Real y = 0.0;
		for (int n = 0; n<5; n++)
		{
			y += 5 * ratio * cord * a[n];
		}
		return y;
	}

	std::vector<Vecd> CreatFlapShape(Real flap_center_x, Real flap_center_y, Real flap_cord, Real flap_thickness, Real flap_resolution)
	{
		Real headtip = flap_center_x;		//head position(cx, cy)
		//Real tailtip = flap_center_x + flap_cord;	//tail position(cx + length, cy)

		//Real h = flap_thickness;				/**< The maximum flap thickness. */
		Real L = flap_cord;                 	/**< The flap length. */

		int Nh = 100;
		Real Lstep = flap_cord / Nh;

		int start_index = 0;
		int end_index = Nh;

		std::vector<Vecd> pnts;
		for (int n = 0; n <= Nh; n++)
		{
			Real t = L - n * Lstep;
			// Real t =n * Lstep;
			// Real x = tailtip - t;
			Real x = headtip + t;
			Real y = flapoutline(t, flap_cord, flap_thickness);
			
			if (y >= flap_resolution)
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
			pnts1.push_back(Vecd(pnts[n][0], pnts[n][1] + flap_center_y));		
		}
		//lower camber line 
		for (int n = N; n >= 0; n--)
		{
			pnts1.push_back(Vecd(pnts1[n][0], -pnts1[n][1] + 2.0*flap_center_y));
		}
		pnts1.push_back(Vecd(pnts[0][0], pnts[0][1] + flap_center_y));

		return pnts1;
	}
};


#endif //FISH_AND_BONES_H