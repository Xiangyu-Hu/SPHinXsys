# ![SPHinXsys Logo](assets/logo.png) SPHinXsys

**Project Status**  
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Linux](https://img.shields.io/badge/os-Linux-green.svg)](https://shields.io/)
[![Windows](https://img.shields.io/badge/os-Windows-green.svg)](https://shields.io/)
[![macOS](https://img.shields.io/badge/os-macOs-green.svg)](https://shields.io/)
![ci workflow](https://github.com/Xiangyu-Hu/SPHinXsys/actions/workflows/ci.yml/badge.svg?event=push)

**Project Communication**  
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/sphinxsys.svg?style=social&label=Follow%20%40sphinxsys)](https://twitter.com/sphinxsys)
[![YouTube](https://img.shields.io/badge/YouTube-FF0000.svg?style=flat&logo=YouTube&logoColor=white)](https://www.youtube.com/channel/UCexdJbxOn9dvim6Jg1dnCFQ)
[![Bilibili](https://img.shields.io/badge/bilibili-%E5%93%94%E5%93%A9%E5%93%94%E5%93%A9-critical)](https://space.bilibili.com/1761273682/video)

## Description

SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle Hydrodynamics for industrial compleX systems.
It provides C++ APIs for physical accurate simulation and aims to model coupled industrial dynamic systems including fluid, solid, multi-body dynamics
and beyond with SPH (smoothed particle hydrodynamics), a meshless computational method using particle discretization.
For more information on the SPHinXsys project, please check the project website: <https://www.sphinxsys.org>.

SPHinXsys is a multi-physics, multi-resolution SPH library based on a unified computational framework in which strong couplings are achieved for all involved physics.
Although it is not a standalone application itself,
many examples designated for the specific type of applications are provided.

## Examples at a glance

Using SPHinXsys library, straightforward and fast multi-physics modeling can be achieved.
Here, we present several short examples in flow, solid dynamics, fluid structure interactions (FSI) and dynamic solid contact.

<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/2d_examples/test_2d_dambreak/Dambreak.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/dambreak.gif" height="192px"></a>
<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/2d_examples/test_2d_fsi2/fsi2.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/fsi-2d.gif" height="192px"></a>
<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/3d_examples/test_3d_elasticSolid_shell_collision/3d_elasticSolid_shell_collision.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/ball-shell.gif" height="192px"></a>
<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/3d_examples/test_3d_shell_stability_half_sphere/test_3d_shell_stability_half_sphere.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/half-sphere.gif" height="192px"></a>
<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/3d_examples/test_3d_twisting_column/twisting_column.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/twisting.gif" height="168px"></a>
<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/user_examples/test_2d_flow_stream_around_fish/2d_flow_stream_around_fish.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/fish-swimming.gif" height="168px"></a>

### Fully compatible to classical FVM method

Through the unified computational framework in SPHinXsys,
the algorithms for particle methods are full compatible to those in the classical finite volume method (FVM).
The following gives an example of the flow around cylinder problem solved by FVM in SPHinXsys.

<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/user_examples/test_2d_FVM_flow_around_cylinder/2d_FVM_flow_around_cylinder.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/fvm-sphinxsys-flow-around-cylinder.gif" height="168px"></a>

Note that the code for FVM algorithm is exact the same one for particle interaction in SPHinXsys.
The only difference is that SPHinXsys reads a predefined mesh, other than generate particles, before the computation.

## Journal publications

Main Reference:

1. C. Zhang, M. Rezavand, Y. Zhu, Y. Yu, D. Wu, W. Zhang, J. Wang, X. Hu,
"SPHinXsys: an open-source multi-physics and multi-resolution library based on smoothed particle hydrodynamics",
Computer Physics Communications, 267, 108066, 2021.  
[![Main Publication](https://img.shields.io/badge/doi-10.1016%2Fj.cpc.2021.108066-d45815.svg)](https://doi.org/10.1016/j.cpc.2021.108066)

The algorithms in SPHinXsys are based on the following publications:

1. Chi Zhang and Yujie Zhu and Dong Wu and Nikolaus A Adams and Xiagnyu Hu,
"Smoothed particle hydrodynamics: Methodology development and recent achievement",
Journal of Hydrodynamics 34(5), 767--805, 2022  
[![Review Publication](https://img.shields.io/badge/doi-10.1007%2Fs42241.022.0052.1-d45815.svg)](https://doi.org/10.1007/s42241-022-0052-1)

2. Luhui Han and Xiangyu Hu,
"SPH modeling of fluid-structure interaction",
Journal of Hydrodynamics, 2018: 30(1):62-69.  
[![FSI Publication](https://img.shields.io/badge/doi-10.1007%2Fs42241.018.0006.9-d45815.svg)](https://doi.org/10.1007/s42241-018-0006-9)

3. Chi Zhang and Massoud Rezavand and Xiangyu Hu,
"Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics",
Journal of Computational Physics 404 (2020) 109135  
[![Dual Time Criteria](https://img.shields.io/badge/doi-10.1016%2Fj.jcp.2019.109135-d45815.svg)](https://doi.org/10.1016/j.jcp.2019.109135)

4. Chi Zhang et al.
"SPHinXsys: An open-source meshless, multi-resolution and multi-physics library",
Software Impacts, 6 (2020) 100033  
[![Software Impact](https://img.shields.io/badge/doi-10.1016%2Fj.simpa.2020.100033-d45815.svg)](https://doi.org/10.1016/j.simpa.2020.100033)

5. Chi Zhang, Massoud Rezavand, Xiangyu Hu,
"A multi-resolution SPH method for fluid-structure interactions",
Journal of Computational Physics,  Journal of Computational Physics, 429 (2021) 110028.  
[![Software Impact](https://img.shields.io/badge/doi-10.1016%2Fj.jcp.2020.110028-d45815.svg)](https://doi.org/10.1016/j.jcp.2020.110028)

6. Chi Zhang, Yanji Wei, Frederic Dias, Xiangyu Hu,
"An efficient fully Lagrangian solver for modeling wave interaction with oscillating wave energy converter",
Ocean Engineering,
Volume 236, (2021) 109540  
[![OWSC](https://img.shields.io/badge/doi-10.1016%2Fj.oceaneng.2021.109540-d45815.svg)](https://doi.org/10.1016/j.oceaneng.2021.109540)

7. Chi Zhang, Jianhang Wang, Massoud Rezavand, Dong Wu, Xiangyu Hu,
"An integrative smoothed particle hydrodynamics framework for modeling cardiac function",
 Computer Methods in Applied Mechanics and Engineering, 381, 113847, 2021.  
 [![Cardiac Function](https://img.shields.io/badge/doi-10.1016%2Fj.cma.2021.113847-d45815.svg)](https://doi.org/10.1016/j.cma.2021.113847)

8. Yujie Zhu, Chi Zhang, Yongchuan Yu, Xiangyu Hu, "A CAD-compatible body-fitted particle generator for arbitrarily complex geometry and its application to wave-structure interaction", Journal of Hydrodynamics, 33(2), 195-206, 2021.  
[![CAD Particle](https://img.shields.io/badge/doi-10.1007%2Fs42241.021.0031.y-d45815.svg)](https://doi.org/10.1007/s42241-021-0031-y)

9. Chi Zhang, Yujie Zhu, Xiuxiu Lyu, Xiangyu Hu, "An efficient and generalized solid boundary condition for SPH: Applications to multi-phase flow and fluid–structure interaction", European Journal of Mechanics - B/Fluids, 94, 276-292, 2022.  
 [![Boundary FSI](https://img.shields.io/badge/doi-10.1016%2Fj.euromechflu.2022.03.011-d45815.svg)](https://doi.org/10.1016/j.euromechflu.2022.03.011)

10. Yujie Zhu, Chi Zhang, Xiangyu Hu, "A dynamic relaxation method with operator splitting and random-choice strategy for SPH", Journal of Computational Physics, 458,
111105, 2022.  
 [![Dynamic Relaxation](https://img.shields.io/badge/doi-10.1016%2Fj.jcp.2022.111105-d45815.svg)](https://doi.org/10.1016/j.jcp.2022.111105)

11. Dong Wu, Chi Zhang, Xiaojing Tang, Xiangyu Hu, "An essentially non-hourglass formulation for total Lagrangian smoothed particle hydrodynamics", Computer Methods in Applied Mechanics and Engineering, 407, 115915, 2023.  
[![Hourglass](https://img.shields.io/badge/doi-10.1016%2Fj.jcp.2022.111105-d45815.svg)](https://doi.org/10.1016/j.jcp.2022.111105)

12. Chi Zhang, Hao Gao, Xiangyu Hu, "A multi-order smoothed particle hydrodynamics method for cardiac electromechanics with the Purkinje network", Computer Methods in Applied Mechanics and Engineering, 407, 115885, 2023.  
[![Purkinje Network](https://img.shields.io/badge/doi-10.1016%2Fj.cma.2023.115885-d45815.svg)](https://doi.org/10.1016/j.cma.2023.115885)

13. Yongchuan Yu, Yujie Zhu, Chi Zhang, Oskar J. Haidn, Xiangyu Hu, "Level-set based pre-processing techniques for particle method", Computer Physics Communications 289, 108744, 2023.  
[![Level Set](https://img.shields.io/badge/doi-10.1016%2Fj.cpc.2023.108744-d45815.svg)](https://doi.org/10.1016/j.cpc.2023.108744)

14. Shuoguo Zhang, Wenbin Zhang, Chi Zhang, Xiangyu Hu, "A Lagrangian free-stream boundary condition for weakly compressible smoothed particle hydrodynamics", Journal of Computational Physics, 490, 112303, 2023.  
[![Free Stream](https://img.shields.io/badge/doi-10.1016%2Fj.jcp.2023.112303-d45815.svg)](https://doi.org/10.1016/j.jcp.2023.112303)

## Software Architecture

SPHinXsys is cross-platform can be compiled and used in Windows, Linux and McOS systems.

## Installation, tutorial and documentation

For installation, program manual and tutorials, please check <https://www.sphinxsys.org/html/sphinx_index.html>.
Please check the documentation of the code at <https://xiangyu-hu.github.io/SPHinXsys/>.

## Get involved to SPHinXsys

You are welcomed to use and get involved in SPHinXsys.

As the code is on git-hub, you can register an account there (if you do not have a github account yet)
and fork out the SPHinXsys repository.
You can work on the forked repository and add new features, and then commit them.

Besides forking the repository and begin to develop by your own,
there are many other ways to make SPHinXsys better for every one.
For example, you can initiate issues on any thing relevant to SPHinXsys, not only bugs and installation issues,
check the pull requests for the current development status,
ask questions, give suggestion and comment on SPHinXsys on the discussion page, etc.

You are also welcomed to join the main repository as a collaborator,
by which you are able to branch directly in the main repository,
and review the pull request.

If you have any further question, you are also welcomed to contact <xiangyu.hu@tum.de>.
