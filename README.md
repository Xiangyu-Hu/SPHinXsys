# ![SPHinXsys Logo](assets/logo.png) SPHinXsys

**Project status**
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Linux](https://img.shields.io/badge/os-Linux-green.svg)](https://shields.io/)
[![Windows](https://img.shields.io/badge/os-Windows-green.svg)](https://shields.io/)
[![macOS](https://img.shields.io/badge/os-macOs-green.svg)](https://shields.io/)
![ci workflow](https://github.com/Xiangyu-Hu/SPHinXsys/actions/workflows/ci.yml/badge.svg?event=push)

**Project communication**  
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/sphinxsys.svg?style=social&label=Follow%20%40sphinxsys)](https://twitter.com/sphinxsys)
[![YouTube](https://img.shields.io/badge/YouTube-FF0000.svg?style=flat&logo=YouTube&logoColor=white)](https://www.youtube.com/channel/UCexdJbxOn9dvim6Jg1dnCFQ)
[![Bilibili](https://img.shields.io/badge/bilibili-%E5%93%94%E5%93%A9%E5%93%94%E5%93%A9-critical)](https://space.bilibili.com/1761273682/video)
[![QQ](https://img.shields.io/badge/QQ_Group-blue?logo=tencentqq&logoColor=white)](https://qm.qq.com/q/BZDAqz70Iw)

## Repository Description

SPHinXsys (pronunciation: s'fink-sis) is an acronym from **S**moothed **P**article **H**ydrodynamics for **in**dustrial comple**X** **sys**tems.
The multi-physics library uses SPH (smoothed particle hydrodynamics) as the underlying numerical method
for both particle-based and mesh-based discretization.
Due to the unified computational framework, SPHinXsys is able to carry out simulation and optimization at the same time.
For more information on the SPHinXsys project, please check the project website: <https://www.sphinxsys.org>.

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
<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/2d_examples/test_2d_flow_stream_around_fish/2d_flow_stream_around_fish.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/fish-swimming.gif" height="168px"></a>
<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/2d_examples/test_2d_column_collapse/column_collapse.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/2d_column_collapse.gif" height="168px"></a>
<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/extra_source_and_tests/test_2d_T_pipe_VIPO_shell/T_pipe_VIPO_shell.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/fluid-shell-interaction.gif" height="168px"></a>

## Fully compatible to classical FVM method

Through the unified computational framework in SPHinXsys,
the algorithms for particle methods are full compatible to those in the classical finite volume method (FVM).
The following gives an example of the flow around cylinder problem solved by FVM in SPHinXsys.

<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/2d_examples/test_2d_FVM_flow_around_cylinder/2d_FVM_flow_around_cylinder.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/fvm-sphinxsys-flow-around-cylinder.gif" height="168px"></a>

Note that the code for FVM algorithm is exact the same one for particle interaction in SPHinXsys.
The only difference is that SPHinXsys reads a predefined mesh, other than generate particles, before the computation.

## Target-driven optimization

The unique target-driven optimization is able to achieve the optimization target and physical solution all-in-once,
which is able to accelerate optimization process greatly.
The following gives an example of optimizing the conductivity distribution
for a thermal domain problem targeting minimum average temperature.

<a href="https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/optimization/test_2d_VP_heat_flux_optimization/VP_heat_flux_optimization.cpp">
<img src="https://github.com/Xiangyu-Hu/SPHinXsys-public-files/blob/master/videos/optimization.gif" height="192px"></a>

Note that the physical solution of the thermal domain (right) and the optimal distribution of conductivity (left)
are obtained at the same time when optimization is finished.
Also note that the entire optimization process is very fast and
only several times slower than that for a single physical solution with given conductivity distribution.  

## Python interface

While SPHinXsys is written in C++, it provides a python interface for users to write python scripts to control the simulation,
including carry out regression tests for continuous integration (CI) and other tasks.
One example is given below for the dambreak case.
Please check the source code of
[2D Dambreak case with python interface](https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/tests/test_python_interface/test_2d_dambreak_python/dambreak_python.cpp)
for the usage.

## Publications

Main publication on the library:

1. C. Zhang, M. Rezavand, Y. Zhu, Y. Yu, D. Wu, W. Zhang, J. Wang, X. Hu,
"SPHinXsys: an open-source multi-physics and multi-resolution library based on smoothed particle hydrodynamics",
Computer Physics Communications, 267, 108066, 2021.  
[![Main Publication](https://img.shields.io/badge/doi-10.1016%2Fj.cpc.2021.108066-d45815.svg)](https://doi.org/10.1016/j.cpc.2021.108066)  
[Google Scholar citations](https://scholar.google.com/scholar?cites=696006064513647619&as_sdt=2005&sciodt=0,5&hl=en)

The numerical methods and computational algorithms in SPHinXsys are based on the following [publications](assets/publication.md).

## Software architecture

SPHinXsys is cross-platform can be compiled and used in Windows, Linux and McOS systems.

## Installation, tutorial and documentation

For installation, program manual and tutorials, please check <https://www.sphinxsys.org/html/sphinx_index.html>.
Please check the documentation of the code at <https://xiangyu-hu.github.io/SPHinXsys/>.
For a Docker image, check <https://hub.docker.com/r/toshev/sphinxsys>.

## Interaction with SPHinXsys and the team

Thank you for using and supporting our open-source project!
We value each feedback.

#### For SPHinXsys users

Your input is crucial to us. We encourage you to report any issues you encounter with the library, including:

* Bug reports
* Poorly written code or algorithm designs
* Benchmark test issues, whether within the library or from literature, especially those highlighting potential deficiencies
* Other issues

We particularly appreciate feedback stemming from practical simulations or projects, as these insights are essential for improving SPHinXsys.

#### For SPHinXsys developers

If you don't have a GitHub account yet, please register for one. Fork the SPHinXsys repository to add new features or improve existing ones. Once your changes are ready, commit them and initiate a pull request to have your contributions merged into the main repository.

To ensure efficient and effective development, we prioritize addressing issues raised by active contributors—whether through code, documentation, or other means. We welcome any interaction with SPHinXsys and our team.

You can also join us as a collaborator, enabling you to branch directly within the main repository and review pull requests.

Together, we can build a leading-edge multi-physics library open for all!

If you have any further question, please contact <xiangyu.hu@tum.de>.
