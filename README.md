# ![](SPHINXsys/logo.png) SPHinXsys

## Description

SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle Hydrodynamics for industrial compleX systems. 
It provides C++ APIs for physical accurate simulation and aims to model coupled industrial dynamic systems including fluid, solid, multi-body dynamics 
and beyond with SPH (smoothed particle hydrodynamics), a meshless computational method using particle discretization. 
For more information on the SPHinXsys project, please check the project website: https://www.sphinxsys.org.

SPHinXsys is a multi-physics, multi-resolution SPH library. 
Although it is not a standalone application itself, 
many examples designated for the specific type of applications are provided.

## Journal publications

The algorithms in SPHinXsys are based on the following publications:

1. Luhui Han and Xiangyu Hu, 
"SPH modeling of fluid-structure interaction", 
Journal of Hydrodynamics, 2018: 30(1):62-69.

2. Chi Zhang and Massoud Rezavand and Xiangyu Hu, 
"Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics", 
Journal of Computational Physics 404 (2020) 109135

3. Chi Zhang et al. 
"SPHinXsys: An open-source meshless, multi-resolution and multi-physics library",
Software Impacts, 6 (2020) 100033

4. Chi Zhang, Massoud Rezavand, Xiangyu Hu,
"A multi-resolution SPH method for fluid-structure interactions",
Journal of Computational Physics,  Journal of Computational Physics, 429 (2021) 110028.

5. Chi Zhang, Yanji Wei, Frederic Dias, Xiangyu Hu,
"An efficient fully Lagrangian solver for modeling wave interaction with oscillating wave energy converter",
Ocean Engineering, 
Volume 236, (2021) 109540

6. Chi Zhang, Jianhang Wang, Massoud Rezavand, Dong Wu, Xiangyu Hu,
"An integrative smoothed particle hydrodynamics framework for modeling cardiac function",
 Computer Methods in Applied Mechanics and Engineering, 381, 113847, 2021.

 7. C. Zhang, M. Rezavand, Y. Zhu, Y. Yu, D. Wu, W. Zhang, J. Wang, X. Hu, "SPHinXsys: an open-source multi-physics and multi-resolution library based on smoothed particle hydrodynamics", Computer Physics Communications, 267, 108066, 2021.

 8. Yujie Zhu, Chi Zhang, Yongchuan Yu, Xiangyu Hu, "A CAD-compatible body-fitted particle generator for arbitrarily complex geometry and its application to wave-structure interaction", Journal of Hydrodynamics, 33(2), 195-206, 2021.

## Software Architecture

SPHinXsys is cross-platform can be compiled and used in Windows, Linux and Apple systems.

## Installation, tutorial and documentation

For installation, program manual and tutorials, please check https://www.sphinxsys.org/html/sphinx_index.html. 
Please check the documentation of the code at https://xiangyu-hu.github.io/SPHinXsys/.

## Contribute to SPHinXsys

You are welcomed to contribute to SPHinXsys. 
As the code is on git-hub, you can register an account there (if you do not have a github account yet) 
and fork out the SPHinXsys repository.
You can work on the forked repository and add new features, and then commit them. 

### Contribute new features

After you have finished a new feature, 
which is a commit together with a new test case which uses the new feature,
you can send a pull request to the SPHinXsys repository. 
The SPHinXsys team will review the new code and give suggestions for revision, if there is.
After the revision is done, your contribution will be merged into the master branch of the original SPHinXsys repository.
When you write the code, please remember to add your name as an author in the header files that you have contributed.

### Contribute new test cases

Test cases are very important for testing the features and showing the ability of SPHinXsys.
If you think that a cases which implemented by you are very useful for the above propose,
you can move this test into the folder cases_test in your local commit 
which is directly branched from the master branch in the original SPHinXsys repository.
After test it, you can carry out a pull request.
Then, a standard review process will be carried out before merge it to the original SPHinXsys repository.
If a new test cases is added, it will be maintained together with code by SPHinXsys team. 

### Contribute challenging benchmarks

Improving numerical algorithms is one of our main drives 
and we are happy to test challenging cases and to be informed 
if there is any other code that does better on any benchmark case.
For this, you can inform us by the contact listed below. 
If you have the case implemented, you can do just as contribute a new test case, 
except move this case to cases_user folder. 
Then, after a standard pull request and review process, 
we will merge you code to master branch and discuss with you on the issues of the case 
and work together to improve the numerical algorithm for that case. 

If you have any further question, you are also welcomed to contact xiangyu.hu@tum.de.