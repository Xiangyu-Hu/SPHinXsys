# Contributing to SPHinXsys

The following is a set of guidelines for contributing to SPHinXsys and its packages,
which are hosted on GitHub. These are mostly guidelines, not rules.
Use your best judgment, and feel free to propose changes to this document in a pull request.

#### Table Of Contents

[Code of Conduct](#code-of-conduct)

[What should I know before I get started?](#what-should-i-know-before-i-get-started)

* [SPHinXsys Libraries](#sphinxsys-libraries)

[How Can I Contribute?](#how-can-i-contribute)

* [Reporting Bugs](#reporting-bugs)
* [Suggesting Enhancements](#suggesting-enhancements)
* [Your First Code Contribution](#your-first-code-contribution)
* [Adding New Feature](#adding-new-feature)
* [Pull Requests](#pull-requests)

[Styleguides](#styleguides)

* [Git Commit Messages](#git-commit-messages)
* [C++ Styleguide](#cpp-styleguide)
* [Documentation Styleguide](#documentation-styleguide)

[Additional Notes](#additional-notes)

* [Issue and Pull Request Labels](#issue-and-pull-request-labels)

## Code of Conduct

This project and everyone participating in it is governed by the [Contributor Covenant][homepage],
version 1.4, available at [https://contributor-covenant.org/version/1/4][version]

## What should I know before I get started?

### SPHinXsys Libraries

SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle Hydrodynamics
for industrial compleX systems. It provides C++ APIs for physical accurate simulation and
aims to model coupled industrial dynamic systems including fluid, solid, multi-body dynamics
and beyond with SPH (smoothed particle hydrodynamics),
a meshless computational method using particle discretization.

SPHinXsys has two major components: one is for the modelling classes,
such as SPHBodies, Materials, Particles, BodyRelations;
the other is for simulation methods for the physical dynamics which are derived from the base class of ParticleDynamics.
Both the models and simulation methods can be extended.

#### SPHinXsys Conventions

There are a few conventions that have developed over time around SPHinXsys:

* The main source is in the folder *src* and the examples and test problems are in the folder tests.
* In *src*, the sources which will be used for both 2D and 3D are in the folder shared.
* In shared, each folder is specific model or physical dynamics.
 For example, the folder materials is for modeling fluid, solid and more complex material properties.
 Another example is the folder particle dynamics. In it, for each discipline of physical dynamics, a separated folder is assigned.
* When different codes are required for 2D of 3D implementations.
 One need to create the folders with the same names in both the folders for 2D and 3D build.
* Each test case has been assigned with a unique folder in tests.

## How Can I Contribute?

### Reporting Bugs

When you are creating a bug or issue report, please include as many details as possible.

### Suggesting Enhancements

You are welcomed to submit an enhancement suggestion for SPHinXsys in the issues,
including completely new features and minor improvements to existing functionality.

### Your First Code Contribution

Unsure where to begin contributing to SPHinXsys?
You can start by looking through these `good first` and `help-wanted` issues:

* [Good first issues] - issues which should only require a few lines of code, and a test or two.
* [Help wanted issues] - issues which should be a bit more involved than `beginner` issues.

Both issue lists are sorted by total number of comments.

### Adding New Feature

In order to be more effective, we may need to follow some proved good practice on contributing a new feature. One thing that discussing or issue first then change the code after the discussion settled. Since SPHinXsys is already quite finished and big, the most concern on adding new feature is to avoid increasing technical debts, which troubles the maintenance greatly.
So, adding new features would usually follow the pipeline of generalization first and new feature second for minimum impact and easy testing. That also means that we first think how the new feature is related to the old code.  As SPHinXsys has a lot of classes and methods defined already, in most cases, the new feature and some old features can be generalized under a same new concept.

Therefore, the suggested practice is first doing a conservative modification, that is, reproduce the old functionality of the code by introduce the new concept using template. Then test the entire library, for this you can use draft pull request, Github will test for you automatically. After this is done, you can add you new feature as an incremental change, and then test your new feature with new test case.

In this way, we introduce minimum technical debt for more sustained development.
Certainly, this will be less straightforward as the simple feature adding approach, but much better in long term.

In summarize, you first introduce a new concept, implement the new concept to the old code and then add new feature using the same concept. In this way, you are contributing the library in design concept level other than simply add something.

### Pull Requests

The process described here has several goals:

* Maintain SPHinXsys's quality
* Fix problems that are important to users
* Engage the community in working toward the best possible SPHinXsys
* Enable a sustainable system for SPHinXsys's maintainers to review contributions

Please follow these steps to have your contribution considered by the maintainers:

1. Follow the [styleguides](#styleguides)
2. After you submit your pull request,
verify that all [status checks](https://help.github.com/articles/about-status-checks/)
are passing <details><summary>What if the status checks are failing?</summary>If a status check is failing,
and you believe that the failure is unrelated to your change,
please leave a comment on the pull request explaining why you believe the failure is unrelated.
A maintainer will re-run the status check for you.
If we conclude that the failure was a false positive,
then we will open an issue to track that problem with our status check suite.</details>

While the prerequisites above must be satisfied prior to having your pull request reviewed,
the reviewer(s) may ask you to complete additional design work, tests,
or other changes before your pull request can be ultimately accepted.

## Styleguides

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* When only changing documentation, include `[ci skip]` in the commit title
* Consider starting the commit message with an applicable emoji:
  * :art: `:art:` when improving the format/structure of the code
  * :racehorse: `:racehorse:` when improving performance
  * :non-potable_water: `:non-potable_water:` when plugging memory leaks
  * :memo: `:memo:` when writing docs
  * :penguin: `:penguin:` when fixing something on Linux
  * :apple: `:apple:` when fixing something on macOS
  * :checkered_flag: `:checkered_flag:` when fixing something on Windows
  * :bug: `:bug:` when fixing a bug
  * :fire: `:fire:` when removing code or files
  * :green_heart: `:green_heart:` when fixing the CI build
  * :white_check_mark: `:white_check_mark:` when adding tests
  * :lock: `:lock:` when dealing with security
  * :arrow_up: `:arrow_up:` when upgrading dependencies
  * :arrow_down: `:arrow_down:` when downgrading dependencies
  * :shirt: `:shirt:` when removing linter warnings

### C++ Styleguide

The code must adhere to [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).

### Documentation Styleguide

* Use [Markdown](https://daringfireball.net/projects/markdown).
