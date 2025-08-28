========================
How to train deep reinforcement learning model in SPHinXsys
========================

By leveraging the Python interface of SPHinXsys, one can rapidly and modularly construct SPHinXsys simulation environments required for deep reinforcement learning (DRL). Effective training can be conducted using various algorithms from the mainstream DRL platform Tianshou. Please refer to the source code of
[2D owsc case with python interface](https://github.com/Xiangyu-Hu/SPHinXsys/tree/master/tests/extra_source_and_tests/test_2d_owsc_python). 

Step-by-Step instructions
---------------------------------------

Use conda or Python to create a virtual Python environment.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Python 3.10 is recommended for compatibility.

..  code-block:: bash

        conda create -n drl_env python=3.10

Activate the python environment and install required libraries.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
..  code-block:: bash

        conda activate drl_env
        pip install tianshou

Set up the SPHinXsys-based reinforcement learning environment in the "drl_gym_environments" directory.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

..  code-block:: bash

        cd path_to_deep_reinforcement_learning_tool/drl_gym_environments
        pip install -e .

Navigate to the `drl_tianshou_training` directory and start the training process.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

..  code-block:: bash

        cd path_to_deep_reinforcement_learning_tool/drl_tianshou_training
        python sac.py

Test the trained model.
^^^^^^^^^^^^^^^^^^^^^^^^^^

Modify three arguments:

'resume-path': the path where your policy was saved,

'test-number': the number of test episodes you want to run,

'watch': True means testing, False means training.

..  code-block:: bash

        python sac.py


Detailed explanation of "drl_gym_environments"
--------------------------------------------------

The main package directory that includes the environment registration logic and all the custom environment implementations.

**drl_gym_environments/setup.py**:  

Provides the configuration for installing the package.  

**drl_gym_environments/gym_env_owsc**:  

This folder contains the OWSC custom environment.

**gym_env_owsc/__init__.py**:  

This file contains the registration logic for all custom environments. Each environment is registered with a unique ID, which can be used with Gymnasium’s `gym.make()` function.

**gym_env_owsc/envs**:  

This sub-directory contains the actual implementation of all the custom environments.

**envs/__init__.py**:  

Imports all custom environments so they can be properly registered.

**envs/owsc.py**:  

Implements the OWSC environment, following the standard Gymnasium `Env` interface. The environment defines unique observation and action spaces and includes specific environment dynamics in the `reset()`, `step()`, and `render()` methods.