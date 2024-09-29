### Detailed Explanation

- **`drl_gym_environments/`**:  
  The main package directory that includes the environment registration logic and all the custom environment implementations.

  - **`drl_gym_environments/setup.py`**:  
    Provides the configuration for installing the package.  
    To install the environments in editable mode, run the following command from the project’s root directory:

    ```bash
    pip install -e .
    ```

  - **`drl_gym_environments/gym_env_owsc/`**:  
    This folder contains the OWSC custom environment.

    - **`gym_env_owsc/__init__.py`**:  
      This file contains the registration logic for all custom environments.  
      Each environment is registered with a unique ID, which can be used with Gymnasium’s `gym.make()` function.

    - **`gym_env_owsc/envs/`**:  
      This sub-directory contains the actual implementation of all the custom environments.

      - **`envs/__init__.py`**:  
        Imports all custom environments so they can be properly registered.

      - **`envs/owsc.py`**:  
        Implements the OWSC environment, following the standard Gymnasium `Env` interface.  
        The environment defines unique observation and action spaces and includes specific environment dynamics in the `reset()`, `step()`, and `render()` methods.
