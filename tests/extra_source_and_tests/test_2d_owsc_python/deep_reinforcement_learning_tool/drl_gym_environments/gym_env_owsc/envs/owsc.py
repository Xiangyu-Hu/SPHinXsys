import sys
import math
import numpy as np
import gymnasium as gym
from gymnasium import spaces
# add dynamic link library or shared object to python env
sys.path.append('/path/to/SPHinXsys/case/lib/dynamic link library or shared object')
import test_2d_owsc_python as test_2d


class OWSCEnv(gym.Env):
    """Custom Environment without rendering."""
    # metadata = {"render_modes": ["human", "rgb_array"], "render_fps": 30}

    def __init__(self, render_mode=None, parallel_envs=0):
        # Initialize environment parameters
        self.parallel_envs = parallel_envs  # Identifier for parallel simulation environments
        self.episode = 1                    # Current episode number
        self.time_per_action = 0.1          # Time interval per action step
        self.low_action = -1.0              # Minimum action value
        self.max_action = 1.0               # Maximum action value
        self.update_per_action = 10         # The action's effect is applied in smaller iterations within one action time step
        self.low_obs = -10.0                # Minimum observation value
        self.high_obs = 10.0                # Maximum observation value
        self.obs_numbers = 16               # Number of observation variables

        # Define action and observation spaces for Gym
        low_action = np.array([self.low_action]).astype(np.float32)
        high_action = np.array([self.max_action]).astype(np.float32)
        low_obs = np.full(self.obs_numbers, self.low_obs).astype(np.float32)
        high_obs = np.full(self.obs_numbers, self.high_obs).astype(np.float32)

        self.action_space = spaces.Box(low_action, high_action)  # Continuous action space
        self.observation_space = spaces.Box(low_obs, high_obs)  # Continuous observation space

    # Reset the environment at the beginning of each episode
    def reset(self, seed=None, options=None):
        super().reset(seed=seed)

        # Initialize the OWSC simulation with the given episode and environment setup
        self.owsc = test_2d.owsc_from_sph_cpp(self.parallel_envs, self.episode)
        self.action_time_steps = 0             # Track the number of action steps
        self.action_time = 0.5                 # Initialize action time
        self.damping_coefficient = 50          # Set damping coefficient for the environment
        self.total_reward_per_episode = 0.0    # Track total reward in each episode

        # Start the simulation with the given action time and damping coefficient
        self.owsc.run_case(self.action_time, self.damping_coefficient)
        
        # Initialize observation array with zero values
        self.observation = np.zeros(self.obs_numbers)
        # Fill the observation array with values from the OWSC simulation
        for i in range(0, 2):  
            self.observation[i] = self.owsc.get_wave_height(i)
            self.observation[i + 2] = self.owsc.get_wave_velocity(i, 0)
            self.observation[i + 4] = self.owsc.get_wave_velocity(i, 1)
            self.observation[i + 6] = self.owsc.get_wave_velocity_on_flap(i, 0)
            self.observation[i + 8] = self.owsc.get_wave_velocity_on_flap(i, 1)
            self.observation[i + 10] = self.owsc.get_flap_position(i, 0)
            self.observation[i + 12] = self.owsc.get_flap_position(i, 1)
        self.observation[14] = self.owsc.get_flap_angle()
        self.observation[15] = self.owsc.get_flap_angle_rate()
        
        self._get_obs = self.observation.astype(np.float32)

        return self._get_obs, {}

    def step(self, action):
        self.action_time_steps += 1
         # Apply the action to change the damping coefficient
        self.damping_change = 5.0 * action[0]
        # Penalty for invalid actions
        penality_0 = 0.0
        # Ensure the damping coefficient stays within valid bounds
        if self.damping_coefficient + self.damping_change < 0.01:
            self.damping_change = 0.01 - self.damping_coefficient
            penality_0 = - 1.0
        if self.damping_coefficient + self.damping_change > 100:
            self.damping_change = 100 - self.damping_coefficient
            penality_0 = - 1.0

        reward_0 = 0.0
        for i in range(self.update_per_action):
            self.flap_angle_rate_previous = self.owsc.get_flap_angle_rate()
            self.damping_coefficient += self.damping_change / self.update_per_action
            self.action_time += self.time_per_action / self.update_per_action
            self.owsc.run_case(self.action_time, self.damping_coefficient)
            self.flap_angle_rate_now = self.owsc.get_flap_angle_rate()
            # Calculate reward based on energy (flap angle rate)
            reward_0 += self.damping_coefficient * math.pow(0.5 * (self.flap_angle_rate_now + self.flap_angle_rate_previous), 2) * self.time_per_action / self.update_per_action
        # Add any penalties to the reward
        reward = reward_0 + penality_0
        self.total_reward_per_episode += reward

        # Update observations from the OWSC simulation
        for i in range(0, 2):  
            self.observation[i] = self.owsc.get_wave_height(i)
            self.observation[i + 2] = self.owsc.get_wave_velocity(i, 0)
            self.observation[i + 4] = self.owsc.get_wave_velocity(i, 1)
            self.observation[i + 6] = self.owsc.get_wave_velocity_on_flap(i, 0)
            self.observation[i + 8] = self.owsc.get_wave_velocity_on_flap(i, 1)
            self.observation[i + 10] = self.owsc.get_flap_position(i, 0)
            self.observation[i + 12] = self.owsc.get_flap_position(i, 1)
        self.observation[14] = self.owsc.get_flap_angle()
        self.observation[15] = self.owsc.get_flap_angle_rate()

        self._get_obs = self.observation.astype(np.float32)

        # Log action and reward information to files
        with open(f'action_env{self.parallel_envs}_epi{self.episode}.txt', 'a') as file:
            file.write(f'action_time:  {self.action_time}  action:  {self.damping_coefficient}\n')

        with open(f'reward_env{self.parallel_envs}_epi{self.episode}.txt', 'a') as file:
            file.write(f'action_time:  {self.action_time}  reward:  {reward}\n')

        # Check if the episode is done after 200 steps
        if self.action_time_steps > 99:
            done = True
            with open(f'reward_env{self.parallel_envs}.txt', 'a') as file:
                file.write(f'episode:  {self.episode}  total_reward:  {self.total_reward_per_episode}\n')
            self.episode += 1
        else:
            done = False

        # Return the updated observation, reward, done flag, and additional info
        return self._get_obs, reward, done, False, {}

    # Render method (optional, no rendering in this case)
    def render(self):
        return 0
    
    # Additional render frame logic (not implemented)
    def _render_frame(self):
        return 0

    # Close the environment and cleanup (optional)
    def close(self):
        return 0