import math
import re
import time
import torch
import io
import numpy as np
import test_2d_free_stream_around_fish_pybind as fish
import gymnasium as gym
from gymnasium import spaces

def checknan(reward):
    if math.isnan(reward):
        return True
    return False


class FISHEnv(gym.Env):
    metadata = {"render_modes": ["human", "rgb_array"], "render_fps": 30}

    def __init__(self, render_mode=None):
        self.alpha = 0.1
        self.episode = 1
        self.time_per_action = 0.05
        low_action = np.array([-1]).astype(np.float32)
        high_action = np.array([1]).astype(np.float32)
        # 4 * 3 pressure point, pos_x, pos_y
        low_obs = np.full(12, -10).astype(np.float32)
        high_obs = np.full(12, 10).astype(np.float32)
        self.obs = np.zeros(12)

        # self.action_space = spaces.Box(low_action, high_action)
        self.action_space = spaces.Discrete(5) # from 0 to 4 [0,1,2,3,4]
        # discrete from

        # self.action_space = spaces.Discrete(2)
        self.observation_space = spaces.Box(low_obs, high_obs)

        self.train_fish = fish.from_sph_relaxation(self.episode)  # generate configuration xml fish file
        self.center_x = 0.4
        self.center_y = 0.4


    def reset(self, seed=None, options=None):
        super().reset(seed=seed)
        self.total_reward = 0.0
        self.action_time = 0.0
        self.action_time_steps = 0

        self.am = 0.12
        self.frequency = 4

        self.train_fish = fish.from_sph_reload_and_train(self.episode)
        self.train_fish.SetFreq(self.frequency)
        self.train_fish.RunCase(self.episode, self.action_time)
        # total 4 points for corresponding point
        for i in range(4):
            self.obs[2 * i] = self.train_fish.GetFishPositionX(i)
            self.obs[2 * i + 1] = self.train_fish.GetFishPositionY(i)
        for i in range(4):
            self.obs[i] = self.train_fish.GetFishPressurePoint(i)

        self._get_obs = self.obs.astype(np.float32)
        return self._get_obs, {}

    def step(self, action):
        self.action_time_steps += 1
        self.freq = action + 1
        print(' frequency:', self.freq,  ' action_time_steps: ', self.action_time_steps)
        file = open(f'action_{self.episode}.txt', 'a')
        file.write('action_time:  ')
        file.write(str(self.action_time))
        file.write('  frequency:  ')
        file.write(str(self.freq)) #?
        file.write('\n')

        for i in range(10):
            self.train_fish.SetFreq(self.freq)
            self.action_time += self.time_per_action * self.alpha
            self.train_fish.RunCase(self.episode, self.action_time)

        self.fish_dx = 0
        self.fish_dy = 0
        for i in range(4):
            self.fish_dx += self.train_fish.GetFishPositionX(i)
            self.fish_dy += self.train_fish.GetFishPositionY(i)
        self.fish_dx /= 4
        self.fish_dy /= 4

        reward = 0
        dist = pow(pow(self.fish_dx - self.center_x, 2) + pow(self.fish_dy - self.center_y, 2), 0.5)
        if dist <= 0.31 and dist >= 0.29:
            reward += 1
        elif dist < 0.29:
            reward += 1 - abs(0.3 - dist) / 0.3
        elif dist > 0.31:
            reward -= 2 * abs(0.3 - dist) / 0.3
        reward_nan = checknan(reward)
        if reward_nan:  # check reward nan
            reward = 0

        for i in range(4):
            self.obs[2 * i] = self.train_fish.GetFishPositionX(i)
            self.obs[2 * i + 1] = self.train_fish.GetFishPositionY(i)
        for i in range(4):
            self.obs[i] = self.train_fish.GetFishPressurePoint(i)

        self._get_obs = self.obs.astype(np.float32)

        done = False
        for i in range(4):
            if self.train_fish.GetFishPositionX(i) < 0.01 or self.train_fish.GetFishPositionX(i) > 0.8:
                reward -= 50
                done = True
                break
            if self.train_fish.GetFishPositionY(i) < 0 or self.train_fish.GetFishPositionY(i) > 0.8:
                reward -= 50
                done = True
                break
        if self.action_time_steps > 1:
            done = True
        if done == False:
            if reward_nan:
                done = True
                reward = 0
        print('reward: ', reward, ' fish_x:', self.fish_dx, ' fish_y: ', self.fish_dy, ' dist: ', dist)

        file = open(f'reward_{self.episode}.txt', 'a')
        file.write('action_time:  ')
        file.write(str(self.action_time))
        file.write('  reward:  ')
        file.write(str(reward))
        file.write('  fish_x  ')
        file.write(str(self.fish_dx))
        file.write('  fish_y:  ')
        file.write(str(self.fish_dy))
        file.write('\n')
        file.close()
        self.total_reward += reward

        if done == True:
            file = open('reward.txt', 'a')
            file.write('episode:  ')
            file.write(str(self.episode))
            file.write('  reward:  ')
            file.write(str(self.total_reward))
            file.write('\n')
            file.close()
            self.episode += 1

        return self._get_obs, reward, done, False, {}

    def render(self):
        return 0

    def _render_frame(self):
        return 0

    def close(self):
        return 0


if __name__ == "__main__":
    env = FISHEnv()
    state_shape = env.observation_space.shape
    print('state_shape:', state_shape)
    action_shape = env.action_space.shape
    print('action_shape: ', action_shape)
    action_space = env.action_space
    print('action_space', action_space)
    init_state = env.reset()

    done = False
    while not done:
        act = env.action_space.sample()
        # print("this is sample action",act)
        s, r, done, _, _ = env.step(act)
        print("Action", act)
        print("Module State", s)
        print("Reward", r)
        print("Done", done)
