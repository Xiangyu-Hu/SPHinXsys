import gym_fish
import gym_fish_test

# import gymnasium as gym
# import numpy as np
# import torch
# import os
# import datetime
# import warnings
# import pprint

# from tianshou.data import Collector, ReplayBuffer, VectorReplayBuffer
# from tianshou.exploration import GaussianNoise
# from tianshou.policy import SACPolicy
# from tianshou.trainer import offpolicy_trainer
# from tianshou.utils.net.common import ActorCritic, Net
# from tianshou.utils.net.continuous import ActorProb, Critic
# from tianshou.utils import TensorboardLogger, WandbLogger

# from torch import nn
# from torch.distributions import Independent, Normal
# from torch.optim.lr_scheduler import LambdaLR
# from torch.utils.tensorboard import SummaryWriter


import argparse
import os
import torch
import numpy as np
from torch import nn
from typing import Any, Dict, Tuple, Union, Optional, Sequence

import gymnasium as gym
import numpy as np
import torch
from torch.utils.tensorboard import SummaryWriter

from tianshou.data import (
    Collector,
    PrioritizedVectorReplayBuffer,
    ReplayBuffer,
    VectorReplayBuffer,
)
from tianshou.env import DummyVectorEnv
from tianshou.policy import DQNPolicy
from tianshou.policy.base import BasePolicy
from tianshou.trainer import OffpolicyTrainer
from tianshou.utils import TensorboardLogger
from tianshou.utils.net.common import Net

task = "FISH"
seed = 0
logdir = "log"
watch = False
eps_test = 0.005
eps_train = 1.
eps_train_final = 0.05
buffer_size = 100000
lr = 0.0001
gamma = 0.99
n_step = 3
target_update_freq = 100
epoch = 3
step_per_epoch = 10
step_per_collect = 10
update_per_step = 0.1
batch_size = 32
hidden_sizes = [128, 128]
dueling_q_hidden_sizes =  [128, 128]
dueling_v_hidden_sizes = [128, 128]

training_num = 1
test_num = 10
render = 0.
device = 'cuda' if torch.cuda.is_available() else 'cpu'
frames_stack = 4
resume_path = None
watch = False
save_buffer_name = None
n_step = 1
resume_path = None
training_num = 1
buffer_size = 650
start_time_steps = 1
epoch = 1
step_per_epoch = 1
step_per_collect = 1
batch_size = 64
update_per_step = 1
test_num = 1
render = 0.


def save_best_fn(policy):
    torch.save(policy.state_dict(), os.path.join(log_path, "policy.pth"))

if __name__ == "__main__":
    train_envs = gym.make('FISH-v0')
    test_envs = gym.make('FISH-TEST-v0')

    state_shape = train_envs.observation_space.shape
    action_shape = train_envs.action_space.shape
    action_space = train_envs.action_space
    # should be N_FRAMES x H x W
    print("Observations shape:", state_shape)
    print("Actions shape:", action_shape)
    # make environments
    # train_envs = SubprocVectorEnv([lambda: make_atari_env(args)
    #                                for _ in range(args.training_num)])
    # test_envs = SubprocVectorEnv([lambda: make_atari_env_watch(args)
    #                               for _ in range(args.test_num)])
    # seed
    np.random.seed(seed)
    torch.manual_seed(seed)

    # train_envs.seed(seed)
    # test_envs.seed(seed)
    # define model

    Q_param = {"hidden_sizes": dueling_q_hidden_sizes}
    V_param = {"hidden_sizes": dueling_v_hidden_sizes}
    net = Net(state_shape, action_shape,
              hidden_sizes=hidden_sizes, device=device,
              dueling_param=(Q_param, V_param)).to(device)
    optim = torch.optim.Adam(net.parameters(), lr=lr)
    # define policy
    policy = DQNPolicy(net, optim, gamma, n_step,
                       target_update_freq=target_update_freq)
    # collector
    # buffer = ReplayBuffer(buffer_size)
    train_collector = Collector(
        policy,
        train_envs,
        ReplayBuffer(buffer_size),
        exploration_noise=True
    )
    test_collector = Collector(policy, test_envs, exploration_noise=True)
    # policy.set_eps(1)
    train_collector.collect(n_step=batch_size * training_num)
    # log
    log_path = os.path.join(logdir, task, 'dqn')
    writer = SummaryWriter(log_path)
    logger = TensorboardLogger(writer)


    def save_fn(policy):
        torch.save(policy.state_dict(), os.path.join(log_path, 'policy.pth'))

    def stop_fn(mean_rewards):
        if env.env.spec.reward_threshold:
            return mean_rewards >= env.spec.reward_threshold
        elif 'Pong' in task:
            return mean_rewards >= 20
        else:
            return False

    def train_fn(epoch, env_step):
        # nature DQN setting, linear decay in the first 1M steps
        if env_step <= 1e6:
            eps = eps_train - env_step / 1e6 * \
                (eps_train - eps_train_final)
        else:
            eps = eps_train_final
        policy.set_eps(eps)
        logger.write('train/eps', env_step, eps)

    def test_fn(epoch, env_step):
        policy.set_eps(eps_test)

        # trainer
    result = offpolicy_trainer(
        policy,
        train_collector,
        test_collector,
        epoch,
        step_per_epoch,
        step_per_collect,
        test_num,
        batch_size,
        update_per_step=update_per_step,
        stop_fn=stop_fn,
        train_fn=train_fn,
        test_fn=test_fn,
        save_fn=save_fn,
        logger=logger
    )


    pprint.pprint(result)
    print(f'success')