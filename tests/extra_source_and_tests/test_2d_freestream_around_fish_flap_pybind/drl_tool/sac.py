import gym_fish
import gym_fish_test

import gymnasium as gym
import numpy as np
import torch
import os
import datetime
import warnings
import pprint

from tianshou.data import Collector, ReplayBuffer, VectorReplayBuffer
from tianshou.exploration import GaussianNoise
from tianshou.policy import SACPolicy
from tianshou.trainer import offpolicy_trainer
from tianshou.utils.net.common import ActorCritic, Net
from tianshou.utils.net.continuous import ActorProb, Critic
from tianshou.utils import TensorboardLogger, WandbLogger

from torch import nn
from torch.distributions import Independent, Normal
from torch.optim.lr_scheduler import LambdaLR
from torch.utils.tensorboard import SummaryWriter

task = "FISH"
seed = 0
logdir = "log"
watch = False
exploration_noise = 0.1
policy_noise = 0.2
noise_clip = 0.5
hidden_sizes = [64, 64]
actor_lr = 0.0001
critic_lr = 0.0001
tau = 0.005
gamma = 0.99
alpha = 0.2
auto_alpha = False
alpha_lr = 0.0003

update_actor_freq = 2
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
    warnings.filterwarnings('ignore')

    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    train_envs = gym.make('FISH-v0')
    test_envs = gym.make('FISH-TEST-v0')

    state_shape = train_envs.observation_space.shape
    action_shape = train_envs.action_space.shape
    action_space = train_envs.action_space
    max_action = train_envs.action_space.high[0]
    exploration_noise = exploration_noise * max_action
    noise_clip = noise_clip * max_action
    # seed
    np.random.seed(seed)
    torch.manual_seed(seed)
    # model
    net_a = Net(state_shape, hidden_sizes=hidden_sizes, device=device)
    actor = ActorProb(net_a, action_shape, max_action=max_action, device=device).to(device)
    actor_optim = torch.optim.Adam(actor.parameters(), lr=actor_lr)
    net_c1 = Net(
        state_shape,
        action_shape,
        hidden_sizes=hidden_sizes,
        concat=True,
        device=device,
    )
    net_c2 = Net(
        state_shape,
        action_shape,
        hidden_sizes=hidden_sizes,
        concat=True,
        device=device,
    )
    critic1 = Critic(net_c1, device=device).to(device)
    critic1_optim = torch.optim.Adam(critic1.parameters(), lr=critic_lr)
    critic2 = Critic(net_c2, device=device).to(device)
    critic2_optim = torch.optim.Adam(critic2.parameters(), lr=critic_lr)

    if auto_alpha:
        target_entropy = -np.prod(env.action_space.shape)
        log_alpha = torch.zeros(1, requires_grad=True, device=args.device)
        alpha_optim = torch.optim.Adam([log_alpha], lr=args.alpha_lr)
        alpha = (target_entropy, log_alpha, alpha_optim)


    policy = SACPolicy(
        actor=actor,
        actor_optim=actor_optim,
        critic1=critic1,
        critic1_optim=critic1_optim,
        critic2=critic2,
        critic2_optim=critic2_optim,
        tau=tau,
        gamma=gamma,
        alpha=alpha,
        estimation_step=n_step,
        action_space=action_space,
    )

    # policy = TD3Policy(
    #     actor,
    #     actor_optim,
    #     critic1,
    #     critic1_optim,
    #     critic2,
    #     critic2_optim,
    #     tau=tau,
    #     gamma=gamma,
    #     exploration_noise=GaussianNoise(sigma=exploration_noise),
    #     policy_noise=policy_noise,
    #     update_actor_freq=update_actor_freq,
    #     noise_clip=noise_clip,
    #     estimation_step=n_step,
    #     action_space=action_space,
    # )

    # load a previous policy
    if resume_path:
        policy.load_state_dict(torch.load(resume_path, map_location=device))
        print("Loaded agent from: ", resume_path)

    # collector
    if training_num > 1:
        buffer = VectorReplayBuffer(buffer_size, len(train_envs))
    else:
        buffer = ReplayBuffer(buffer_size)
    train_collector = Collector(policy, train_envs, buffer, exploration_noise=True)
    test_collector = Collector(policy, test_envs)
    train_collector.collect(n_step=start_time_steps, random=True)

    # log
    now = datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    algo_name = "sac"
    log_name = os.path.join(task, algo_name, str(seed), now)
    log_path = os.path.join(logdir, log_name)

    # logger
    writer = SummaryWriter(log_path)
    logger = TensorboardLogger(writer)

    if not watch:
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
            save_best_fn=save_best_fn,
            logger=logger,
            update_per_step=update_per_step,
            test_in_train=False,
        )
        pprint.pprint(result)

    # Let's watch its performance!
    policy.eval()
    # test_envs.seed()
    test_collector.reset()
    result = test_collector.collect(n_episode=test_num)
    print(f'success')