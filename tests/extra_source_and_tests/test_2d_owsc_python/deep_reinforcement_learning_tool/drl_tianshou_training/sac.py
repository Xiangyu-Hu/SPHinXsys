#!/usr/bin/env python3
import gymnasium as gym
# Custom OWSC environment
import gym_env_owsc

import argparse
import datetime
import os
import pprint

import numpy as np
import torch

from torch.utils.tensorboard import SummaryWriter

from tianshou.data import Collector, ReplayBuffer, VectorReplayBuffer
from tianshou.policy import SACPolicy
from tianshou.trainer import OffpolicyTrainer
from tianshou.env import SubprocVectorEnv
from tianshou.utils import TensorboardLogger, WandbLogger
from tianshou.utils.net.common import Net
from tianshou.utils.net.continuous import ActorProb, Critic


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=str, default="OWSC-v0")  # Environment ID
    parser.add_argument("--seed", type=int, default=0) 
    parser.add_argument("--buffer-size", type=int, default=10100)  
    parser.add_argument("--hidden-sizes", type=int, nargs="*", default=[64, 64])  
    parser.add_argument("--actor-lr", type=float, default=1e-4)  
    parser.add_argument("--critic-lr", type=float, default=1e-4)  
    parser.add_argument("--gamma", type=float, default=0.95)  
    parser.add_argument("--tau", type=float, default=0.005)  
    parser.add_argument("--alpha", type=float, default=0.25)  
    parser.add_argument("--auto-alpha", default=False, action="store_true") 
    parser.add_argument("--alpha-lr", type=float, default=5e-4)  
    parser.add_argument("--start-timesteps", type=int, default=100)  
    parser.add_argument("--epoch", type=int, default=100)  
    parser.add_argument("--step-per-epoch", type=int, default=100)  
    parser.add_argument("--step-per-collect", type=int, default=1)  
    parser.add_argument("--update-per-step", type=int, default=1) 
    parser.add_argument("--n-step", type=int, default=1) 
    parser.add_argument("--batch-size", type=int, default=64) 
    parser.add_argument("--training-num", type=int, default=1) 
    parser.add_argument('--test-num', type=int, default=0)
    parser.add_argument("--logdir", type=str, default="log") 
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu") 
    parser.add_argument("--resume-path", type=str, default=None) 
    parser.add_argument("--logger", type=str, default="tensorboard", choices=["tensorboard", "wandb"]) 
    parser.add_argument("--wandb-project", type=str, default="mujoco.benchmark")  
    parser.add_argument("--watch", default=False, action="store_true", help="Watch the play of pre-trained policy only")
    return parser.parse_args()


def training_sac(args=get_args()):
    """Main function for setting up and training a SAC agent in the OWSC environment."""
    
    envs = gym.make(args.task)
    # Create vectorized environments for parallel training
    # envs = SubprocVectorEnv([
    #     lambda i=i: gym.make(args.task, parallel_envs=i)  # Create environments with different parallel_envs IDs
    #     for i in range(args.training_num)])

    # Retrieve environment observation and action space details
    args.state_shape = envs.observation_space.shape or envs.observation_space.n
    args.action_shape = envs.action_space.shape or envs.action_space.n
    args.max_action = envs.action_space.high[0]

    # seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    # Define the actor and critic neural networks for SAC
    net_a = Net(args.state_shape, hidden_sizes=args.hidden_sizes, device=args.device)
    actor = ActorProb(net_a, args.action_shape, device=args.device, max_action=args.max_action, conditioned_sigma=True).to(args.device)
    actor_optim = torch.optim.Adam(actor.parameters(), lr=args.actor_lr)

    net_c1 = Net(args.state_shape, args.action_shape, hidden_sizes=args.hidden_sizes, concat=True, device=args.device)
    critic1 = Critic(net_c1, device=args.device).to(args.device)
    critic1_optim = torch.optim.Adam(critic1.parameters(), lr=args.critic_lr)

    net_c2 = Net(args.state_shape, args.action_shape, hidden_sizes=args.hidden_sizes, concat=True, device=args.device)
    critic2 = Critic(net_c2, device=args.device).to(args.device)
    critic2_optim = torch.optim.Adam(critic2.parameters(), lr=args.critic_lr)

    # Optionally, use automatic tuning of the temperature parameter (alpha)
    if args.auto_alpha:
        target_entropy = -np.prod(envs.action_space.shape)
        log_alpha = torch.zeros(1, requires_grad=True, device=args.device)
        alpha_optim = torch.optim.Adam([log_alpha], lr=args.alpha_lr)
        args.alpha = (target_entropy, log_alpha, alpha_optim)

    # Initialize the SAC policy
    policy = SACPolicy(
        actor=actor,
        actor_optim=actor_optim,
        critic1=critic1,
        critic1_optim=critic1_optim,
        critic2=critic2,
        critic2_optim=critic2_optim,
        tau=args.tau,
        gamma=args.gamma,
        alpha=args.alpha,
        estimation_step=args.n_step,
        action_space=envs.action_space,
    )

    # load a previous policy
    if args.resume_path:
        policy.load_state_dict(torch.load(args.resume_path, map_location=args.device))
        print("Loaded agent from: ", args.resume_path)

    # Setup replay buffer and collector for training
    buffer = VectorReplayBuffer(args.buffer_size, len(envs)) if args.training_num > 1 else ReplayBuffer(args.buffer_size)
    train_collector = Collector(policy, envs, buffer, exploration_noise=True)
    test_collector=None
    train_collector.collect(n_step=args.start_timesteps, random=True)

    # Setup logging
    now = datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    log_path = os.path.join(args.logdir, os.path.join(args.task, "sac", str(args.seed), now))
    writer = SummaryWriter(log_path)
    writer.add_text("args", str(args))
    logger = TensorboardLogger(writer) if args.logger == "tensorboard" else None  # Only Tensorboard logging is used

    # Save functions for best policy and checkpoints
    def save_best_fn(policy):
        torch.save(policy.state_dict(), os.path.join(log_path, "policy.pth"))

    # Training loop using OffpolicyTrainer
    if not args.watch:
        result = OffpolicyTrainer(
            policy=policy,
            train_collector=train_collector,
            test_collector=test_collector,
            episode_per_test=args.test_num,
            max_epoch=args.epoch,
            step_per_epoch=args.step_per_epoch,
            step_per_collect=args.step_per_collect,
            batch_size=args.batch_size,
            save_best_fn=save_best_fn,
            logger=logger,
            update_per_step=args.update_per_step,
        ).run()
        pprint.pprint(result)


if __name__ == "__main__":
    training_sac()
