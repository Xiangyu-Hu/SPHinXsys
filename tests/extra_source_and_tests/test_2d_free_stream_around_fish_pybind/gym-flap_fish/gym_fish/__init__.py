from gymnasium.envs.registration import register

register(
    id="FISH-v0",
    entry_point="gym_fish.envs:FISHEnv",
    max_episode_steps=5000,
    reward_threshold=5000.0,
)
