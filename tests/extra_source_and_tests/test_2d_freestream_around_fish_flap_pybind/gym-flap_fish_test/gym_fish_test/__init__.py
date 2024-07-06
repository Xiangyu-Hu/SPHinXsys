from gymnasium.envs.registration import register

register(
    id="FISH-TEST-v0",
    entry_point="gym_fish_test.envs:FISHTestEnv",
    max_episode_steps=5000,
    reward_threshold=5000.0,
)
