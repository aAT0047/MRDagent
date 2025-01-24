import random
import gym
import numpy as np
from collections import deque
import torch
import torch.nn as nn
import torch.optim as optim
from Env import CustomEnv  # 假设你的自定义环境

class DQN(nn.Module):
    def __init__(self, state_size, action_size):
        super(DQN, self).__init__()
        self.fc1 = nn.Linear(state_size, 64)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(64, 64)
        self.fc3 = nn.Linear(64, action_size)

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.relu(self.fc2(x))
        return self.fc3(x)

class Agent:
    def __init__(self, state_size, action_size, num_discrete_bins=10):
        self.state_size = state_size
        self.action_size = action_size
        self.num_discrete_bins = num_discrete_bins
        self.memory = deque(maxlen=2000)
        self.gamma = 0.95
        self.epsilon = 1.0
        self.epsilon_min = 0.01
        self.epsilon_decay = 0.995
        self.model = DQN(state_size, action_size * num_discrete_bins)
        self.target_model = DQN(state_size, action_size * num_discrete_bins)
        self.optimizer = optim.Adam(self.model.parameters(), lr=0.001)
        self.update_target_model()

    def discretize_action(self, continuous_action):
        bins = np.linspace(-1, 1, self.num_discrete_bins)
        return np.array([np.digitize(action, bins) - 1 for action in continuous_action])

    def map_discrete_to_continuous(self, discrete_action):
        bins = np.linspace(-1, 1, self.num_discrete_bins)
        return np.array([bins[action] for action in discrete_action])

    def act(self, state):
        if random.random() <= self.epsilon:
            discrete_action = np.random.randint(0, self.num_discrete_bins, self.action_size)
            return self.map_discrete_to_continuous(discrete_action)
        
        state = torch.FloatTensor(state).unsqueeze(0)
        q_values = self.model(state).view(self.action_size, self.num_discrete_bins)
        discrete_action = torch.argmax(q_values, dim=1).cpu().numpy()
        return self.map_discrete_to_continuous(discrete_action)

    def remember(self, state, action, reward, next_state, done):
        self.memory.append((state, action, reward, next_state, done))

    def process_states(self, states):
        processed_states = [state[0] if isinstance(state, tuple) else state for state in states]
        states_array = np.array(processed_states, dtype=np.float32)
        return torch.tensor(states_array)

    def replay(self, batch_size):
        if len(self.memory) < batch_size:
            return

        minibatch = random.sample(self.memory, batch_size)
        states, actions, rewards, next_states, dones = zip(*minibatch)
        states = self.process_states(states)
        next_states = self.process_states(next_states)
        actions = np.array([self.discretize_action(action) for action in actions])
        actions = torch.LongTensor(actions)

        rewards = torch.FloatTensor(rewards)
        dones = torch.tensor(dones, dtype=torch.bool)

        current_qs = self.model(states).view(-1, self.action_size, self.num_discrete_bins)
        next_qs = self.target_model(next_states).view(-1, self.action_size, self.num_discrete_bins)

        max_next_qs = torch.max(next_qs, dim=2)[0]
        max_next_qs[dones] = 0
        targets = rewards + self.gamma * max_next_qs.sum(dim=1)

        targets_f = current_qs.clone()
        for i in range(batch_size):
            for j in range(self.action_size):
                targets_f[i, j, actions[i, j]] = targets[i]

        loss = nn.MSELoss()(current_qs, targets_f)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()
        self.update_epsilon()

    def update_target_model(self):
        self.target_model.load_state_dict(self.model.state_dict())

    def update_epsilon(self):
        if self.epsilon > self.epsilon_min:
            self.epsilon *= self.epsilon_decay

def train_dqn(episodes, env_kwargs):
    env = gym.make('CustomEnv-v0', **env_kwargs)
    state_size = env.observation_space.shape[0]
    action_size = env.action_space.shape[0]
    agent = Agent(state_size, action_size)
    best_score = -float('inf')
    best_model_state = None
    batch_size = 4
    UPDATE_TARGET_EVERY = 5

    for e in range(episodes):
        state, _ = env.reset()
        total_reward = 0
        terminated = truncated = False

        while not (terminated or truncated):
            action = agent.act(state)
            next_state, reward, terminated, truncated, _ = env.step(action)

            total_reward += reward
            agent.remember(state, action, reward, next_state, terminated)
            state = next_state

            if terminated or truncated:
                print(f"Episode: {e+1}/{episodes} | Score: {total_reward}")
                if total_reward > best_score:
                    best_score = total_reward
                    best_model_state = agent.model.state_dict()
                break

            if len(agent.memory) > batch_size:
                agent.replay(batch_size)

        if e % UPDATE_TARGET_EVERY == 0:
            agent.update_target_model()

    env.close()
    print(f"Best Score: {best_score}")
    return best_model_state


# env_kwargs = {
#     'ref_fasta': 'ref_fasta',
#     'input_bam': 'input_bam',
#     'output_path': 'output_path',
#     'base_tar_gz': 'base_tar_gz',
#     'base_vcf': 'base_vcf',
#     'truth_vcf_path': 'truth_vcf_path'
# }

# best_model_state = train_dqn(100, env_kwargs)
