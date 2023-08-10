import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset
class Dataset(Dataset):

    def __init__(self, keys , all_torch):
        self.keys = keys
        self.all_torch = all_torch

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, index):

        key = self.keys[index]
        data = self.all_torch[key].to(torch.float32)
        sample = {'data': data,
                  'key': key \
                  }
        return sample

class MLP1Layer(nn.Module):
    def __init__(self, n_input, n_hidden, n_output):
        super(MLP1Layer, self).__init__()
        self.fc1 = nn.Linear(n_input, n_hidden)
        self.fc2 = nn.Linear(n_hidden, n_output)

    def forward(self, x):
        out = self.fc1(x)
        out = F.relu(out)
        out = self.fc2(out)
        return out
class MLP2Layer(nn.Module):
    def __init__(self, n_input, n_hidden_1, n_hidden_2, n_output):
        super(MLP2Layer, self).__init__()
        self.fc1 = nn.Linear(n_input, n_hidden_1)
        self.fc2 = nn.Linear(n_hidden_1, n_hidden_2)
        self.fc3 = nn.Linear(n_hidden_2, n_output)

    def forward(self, x):
        out = self.fc1(x)
        out = F.relu(out)
        out = self.fc2(out)
        out = F.relu(out)
        out = self.fc3(out)
        return out
class MLP3Layer(nn.Module):
    def __init__(self, n_input, n_hidden_1, n_hidden_2,n_hidden_3, n_output):
        super(MLP3Layer, self).__init__()
        self.fc1 = nn.Linear(n_input, n_hidden_1)
        self.fc2 = nn.Linear(n_hidden_1, n_hidden_2)
        self.fc3 = nn.Linear(n_hidden_2, n_hidden_3)
        self.fc4 = nn.Linear(n_hidden_3, n_output)

    def forward(self, x):
        out = self.fc1(x)
        out = F.relu(out)
        out = self.fc2(out)
        out = F.relu(out)
        out = self.fc3(out)
        out = F.relu(out)
        out = self.fc4(out)
        return out
class MLP4Layer(nn.Module):
    def __init__(self, n_input, n_hidden_1, n_hidden_2,n_hidden_3,n_hidden_4,n_output):
        super(MLP4Layer, self).__init__()
        self.fc1 = nn.Linear(n_input, n_hidden_1)
        self.fc2 = nn.Linear(n_hidden_1, n_hidden_2)
        self.fc3 = nn.Linear(n_hidden_2, n_hidden_3)
        self.fc4 = nn.Linear(n_hidden_3, n_hidden_4)
        self.fc5 = nn.Linear(n_hidden_4, n_output)

    def forward(self, x):
        out = self.fc1(x)
        out = F.relu(out)
        out = self.fc2(out)
        out = F.relu(out)
        out = self.fc3(out)
        out = F.relu(out)
        out = self.fc4(out)
        out = F.relu(out)
        out = self.fc5(out)
        return out
class MLP5Layer(nn.Module):
    def __init__(self, n_input, n_hidden_1, n_hidden_2,n_hidden_3,n_hidden_4,n_hidden_5,n_output):
        super(MLP5Layer, self).__init__()
        self.fc1 = nn.Linear(n_input, n_hidden_1)
        self.fc2 = nn.Linear(n_hidden_1, n_hidden_2)
        self.fc3 = nn.Linear(n_hidden_2, n_hidden_3)
        self.fc4 = nn.Linear(n_hidden_3, n_hidden_4)
        self.fc5 = nn.Linear(n_hidden_4, n_hidden_5)
        self.fc6 = nn.Linear(n_hidden_5, n_output)

    def forward(self, x):
        out = self.fc1(x)
        out = F.relu(out)
        out = self.fc2(out)
        out = F.relu(out)
        out = self.fc3(out)
        out = F.relu(out)
        out = self.fc4(out)
        out = F.relu(out)
        out = self.fc5(out)
        out = F.relu(out)
        out = self.fc6(out)
        return out
