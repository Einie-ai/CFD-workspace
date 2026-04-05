"""PINN Neural Network Architecture"""
import torch
import torch.nn as nn

class PINN(nn.Module):
    def __init__(self, input_dim=2, output_dim=1, hidden_layers=4, hidden_units=64):
        super().__init__()
        layers = [nn.Linear(input_dim, hidden_units), nn.Tanh()]
        for _ in range(hidden_layers - 1):
            layers += [nn.Linear(hidden_units, hidden_units), nn.Tanh()]
        layers.append(nn.Linear(hidden_units, output_dim))
        self.net = nn.Sequential(*layers)
        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x):
        return self.net(x)
