import torch
from model import PINN

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
model = PINN(input_dim=2, output_dim=1).to(DEVICE)
optimiser = torch.optim.Adam(model.parameters(), lr=1e-3)
print(f'Training on: {DEVICE}')
