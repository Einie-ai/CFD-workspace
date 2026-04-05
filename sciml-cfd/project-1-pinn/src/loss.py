import torch

def boundary_loss(model, bc_pts, bc_vals):
    return torch.mean((model(bc_pts) - bc_vals) ** 2)
