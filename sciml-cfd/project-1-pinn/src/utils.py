import numpy as np
import matplotlib.pyplot as plt

def sample_collocation(n, x_range=(0,1), t_range=(0,1)):
    x = np.random.uniform(*x_range, (n, 1))
    t = np.random.uniform(*t_range, (n, 1))
    return np.hstack([x, t])

def plot_loss(history, save_path=None):
    plt.figure(figsize=(8, 4))
    plt.semilogy(history)
    plt.xlabel('Epoch'); plt.ylabel('Loss')
    plt.title('PINN Training Loss')
    plt.tight_layout()
    if save_path: plt.savefig(save_path, dpi=150)
    plt.show()
