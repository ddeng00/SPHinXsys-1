import argparse
from .utils import dir_path_validator
from pathlib import Path
import numpy as np
import pandas as pd
import wandb
from sklearn.model_selection import train_test_split

import torch
from torch.autograd import Variable
from torch.nn import Linear, ReLU, CrossEntropyLoss, Sequential, Conv2d, MaxPool2d, Module, Softmax, BatchNorm2d, Dropout
from torch.optim import Adam, SGD
from torch.utils.data import TensorDataset, DataLoader


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'dir', nargs='?', type=dir_path_validator, default=None)
    args, _ = parser.parse_known_args()

    # Set up wandb
    wandb.init(project="SPH-CNN", entity="ddeng00")
    wandb.config = {
        "epochs": 25,
        "batch_size": 32,
        "validation_split": 0.05
    }

    # Load data
    data = np.load(args.dir/'dataset.npz', mmap_mode='r')
    X = torch.from_numpy(data['X'])
    Y = torch.from_numpy(data['Y'])
    dataloader = DataLoader(TensorDataset(X, Y), batch_size=wandb.config.batch_size, shuffle=True)


if __name__ == '__main__':
    main()
