#!/usr/bin/env python
""" run classifier on the emedding training """

import os
import time
import argparse
import logging
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.utils.data.dataset import TensorDataset
from sklearn.model_selection import KFold
import logging

class LinearClassifier(nn.Module):
    """ Linear classifier """
    def __init__(self, input_dim, num_classes):
        super(LinearClassifier, self).__init__()
        self.fc = nn.Linear(input_dim, num_classes)

    def forward(self, x):
        """ linear layer forward """
        return self.fc(x)

def train_with_kfold(indim, dataset, num_classes, logger, outdir, names, k=5):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    kf = KFold(n_splits=k, shuffle=True, random_state=42)  # Define k-fold
    fold_accuracies = []

    for fold, (train_indices, test_indices) in enumerate(kf.split(dataset)):
        logger.info(f"Training on fold {fold + 1}/{k}")
        
        # Prepare train and test datasets
        train_subset = torch.utils.data.Subset(dataset, train_indices)
        test_subset = torch.utils.data.Subset(dataset, test_indices)
        
        train_loader = DataLoader(train_subset, batch_size=4096, shuffle=True, num_workers=4, pin_memory=True)
        test_loader = DataLoader(test_subset, batch_size=4096, shuffle=False, num_workers=4, pin_memory=True)
        
        # Initialize the model, criterion, and optimizer
        classifier = LinearClassifier(input_dim=indim, num_classes=num_classes + 1).to(device)
        criterion = nn.CrossEntropyLoss()
        optimizer_lc = torch.optim.Adam(classifier.parameters(), lr=0.001)
        
        # Train the classifier
        for epoch in range(300):  # Number of epochs
            classifier.train()
            for embedding, labels in train_loader:
                embedding, labels = embedding.to(device), labels.to(device)
                outputs = classifier(embedding)
                loss = criterion(outputs, labels)
                
                optimizer_lc.zero_grad()
                loss.backward()
                optimizer_lc.step()
        
            logger.info(f"Fold {fold + 1}, Epoch {epoch + 1}, Loss: {loss.detach().item()}")
        
        # Evaluate on the test set
        classifier.eval()
        correct = 0
        total = 0
        with torch.no_grad():
            for embedding, labels in test_loader:
                embedding, labels = embedding.to(device), labels.to(device)
                outputs = classifier(embedding)
                _, predicted = torch.max(outputs.data, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()
        
        accuracy = 100 * correct / total
        fold_accuracies.append(accuracy)
        logger.info(f"Fold {fold + 1} Accuracy: {accuracy:.2f}%")
    
    # Log the average accuracy across folds
    avg_accuracy = sum(fold_accuracies) / len(fold_accuracies)
    logger.info(f"Average Accuracy across {k} folds: {avg_accuracy:.2f}%")

def main() -> None:

    """ Assess embedding """
    start = time.time()
    parser = argparse.ArgumentParser(
        prog="mcdevol",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s \
        --latent --otuids --names --outdir [options]",
        add_help=False,
    )

    parser.add_argument("--latent", type=str, \
        help="latent space embedding", required=True)
    parser.add_argument("--otuids", type=str, \
        help="otuids of contigs", required=True)
    parser.add_argument("--names", type=str, \
        help="ids of contigs", required=True)
    parser.add_argument("--outdir", type=str, \
        help="output directory", required=True)

    args = parser.parse_args()
    args.latent = np.load(args.latent, allow_pickle=True)
    args.latent = args.latent.astype(np.float32)
    args.names = np.load(args.names, allow_pickle=True)

    args.outdir = os.path.join(args.outdir, '')
    print(args.outdir ,'output directory')

    try:
        if not os.path.exists(args.outdir):
            print('create output folder')
            os.makedirs(args.outdir)
    except RuntimeError as e:
        print(f'output directory already exist. Using it {e}')

    logging.basicConfig(format='%(asctime)s - %(message)s', \
    level=logging.INFO, datefmt='%d-%b-%y %H:%M:%S',
    filename=args.outdir + '/kfold_classifier.log', filemode='w')
    args.logger = logging.getLogger()

    args.otuids = pd.read_csv(args.otuids, header=None)
    unique_otu_ids = args.otuids[0].unique()
    otu_mapping = {otu_id: idx for idx, otu_id in enumerate(unique_otu_ids)}
    args.otuids[1] = args.otuids[0].map(otu_mapping)
    labels = args.otuids[1].to_numpy()

    dataset = TensorDataset(torch.from_numpy(args.latent), torch.from_numpy(labels))

    train_with_kfold(args.latent.shape[1], dataset, np.max(labels), args.logger, args.outdir, args.names, k=5)
    
    args.logger.info(f'{time.time()-start}, seconds to complete')
if __name__ == "__main__" :
    main()


### Note: COMEBIN embedding is not sorted by name. Hence, first sort by name and then use it for linear classification

