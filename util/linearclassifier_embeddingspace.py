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

import logging

class LinearClassifier(nn.Module):
    """ Linear classifier """
    def __init__(self, input_dim, num_classes):
        super(LinearClassifier, self).__init__()
        self.fc = nn.Linear(input_dim, num_classes)

    def forward(self, x):
        """ linear layer forward """
        return self.fc(x)


def train_linear_classifier(indim, whole_dataloader, train_loader, test_loader, num_classes, logger, outdir, names):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    classifier = LinearClassifier(input_dim=indim, \
        num_classes=num_classes+1).to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer_lc = torch.optim.Adam(classifier.parameters(), lr=0.001)

    # Train classifier
    for epoch in range(300):  # Number of epochs
        classifier.train()
        for embedding, labels in train_loader:
            embedding, labels = embedding.to(device), labels.to(device)
            # Train classifier on the representations
            outputs = classifier(embedding)
            loss = criterion(outputs, labels)

            optimizer_lc.zero_grad()
            loss.backward()
            optimizer_lc.step()
        # print(epoch+1, "=", loss.detach())
        logger.info(f"Epoch {epoch+1}, Loss: {loss.detach().item()}")

    # Evaluate classifier
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

    logger.info(f'Accuracy of the network on the \
        {len(test_loader)*4096} test data: {100 * correct / total} %')

    classifier.eval()
    correct = 0
    total = 0
    predicted_labels = []
    probabilities = []
    with torch.no_grad():
        for embedding, labels in whole_dataloader:
            embedding, labels = embedding.to(device), labels.to(device)
            outputs = classifier(embedding)
            _, predicted = torch.max(outputs.data, 1)
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
            prob = nn.functional.softmax(outputs.data, dim=1)
            prob, _ = torch.max(prob, 1)
            predicted_labels.extend(predicted.detach().cpu().numpy())
            probabilities.extend(prob.detach().cpu().numpy())
    # print(predicted_labels, 'predicted labels')
    assignment = np.vstack(predicted_labels).flatten()
    np.save(os.path.join(outdir, "assignment.npy"),assignment)
    binids = pd.DataFrame(np.vstack((names, assignment)).T)
    binids.to_csv(os.path.join(outdir, 'linearclassifier_binids'), header=None,sep='\t', index=False)
    probabilities = pd.DataFrame(np.vstack(probabilities).flatten())
    probabilities.to_csv(os.path.join(outdir,'probabilities'), header=None, sep='\t')
    logger.info(f'Accuracy of the network on the \
        {len(whole_dataloader)*4096} whole data: {100 * correct / total} %')

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
    filename=args.outdir + '/classifier.log', filemode='w')
    args.logger = logging.getLogger()

    args.otuids = pd.read_csv(args.otuids, header=None)
    unique_otu_ids = args.otuids[0].unique()
    otu_mapping = {otu_id: idx for idx, otu_id in enumerate(unique_otu_ids)}
    args.otuids[1] = args.otuids[0].map(otu_mapping)
    labels = args.otuids[1].to_numpy()

    dataset = TensorDataset(torch.from_numpy(args.latent), torch.from_numpy(labels))

    train_size = int(args.latent.shape[0] * 0.8)
    test_size = int(args.latent.shape[0] - train_size)
    dataset_train, dataset_test = torch.utils.data.random_split(dataset,[train_size,test_size])
    train_loader = DataLoader(dataset_train, batch_size=4096, shuffle=True,\
            num_workers=4, pin_memory=True)
    test_loader = DataLoader(dataset_test,batch_size=4096, shuffle=False, \
            num_workers=4, pin_memory=True)

    whole_dataloader = DataLoader(dataset, batch_size=4096, shuffle=False,\
            num_workers=4, pin_memory=True)
    train_linear_classifier(args.latent.shape[1],whole_dataloader, \
        train_loader, test_loader, np.max(labels), args.logger, args.outdir, args.names)
    args.logger.info(f'{time.time()-start}, seconds to complete')
if __name__ == "__main__" :
    main()


### Note: COMEBIN embedding is not sorted by name. Hence, first sort by name and then use it for linear classification