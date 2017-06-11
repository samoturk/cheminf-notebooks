#!/usr/bin/env python

__author__ = "Samo Turk"
__license__ = "BSD 3-clause"

import gzip
from random import sample
import sys
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from joblib import Parallel, delayed
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans

def get_fp_from_smi(smi, r, n=1024):
    smi, name = smi.split()[:2]
    m = Chem.MolFromSmiles(smi)
    if m is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(m,r,n).ToBitString()
        return fp + "\t" + str(Chem.MolToSmiles(m, isomericSmiles=True)) +"\t" + str(name)
    else:
        return None

def arg_parser():
    parser = argparse.ArgumentParser(description='K-means clustering of large sets of compounds with intermediate PCA step. \
                                                 Creates 3 files: infile.fp.gz, infile.pca_train_sample.fp.gz and results file \
                                                 infile.cluster.csv')
    parser.add_argument('-i', '--infile', help = "Input smi.gz file to cluster.")
    parser.add_argument('-n', '--n_clusters', help = "Number of clusters.", type=int)
    parser.add_argument('-p', '--pca_train_size', help = "Size of training set for PCA (randomly selected from input).", type=int)
    parser.add_argument('-j', '--jobs', default = 1, help = "Number of cores to use for FP calculation. Default: 1", type=int)
    parser.add_argument('-c', '--pca_components', default = 10, help = "Number of PCA components to calculate. Default: 10", type=int)
    parser.add_argument('-r', '--radius', default = 1, help = "Morgan fingerprint bit radius. Default: 1", type=int)
    parser.add_argument('-b', '--bit_length', default = 1024, help = "Morgan fingerprint bit length. Default: 1024", type=int)
    return parser

if __name__ == "__main__":
    parser = arg_parser()
    if len(sys.argv) == 1:
        argv = ['-h']
    else:
        argv = sys.argv[1:]
    args = parser.parse_args(argv)
    if not all([args.infile, args.n_clusters, args.pca_train_size]):
        print('Not enough arguments. -i, -n and -p have to be provided')
        sys.exit()

    infile = args.infile
    fp_file = '.'.join(infile.split('.')[:-1] + ['fp.gz'])
    pca_train_fp_file = '.'.join(infile.split('.')[:-1] + ['pca_train_sample.fp.gz'])
    out_file = '.'.join(infile.split('.')[:-1] + ['clusters.csv'])
    b_size = 100000  # 100k chunks fot pca.transform(), otherwise we get MemoryError

    # Open File
    print('Reading input smiles.')
    with gzip.open(infile, 'rt') as f:
        lines = f.readlines()

    # Calculate FPs
    print('Calculating FPs.')
    result = Parallel(n_jobs=args.jobs, verbose=1)(delayed(get_fp_from_smi)(x, args.radius, args.bit_length) for x in lines)

    # Save FPs including a random sample to train the PCA
    print('Saving FPs.')
    random_sample = set(sample(range(len(result)), args.pca_train_size))
    f_fp = gzip.open(fp_file, "wt")
    f_random_sample = gzip.open(pca_train_fp_file, "wt")
    for i, x in enumerate(result):
        if x is not None:
            f_fp.write(x + "\n")
            if i in random_sample:
                f_random_sample.write(x + "\n")
    f_fp.close()
    f_random_sample.close()

    # Read PCA training data
    print('Reading PCA training data.')
    with gzip.open(pca_train_fp_file, 'rt') as f:
        lines = f.readlines()

    # Fit PCA model
    print('Fitting PCA model on training data: %i samples.' % len(lines))
    X_pca = np.array([np.fromstring(x.split()[0], dtype=np.uint8) - ord('0') for x in lines])
    pca = PCA(n_components=args.pca_components)
    pca.fit(X_pca)

    # Read all FP data and project it in our PCA space
    print('Reading all FPs.')
    with gzip.open(fp_file, 'rt') as f:
        lines = f.readlines()

    print('Projecting FPs to our PCA space. %i samples.' % len(lines))
    X_clustering = []
    processed = 0
    # pca.transform() has to be done in chunks to avoid MemoryError
    while processed < len(lines):
        X_pca = np.array(
            [np.fromstring(x.split()[0], dtype=np.uint8) - ord('0') for x in lines[processed:processed + b_size]])
        X_clustering += list(pca.transform(X_pca))
        processed += b_size

    # Perform clustering
    print('Performing clustering.')
    clust_pc = MiniBatchKMeans(n_clusters=args.n_clusters, batch_size=args.n_clusters)
    clusters = clust_pc.fit_predict(X_clustering)

    # Save results
    print('Saving final results.')
    out_f = open(out_file, "w")
    out_f.write(','.join(['Smiles', 'ID', 'Cluster_ID']) + '\n')
    for l, c in zip(lines, clusters):
        fp, smi, name = l.split()
        out_f.write(','.join([smi, name, str(c)]) + '\n')
    out_f.close()