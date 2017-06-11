# Collection of useful cheminformatics scripts
### anymol2vina.py
Script that automatically docks molecules with Autodock Vina.

### PCA-k-means.py
Script for k-means clustering of large sets of compounds (millions).  
How it works?:  

* Calculation of Morgan fingerprints
* Fitting of PCA model on smaller randomly selected set of compounds
* Transforming all FPs to fitted PCA space
* [Mini Batch K-Means clustering](http://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html)

Notes:

* Memory usage and compute time grows with number of PCA components
* 10 components and 4M compounds takes ~1 hour on i5 CPU
* Can work with 16GB of RAM, 32GB is better

##### Usage
```
K-means clustering of large sets of compounds with intermediate PCA step.
Creates 3 files: infile.fp.gz, infile.pca_train_sample.fp.gz and results file
infile.cluster.csv

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input smi.gz file to cluster.
  -n N_CLUSTERS, --n_clusters N_CLUSTERS
                        Number of clusters.
  -p PCA_TRAIN_SIZE, --pca_train_size PCA_TRAIN_SIZE
                        Size of training set for PCA (randomly selected from
                        input).
  -j JOBS, --jobs JOBS  Number of cores to use for FP calculation. Default: 1
  -c PCA_COMPONENTS, --pca_components PCA_COMPONENTS
                        Number of PCA components to calculate. Default: 10
  -r RADIUS, --radius RADIUS
                        Morgan fingerprint bit radius. Default: 1
  -b BIT_LENGTH, --bit_length BIT_LENGTH
                        Morgan fingerprint bit length. Default: 1024
```