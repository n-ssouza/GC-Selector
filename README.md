# GC-Selector
Experimental python script to select extragalactic globular cluster candidates from multi-dimensional (usually photometric) data catalogs.
The code performs dimensionality reduction upon the input dataset, then statistical matching is applied within the reduced space to select candidates.

**This is by no means a definite or universal procedure to select extragalactic globular cluster candidates.** It is one of the ways to approach this problem.  

# Requirements
Both Python and R languages are required because the MatchIt R library is used for statistical matching.

**Python libraries required:**
- os
- numpy
- scipy
- pandas
- argparse
- configparser
- scikit-learn (to use PCA)
- umap (to use UMAP)
- joblib (to save and load umap objects)
- rpy2 (which requires R to be installed)

**R libraries required:**
- MatchIt

# Assumptions about the input dataset
This code assumes the input dataset to be a CSV file.
It also assumes the input CSV to have:
- columns with errors for each of the quantities of interest (in case you want to use maximum likelihood PCA),
- a column of integers to label the confirmed globular clusters (the training sample, the treatment units):
  - 0: not labeled
  - 1 (or any positive integer): confirmed globular cluster 

# Usage
The input-output (io) configurations, the method used to reduce dimensionality and the matching parameters are to be set in the text file `params.ini` (default name).
The code is then run:
```
python3 gc-selector.py -i <name_of_the_parameters_file>.ini
```

# Important notes
It is expected that the user is familiar with how MatchIt and UMAP (in case it's used) work.
Make sure to check the their documentations: [MatchIt¹](https://cran.r-project.org/web/packages/MatchIt/MatchIt.pdf), [MatchIt²](https://imai.fas.harvard.edu/research/files/matchit.pdf), [UMAP](https://umap-learn.readthedocs.io/en/latest/).

There are many other matching possibilities available via MatchIt whose respective parameters are not present in the parameters file of this code. If you want use some of these other possibilities, I encourage you to edit the code yourself and send pull requests with their implementations.

# Acknowledgement
This selection procedure was mainly designed by [Rafael S. de Souza](https://scholar.google.com/citations?user=ozs9GcgAAAAJ).
