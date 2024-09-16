# GC-Selector
Python script to select extragalactic globular cluster candidates from multi-band photometry catalogs.
The code performs dimensionality reduction upon the input dataset, then statistical matching is applied within the reduced space.

# Required python libraries
- scipy
- numpy
- pandas
- configparser
- scikit-learn (to use PCA)
- umap (to use UMAP)
- joblib (to save and load umap objects)
- psmpy (to use propensity score matching)

# Usage
The input and ouput paths, the method used to reduce dimensionality and the matching parameters are to be set in the text file `params.ini` (default name).
The code is then run:
```
python3 PyFNB.py -i <name_of_the_parameters_file>.ini
```

# Acknowledgement
This selection procedure was mainly developed by [Rafael S. de Souza](https://scholar.google.com/citations?user=ozs9GcgAAAAJ).
