# GC-Selector
Python script to select extragalactic globular cluster candidates from multi-band photometry catalogs.
The code performs dimensionality reduction upon the input dataset, then statistical matching is applied within the reduced space.

# Requirements
Both python and R languages are required because the MatchIt R library is used for statistical matching.

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

# Usage
The input and ouput (io) configurations, the method used to reduce dimensionality and the matching parameters are to be set in the text file `params.ini` (default name).
The code is then run:
```
python3 gc-selector.py -i <name_of_the_parameters_file>.ini
```

# Acknowledgement
This selection procedure was mainly designed by [Rafael S. de Souza](https://scholar.google.com/citations?user=ozs9GcgAAAAJ).
