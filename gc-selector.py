import os
import numpy as np
from numpy.linalg import inv 
import scipy.sparse.linalg as sp
import pandas as pd
from sklearn.decomposition import PCA

from argparse import ArgumentParser as parser
from configparser import ConfigParser

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from rpy2.robjects import pandas2ri

from mlpca import MLPCA



def init():
    # io 
    global inputdir, inputfile, outputdir, outputfile, separator, na_flag, copy_initial_data   
    # dataset description
    global filters, filters_error, GC_flag_column_name
    # dimensionality reduction
    global method_dim_reduction, N_dim_reduction
    # matching
    global N_dim_matching, distance, discard, method_matching, ratio, replace, weight_cut
    
    flags = parser(description="This script selects extragalactic globular cluster candidates from an input\
                                multi-band photometry catalog. It uses dimensionaility reduction and statistical\
                                matching to do so.")

    flags.add_argument('-i', help='The name of the .ini file.',
                       metavar="params.ini", default="params.ini")

    args = flags.parse_args()
    params_input = args.i

    if not os.path.exists(params_input):
        print("Parameters file not found:", params_input)
        exit(0)

    # loading parameters
    config = ConfigParser()
    config.read(params_input)

    inputdir = config.get('io', 'input-directory')          
    outputdir = config.get('io', 'output-directory')          
    inputfile = config.get('io', 'input-filename')          
    outputfile = config.get('io', 'output-filename')          
    separator = config.get('io', 'separator')            
    na_flag = config.get('io', 'na-flag')            
    copy_initial_data = config.getint('io', 'copy-initial-data')

    filters = config.get('data-description', 'filters').split(' ')
    filters_error = config.get('data-description', 'filters-error').split(' ')
    GC_flag_column_name = config.get('data-description', 'GC-flag-column-name')

    method_dim_reduction = config.get('auxiliary-space', 'method-dim-reduction')
    N_dim_reduction = config.getint('auxiliary-space', 'N-dim-reduction')

    N_dim_matching = config.getint('matching', 'N-dim-matching')
    distance = config.get('matching', 'distance')
    discard = config.get('matching', 'discard')
    method_matching = config.get('matching', 'method-matching')
    ratio = config.getint('matching', 'ratio')
    replace = config.getint('matching', 'replace')
    weight_cut = config.getint('matching', 'weight-cut')

    if method_dim_reduction == 'UMAP':
        # initialize umap main parameters
        global n_neighbors, min_dist, metric_umap
        global save_umap, umap_obj_filename
        
        n_neighbors = config.getint('auxiliary-space', 'n-neighbors')
        min_dist = config.getfloat('auxiliary-space', 'min-dist')
        metric_umap = config.get('auxiliary-space', 'metric-umap')
        
        save_umap = config.getint('io', 'save-umap-object')
        umap_obj_filename = config.get('io', 'umap-object-filename')

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)


'''
not being used for now
def matching():
    intended to define matching procedures not available via MatchIt
'''

def main():
    
    # INITIALIZE PARAMETERS
    print('Importing parameters')
    init()    

    # IMPORT DATA
    print('Importing data')
    data = pd.read_csv(inputdir + inputfile, sep=separator, na_values=na_flag) 
    
    print('Preparing arrays')
    X_ready = data[filters].to_numpy()
    Xsd_ready = data[filters_error].to_numpy() 
    
    # AUXILIARY SPACE
    # creating the auxiliary space using dimensionality reduction:
    if method_dim_reduction == 'MLPCA':
        # Data de-noising using maximum likelihood PCA:
        print('Computing MLPCs')
        U, S, V = MLPCA(X=X_ready, Xsd=Xsd_ready, p=len(filters)-1)
        cleanData = U @ S @ V

        cleanData = pd.DataFrame(data=cleanData, columns=filters)

        # Principal component analysis
        pca = PCA()
        new_dim = pca.fit_transform(cleanData)

        list_new_dim_str = [f'MLPC{i+1}' for i in range(new_dim.shape[1])]

    elif method_dim_reduction == 'PCA':
        # Principal component analysis
        print('Computing PCs')
        pca = PCA()
        new_dim = pca.fit_transform(X_ready)

        list_new_dim_str = [f'PC{i+1}' for i in range(new_dim.shape[1])]
    
    elif method_dim_reduction == 'UMAP':
        import umap 
        import joblib

        print('Computing UMAP space')
        new_space = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=N_dim_reduction, metric=metric_umap).fit(X_ready)
        new_dim = new_space.embedding_

        list_new_dim_str = [f'UMAP{i+1}' for i in range(new_dim.shape[1])]
        
        if save_umap == True: 
            print('Saving UMAP model object')
            joblib.dump(new_space, outputdir+umap_obj_filename)

    else:
        exit("No valid method to reduce dimensionality was selected!")
   

    # MATCHING 
    print('Importing MatchIt library from R') 
    MatchIt = importr('MatchIt')

    r_formula_str = f"{GC_flag_column_name} ~"
    for i in range(N_dim_matching):
       r_formula_str = r_formula_str + ' ' + list_new_dim_str[i]
       if i < (N_dim_matching - 1):
           r_formula_str = r_formula_str + ' +'

    df_for_matching = pd.DataFrame(data=new_dim[:,:N_dim_matching], columns=list_new_dim_str[:N_dim_matching])
    df_for_matching[GC_flag_column_name] = np.where(data[GC_flag_column_name] > 0, 1, 0)

    with (ro.default_converter + pandas2ri.converter).context():
        r_df_for_matching = ro.conversion.get_conversion().py2rpy(df_for_matching)
    
    print('Performing matching via MatchIt')
    matchItOut = MatchIt.matchit(Formula(r_formula_str), data=r_df_for_matching, distance=distance, discard=discard, method=method_matching, replace=bool(replace), ratio=ratio)
    
    # index 1 corresponds to the matching weights vector computed by MatchIt
    weights = np.array(matchItOut[1])

    new_labels = []
    for i in range(len(data)):
        if i in data[data['GCs'] > 0].index:
            new_labels.append('Confirmed')
        elif weights[i] > weight_cut:
            new_labels.append('Candidate')
        else:
            new_labels.append('No label')
    
    print('\n' + str(len(np.where(weights > 1)[0])) + ' candidates were selected!\n')
   

    # OUTPUT
    print('Preparing output')
    if copy_initial_data == 1:
        # add columns with coordinates of the objects in the new dimensions to the original dataset
        data[list_new_dim_str] = new_dim
        # add matching informations
        data['weights'] = weights
        # add new_labels to the original dataset
        data['labels'] = new_labels
        # write output 
        data.to_csv(outputdir + outputfile, sep=separator, na_rep=na_flag, index=False)

    elif copy_initial_data == 0:
        output_dataframe = pd.DataFrame(data=new_dim, columns=list_new_dim_str)
        output_dataframe['weights'] = weights 
        output_dataframe['labels'] = new_labels
        
        output_dataframe.to_csv(outputdir + outputfile, sep=separator, na_rep=na_flag, index=False)
       
    else:
        exit('No valid output setting was selected!')


    print('Finished!')



if __name__ == '__main__':
    main()
    
    
