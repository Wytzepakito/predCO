#!/usr/bin/python3

import gzip
import os
import sys
import re
import numpy as np
import prediction_v4_module as pr 
import pandas as pd 
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn import linear_model
from sklearn import tree


def read_features(f_handle, label):
    
    df = pd.read_csv(f_handle, sep="\t", header=0, index_col=0)
    
    X = np.array(df)
    
    Y = np.repeat(label, X.shape[0])
    
    names = list(df.columns)
    
    return X, Y, names


if __name__=="__main__":

    f_path = sys.argv[1]
    
    #path = "/mnt/scratch/demir004/prediction/data/"
    #sp_path ="rice/sampled_gene_bw_opt/prediction_with_motifs"

    Xpos, Ypos, names1 = read_features(f_path + "positive.features.txt", 1) # sys.argv[1]

    Xneg, Yneg, names2 = read_features(f_path + "negative.features.txt", 0) # sys.argv[2]

    if names1 == names2:
        
        names = names1
        
    X_all = np.r_[Xpos, Xneg]

    Y = np.r_[Ypos, Yneg]


    RANDOM_SEED = 123
    PRNG = np.random.RandomState(RANDOM_SEED)

    outfile_root = f_path + "D_"  # sys.argv[3]

    pr.plot_ind_histograms(X_all, names, outfile_root + 'ind_histograms.pdf')

    X_all, names, adjpvals = pr.apply_ttest(X_all, names, outfile_root + "adj_pvals_features.txt")

    X_scaled = preprocessing.scale(X_all)

    pr.make_scatter_plots(X_scaled, names, outfile_root + 'scatter_plots.pdf')

    data = np.zeros((3,4))

    dfper = pd.DataFrame(data, columns=['accuracy', 'precision', 'recall', 'roc_auc'], index=['DT', 'LR', 'RF']) # performance table

    clf = linear_model.LogisticRegressionCV(refit=True, random_state=PRNG)

    LR_imp, dfper = pr.build_model(clf, X_scaled, Y, names, "LR", dfper, outfile_root  + "LR_features.txt")

    clf = tree.DecisionTreeClassifier(class_weight=None, criterion='gini', max_depth=None,
                max_features=None, max_leaf_nodes=None,
                min_impurity_split=1e-07, min_samples_leaf=5,
                min_samples_split=2, min_weight_fraction_leaf=0.0,
                presort=False, random_state=PRNG, splitter='best')

    DT_imp, dfper = pr.build_model(clf, X_scaled, Y, names, "DT", dfper, outfile_root + "DT_features.txt")

    clf = RandomForestClassifier(n_estimators=1000, random_state=PRNG)

    #RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
    #			max_depth=None, max_features='auto', max_leaf_nodes=None,
    #			min_impurity_split=1e-07, min_samples_leaf=1,
    #			min_samples_split=2, min_weight_fraction_leaf=0.0,
    #			n_estimators=10, n_jobs=1, oob_score=False, random_state=None,
    #			verbose=0, warm_start=False)
    ## maximum features to look for a split is sqrt(n_features). 

    RF_imp, dfper = pr.build_model(clf, X_scaled, Y, names, "RF", dfper, outfile_root + "RF_features.txt")

    pr.plot_ROC(X_scaled, Y, outfile_root + 'ROC.png', PRNG)

    dfper.to_csv(path_or_buf= outfile_root + "performance.txt", sep=',')

    direction = pr.find_direction(X_all)

    ## feature importance table for different methods
    dffeat = pd.DataFrame({'-log(p)' : -np.log(adjpvals), 
                        'LR' : LR_imp,
                        'DT' : DT_imp,
                        'RF' : RF_imp,
                        'effect' : direction}, index= names)

    dffeat.to_csv(path_or_buf= outfile_root + "features_report.txt" , sep=',', header=True, index=True)