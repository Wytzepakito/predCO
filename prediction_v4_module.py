#!/usr/bin/python3

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pydotplus
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn import model_selection 
from sklearn import svm
from sklearn import tree
from sklearn import preprocessing
from sklearn import metrics
from sklearn import linear_model
from IPython.display import Image  
from sklearn import metrics
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from statsmodels.sandbox.stats.multicomp import multipletests
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interp
from itertools import cycle
from sklearn.model_selection import StratifiedKFold

## cross_val_score uses model_selection.StratifiedKFold to split the data 
## and it also uses several scoring method from here: 
## http://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter


def read_data(path, datatype, feat_name_file=""):

	" This function reads from different files depend on the datatype and returns a matrix with data (X_all), their classes (Y) and the feature names (names)"

	# possible values for datatype is "A_" "B_" or "C_" or "Bg_"

	f_seq_pos = path + "positive.all.4000.maxoverlap1kb.fasta.fromSeq.txt"
	f_seq_neg = path + "negative.all.4000.maxoverlap1kb.fasta.fromSeq.txt"

	with open(f_seq_pos, 'r') as f:
		num_cols = len(f.readline().split())
		f.seek(0)
		X_seq_pos = np.loadtxt(f, delimiter="\t", usecols=range(1,num_cols))
	
	with open(f_seq_neg, 'r') as f:
		num_cols = len(f.readline().split())
		f.seek(0)
		X_seq_neg = np.loadtxt(f, delimiter="\t", usecols=range(1,num_cols))

	#X_seq_pos = np.loadtxt(f_seq_pos, delimiter="\t", usecols=range(1,29)) # 0th row:names  
	#X_seq_neg = np.loadtxt(f_seq_neg, delimiter="\t", usecols=range(1,29)) # 0th row:names,16 nuc + 12 Dnastruct


	if datatype == "B_" or datatype == "C_" or datatype == "Bg_": 

		f_ann_pos = path + "positive.all.4000.maxoverlap1kb.bed.fromAnn.txt"
		f_ann_neg = path + "negative.all.4000.maxoverlap1kb.bed.fromAnn.txt"
			   
		with open(f_ann_pos, 'r') as f:
			num_cols = len(f.readline().split())
			f.seek(0)
			X_ann_pos = np.loadtxt(f, delimiter="\t", usecols=range(1,num_cols))
	
		with open(f_ann_neg, 'r') as f:
			num_cols = len(f.readline().split())
			f.seek(0)
			X_ann_neg = np.loadtxt(f, delimiter="\t", usecols=range(1,num_cols))
	
	#	X_ann_pos = np.loadtxt(f_ann_pos, delimiter="\t", usecols=range(1,n)) # 0th row:names  
	#	X_ann_neg = np.loadtxt(f_ann_neg, delimiter="\t", usecols=range(1,n)) # 0th row:names, exon + some repeats


	if datatype == "C_" : 

		f_snpindel_pos = path + "positive.all.4000.maxoverlap1kb.bed.snpindel.txt"
		f_snpindel_neg = path + "negative.all.4000.maxoverlap1kb.bed.snpindel.txt"

		with open(f_snpindel_pos, 'r') as f:
			num_cols = len(f.readline().split())
			f.seek(0)
			X_snpindel_pos = np.loadtxt(f, delimiter="\t", usecols=range(1,num_cols))
	
		with open(f_snpindel_neg, 'r') as f:
			num_cols = len(f.readline().split())
			f.seek(0)
			X_snpindel_neg = np.loadtxt(f, delimiter="\t", usecols=range(1,num_cols))

	#	X_snpindel_pos = np.loadtxt(f_snpindel_pos, delimiter="\t", usecols=range(1,6)) # 0th row:names
	#	X_snpindel_neg = np.loadtxt(f_snpindel_neg, delimiter="\t", usecols=range(1,6)) # 0th row:names, 5 feature




	#!!!! ORDER IS IMPORTANT

	if datatype == "A_" :
	
		f_list = [f_seq_pos]

	if datatype == "B_" or datatype == "Bg_": 
	
		f_list = [f_seq_pos, f_ann_pos]
	
	if datatype == "C_" : 
	
		f_list = [f_seq_pos, f_ann_pos, f_snpindel_pos]
	

	names = []

	for fi in f_list:

		f = open(fi, 'r')

		for line in f: 

			line = line.strip();

			if line.startswith('#'):
			
				names.extend(line.split('\t')[1:])	  

			else:

				break

		f.close()


	if datatype == "A_" :
	
		X_all_pos = X_seq_pos

		X_all_neg = X_seq_neg

	if datatype == "B_" or datatype == "Bg_": 
	
		X_all_pos = np.c_[X_seq_pos, X_ann_pos]

		X_all_neg = np.c_[X_seq_neg, X_ann_neg]
	
	if datatype == "C_" : 
	
		X_all_pos = np.c_[X_seq_pos, X_ann_pos, X_snpindel_pos]

		X_all_neg = np.c_[X_seq_neg, X_ann_neg, X_snpindel_neg]


	X_all = np.r_[X_all_pos, X_all_neg]


	a = np.repeat(1, X_all_pos.shape[0]) ## in F1 score; pos_label : str or int, 1 by default

	b = np.repeat(0, X_all_neg.shape[0])

	Y = np.r_[a,b]


	## specify the features ##########

	if datatype == "Bg_": 

		general_feats = list(pd.read_csv(filepath_or_buffer= feat_name_file , sep=',', header=None).iloc[:,0])

		general_feats_id = [names.index(i) for i in general_feats ]

		X_all = X_all[:,general_feats_id]

		names = general_feats

	
	return X_all, Y, names


def plot_ind_histograms(X_all, names, outfile):

	"plots many histograms and writes to path"

	with PdfPages(outfile) as pdf:

		for i in range(X_all.shape[1]):

			buffer = (X_all[:, i].max() - X_all[:, i].min()) / 100

			r_min, r_max = X_all[:, i].min() - buffer , X_all[:, i].max() + buffer

			bins = np.linspace(r_min, r_max, 100)

			m = int(X_all.shape[0] / 2)
			plt.hist(X_all[:(m-1), i], bins, alpha=0.5, label='CO')
			plt.hist(X_all[m:, i], bins, alpha=0.5, label='random')
			plt.axvline(X_all[m:, i].mean(), color='green', linestyle='dashed', linewidth=2)
			plt.axvline(X_all[:(m-1), i].mean(), color='b', linestyle='dashed', linewidth=2)
			plt.legend(loc='upper right')
			plt.title(names[i])

			pdf.savefig()
			plt.close()



def apply_ttest(X_all, names, outfile):

	"performs t-test and updates the data. Returns (X_all, names) and the probabilities (adjpvals)"

	uselessFeatures =list()
	idX = list()


	for i in range(X_all.shape[1]):
	
		if np.unique(X_all[:,i]).shape[0] == 1:
		
			uselessFeatures.append(names[i])
			idX.append(i)

	### UPDATE X_all columns by removing uselessFeatures. 

	X_all = np.delete(X_all, idX, 1) 

	for j in uselessFeatures:
		names.remove(j)
	
	## Individual t-test for features

	uselessFeatures= list()

	tprob = list()
	t2prob = list()

	################## Scaled data has different pvals 

	for i in range(X_all.shape[1]):
	
		m = int(X_all.shape[0] /2 )
	
		pos = X_all[:(m-1), i]
	
		neg = X_all[m:, i]

		t, prob = stats.ttest_ind(pos, neg)
		  
		if np.isnan(prob):
		
			print(names[i])

	#		tprob.append(1)
		
		else:
		
			tprob.append(prob)
	
		t, prob = stats.ttest_ind(pos, neg, equal_var = False)
	
		if np.isfinite(prob):
		
			t2prob.append(prob)
	#	else:
	#		t2prob.append(1)


	## multiple-testing for tprob


	a, adjpvals, c, d = multipletests(tprob, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

	tprob = adjpvals

	imp = np.c_[names, tprob] 

	df = pd.DataFrame(imp, columns=['names', 'tprob'])

	df[['tprob']] = df[['tprob']].apply(pd.to_numeric)


	dfs = df.sort_values('tprob', axis=0, ascending=False)

	dfs = dfs.set_index(['names'])

	dfs = dfs.sort_values('tprob', axis=0, ascending=True)

	dfs.to_csv(path_or_buf= outfile, sep=',')


	return X_all, names, adjpvals



def find_direction(X_all):

	direction = list()

	for i in range(X_all.shape[1]):
	
		m = int(X_all.shape[0] /2 )
	
		pos = X_all[:(m-1), i]
	
		neg = X_all[m:, i] 
	
		di =  (np.mean(pos) - np.mean(neg)) / np.std(X_all[:, i])
	
		direction.append(di)

	return direction



def make_scatter_plots(X_scaled, names, outfile):

	" pairwise scatters of features and write it to path"

	m = int(X_scaled.shape[0] / 2)

	cs1= ["red"] *m   # CO 

	cs2 = ["blue"] *m	# random

	cs = cs1 + cs2

	n = X_scaled.shape[1]


	with PdfPages(outfile) as pdf:

		for i in range(n-1):

			j = i+1

			while j <= n-1: 

				plt.scatter(X_scaled[:, i], X_scaled[:, j], color=cs)  

				plt.xlabel(names[i])

				plt.ylabel(names[j])

				pdf.savefig()

				plt.close()

				j = j + 1



def plot_tree_model_DT(X_all, Y, names, path, datatype, PRNG):

	" Decision Tree 3 layer tree on ALL data - non-scaled"

	clf = tree.DecisionTreeClassifier(class_weight=None, criterion='gini', max_depth=3,
				max_features=None, max_leaf_nodes=None,
				min_impurity_split=1e-07, min_samples_leaf=5,
				min_samples_split=2, min_weight_fraction_leaf=0.0,
				presort=False, random_state=PRNG, splitter='best') # min_samples_split is common in literature

	clf = clf.fit(X_all, Y) 

	dot_data = tree.export_graphviz(clf, out_file=None, 
							 feature_names=names,  
							 class_names=(["random", "CO"]),  
							 filled=True, rounded=True,  
							 special_characters=True)  
	graph = pydotplus.graph_from_dot_data(dot_data) 

	graph.write_png(path + datatype + "decision_tree.png")


def build_model(clf, X_scaled, Y, names, model, dfper, outfile):

	
	for smethod in ['accuracy', 'precision', 'recall', "f1", 'roc_auc']:

		scores = cross_val_score(clf, X_scaled, Y, cv=10, scoring= smethod)

		dfper.loc[model, smethod] = "%0.2f (s= %0.2f)" % (scores.mean(), scores.std())
	

	clf = clf.fit(X_scaled, Y)

	if model == "LR":

		imp = abs(clf.coef_[0,:])

	else:

		imp = clf.feature_importances_

	df= pd.DataFrame(np.c_[names, imp])

	df[[1]] = df[[1]].apply(pd.to_numeric)

	df = df.sort_values(1, axis=0, ascending=True)

	df.columns=['Features', 'importances']

	df = df.set_index(['Features'])

	dfs = df.sort_values('importances', axis=0, ascending=True)

	dfs.to_csv(path_or_buf= outfile, sep=',')

	return imp , dfper



def plot_ROC(X, y, outfile_name, PRNG):

	"ROC curve summary for the classifiers - with cross-validation"

	cv = StratifiedKFold(n_splits=10)

	clf1 = tree.DecisionTreeClassifier(class_weight=None, criterion='gini', max_depth=None,
				max_features=None, max_leaf_nodes=None,
				min_impurity_split=1e-07, min_samples_leaf=5,   ## min 5 in each leaf
				min_samples_split=2, min_weight_fraction_leaf=0.0,
				presort=False, random_state=PRNG, splitter='best') # min_samples_split is common in literature

	clf2 = linear_model.LogisticRegressionCV(refit=True, random_state=PRNG)

	clf3 = RandomForestClassifier(n_estimators=1000, random_state=PRNG)

	clfnames = ['DT', 'LR', 'RF']

	colors = cycle(['blue', 'green', 'darkorange'])

	plt.figure(figsize=(4,3), dpi=300)

	for clf, cnames, color in zip([clf1, clf2, clf3], clfnames, colors):

		mean_tpr = 0.0
		mean_fpr = np.linspace(0, 1, 100)
	
		lw = 1.5

		i = 0
		for (train, test), color in zip(cv.split(X, y), colors):
		
			probas_ = clf.fit(X[train], y[train]).predict_proba(X[test])
		
			# Compute ROC curve and area the curve
			fpr, tpr, thresholds = metrics.roc_curve(y[test], probas_[:, 1])
			mean_tpr += interp(mean_fpr, fpr, tpr)
			mean_tpr[0] = 0.0
			roc_auc = metrics.auc(fpr, tpr)

			i += 1
	
		mean_tpr /= cv.get_n_splits(X, y)
		mean_tpr[-1] = 1.0
		mean_auc = metrics.auc(mean_fpr, mean_tpr)
		plt.plot(mean_fpr, mean_tpr, color=color, 
				 label='%s (%0.2f)' % (cnames, mean_auc), lw=lw)


	font = {'size'   : 12,
			'weight' : 'normal',
			'family' : 'sans-serif'}

	plt.rc('font', **font)


	plt.plot([0, 1], [0, 1], linestyle='--', lw=lw, color='k')
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('FPR')
	plt.ylabel('TPR')

	plt.title("")
	plt.legend(loc="lower right", prop={'size':12})
	plt.savefig(outfile_name)



if __name__=="__main__":
  
	RANDOM_SEED = 123
	PRNG = np.random.RandomState(RANDOM_SEED)
	
	path = sys.argv[1]
	datatype = sys.argv[2]
	if datatype == "Bg_":
		featfile = sys.argv[3]
	else:
		featfile = ""

	X_all, Y, names = read_data(path, datatype, featfile)

	outfile_root = path + datatype

	plot_ind_histograms(X_all, names, outfile_root + 'ind_histograms.pdf')

	X_all, names, adjpvals = apply_ttest(X_all, names, outfile_root + "adj_pvals_features.txt")

	X_scaled = preprocessing.scale(X_all)

	make_scatter_plots(X_scaled, names, outfile_root + 'scatter_plots.pdf')

	data = np.zeros((3,4))

	dfper = pd.DataFrame(data, columns=['accuracy', 'precision', 'recall', 'roc_auc'], index=['DT', 'LR', 'RF']) # performance table


	clf = linear_model.LogisticRegressionCV(refit=True, random_state=PRNG)

	LR_imp, dfper = build_model(clf, X_scaled, Y, names, "LR", dfper, outfile_root  + "LR_features.txt")

	clf = tree.DecisionTreeClassifier(class_weight=None, criterion='gini', max_depth=None,
				max_features=None, max_leaf_nodes=None,
				min_impurity_split=1e-07, min_samples_leaf=5,
				min_samples_split=2, min_weight_fraction_leaf=0.0,
				presort=False, random_state=PRNG, splitter='best')

	DT_imp, dfper = build_model(clf, X_scaled, Y, names, "DT", dfper, outfile_root + "DT_features.txt")

	clf = RandomForestClassifier(n_estimators=1000, random_state=PRNG)

	#RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
	#			max_depth=None, max_features='auto', max_leaf_nodes=None,
	#			min_impurity_split=1e-07, min_samples_leaf=1,
	#			min_samples_split=2, min_weight_fraction_leaf=0.0,
	#			n_estimators=10, n_jobs=1, oob_score=False, random_state=None,
	#			verbose=0, warm_start=False)
	## maximum features to look for a split is sqrt(n_features). 

	RF_imp, dfper = build_model(clf, X_scaled, Y, names, "RF", dfper, outfile_root + "RF_features.txt")

	plot_ROC(X_scaled, Y, outfile_root + 'ROC.png', PRNG)

	dfper.to_csv(path_or_buf= outfile_root + "performance.txt", sep=',')

	direction = find_direction(X_all)

	## feature importance table for different methods
	dffeat = pd.DataFrame({'-log(p)' : -np.log(adjpvals), 
						  'LR' : LR_imp,
						  'DT' : DT_imp,
						  'RF' : RF_imp,
						  'effect' : direction}, index= names)

	dffeat.to_csv(path_or_buf= outfile_root + "features_report.txt" , sep=',', header=True, index=True)

