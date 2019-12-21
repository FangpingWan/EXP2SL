import numpy as np
import pickle
import sys
from sklearn.model_selection import KFold, StratifiedKFold
from scipy.spatial.distance import jaccard
from github_model import *
from sklearn.metrics import roc_auc_score, average_precision_score
#tissue_set = {'A375', 'A549', 'HEK', 'HT29'}

def get_std(scores):
    fold_avg = []
    for i in range(5):
        fold_score = scores[i*5:(i+1)*5]
        fold_avg.append(np.mean(fold_score, axis=0))
    fold_avg = np.array(fold_avg)
    return np.std(fold_avg, axis=0)
    
def load_label(tissue, feature_dict):
	symbolA_list, symbolB_list, label_list = [], [], []
	selected = set()
	pos_gene = set()
	f = open('./GEMINI/gemini_'+tissue+'_labels.tsv')
	for line in f.readlines():
		symbolA, geneA, symbolB, geneB, label = line.strip().split('\t')
		if symbolA not in feature_dict or symbolB not in feature_dict:
			continue
		if label == '1':
			pos_gene.add(geneA)
			pos_gene.add(geneB)
		else:
			pass
		if symbolA+' ' +symbolB not in selected:
			symbolA_list.append(symbolA)
			symbolB_list.append(symbolB)
			label_list.append(int(label))
			selected.add(symbolA+' ' +symbolB)
		assert symbolB + ' ' + symbolA not in selected
	f.close()
	print ('number of samples', len(symbolA_list), len(symbolB_list), len(label_list))
	print ('positive', np.sum(label_list), 'negative', len(label_list)-np.sum(label_list))

	return np.array(symbolA_list), np.array(symbolB_list), np.array(label_list)


def load_feature_list(tissue, gene_list):
	if tissue == 'HEK':
		tissue = 'HEK293T'
	with open('./L1000/shRNA_cgs') as f:
		shrna_dict = pickle.load(f)
	feature_dict = {}
	for symbol, t in shrna_dict.keys():
		if t == tissue:
			feature_dict[symbol] = shrna_dict[(symbol,t)]
	feature_list = [feature_dict[gene] for gene in gene_list]
	return np.array(feature_list), feature_dict


def get_sl_mat(symbolA_list, symbolB_list, labels, gene_list, feature_dict):
	gene_to_idx = {gene_list[i]:i for i in range(len(gene_list))}
	sl_mat = np.zeros((len(gene_list), len(gene_list)))
	sl_mask = np.zeros((len(gene_list), len(gene_list)))
	sl_gene_set = set()
	for i in range(len(symbolA_list)):
		geneA = symbolA_list[i]
		geneB = symbolB_list[i]
		if geneA not in gene_to_idx or geneB not in gene_to_idx:
			continue
		sl_mat[gene_to_idx[geneA], gene_to_idx[geneB]] = labels[i]
		sl_mask[gene_to_idx[geneA], gene_to_idx[geneB]] = 1
		sl_gene_set.add(geneA)
		sl_gene_set.add(geneB)
	return sl_mat, sl_mask, sl_gene_set


def train_test_split_fixseed(sl_mat, sl_mask, fold, start_seed):
	seed = start_seed
	all_idx = np.where(sl_mask==1)
	idx = np.array(range(len(all_idx[0])))
	
	success = False
	while not success:
		success = True

		np.random.seed(seed)
		np.random.shuffle(idx)
		np.random.seed(seed)
		kf = KFold(n_splits=fold, shuffle=True, random_state = seed)
		train_mat_list, train_mask_list, test_mat_list, test_mask_list = [], [], [], []
		for train_i, test_i in kf.split(idx):
			train_mat = np.zeros(sl_mat.shape, dtype=sl_mat.dtype)
			train_mask = np.zeros(sl_mat.shape, dtype=sl_mat.dtype)
			test_mat = np.zeros(sl_mat.shape, dtype=sl_mat.dtype)
			test_mask = np.zeros(sl_mat.shape, dtype=sl_mat.dtype)
			train_mat[all_idx[0][train_i], all_idx[1][train_i] ] = sl_mat[all_idx[0][train_i], all_idx[1][train_i] ]
			train_mask[all_idx[0][train_i], all_idx[1][train_i] ] = 1
			test_mat[all_idx[0][test_i], all_idx[1][test_i] ] = sl_mat[all_idx[0][test_i], all_idx[1][test_i] ]
			test_mask[all_idx[0][test_i], all_idx[1][test_i] ] = 1
			print ('train-test split', 'train', train_mat.shape, train_mat.sum(), train_mask.sum(), 'test', test_mat.shape, test_mat.sum(), test_mask.sum())
			if train_mat.sum() == 0 or train_mat.sum() == train_mask.sum() or test_mat.sum() == 0 or test_mat.sum() == test_mask.sum():
				success = False
				seed += 1
				print ('re-splitting....')
				break
			train_mat_list.append(train_mat+train_mat.transpose())
			train_mask_list.append(train_mask+train_mask.transpose())
			test_mat_list.append(test_mat+test_mat.transpose())
			test_mask_list.append(test_mask+test_mask.transpose())
	seed += 1
	return train_mat_list, train_mask_list, test_mat_list, test_mask_list, seed


def leave_gene_split_fixseed(sl_mat, sl_mask, fold, start_seed):
	seed = start_seed

	all_idx = np.where(sl_mask==1)
	idx = np.array(range(sl_mat.shape[0]))
	
	success = False
	while not success:
		success = True
		np.random.seed(seed)
		np.random.shuffle(idx)
		np.random.seed(seed)
		kf = KFold(n_splits=fold, shuffle=True, random_state = seed)
		train_mat_list, train_mask_list, test_mat_list, test_mask_list = [], [], [], []
		for train_i, test_i in kf.split(idx):
			train_i_ = idx[train_i]
			test_i_ = idx[test_i]
			train_mat = np.zeros(sl_mat.shape, dtype=sl_mat.dtype)
			train_mask = np.zeros(sl_mat.shape, dtype=sl_mat.dtype)
			test_mat = np.zeros(sl_mat.shape, dtype=sl_mat.dtype)
			test_mask = np.zeros(sl_mat.shape, dtype=sl_mat.dtype)
			
			for i in range(len(all_idx[0])):
				x_i = all_idx[0][i]
				y_i = all_idx[1][i]
				
				if x_i in train_i_ and y_i in train_i_:
					train_mat[x_i, y_i] = sl_mat[x_i, y_i]
					train_mask[x_i, y_i] = 1
				elif x_i in test_i_ or y_i in test_i_:
					test_mat[x_i, y_i] = sl_mat[x_i, y_i]
					test_mask[x_i, y_i] = 1
				else:
					print ('error')
			print ('train-test split', 'train', train_mat.shape, train_mat.sum(), train_mask.sum(), 'test', test_mat.shape, test_mat.sum(), test_mask.sum())
			if train_mat.sum() == 0 or train_mat.sum() == train_mask.sum() or test_mat.sum() == 0 or test_mat.sum() == test_mask.sum():
				success = False
				seed += 1
				print ('re-splitting....')
				break
			train_mat_list.append(train_mat+train_mat.transpose())
			train_mask_list.append(train_mask+train_mask.transpose())
			test_mat_list.append(test_mat+test_mat.transpose())
			test_mask_list.append(test_mask+test_mask.transpose())
	seed += 1
	return train_mat_list, train_mask_list, test_mat_list, test_mask_list, seed

def masked_auc(test_mat, test_mask, test_pred):
	test_idx = np.where(test_mask==1)
	test_label = test_mat[test_idx]
	pred =test_pred[test_idx]
	#print 'test_label', test_label.shape, 'pred', pred.shape
	return roc_auc_score(test_label, pred), average_precision_score(test_label, pred)

# main


tissue = sys.argv[1]
split = sys.argv[2]
print ('tissue', tissue, 'split', split)

n_fold = int(sys.argv[7])
n_rep = int(sys.argv[8])

gene_list = np.load('./L1000/shrna_gene_list_'+tissue+'.npy')
feature_list, feature_dict = load_feature_list(tissue, gene_list)

symbolA_list, symbolB_list, labels = load_label(tissue, feature_dict)
symbol_set = set(symbolA_list).union(set(symbolB_list))
print ('number of unique genes with SL labels', len(symbol_set))

sl_mat, sl_mask, sl_gene_set = get_sl_mat(symbolA_list, symbolB_list, labels, gene_list, feature_dict)
print ('sl_mat', sl_mat.shape, sl_mat.sum(), sl_mask.sum())

initial_gene_features = np.array([feature_dict[gene] for gene in gene_list])
print ('initial_gene_features', initial_gene_features.shape)

sw = sys.argv[3]
dnn_layer = sys.argv[4]
dim = sys.argv[5]
l2 = sys.argv[6]



print ('='*30)
print ('tissue', tissue, 'split', split)
				
auc_list = []
seed = 0
for rep_i in range(n_rep):
	if split == 'edge':
		train_mat_list, train_mask_list, test_mat_list, test_mask_list, seed  = train_test_split_fixseed(sl_mat, sl_mask, n_fold, seed)  # random split
	elif split == 'gene':
		train_mat_list, train_mask_list, test_mat_list, test_mask_list, seed = leave_gene_split_fixseed(sl_mat, sl_mask, n_fold, seed)
	
	for fold_i in range(n_fold):
		print ('repeat', rep_i, 'fold', fold_i)
		train_mat, train_mask, test_mat, test_mask = train_mat_list[fold_i], train_mask_list[fold_i], test_mat_list[fold_i], test_mask_list[fold_i]
		print ('train_mat', train_mat.shape, train_mask.sum(), 'test_mat', test_mat.shape, test_mask.sum())
		model = DeepModel(initial_gene_features, float(sw), int(dnn_layer), int(dim), float(l2), n_epoch=1000, n_ensemble=1)
		model.fit(train_mat, train_mask, test_mat, test_mask)
		
		test_pred = model.predict(train_mat, train_mask, test_mask)
		print ('test_pred', test_pred.shape)
		test_label = test_mat[np.where(test_mask==1)]
		auc = roc_auc_score(test_label.reshape(-1), test_pred.reshape(-1))
		aupr = average_precision_score(test_label.reshape(-1), test_pred.reshape(-1))
		#auc, aupr = masked_auc(test_mat, test_mask, test_pred)
		print ('model auc, aupr', auc, aupr)
		auc_list.append([auc, aupr])

auc_list = np.array(auc_list)
print ('repeat avg', auc_list.shape, 'mean', np.mean(auc_list, axis=0), get_std(auc_list))
np.save('./EXP2SL_'+tissue+'_'+split+'_'+str(sw)+'_'+str(dnn_layer)+'_'+str(dim)+'_'+str(l2), auc_list)
