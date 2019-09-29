import numpy as np
import pickle
import time

import torch
from torch import nn
import torch.optim as optim
import torch.nn.functional as F
from torch.autograd import Variable
from torch.nn.utils import weight_norm
from sklearn.metrics import roc_auc_score, average_precision_score

class Rank_Loss2(nn.Module):
	def __init__(self):
		super(Rank_Loss2, self).__init__()
		self.criterion = nn.BCELoss(reduce=False)
		self.MSEcriterion = nn.MSELoss(reduce=False)
	def forward(self, output_pos, output_neg, output_neg_sampled, output_unknown, semi_weight):
		label_pos = 1*torch.ones(output_pos.size()).cuda()
		loss_pos = self.MSEcriterion((output_pos), label_pos)
		label_neg = -1*torch.ones(output_neg.size()).cuda()
		loss_neg = self.MSEcriterion((output_neg), label_neg)
		loss = torch.sum(loss_pos) + torch.sum(loss_neg) 

		un_score = torch.sigmoid(output_unknown - output_neg_sampled)
		pu_score = torch.sigmoid(output_pos - output_unknown)
		label = torch.ones(un_score.size()).cuda()
		loss += semi_weight*(torch.sum(self.criterion(un_score, label)) + torch.sum(self.criterion(pu_score, label)))
		return loss


class Net(nn.Module):
	def __init__(self, init_features, dnn_layers, dim):
		super(Net, self).__init__()
		self.init_features = torch.Tensor(init_features).cuda()
		self.dim = dim #hidden_size
		self.dnn_layers = dnn_layers #dnn_layers
		self.MLP1 = nn.Linear(978, self.dim)
		self.MLP2 = nn.ModuleList([nn.Linear(self.dim, self.dim) for i in range(self.dnn_layers)])
		self.LR = nn.Linear(2*self.dim, 1)

	def forward(self, train_mat, pos_idx, neg_all_idx, unknown_all_idx, test_idx, epoch):
		#n_gene = self.init_features.size(0)
		time0 = time.time()

		gene_features = self.init_features
		
		gene_features = F.relu(self.MLP1(gene_features))
		for i in range(self.dnn_layers):
			gene_features = F.relu(self.MLP2[i](gene_features))

		# get_index
		time1 = time.time()
		n_pos = len(pos_idx[0])
		neg_sampled = np.random.choice(range(neg_all_idx.shape[1]), n_pos, replace=True)
		neg_all = np.random.choice(range(neg_all_idx.shape[1]), neg_all_idx.shape[1], replace=False)
		#neg_idx = neg_all_idx
		neg_idx_sampled = (neg_all_idx[0][neg_sampled], neg_all_idx[1][neg_sampled])
		neg_idx = (neg_all_idx[0][neg_all], neg_all_idx[1][neg_all])
		unknown_sampled = range(epoch*n_pos, (epoch+1)*n_pos)  #np.random.choice(range(unknown_all_idx.shape[1]), len(pos_idx[0]), replace=True)
		unknown_idx = (unknown_all_idx[0][unknown_sampled], unknown_all_idx[1][unknown_sampled])
		
		time2 = time.time()
		
		sl_pos_left = gene_features[pos_idx[0].tolist()+pos_idx[1].tolist()]
		sl_pos_right = gene_features[pos_idx[1].tolist()+pos_idx[0].tolist()]
		sl_neg_left = gene_features[neg_idx[0].tolist()+neg_idx[1].tolist()]
		sl_neg_right = gene_features[neg_idx[1].tolist()+neg_idx[0].tolist()]
		sl_neg_left_sampled = gene_features[neg_idx_sampled[0].tolist()+neg_idx_sampled[1].tolist()]
		sl_neg_right_sampled = gene_features[neg_idx_sampled[1].tolist()+neg_idx_sampled[0].tolist()]


		sl_unknown_left = gene_features[unknown_idx[0].tolist()+unknown_idx[1].tolist()]
		sl_unknown_right = gene_features[unknown_idx[1].tolist()+unknown_idx[0].tolist()]
		sl_test_left = gene_features[test_idx[0].tolist()+test_idx[1].tolist()]
		sl_test_right = gene_features[test_idx[1].tolist()+test_idx[0].tolist()]
		
		time3 = time.time()
		

		concat_pos = torch.cat((sl_pos_left, sl_pos_right), dim=1)
		concat_neg = torch.cat((sl_neg_left, sl_neg_right), dim=1)
		concat_neg_sampled = torch.cat((sl_neg_left_sampled, sl_neg_right_sampled), dim=1)
		concat_unknown = torch.cat((sl_unknown_left, sl_unknown_right),dim=1)
		concat_test = torch.cat((sl_test_left, sl_test_right), dim=1)

		output_pos = self.LR(concat_pos)
		output_neg = self.LR(concat_neg)
		output_neg_sampled = self.LR(concat_neg_sampled)
		output_unknown = self.LR(concat_unknown)
		output_test = self.LR(concat_test)

		output_test_size = output_pos.size(0)/2
		output_pos = 0.5*(output_pos[:output_test_size] + output_pos[output_test_size:])

		output_test_size = output_neg.size(0)/2
		output_neg = 0.5*(output_neg[:output_test_size] + output_neg[output_test_size:])

		output_test_size = output_neg_sampled.size(0)/2
		output_neg_sampled = 0.5*(output_neg_sampled[:output_test_size] + output_neg_sampled[output_test_size:])
		
		output_test_size = output_unknown.size(0)/2
		output_unknown = 0.5*(output_unknown[:output_test_size] + output_unknown[output_test_size:])


		output_test_size = output_test.size(0)/2
		output_test = 0.5*(output_test[:output_test_size] + output_test[output_test_size:])

		time4 = time.time()
		return output_pos, output_neg, output_neg_sampled, output_unknown, output_test


class DeepModel():
	def __init__(self, initial_gene_features, semi_weight, dnn_layers, dim, l2, n_epoch, n_ensemble=1):
		#self.batch_size = batch_size
		self.n_ensemble = n_ensemble
		self.n_epoch = n_epoch
		self.model_list = []
		self.l2 = l2
		self.semi_weight = semi_weight
		for ensemble in range(self.n_ensemble):
			model = Net(initial_gene_features, dnn_layers, dim)
			model.cuda()
			self.model_list.append(model)
	
	def fit(self, train_mat, train_mask, test_mat, test_mask):
		i = 0
		# get_index
		unknown_mask = np.ones(train_mask.shape) - train_mask - test_mask
		
		pos_idx = np.array(np.where(train_mat != 0))
		neg_all_idx = np.array(np.where( (train_mask-train_mat) != 0))
		unknown_all_idx = np.array(np.where(unknown_mask != 0))
		shuffle_unknown = np.array(range(unknown_all_idx.shape[1]))
		np.random.shuffle(shuffle_unknown)
		unknown_all_idx = unknown_all_idx[:, shuffle_unknown]
		unknown_all_idx = np.hstack([unknown_all_idx, unknown_all_idx, unknown_all_idx, unknown_all_idx, unknown_all_idx])
		test_idx = np.array(np.where(test_mask != 0))
		print 'pos_idx', pos_idx.shape, 'neg_all_idx', neg_all_idx.shape, 'unknown_all_idx', unknown_all_idx.shape, 'test_idx',test_idx.shape
		
		train_mat = torch.Tensor(train_mat).cuda()
		train_mask = torch.Tensor(train_mask).cuda()
		for model in self.model_list:
			i += 1
			print 'ensemble', i
			self.fit_single_mode(train_mat, train_mask, model, test_mat, test_mask, pos_idx, neg_all_idx, unknown_all_idx, test_idx)
	
	def fit_single_mode(self, train_mat, train_mask, model, test_mat, test_mask, pos_idx, neg_all_idx, unknown_all_idx, test_idx):
		optimizer = optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=0.001, weight_decay=self.l2, amsgrad=False)
		scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=3000, gamma=0.5)
		criterion1 = Rank_Loss2()
		for epoch in range(self.n_epoch):
			#scheduler.step()
			optimizer.zero_grad()
			output_pos, output_neg, output_neg_sampled, output_unknown, output_test = model(train_mat, pos_idx, neg_all_idx, unknown_all_idx, test_idx, epoch)
	
			loss = criterion1(output_pos, output_neg, output_neg_sampled, output_unknown, self.semi_weight)
			total_loss = float(loss.data)
			loss.backward()
			nn.utils.clip_grad_norm_(model.parameters(), 5)
			optimizer.step()
			if epoch % 50 == 0:
				test_label = test_mat[np.where(test_mask==1)]
				auc = roc_auc_score(test_label.reshape(-1), output_test.cpu().detach().numpy().reshape(-1))
				aupr = average_precision_score(test_label.reshape(-1), output_test.cpu().detach().numpy().reshape(-1))
				print 'epoch', epoch, 'total_loss', total_loss, 'test auc', auc, 'test aupr', aupr

	def predict(self, train_mat, train_mask, test_mask):		
		pos_idx = np.array(np.where(train_mat != 0))
		neg_all_idx = np.array(np.where( (train_mask-train_mat) != 0))
		unknown_all_idx = pos_idx#neg_all_idx #np.array(range(len(pos_idx)))#
		test_idx = np.array(np.where(test_mask != 0))

		train_mat = torch.Tensor(train_mat).cuda()

		ensemble_sl_pred = None
		for model in self.model_list:
			output_pos, output_neg, output_neg_sampled, output_unknown, output_test = \
			model(train_mat, pos_idx, neg_all_idx, unknown_all_idx, test_idx, 0)
			if ensemble_sl_pred is None:
				ensemble_sl_pred = output_test.cpu().detach().numpy()
			else:
				ensemble_sl_pred += output_test.cpu().detach().numpy()
		return  ensemble_sl_pred/self.n_ensemble

def masked_auc(test_mat, test_mask, test_pred):
	test_idx = np.where(test_mask==1)
	test_label = test_mat[test_idx]
	pred = test_pred[test_idx]
	print 'test_label', test_label.shape, 'pred', pred.shape
	return roc_auc_score(test_label, pred)
	
	
