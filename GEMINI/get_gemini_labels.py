import numpy as np
import mygene
"""
# ============================ A549 ===============================
# Step1 
pair_list = []
score_list = []
f = open('Big-Papi.strong')
for line in f.readlines()[2:]:
	lines = line.strip().split('\t')
	symbolA, symbolB = lines[0].split(';')
	A549_score = float(lines[3])
	pair_list.append([symbolA, symbolB])
	score_list.append(A549_score)
f.close()
print 'pair_list', len(pair_list)
print 'score_list', len(score_list)

pair_list = np.array(pair_list)
score_list = np.array(score_list)

select_idx = np.argsort(-score_list)[:len(pair_list)/20]
assert score_list[np.argsort(-score_list)[len(pair_list)/20]] > 0
print 'select_idx', select_idx.shape
neg_idx = np.where(score_list<0)[0]
assert len(neg_idx)*2 < len(score_list)
print 'neg_idx', neg_idx.shape

select_pairs_bigpapi = pair_list[select_idx]
neg_pairs_bigpapi = pair_list[neg_idx]
print 'select_pairs_bigpapi', len(select_pairs_bigpapi), 'neg_pairs_bigpapi', len(neg_pairs_bigpapi)

# Step2 
pair_list = []
score_list = []
f = open('Shen-Mali.strong')
for line in f.readlines()[2:]:
	lines = line.strip().split('\t')
	symbolA, symbolB = lines[0].split(';')
	A549_score = float(lines[2])
	pair_list.append([symbolA, symbolB])
	score_list.append(A549_score)
f.close()
print 'pair_list', len(pair_list)
print 'score_list', len(score_list)

pair_list = np.array(pair_list)
score_list = np.array(score_list)

select_idx = np.argsort(-score_list)[:len(pair_list)/20]
assert score_list[np.argsort(-score_list)[len(pair_list)/20]] > 0
print 'select_idx', select_idx.shape
neg_idx = np.where(score_list<0)[0]
assert len(neg_idx)*2 < len(score_list)
print 'neg_idx', neg_idx.shape

select_pairs_shenmali = pair_list[select_idx]
neg_pairs_shenmali = pair_list[neg_idx]
print 'select_pairs_shenmali', len(select_pairs_shenmali),  'neg_pairs_shenmali', len(neg_pairs_shenmali)

# Step3 
pair_list = []
score_list = []
f = open('Zhao-Mali.strong')
for line in f.readlines()[2:]:
	lines = line.strip().split('\t')
	symbolA, symbolB = lines[0].split(';')
	A549_score = float(lines[2])
	pair_list.append([symbolA, symbolB])
	score_list.append(A549_score)
f.close()
print 'pair_list', len(pair_list)
print 'score_list', len(score_list)

pair_list = np.array(pair_list)
score_list = np.array(score_list)

select_idx = np.argsort(-score_list)[:len(pair_list)/20]
assert score_list[np.argsort(-score_list)[len(pair_list)/20]] > 0
print 'select_idx', select_idx.shape
neg_idx = np.where(score_list<0)[0]
assert len(neg_idx)*2 < len(score_list)
print 'neg_idx', neg_idx.shape

select_pairs_zhaomali = pair_list[select_idx]
neg_pairs_zhaomali = pair_list[neg_idx]
print 'select_pairs_zhaomali', len(select_pairs_zhaomali), 'neg_pairs_zhaomali', len(neg_pairs_zhaomali)

A549_pairs = select_pairs_bigpapi.tolist() + select_pairs_shenmali.tolist() + select_pairs_zhaomali.tolist()
A549_negs = neg_pairs_bigpapi.tolist() + neg_pairs_shenmali.tolist() + neg_pairs_zhaomali.tolist()
print 'A549_pairs', len(A549_pairs), 'A549_negs', len(A549_negs)

mg = mygene.MyGeneInfo()

symbol_to_gene = {}
fw = open('gemini_A549_labels.tsv', 'w')
for symbolA, symbolB in A549_pairs:
	if symbolA not in symbol_to_gene:
		geneA_hits = mg.query('symbol:'+symbolA, species='human')['hits']
		assert len(geneA_hits) == 1
		geneA = geneA_hits[0]['_id']
		symbol_to_gene[symbolA] = geneA
	else:
		geneA = symbol_to_gene[symbolA]
	if symbolB not in symbol_to_gene:
		geneB_hits = mg.query('symbol:'+symbolB, species='human')['hits']
		assert len(geneB_hits) == 1
		geneB = geneB_hits[0]['_id']
		symbol_to_gene[symbolB] = geneB
	else:
		geneB = symbol_to_gene[symbolB]
	fw.write(symbolA + '\t' + geneA + '\t' + symbolB + '\t' + geneB + '\t1\n')
for symbolA, symbolB in A549_negs:
	if symbolA not in symbol_to_gene:
		geneA_hits = mg.query('symbol:'+symbolA, species='human')['hits']
		if not len(geneA_hits) == 1:
			continue
		geneA = geneA_hits[0]['_id']
		symbol_to_gene[symbolA] = geneA
	else:
		geneA = symbol_to_gene[symbolA]
	if symbolB not in symbol_to_gene:
		geneB_hits = mg.query('symbol:'+symbolB, species='human')['hits']
		if not len(geneB_hits) == 1:
			continue
		geneB = geneB_hits[0]['_id']
		symbol_to_gene[symbolB] = geneB
	else:
		geneB = symbol_to_gene[symbolB]
	fw.write(symbolA + '\t' + geneA + '\t' + symbolB + '\t' + geneB + '\t0\n')
fw.close()
print 'symbol_to_gene', len(symbol_to_gene)
"""
"""
# ============================ A375 ===============================
# Step1 
pair_list = []
score_list = []
f = open('Big-Papi.strong')
for line in f.readlines()[2:]:
	lines = line.strip().split('\t')
	symbolA, symbolB = lines[0].split(';')
	A375_score = float(lines[2])
	pair_list.append([symbolA, symbolB])
	score_list.append(A375_score)
f.close()
print 'pair_list', len(pair_list)
print 'score_list', len(score_list)

pair_list = np.array(pair_list)
score_list = np.array(score_list)

select_idx = np.argsort(-score_list)[:len(pair_list)/20]
assert score_list[np.argsort(-score_list)[len(pair_list)/20]] > 0
print 'select_idx', select_idx.shape
neg_idx = np.where(score_list<0)[0]
assert len(neg_idx)*2 < len(score_list)
print 'neg_idx', neg_idx.shape

select_pairs_bigpapi = pair_list[select_idx]
neg_pairs_bigpapi = pair_list[neg_idx]
print 'select_pairs_bigpapi', len(select_pairs_bigpapi), 'neg_pairs_bigpapi', len(neg_pairs_bigpapi)

A375_pairs = select_pairs_bigpapi.tolist() #+ select_pairs_shenmali.tolist() + select_pairs_zhaomali.tolist()
A375_negs = neg_pairs_bigpapi.tolist() #+ neg_pairs_shenmali.tolist() + neg_pairs_zhaomali.tolist()
print 'A375 pos', len(A375_pairs), 'A375 neg', len(A375_negs)

mg = mygene.MyGeneInfo()

symbol_to_gene = {}
fw = open('gemini_A375_labels.tsv', 'w')
for symbolA, symbolB in A375_pairs:
	if symbolA not in symbol_to_gene:
		geneA_hits = mg.query('symbol:'+symbolA, species='human')['hits']
		assert len(geneA_hits) == 1
		geneA = geneA_hits[0]['_id']
		symbol_to_gene[symbolA] = geneA
	else:
		geneA = symbol_to_gene[symbolA]
	if symbolB not in symbol_to_gene:
		geneB_hits = mg.query('symbol:'+symbolB, species='human')['hits']
		assert len(geneB_hits) == 1
		geneB = geneB_hits[0]['_id']
		symbol_to_gene[symbolB] = geneB
	else:
		geneB = symbol_to_gene[symbolB]
	fw.write(symbolA + '\t' + geneA + '\t' + symbolB + '\t' + geneB + '\t1\n')
for symbolA, symbolB in A375_negs:
	if symbolA not in symbol_to_gene:
		geneA_hits = mg.query('symbol:'+symbolA, species='human')['hits']
		if not len(geneA_hits) == 1:
			continue
		geneA = geneA_hits[0]['_id']
		symbol_to_gene[symbolA] = geneA
	else:
		geneA = symbol_to_gene[symbolA]
	if symbolB not in symbol_to_gene:
		geneB_hits = mg.query('symbol:'+symbolB, species='human')['hits']
		if not len(geneB_hits) == 1:
			continue
		geneB = geneB_hits[0]['_id']
		symbol_to_gene[symbolB] = geneB
	else:
		geneB = symbol_to_gene[symbolB]
	fw.write(symbolA + '\t' + geneA + '\t' + symbolB + '\t' + geneB + '\t0\n')
fw.close()
print 'symbol_to_gene', len(symbol_to_gene)
"""

# ============================ HEK293T ===============================
# Step2 
pair_list = []
score_list = []
f = open('Shen-Mali.strong')
for line in f.readlines()[2:]:
	lines = line.strip().split('\t')
	symbolA, symbolB = lines[0].split(';')
	HEK_score = float(lines[1])
	pair_list.append([symbolA, symbolB])
	score_list.append(HEK_score)
f.close()
print 'pair_list', len(pair_list)
print 'score_list', len(score_list)

pair_list = np.array(pair_list)
score_list = np.array(score_list)

select_idx = np.argsort(-score_list)[:len(pair_list)/10]
assert score_list[np.argsort(-score_list)[len(pair_list)/10]] > 0
print 'select_idx', select_idx.shape
neg_idx = np.argsort(score_list)[:len(pair_list)/2]
assert len(set(neg_idx.tolist()).intersection(select_idx.tolist())) == 0
#neg_idx = np.where(score_list<0)[0]
#assert len(neg_idx)*2 < len(score_list)
print 'neg_idx', neg_idx.shape

select_pairs_shenmali = pair_list[select_idx]
neg_pairs_shenmali = pair_list[neg_idx]
print 'select_pairs_shenmali', len(select_pairs_shenmali),  'neg_pairs_shenmali', len(neg_pairs_shenmali)


HEK_pairs = select_pairs_shenmali.tolist() 
HEK_negs = neg_pairs_shenmali.tolist()
print 'HEK pos', len(HEK_pairs), 'HEK neg', len(HEK_negs)


mg = mygene.MyGeneInfo()

symbol_to_gene = {}
fw = open('gemini_HEK_labels.tsv', 'w')
for symbolA, symbolB in HEK_pairs:
	if symbolA not in symbol_to_gene:
		geneA_hits = mg.query('symbol:'+symbolA, species='human')['hits']
		assert len(geneA_hits) == 1
		geneA = geneA_hits[0]['_id']
		symbol_to_gene[symbolA] = geneA
	else:
		geneA = symbol_to_gene[symbolA]
	if symbolB not in symbol_to_gene:
		geneB_hits = mg.query('symbol:'+symbolB, species='human')['hits']
		assert len(geneB_hits) == 1
		geneB = geneB_hits[0]['_id']
		symbol_to_gene[symbolB] = geneB
	else:
		geneB = symbol_to_gene[symbolB]
	fw.write(symbolA + '\t' + geneA + '\t' + symbolB + '\t' + geneB + '\t1\n')
for symbolA, symbolB in HEK_negs:
	if symbolA not in symbol_to_gene:
		geneA_hits = mg.query('symbol:'+symbolA, species='human')['hits']
		if not len(geneA_hits) == 1:
			continue
		geneA = geneA_hits[0]['_id']
		symbol_to_gene[symbolA] = geneA
	else:
		geneA = symbol_to_gene[symbolA]
	if symbolB not in symbol_to_gene:
		geneB_hits = mg.query('symbol:'+symbolB, species='human')['hits']
		if not len(geneB_hits) == 1:
			continue
		geneB = geneB_hits[0]['_id']
		symbol_to_gene[symbolB] = geneB
	else:
		geneB = symbol_to_gene[symbolB]
	fw.write(symbolA + '\t' + geneA + '\t' + symbolB + '\t' + geneB + '\t0\n')
fw.close()
print 'symbol_to_gene', len(symbol_to_gene)

"""
# ============================ HT29 ===============================
# Step1 
pair_list = []
score_list = []
f = open('Big-Papi.strong')
for line in f.readlines()[2:]:
	lines = line.strip().split('\t')
	symbolA, symbolB = lines[0].split(';')
	HT29_score = float(lines[4])
	pair_list.append([symbolA, symbolB])
	score_list.append(HT29_score)
f.close()
print 'pair_list', len(pair_list)
print 'score_list', len(score_list)

pair_list = np.array(pair_list)
score_list = np.array(score_list)

select_idx = np.argsort(-score_list)[:len(pair_list)/20]
assert score_list[np.argsort(-score_list)[len(pair_list)/20]] > 0
print 'select_idx', select_idx.shape
neg_idx = np.argsort(score_list)[:len(pair_list)/2]
assert len(set(neg_idx.tolist()).intersection(select_idx.tolist())) == 0
#neg_idx = np.where(score_list<0)[0]
#assert len(neg_idx)*2 < len(score_list)
print 'neg_idx', neg_idx.shape

select_pairs_bigpapi = pair_list[select_idx]
neg_pairs_bigpapi = pair_list[neg_idx]
print 'select_pairs_bigpapi', len(select_pairs_bigpapi), 'neg_pairs_bigpapi', len(neg_pairs_bigpapi)

HT29_pairs = select_pairs_bigpapi.tolist() #+ select_pairs_shenmali.tolist() + select_pairs_zhaomali.tolist()
HT29_negs = neg_pairs_bigpapi.tolist() #+ neg_pairs_shenmali.tolist() + neg_pairs_zhaomali.tolist()
print 'HT29 pos', len(HT29_pairs), 'HT29 neg', len(HT29_negs)

mg = mygene.MyGeneInfo()

symbol_to_gene = {}
fw = open('gemini_HT29_labels.tsv', 'w')
for symbolA, symbolB in HT29_pairs:
	if symbolA not in symbol_to_gene:
		geneA_hits = mg.query('symbol:'+symbolA, species='human')['hits']
		assert len(geneA_hits) == 1
		geneA = geneA_hits[0]['_id']
		symbol_to_gene[symbolA] = geneA
	else:
		geneA = symbol_to_gene[symbolA]
	if symbolB not in symbol_to_gene:
		geneB_hits = mg.query('symbol:'+symbolB, species='human')['hits']
		assert len(geneB_hits) == 1
		geneB = geneB_hits[0]['_id']
		symbol_to_gene[symbolB] = geneB
	else:
		geneB = symbol_to_gene[symbolB]
	fw.write(symbolA + '\t' + geneA + '\t' + symbolB + '\t' + geneB + '\t1\n')
for symbolA, symbolB in HT29_negs:
	if symbolA not in symbol_to_gene:
		geneA_hits = mg.query('symbol:'+symbolA, species='human')['hits']
		if not len(geneA_hits) == 1:
			continue
		geneA = geneA_hits[0]['_id']
		symbol_to_gene[symbolA] = geneA
	else:
		geneA = symbol_to_gene[symbolA]
	if symbolB not in symbol_to_gene:
		geneB_hits = mg.query('symbol:'+symbolB, species='human')['hits']
		if not len(geneB_hits) == 1:
			continue
		geneB = geneB_hits[0]['_id']
		symbol_to_gene[symbolB] = geneB
	else:
		geneB = symbol_to_gene[symbolB]
	fw.write(symbolA + '\t' + geneA + '\t' + symbolB + '\t' + geneB + '\t0\n')
fw.close()
print 'symbol_to_gene', len(symbol_to_gene)
"""