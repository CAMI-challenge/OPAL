# NOTE: Something is odd about the diffab vectors for weighted and weighted_flow since different vectors are being returned...
# Using real data, I'm pretty sure the problem is with the weighted_flow version...Now I'm absolutely confident of it. Weighted is ok, weighted_flow has issues...
import numpy as np
import dendropy
import matplotlib.pyplot as plt
import warnings


def parse_tree(tree_str):
	'''
	(Tint,lint,nodes_in_order) = parse_tree(tree_str) 
	This function will parse a newick tree string and return the dictionary of ancestors Tint.
	Tint indexes the nodes by integers, Tint[i] = j means j is the ancestor of i.
	lint is a dictionary returning branch lengths: lint[(i,j)] = w(i,j) the weight of the edge connecting i and j.
	nodes_in_order is a list of the nodes in the input tree_str such that Tint[i]=j means nodes_in_order[j] is an ancestor
	of nodes_in_order[i]. Nodes are labeled from the leaves up.
	'''
	dtree = dendropy.Tree.get(data=tree_str, schema="newick", suppress_internal_node_taxa=False,store_tree_weights=True)
	#Name all the internal nodes
	nodes = dtree.nodes()
	i=0
	for node in nodes:
		if node.taxon == None:
			node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
			i = i+1
	full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]  # i in path from root to j only if i>j
	full_nodes_in_order.reverse()
	nodes_in_order = [item.taxon.label for item in full_nodes_in_order]  # i in path from root to j only if i>j
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = full_nodes_in_order[i]
		parent = node.parent_node
		if parent != None:
			Tint[i] = nodes_to_index[parent.taxon.label]
			lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
	return (Tint,lint,nodes_in_order)



def parse_tree_file(tree_str_file, suppress_internal_node_taxa=False, suppress_leaf_node_taxa=False):
	'''
	(Tint,lint,nodes_in_order) = parse_tree(tree_str_file) 
	This function will parse a newick tree file (in the file given by tree_str_file) and return the dictionary of ancestors Tint.
	Tint indexes the nodes by integers, Tint[i] = j means j is the ancestor of i.
	lint is a dictionary returning branch lengths: lint[i,j] = w(i,j) the weight of the edge connecting i and j.
	nodes_in_order is a list of the nodes in the input tree_str such that T[i]=j means nodes_in_order[j] is an ancestor
	of nodes_in_order[i]. Nodes are labeled from the leaves up.
	'''
	dtree = dendropy.Tree.get(path=tree_str_file, schema="newick",
							suppress_internal_node_taxa=suppress_internal_node_taxa,
							store_tree_weights=True,
							suppress_leaf_node_taxa = suppress_leaf_node_taxa)
	#Name all the internal nodes
	nodes = dtree.nodes()
	i=0
	for node in nodes:
		if node.taxon == None:
			node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
			i = i+1
	full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]  # i in path from root to j only if i>j
	full_nodes_in_order.reverse()
	nodes_in_order = [item.taxon.label for item in full_nodes_in_order]  # i in path from root to j only if i>j
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = full_nodes_in_order[i]
		parent = node.parent_node
		if parent != None:
			Tint[i] = nodes_to_index[parent.taxon.label]
			if isinstance(node.edge.length, float):
				lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
			else:
				lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = 0.0
	return (Tint,lint,nodes_in_order)

def simulate_data(basis):
	'''
	Simulate environments, suitable for feeding directly to FastUnifrac. 
	Input is a list of nodes on which the distribution will be given.
	Will return distribution only on labeled nodes. Returns (envs)
	'''
	weights_sample1 = np.random.exponential(scale=1.0, size=(1,len(basis))).transpose()
	#Multiply since fast_unifrac expects numbers > 1
	weights_sample1 = 1000*weights_sample1 
	weights_sample2 = np.random.exponential(scale=1.0, size=(1,len(basis))).transpose()
	#make it sparse
	weights_sample2 = 1000*weights_sample2
	envs = dict()
	i=0
	for node in basis:
		envs[node] = {'sample1':weights_sample1[i],'sample2':weights_sample2[i]}
		i = i+1
	return (envs)


def parse_envs(envs, nodes_in_order):
	'''
	(envs_prob_dict, samples) = parse_envs(envs, nodes_in_order)
	This function takes an environment envs and the list of nodes nodes_in_order and will return a dictionary envs_prob_dict
	with keys given by samples. envs_prob_dict[samples[i]] is a probability vector on the basis nodes_in_order denoting for sample i.
	'''
	nodes_in_order_dict = dict(zip(nodes_in_order,range(len(nodes_in_order))))
	for node in envs.keys():
		if node not in nodes_in_order_dict:
			print("Warning: environments contain taxa " + node + " not present in given taxonomic tree. Ignoring")
	envs_prob_dict = dict()
	for i in range(len(nodes_in_order)):
		node = nodes_in_order[i]
		if node in envs:
			samples = envs[node].keys()
			for sample in samples:
				if sample not in envs_prob_dict:
					envs_prob_dict[sample] = np.zeros(len(nodes_in_order))
					envs_prob_dict[sample][i] = envs[node][sample]
				else:
					envs_prob_dict[sample][i] = envs[node][sample]
	#Normalize samples
	samples = envs_prob_dict.keys()
	for sample in samples:
		if envs_prob_dict[sample].sum() == 0:
			warnings.warn("Warning: the sample %s has non-zero counts, do not use for Unifrac calculations" % sample)
		envs_prob_dict[sample] = envs_prob_dict[sample]/envs_prob_dict[sample].sum()
	return (envs_prob_dict, samples)


def EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q):
	'''
	(Z, F, diffab) = EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q)
	This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
	and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
	Returns the weighted Unifrac distance Z, the flow F, and the differential abundance vector diffab.
	The flow F is a dictionary with keys of the form (i,j) where F[(i,j)] == num means that in the calculation of the
	Unifrac distance, a total mass of num was moved from the node nodes_in_order[i] to the node nodes_in_order[j].
	The differential abundance vector diffab is	a dictionary with tuple keys using elements of Tint and values
	diffab[(i, j)] equal to the signed difference of abundance between the two samples restricted to the sub-tree
	defined by nodes_in_order(i) and weighted by the edge length lint[(i,j)].
	'''
	num_nodes = len(nodes_in_order)
	F = dict()
	G = dict()
	diffab = dict()
	Z = 0
	w = np.zeros(num_nodes)
	pos = dict()
	neg = dict()
	for i in range(num_nodes):
		pos[i] = set([])
		neg[i] = set([])
	for i in range(num_nodes):
		if P[i] > 0 and Q[i] > 0:
			F[(i, i)] = np.minimum(P[i], Q[i])
		G[(i, i)] = P[i] - Q[i]
		if P[i] > Q[i]:
			pos[i].add(i)
		elif P[i] < Q[i]:
			neg[i].add(i)
		posremove = set()
		negremove = set()
		for j in pos[i]:
			for k in neg[i]:
				if (j not in posremove) and (k not in negremove):
					val = np.minimum(G[(i, j)], -G[(i, k)])
					if val > 0:
						F[(j, k)] = np.minimum(G[(i, j)], -G[(i, k)])
						G[(i, j)] = G[(i, j)] - val
						G[(i, k)] = G[(i, k)] + val
						Z = Z + (w[j] + w[k])*val
					if G[(i, j)] == 0:
						posremove.add(j)
					if G[(i, k)] == 0:
						negremove.add(k)
		pos[i].difference_update(posremove)
		neg[i].difference_update(negremove)
		if i < num_nodes-1:
			for j in pos[i].union(neg[i]):
					if (Tint[i], j) in G:
						G[(Tint[i], j)] = G[(Tint[i], j)] + G[(i, j)]
						diffab[(i, Tint[i])] = diffab[(i, Tint[i])] + G[(i, j)]  # Added to capture 'slack' in subtree JLMC
					else:
						G[(Tint[i], j)] = G[(i, j)]
						diffab[(i, Tint[i])] = G[(i, j)]  # Added to capture 'slack' in subtree JLMC
					w[j] = w[j] + lint[i, Tint[i]]
			if (i, Tint[i]) in diffab:
				#diffab[(i,Tint[i])] = lint[i,Tint[i]]*abs(diffab[(i,Tint[i])])  # Captures DiffAbund, scales by length of edge JLMC
				diffab[(i, Tint[i])] = lint[i, Tint[i]]*diffab[(i, Tint[i])]  # Captures DiffAbund, scales by length of edge JLMC
			pos[Tint[i]] |= pos[i]
			neg[Tint[i]] |= neg[i]
	return (Z, F, diffab)  # The returned flow and diffab are on the basis nodes_in_order and is given in sparse matrix dictionary format. eg {(0,0):.5,(1,2):.5}


# This will return the EMDUnifrac distance only
def EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q):
	'''
	(Z, diffab) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
	This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
	and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
	Returns the weighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
	F[(i,j)] == num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node
	nodes_in_order[i] to the node nodes_in_order[j].
	'''
	num_nodes = len(nodes_in_order)
	Z = 0
	diffab = dict()
	partial_sums = P - Q
	for i in range(num_nodes - 1):
		val = partial_sums[i]
		partial_sums[Tint[i]] += val
		if val != 0:
			diffab[(i, Tint[i])] = lint[i, Tint[i]]*val  # Captures diffab
		Z += lint[i, Tint[i]]*abs(val)
	return (Z, diffab)


# This will return the EMDUnifrac distance only
def EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q):
	'''
	(Z, diffab) = EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
	This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
	and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
	Returns the unweighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
	F[(i,j)] == 1 means that in the calculation of the Unifrac distance, a total mass of 1 was moved from the node
	nodes_in_order[i] to the node nodes_in_order[j].
	'''
	num_nodes = len(nodes_in_order)
	Z = 0
	diffab = dict()
	for i in range(num_nodes):
		if P[i]>0:
			P[i] = 1
		if Q[i]>0:
			Q[i] = 1
	partial_sums = P - Q
	for i in range(num_nodes - 1):
		val = partial_sums[i]
		partial_sums[Tint[i]] += val
		if val != 0:
			diffab[(i, Tint[i])] = lint[i, Tint[i]]*val  # Captures diffab
		Z += lint[i, Tint[i]]*abs(val)
	return Z, diffab


def EMDUnifrac_unweighted_flow(Tint, lint, nodes_in_order, P, Q):
	'''
	(Z, F) = EMDUnifrac_unweighted_flow(Tint, lint, nodes_in_order, P, Q)
	This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
	and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
	Returns the unweighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
	F[(i,j)] == 1 means that in the calculation of the Unifrac distance, a total mass of 1 was moved from the node
	nodes_in_order[i] to the node nodes_in_order[j].
	'''
	num_nodes = len(nodes_in_order)
	F = dict()
	G = dict()
	diffab = dict()
	Z = 0
	w = np.zeros(num_nodes)
	pos = dict()
	neg = dict()
	for i in range(num_nodes):
		pos[i] = set([])
		neg[i] = set([])
	for i in range(num_nodes):
		if P[i] > 0:
			P[i] = 1
		if Q[i] > 0:
			Q[i] = 1
		if P[i]>0 and Q[i]>0:
			F[(i,i)] = np.minimum(P[i],Q[i])
		G[(i,i)] = P[i] - Q[i]
		if P[i] > Q[i]:
			pos[i].add(i)
		elif P[i] < Q[i]:
			neg[i].add(i)
		posremove = set()
		negremove = set()
		for j in pos[i]:
			for k in neg[i]:
				if (j not in posremove) and (k not in negremove):
					val = np.minimum(G[(i,j)], -G[(i,k)])
					if val > 0:
						F[(j,k)] = np.minimum(G[(i,j)], -G[(i,k)])
						G[(i,j)] = G[(i,j)] - val
						G[(i,k)] = G[(i,k)] + val
						Z = Z + (w[j] + w[k])*val
					if G[(i,j)] == 0:
						posremove.add(j)
					if G[(i,k)] == 0:
						negremove.add(k)
		pos[i].difference_update(posremove)
		neg[i].difference_update(negremove)
		if i < num_nodes-1:
			for j in pos[i].union(neg[i]):
					if (Tint[i],j) in G:
						G[(Tint[i],j)] = G[(Tint[i],j)] + G[(i,j)]
						diffab[(i,Tint[i])] = diffab[(i,Tint[i])] + G[(i,j)]  # Added to capture 'slack' in subtree JLMC
					else:
						G[(Tint[i],j)] = G[(i,j)]
						diffab[(i,Tint[i])] = G[(i,j)]  # Added to capture 'slack' in subtree JLMC
					w[j] = w[j] + lint[i,Tint[i]]
			if (i, Tint[i]) in diffab:
				diffab[(i,Tint[i])] = lint[i,Tint[i]]*diffab[(i,Tint[i])]  # Captures DiffAbund, scales by length of edge JLMC
			pos[Tint[i]] |= pos[i]
			neg[Tint[i]] |= neg[i]
	return (Z, F, diffab)  # The returned flow and diffab are on the basis nodes_in_order and is given in sparse matrix dictionary format. eg {(0,0):.5,(1,2):.5}


def plot_diffab(nodes_in_order, diffab, P_label, Q_label, plot_zeros=True, thresh=0):
	'''
	plot_diffab(nodes_in_order, diffab, P_label, Q_label)
	Plots the differential abundance vector.
	:param nodes_in_order: list returned from parse_envs
	:param diffab: differential abundance vector (returned from one flavor of EMDUnifrac)
	:param P_label: label corresponding to the sample name for P (e.g. when calling EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q))
	:param Q_label: label corresponding to the sample name for P (e.g. when calling EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q))
	:param plot_zeros: flag (either True or False) that specifies if the zero locations should be plotted. Warning, if your tree is large and plot_zeros=True, this can cause a crash.
	:param thresh: only plot those parts of the diffab vector that are above thresh, specify everything else as zero
	:return: None (makes plot)
	'''
	x = range(len(nodes_in_order))
	y = np.zeros(len(nodes_in_order))
	keys = diffab.keys()
	#for i in range(len(nodes_in_order)):
	#	for key in keys:
	#		if key[0] == i:
	#			y[i] = diffab[key]
	# Much faster way to compute this
	for key in keys:
		y[key[0]] = diffab[key]

	pos_loc = [x[i] for i in range(len(y)) if y[i] > thresh]
	neg_loc = [x[i] for i in range(len(y)) if y[i] < -thresh]
	zero_loc = [x[i] for i in range(len(y)) if -thresh <= y[i] <= thresh]
	if not pos_loc:
		raise Exception('Threshold too high! Please lower and try again.')
	if not neg_loc:
		raise Exception('Threshold too high! Please lower and try again.')

	pos_val = [y[i] for i in range(len(y)) if y[i] > thresh]
	neg_val = [y[i] for i in range(len(y)) if y[i] < -thresh]
	zero_val = [y[i] for i in range(len(y)) if -thresh <= y[i] <= thresh]

	# The following is to get the indicies in order. Basically, I iterate down both pos_loc and neg_loc simultaneously
	# and create new lists (pos_loc_adj and neg_loc_adj) that are in the same order as pos_loc and neg_loc, but whose
	# union of indicies is equal to range(len(pos_loc + neg_loc)). Simply to make things pretty
	if plot_zeros:
		pos_loc_adj = pos_loc
		neg_loc_adj = neg_loc
		zero_loc_adj = zero_loc
	else:
		pos_loc_adj = []
		neg_loc_adj = []
		tick_names = []

		# rename the indicies so they are increasing by 1
		pos_ind = 0
		neg_ind = 0
		it = 0
		while pos_ind < len(pos_loc) or neg_ind < len(neg_loc):
			if pos_ind >= len(pos_loc):
				neg_loc_adj.append(it)
				tick_names.append(nodes_in_order[neg_loc[neg_ind]])
				it += 1
				neg_ind += 1
			elif neg_ind >= len(neg_loc):
				pos_loc_adj.append(it)
				tick_names.append(nodes_in_order[pos_loc[pos_ind]])
				it += 1
				pos_ind += 1
			elif pos_loc[pos_ind] < neg_loc[neg_ind]:
				pos_loc_adj.append(it)
				tick_names.append(nodes_in_order[pos_loc[pos_ind]])
				it += 1
				pos_ind += 1
			elif pos_loc[pos_ind] > neg_loc[neg_ind]:
				neg_loc_adj.append(it)
				tick_names.append(nodes_in_order[neg_loc[neg_ind]])
				it += 1
				neg_ind +=1
			else:
				print('Something went wrong')
				break


	fig, ax = plt.subplots()

	markerline, stemlines, baseline = ax.stem(neg_loc_adj, neg_val)
	plt.setp(baseline, 'color', 'k', 'linewidth', 1)
	plt.setp(markerline, 'color','r')
	for i in range(len(neg_loc)):
		plt.setp(stemlines[i], 'linewidth', 3)
		plt.setp(stemlines[i], 'color', 'r')

	markerline, stemlines, baseline = ax.stem(pos_loc_adj, pos_val)
	plt.setp(baseline, 'color', 'k', 'linewidth', 1)
	plt.setp(markerline, 'color','b')
	for i in range(len(pos_loc)):
		plt.setp(stemlines[i], 'linewidth', 3)
		plt.setp(stemlines[i], 'color', 'b')

	if plot_zeros:
		markerline, stemlines, baseline = ax.stem(zero_loc, zero_val)
		plt.setp(baseline, 'color', 'k', 'linewidth', 1)
		plt.setp(markerline, 'color','k')
		for i in range(len(zero_loc)):
			plt.setp(stemlines[i], 'linewidth', 3)
			plt.setp(stemlines[i], 'color', 'k')

	plt.ylabel('DiffAbund', fontsize=16)
	plt.gcf().subplots_adjust(right=0.93, left=0.15)
	# If you want the zeros plotted, label EVERYTHING, otherwise just label the things that are there...
	if plot_zeros:
		plt.xticks(x, nodes_in_order, rotation='vertical', fontsize=14)
	else:
		#tick_names = [nodes_in_order[i] for i in pos_loc] + [nodes_in_order[i] for i in neg_loc]  # Don't need this with new code
		plt.xticks(range(len(pos_loc_adj + neg_loc_adj)), tick_names, rotation='vertical', fontsize=14)

	plt.subplots_adjust(bottom=0.3, top=.93)
	plt.text(plt.xticks()[0][-1]+0.1, max(pos_val), P_label, rotation=90, horizontalalignment='center', verticalalignment='top', multialignment='center', color='b', fontsize=14)
	plt.text(plt.xticks()[0][-1]+0.1, min(neg_val), Q_label, rotation=90, horizontalalignment='center', verticalalignment='bottom', multialignment='center', color='r', fontsize=14)
	plt.show()

def EMDUnifrac_weighted_plain(ancestors, edge_lengths, nodes_in_order, P, Q):
	'''
	Z = EMDUnifrac_weighted(ancestors, edge_lengths, nodes_in_order, P, Q)
	This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
	and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
	Returns the weighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
	F[(i,j)] == num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node
	nodes_in_order[i] to the node nodes_in_order[j].
	'''
	num_nodes = len(nodes_in_order)
	Z = 0
	eps = 1e-8
	partial_sums = P - Q
	total_mass = 1
	for i in range(num_nodes - 1):
		val = partial_sums[i]
		if abs(val) > eps:  # if the value is big enough to care about (i.e. ignore zeros)
			#if np.sign(val) * np.sign(partial_sums[ancestors[i]]) == -1:  # if mass is getting destroyed
				#total_mass -= abs(val + partial_sums[ancestors[i]])  # keep track of how much was lost
				#print(total_mass)  # Can use this to break once the change is small enough
			partial_sums[ancestors[i]] += val
			Z += edge_lengths[i, ancestors[i]]*abs(val)
	return Z


def EMDUnifrac_group(ancestors, edge_lengths, nodes_in_order, rel_abund):
	eps = 1e-8
	num_nodes = len(nodes_in_order)
	num_samples = len(rel_abund)
	Z = np.zeros((num_samples, num_samples))
	partial_sums = np.zeros((num_samples, num_samples, num_nodes))
	for i in range(num_samples):
		for j in range(num_samples):
			partial_sums[i, j] = rel_abund[i] - rel_abund[j]
	for n in range(num_nodes - 1):
		for i in range(num_samples):
			for j in range(i, num_samples):
				val = partial_sums[i, j, n]
				if abs(val) > eps:  # only do the work if it's big enough
					partial_sums[i, j, ancestors[n]] += val
					Z[i, j] += edge_lengths[n, ancestors[n]] * abs(val)
	Z = Z + np.transpose(Z)
	return Z

#########################################
# Tests

def test_parse_tree():
	tree_str = '((B:0.1,C:0.2)A:0.3);'
	(Tint, lint, nodes_in_order) = parse_tree(tree_str)
	assert Tint == {0: 2, 1: 2, 2: 3}
	assert lint == {(1, 2): 0.1, (2, 3): 0.3, (0, 2): 0.2}
	assert nodes_in_order == ['C', 'B', 'A', 'temp0']  # temp0 is the root node

def test_simulate_data():
	basis = ['A', 'B', 'C', 'temp0']  # temp0 is the root node
	basis_sim = simulate_data(basis)
	assert basis_sim.keys().sort() == basis.sort()
	assert basis_sim['A'].keys() == ['sample1', 'sample2']
	assert basis_sim['B'].keys() == ['sample1', 'sample2']
	assert basis_sim['C'].keys() == ['sample1', 'sample2']
	assert basis_sim['temp0'].keys() == ['sample1', 'sample2']  # temp0 is the root node
	
def test_parse_envs():
	basis = ['C', 'B', 'A', 'temp0']  # temp0 is the root node
	basis_samples = {
		'C': {'sample1': 1, 'sample2': 0},
		'B': {'sample1': 1, 'sample2': 1},
		'A': {'sample1': 0, 'sample2': 0},
		'temp0': {'sample1': 0, 'sample2': 1}}  # temp0 is the root node
	(basis_weighted, samples) = parse_envs(basis_samples, basis)
	assert [basis_weighted['sample1'][i] for i in range(4)] == [0.5, 0.5, 0.0, 0.0]
	assert [basis_weighted['sample2'][i] for i in range(4)] == [0.0, 0.5, 0.0, 0.5]
	assert samples == ['sample1', 'sample2']
	
def test_EMDUnifrac_weighted_flow():
	tree_str = '((B:0.1,C:0.2)A:0.3);'  # input taxonomic tree
	(Tint, lint, nodes_in_order) = parse_tree(tree_str)  # parse the tree, getting nodes (Tint), edge lengths (lint), and node names (nodes_in_order)
	# Create a toy environment
	nodes_samples = {
		'C': {'sample1': 1, 'sample2': 0},
		'B': {'sample1': 1, 'sample2': 1},
		'A': {'sample1': 0, 'sample2': 0},
		'temp0': {'sample1': 0, 'sample2': 1}}  # temp0 is the root node, not named in Newick format, but included in nodes_in_order
	(nodes_weighted, samples) = parse_envs(nodes_samples, nodes_in_order)  # Parse the environments.
	(Z, F, diffab) = EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, nodes_weighted['sample1'], nodes_weighted['sample2'])  # Run the weighted version of EMDUnifrac that returns the flow
	# Check to make sure results make sense
	assert Z == 0.25  # This is the Unifrac distance
	assert F[(0, 3)] == 0.5  # F is the flow and is in a sparse matrix format: a dictionary with tuple keys using elements of Tint and values T[(i, j)] equal to amount of abundance moved from organism nodes_in_order(i) to nodes_in_order(j)
	assert F[(1, 1)] == 0.5
	assert sum(F.values()) == 1
	assert diffab == {(2, 3): 0.14999999999999999, (0, 2): 0.10000000000000001}  # diffab is the differential abundance vector, also in a sparse matrix format: a dictionary with tuple keys using elements of Tint and values diffab[(i, j)] equal to the signed difference of abundance between the two samples restricted to the sub-tree defined by nodes_in_order(i) and weighted by the edge length lint[(i,j)].
	
def test_EMDUnifrac_weighted():
	tree_str = '((B:0.1,C:0.2)A:0.3);'  # there is an internal node (temp0) here.
	(Tint, lint, nodes_in_order) = parse_tree(tree_str)
	nodes_samples = {
		'C': {'sample1': 1, 'sample2': 0},
		'B': {'sample1': 1, 'sample2': 1},
		'A': {'sample1': 0, 'sample2': 0},
		'temp0': {'sample1': 0, 'sample2': 1}}  # temp0 is the root node
	(nodes_weighted, samples) = parse_envs(nodes_samples, nodes_in_order)
	(Z, diffab) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, nodes_weighted['sample1'], nodes_weighted['sample2'])
	assert Z == 0.25
	assert diffab == {(2, 3): 0.14999999999999999, (0, 2): 0.10000000000000001}
	
def test_EMDUnifrac_unweighted():
	tree_str = '((B:0.1,C:0.2)A:0.3);'
	(Tint, lint, nodes_in_order) = parse_tree(tree_str)
	nodes_samples = {
		'C': {'sample1': 1, 'sample2': 0},
		'B': {'sample1': 1, 'sample2': 1},
		'A': {'sample1': 0, 'sample2': 0},
		'temp0': {'sample1': 0, 'sample2': 1}}
	(nodes_weighted, samples) = parse_envs(nodes_samples, nodes_in_order)
	(Z, diffab) = EMDUnifrac_unweighted(Tint, lint, nodes_in_order, nodes_weighted['sample1'], nodes_weighted['sample2'])
	assert Z == 0.5
	assert diffab == {(2, 3): 0.29999999999999999, (0, 2): 0.20000000000000001}
	
def test_EMDUnifrac_unweighted_flow():
	tree_str = '((B:0.1,C:0.2)A:0.3);'
	(Tint, lint, nodes_in_order) = parse_tree(tree_str)
	nodes_samples = {
		'C': {'sample1': 1, 'sample2': 0},
		'B': {'sample1': 1, 'sample2': 1},
		'A': {'sample1': 0, 'sample2': 0},
		'temp0': {'sample1':0,'sample2': 1}}
	(nodes_weighted, samples) = parse_envs(nodes_samples, nodes_in_order)
	(Z, F, diffab) = EMDUnifrac_unweighted_flow(Tint, lint, nodes_in_order, nodes_weighted['sample1'], nodes_weighted['sample2'])
	assert Z == 0.5
	assert F[(0, 3)] == 1
	assert F[(1, 1)] == 1
	assert sum(F.values()) == 2
	assert diffab == {(2, 3): 0.29999999999999999, (0, 2): 0.20000000000000001}


def run_tests():
	test_parse_tree()
	test_simulate_data()
	test_parse_envs()
	test_EMDUnifrac_weighted_flow()
	test_EMDUnifrac_weighted()
	test_EMDUnifrac_unweighted()
	test_EMDUnifrac_unweighted_flow()
