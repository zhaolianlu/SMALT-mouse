#####This script is used to perform tree pruning in 'Estimating the number of progenitor cells seeding each mouse neoplasm' of methods


import copy
from Bio import Phylo
import sys
import os
sys.setrecursionlimit(5000)  # or a higher number



def get_descendant_tips(clade):
	"""Recursively retrieve all tips under a given clade."""
	tips = []
	if clade.is_terminal():  # If it's a tip
		tips.append(clade.name)
	else:
		for child in clade:
			tips.extend(get_descendant_tips(child))
	return tips


def preorder_traversal(clade, cutoff, normal_tips):
	# If it's an internal node or a tip
	if not clade.is_terminal():
		descendant_tips = get_descendant_tips(clade)
		elements_start_with_N = [item for item in descendant_tips if item.startswith('N_')]
		# Calculate the percentage
		percentage = (len(elements_start_with_N) / len(descendant_tips)) * 100
		if percentage <= cutoff:
			# Extend the empty list with the extracted elements
			normal_tips.extend(elements_start_with_N)
			next
		else:
			for child in clade:
				preorder_traversal(child, cutoff, normal_tips)
	return normal_tips
def is_tumor_tip(tip):
	return not tip.name.startswith('N')
	#return 'T' in tip.name
def name_ancestral_nodes(tree, prefix="node"):
    """
    Name all unnamed ancestral nodes in a given tree.
    The naming will follow a pattern: prefix1, prefix2, ...
    """
    # Counter for the new node names
    count = 1

    # Walk through the tree's clades (pre-order traversal)
    for clade in tree.find_clades(order='preorder'):
        # If the clade (node) doesn't have a name and is not a terminal, name it
        if clade.name is None and not clade.is_terminal():
            clade.name = f"{prefix}{count}"
            count += 1

def is_monophyletic_tumor_group(clade):
	tumor_tips = [tip for tip in clade.get_terminals() if is_tumor_tip(tip)]
	return len(tumor_tips) == len(clade.get_terminals()) and len(tumor_tips) >= 3

def check_descendants(clade,flag):
	# Check the current clade
#	print(clade.name)
	if is_monophyletic_tumor_group(clade):
		return clade
	else:
		for child_clade in clade.clades:
			result = check_descendants(child_clade, flag)
			if result:
				return result
def clade_traversal(clade, tumor_tips):
	for child_clade in clade.clades:
		print(child_clade.name)
		flag = 0
		flag = check_descendants(child_clade,flag)
		if flag !=0 and flag is not None:
			if child_clade.name == flag.name:
				continue
			else:
				clade_traversal(child_clade,tumor_tips)
		else:
			descendant_tips = get_descendant_tips(child_clade)
			elements_start_with_T = [item for item in descendant_tips if item.startswith('T')]
			tumor_tips.extend(elements_start_with_T)
	return tumor_tips

tree = Phylo.read(sys.argv[1], 'newick')
cutoff = float(sys.argv[2])
normal_tips = []
normal_tips = preorder_traversal(tree.clade, cutoff, normal_tips)
# Prune the tree by excluding the specified tips
for tip in normal_tips:
	tree.prune(tip)
path_without_extension, extension = os.path.splitext(sys.argv[1])
Phylo.write(tree, f'{path_without_extension}.pruned_tree.{cutoff}.nwk', 'newick')


original_tree = copy.deepcopy(tree)
name_ancestral_nodes(tree)
tumor_tips = []
tumor_tips = clade_traversal(tree.root, tumor_tips)
# Create a dictionary mapping from original labels to modified labels
modified_list = ["N_" + item.split('_')[1] for item in tumor_tips]
label_mapping = dict(zip(tumor_tips, modified_list))

# Traverse the tree and modify the labels
for tip in original_tree.get_terminals():
	if tip.name in label_mapping:
		tip.name = label_mapping[tip.name]
Phylo.write(original_tree, f'{path_without_extension}.substituted_tumor_treev2.{cutoff}.nwk', 'newick')


