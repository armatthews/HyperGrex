#coding: utf8
import sys
import math
import heapq
from tree import TreeNode, NonTerminalNode, TerminalNode
from collections import defaultdict
from helpers import computeSpans, Alignment, compute_generations
from hypergraph import Hypergraph
from itertools import izip

# Reads the next line from a file stream representing input trees.
# TODO: Allow this to read in Hypergraphs or Berkeley k-best lists
def read_tree_file(stream):
	while True:
		line = stream.readline()
		if not line:
			break
		line = line.decode('utf-8').strip()
		# Sometimes Berkeley Parser likes to given broken trees
		if line == '(())':
			yield None
			continue

		# Otherwise, turn the string into a tree, and then into a hypergraph
		tree = TreeNode.from_string(line)
		computeSpans(tree)
		hg = Hypergraph.from_tree(tree)
		yield hg

def read_kbest_tree_file(stream):
	# Input format is
	# for each sentence:
	# logprob1	( tree1 )
	# logprob2	( tree2 )
	# (empty line)
	# Note the blank line after each seentence's k-best list
	# and the extra parentheses that need stripped

	def combine_trees(trees_to_combine):
		if len(trees_to_combine) == 0:
			return None
		hypergraphs_to_combine = []
		total_scores = sum(score for _, score in trees_to_combine)
		for tree, score in trees_to_combine:
			computeSpans(tree)	
			tree_hg = Hypergraph.from_tree(tree, score/total_scores)
			tree_hg.sanity_check()
			hypergraphs_to_combine.append(tree_hg)

		final_hypergraph = hypergraphs_to_combine[0]
		for hypergraph in hypergraphs_to_combine[1:]:
			final_hypergraph.combine(hypergraph)
		return final_hypergraph

	trees_to_combine = []
	hg = None
	while True:
		line = stream.readline()
		if not line:
			break

		line = line.decode('utf-8').strip()
		if not line:
			yield combine_trees(trees_to_combine)
			trees_to_combine = []
			continue
		# TODO: It's as of yet unknown what Berkeley does in a kbest
		#       list when its best tree is (())
		if line.startswith('Don\'t have a'):
			continue
		score, line = line.split('\t')
		score = math.exp(float(score))

		# Strip the extra parens
		# TODO: What if the trees have different root nodes?
		assert line.startswith('( (')
		assert line.endswith(') )')
		line = line[2:-2]

		tree = TreeNode.from_string(line)
		trees_to_combine.append((tree, score))

	if len(trees_to_combine) > 0:
		yield combine_trees(trees_to_combine)

# Reads the next line from a file stream representing alignments.
# TODO: Allow this to add probabilities to alignment links.
def read_alignment_file(stream):
	while True:
		line = stream.readline()
		if not line:
			break
		line = line.decode('utf-8').strip()
		alignment = Alignment.from_string(line)
		yield alignment

# Determines whether source_node and target_node are node-aligned.
def are_aligned(source_node, target_node, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments):
	source_node_terminals = source_terminals[source_node.span.start : source_node.span.end]
	target_node_terminals = target_terminals[target_node.span.start : target_node.span.end]

	has_alignments = False
	for terminal in source_node_terminals:
		if len(s2t_word_alignments[terminal]) > 0:
			has_alignments = True
			break
	if not has_alignments:
		return False

	for source_terminal in source_node_terminals:
		for target_terminal in s2t_word_alignments[source_terminal]:
			if target_terminal not in target_node_terminals:
				return False

	for target_terminal in target_node_terminals:
		for source_terminal in t2s_word_alignments[target_terminal]:
			if source_terminal not in source_node_terminals:
				return False

	return True

# Extracts all available rules from the node pair source_node and target_node
def extract_rules(source_node, target_node, s2t_node_alignments, t2s_node_alignments, source_root, target_root, max_rule_size, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments, minimal_only):
	rules = list()
	for source_child_edge in source_node.get_child_edges(source_root):
		source_children = source_child_edge.tails
		source_nonterminal_count = len([node for node in source_children if not node.is_terminal(source_root)])
		for target_child_edge in target_node.get_child_edges(target_root):
			target_children = target_child_edge.tails
			target_nonterminal_count = len([node for node in target_children if not node.is_terminal(target_root)])
			if source_nonterminal_count != target_nonterminal_count:
				continue

			# work out which of the source_node's children correspond to which of the target_node's
			child_alignments = []
			for source_child_node in source_children:
				alignments = []
				source_child_node_alignments = s2t_node_alignments[source_child_node]
				for target_child_node in target_children:
					if target_child_node in source_child_node_alignments:
						if source_child_node.is_terminal(source_root) == target_child_node.is_terminal(target_root):
							alignments.append(target_child_node)
				child_alignments.append(alignments)

			# Build up lists of the parts that make up the source and target RHS
			# index is the number of NT-NT pairs that have been added to the rule so far
			# the rule_part_maps give, for each NT, what opposite-side NT it corresponds to along with its index
			index = 1
			source_parts = []
			target_parts = []
			s2t_rule_part_map = {}
			t2s_rule_part_map = {}
			unused_target_children = list(target_children)
			non_minimal = False
			for i in range(len(source_children)):
				assert len(child_alignments[i]) in [0, 1]
				if len(child_alignments[i]) == 1:
					s = source_children[i]
					t = child_alignments[i][0]
					source_parts.append(s)
					target_parts.append(t)
					unused_target_children.remove(child_alignments[i][0])
					if s.is_terminal(source_root):
						pass
					else:
						s2t_rule_part_map[s] = (t, index)
						t2s_rule_part_map[t] = (s, index)
						index += 1
				else:	
					source_parts += source_terminals[source_children[i].span.start : source_children[i].span.end]
					if not source_children[i].is_terminal(source_root):
						non_minimal = True
			for node in unused_target_children:
				target_parts += target_terminals[node.span.start : node.span.end]
				if not node.is_terminal(target_root):
					non_minimal = True

			if minimal_only and non_minimal:
				continue

			if non_minimal:
				continue

			if len(source_parts) > max_rule_size or len(target_parts) > max_rule_size:
				continue
			# sort the RHS parts, as the NTs and Terminals are in there intermixed now
			source_parts = sorted(source_parts, key=lambda node: node.span.start)
			target_parts = sorted(target_parts, key=lambda node: node.span.start)

			# Turn the source and target RHS's into strings
			source_rhs = []
			for part in source_parts:
				if part.is_terminal(source_root):
					source_rhs.append(part.label)
				else:
					target, index = s2t_rule_part_map[part]
					source_rhs.append('[%s::%s,%d]' % (part.label, target.label, index))

			target_rhs = []
			for part in target_parts:
				if part.is_terminal(target_root):
					target_rhs.append(part.label)
				else:
					source, index = t2s_rule_part_map[part]
					target_rhs.append('[%s::%s,%d]' % (source.label, part.label, index))

			# Work out the rule type
			is_abstract_source = True
			is_abstract_target = True
			is_phrase = True
			for node in source_parts:
				if node.is_terminal(source_root):
					is_abstract_source = False
				else:
					is_phrase = False

			for node in target_parts:
				if node.is_terminal(target_root):
					is_abstract_target = False

			node_types = []
			s = 'V' if '-' in source_node.label else 'O'
			t = 'V' if '-' in target_node.label else 'O'
			node_types.append(s + t)
			for source in source_parts:
				if source.is_terminal(source_root):
					continue
				target = s2t_rule_part_map[source]
				s = 'V' if '-' in source.label else 'O'
				t = 'V' if '-' in source.label else 'O'
				node_types.append(s + t)

			# Calculate the node-to-node and word-to-word alignments within this rule
			alignments = []
			for i, source_part in enumerate(source_parts):
				for j, target_part in enumerate(target_parts):
					if target_part in s2t_node_alignments[source_part] or source_part in t2s_node_alignments[target_part] or \
					   target_part in s2t_word_alignments[source_part] or source_part in t2s_word_alignments[target_part]:
						alignments.append((i, j))

			# output !
			if is_phrase:
				rule_type = 'P'
			elif is_abstract_source and is_abstract_target:
				rule_type = 'A'
			elif is_abstract_source:
				rule_type = 'A'
			elif is_abstract_target:
				rule_type = 'A'
			else:
				rule_type = 'G'
			lhs = '[%s::%s]' % (source_node.label, target_node.label)
			weight = source_root.weights[source_child_edge] * target_root.weights[target_child_edge]
			#debug_str = ' ||| '.join([' '.join([node.label for node in source_children]), ' '.join([node.label for node in target_children]), str(source_root.weights[source_child_edge]), str(target_root.weights[target_child_edge])])

			parts = [rule_type, lhs, ' '.join(source_rhs), ' '.join(target_rhs), ' '.join('%d-%d' % link for link in alignments), ' '.join(node_types)]
			if args.kbest_input:
				parts.append(str(weight))
			#parts.append(debug_str)
			rule = ' ||| '.join(parts)
			rules.append(rule)
	for rule in rules:
		print rule.encode('utf-8')
	sys.stdout.flush()

def find_best_minimal_alignment(node, target_nodes, taken_target_nodes, target_generations):
	target_nodes = [target_node for target_node in target_nodes if target_node not in taken_target_nodes]
	if len(target_nodes) == 0:
		return None

	min_span_size = min([target_node.span.end - target_node.span.start for target_node in target_nodes])
	minimal_alignments = [target_node for target_node in target_nodes
                                          if target_node.span.end - target_node.span.start == min_span_size]

	# If there is more than one, such as in a unary chain
	# return the aligned node lowest in the tree, i.e. with max generation
	return max(minimal_alignments, key=lambda target_node: target_generations[target_node])

def minimize_alignments_helper(source_node, s2t, t2s, target_generations, source_root, taken_target_nodes, visited_nodes):
	for edge in source_root.head_index[source_node]:
		for child in edge.tails:
			if child not in visited_nodes:
				minimize_alignments_helper(child, s2t, t2s, target_generations, source_root, taken_target_nodes, visited_nodes)

	best_alignment = find_best_minimal_alignment(source_node, s2t[source_node], taken_target_nodes, target_generations)
	if best_alignment is not None:
		for target_node in s2t[source_node].copy():
			if target_node != best_alignment:
				s2t[source_node].remove(target_node)
				t2s[target_node].remove(source_node)
				if len(t2s[target_node]) == 0:
					del t2s[target_node]
		taken_target_nodes.add(best_alignment)
	visited_nodes.add(source_node)

# Removes non-minimal node alignments from a hypergraph. We define an edge as
# non-minimal if it skips (i.e. composes over) a node that is
# node-aligned to something in the opposite tree.
def minimize_alignments(source_root, target_root, s2t, t2s):
	target_generations = compute_generations(target_root)
	visited_nodes = set()
	taken_target_nodes = set()
	minimize_alignments_helper(source_root.start, s2t, t2s, target_generations, source_root, taken_target_nodes, visited_nodes)

# Takes two hypergraphs representing source and target trees, as well as a word
# alignment, and finds all rules extractable there from.
def handle_sentence(source_tree, target_tree, alignment):
		# Add any necessary edges to the tree structures
		source_tree.add_virtual_nodes(args.virtual_size, False)
		target_tree.add_virtual_nodes(args.virtual_size, False)

		if not args.minimal_rules:
			source_tree.add_composed_edges(args.max_rule_size)
			target_tree.add_composed_edges(args.max_rule_size)

		# Build word alignment maps
		s2t_word_alignments = defaultdict(list)
		t2s_word_alignments = defaultdict(list)

		source_terminals = source_tree.start.find_terminals(source_tree)
		target_terminals = target_tree.start.find_terminals(target_tree)

		for s, t in alignment.links:
			s_node = source_terminals[s]
			t_node = target_terminals[t]
			s2t_word_alignments[s_node].append(t_node)
			t2s_word_alignments[t_node].append(s_node)

		# Build node alignments maps
		s2t_node_alignments = defaultdict(set)
		t2s_node_alignments = defaultdict(set)
		for s_node in source_tree.nodes:
			s2t_node_alignments[s_node] = set()
			for t_node in target_tree.nodes:
					if are_aligned(s_node, t_node, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments):
						s2t_node_alignments[s_node].add(t_node)
						t2s_node_alignments[t_node].add(s_node)

		# The roots of the two trees are always node-aligned, even when there are no alignment links
		s2t_node_alignments[source_tree.start].add(target_tree.start)
		t2s_node_alignments[target_tree.start].add(source_tree.start)

		if args.minimal_rules:
			#minimize_alignments(s2t_node_alignments, t2s_node_alignments, source_tree, target_tree)
			#minimize_alignments(t2s_node_alignments, s2t_node_alignments, target_tree, source_tree)	
			minimize_alignments2(source_tree, target_tree, s2t_node_alignments, t2s_node_alignments)
			minimize_alignments2(target_tree, source_tree, t2s_node_alignments, s2t_node_alignments)

		# Finally extract rules
		for source_node, target_nodes in s2t_node_alignments.copy().iteritems():
			for target_node in target_nodes:
				if not source_node.is_terminal(source_tree) and not target_node.is_terminal(target_tree):
					extract_rules(source_node, target_node, s2t_node_alignments, t2s_node_alignments, source_tree, target_tree, args.max_rule_size, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments, args.minimal_rules)

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('source_trees')
        parser.add_argument('target_trees')
        parser.add_argument('alignments')
        parser.add_argument('--kbest_input', '-k', action='store_true', help='Assume the source tree input file is a kbest list')
        parser.add_argument('--virtual_size', '-v', type=int, default=1, help='Maximum number of components in a virtual node')
        parser.add_argument('--minimal_rules', '-m', action='store_true', help='Only extract minimal rules')
        parser.add_argument('--max_rule_size', '-s', type=int, default=5, help='Maximum number of parts (terminal or non-terminal) in the RHS of a rule')
        args = parser.parse_args()

	source_tree_file = open(args.source_trees)
	target_tree_file = open(args.target_trees)
	alignment_file = open(args.alignments)

	sentence_number = 1
	read_tree_file_func = read_kbest_tree_file if args.kbest_input else read_tree_file
	for source_tree, target_tree, alignment in izip(read_tree_file_func(source_tree_file), read_tree_file_func(target_tree_file), read_alignment_file(alignment_file)):
		print 'Sentence', sentence_number
		# Can happen if Berkeley gives borked trees	
		if source_tree == None or target_tree == None:
			continue
		handle_sentence(source_tree, target_tree, alignment)
		sys.stdout.flush()
		sentence_number += 1
