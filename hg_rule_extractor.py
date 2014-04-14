#coding: utf8
import sys
import math
import heapq
from tree import TreeNode, NonTerminalNode, TerminalNode
from collections import defaultdict
from helpers import computeSpans, Alignment, compute_generations
from hypergraph import Hypergraph
from itertools import izip

# Turns a line of input into a hypergraph.
# Returns a hypergraph, a log weight to be used when combining this HG with others,
# and a boolean that indicates whether this tree is one of a k-best list.
# Returns None on error, including some quirky cases caused by Berkeley parser.
# Input format is one of these for each sentence:
# Format 1 (just 1-best trees):
# tree1
# Format 2:
# log p(tree1 | sent)	( tree1 )
# log p(tree2 | sent)	( tree2 )
# ...
# (empty line)
# Format 3:
# log p(sent)	log p(tree1, sent)	( tree1 )
# log p(sent)	log p(tree2, sent)	( tree2 )
# ...
# (empty line)
# Note the blank line after each seentence's k-best list
# and the extra parentheses that need stripped

def hypergraph_from_line(line):
	line = line.strip()
	# Berkeley parser will output lines like
	# "Don't have a 7-best tree" when you ask it
	# for 10 best, but it only has 6.
	if line.startswith('Don\'t have a'):
		return None

	is_one_of_kbest = True
	parts = line.split('\t')
	if len(parts) == 1:
		score = 0.0
		line, = parts
		is_one_of_kbest = False
	elif len(parts) == 2:
		score, line = parts
		score = float(score)
	else:
		sent_prob, joint_prob, line = parts
		score = float(joint_prob) - float(sent_prob)

	if score == '-Infinity' or line == '(())':
		return None

	score = math.exp(score)
	# Strip the extra parens
	# TODO: What if the trees have different root nodes?
	if line.startswith('( (') and line.endswith(') )'):
		line = line[2:-2]

	tree = TreeNode.from_string(line)
	return tree, score, is_one_of_kbest

# Input is a list of (tree, weight) pairs.
# Weights are normalized to sum to 1, then the trees are combined into a single hypergraph.
def combine_trees(trees_to_combine):
	if len(trees_to_combine) == 0:
		return None
	hypergraphs_to_combine = []
	total_scores = sum(score for _, score in trees_to_combine)
	for tree, score in trees_to_combine:
		if total_scores != 0.0:
			score = score / total_scores
		else:
			score = 1.0 / len(trees_to_combine)
		computeSpans(tree)
		tree_hg = Hypergraph.from_tree(tree, score)
		tree_hg.sanity_check()
		hypergraphs_to_combine.append(tree_hg)

	final_hypergraph = hypergraphs_to_combine[0]
	for hypergraph in hypergraphs_to_combine[1:]:
		final_hypergraph.combine(hypergraph)
	return final_hypergraph


def read_tree_file(stream):
	trees_to_combine = []
	while True:
		line = stream.readline()
		if not line:
			break

		line = line.decode('utf-8').strip()
		if not line:
			yield combine_trees(trees_to_combine)
			trees_to_combine = []
			continue

		parts = hypergraph_from_line(line)
		if parts == None:
			continue
		else:
			tree, score, can_combine = parts
			trees_to_combine.append((tree, score))
		
		if not can_combine:
			assert len(trees_to_combine) == 1
			yield combine_trees(trees_to_combine)
			trees_to_combine = []

	if len(trees_to_combine) > 0:
		yield combine_trees(trees_to_combine)

def read_string_file(stream):
	while True:
		line = stream.readline()
		if not line:
			break
		line = line.decode('utf-8').strip()
		yield Hypergraph.from_surface_string(line)

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
# Two nodes are aligned if the alignment links emanating from their terminals
# align only to terminals of the other, or to NULL.
# The one exception is if both nodes have no alignment links coming from their
# terminals at all. In this case the nodes are not considered to be aligned. 
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

def build_mini_hypergraph(edges):
	hg = Hypergraph(edges[0].head)
	for edge in edges:
		hg.add(edge)	
	return hg

# Extracts all available rules from the node pair source_node and target_node
# Each rule will come from a pair of edges, one with head at source_node and
# the other with head at target_node.
# These two edges must have the same number of non-terminal children.
def extract_rules(source_node, target_node, s2t_node_alignments, t2s_node_alignments, source_root, target_root, max_rule_size, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments, minimal_only):
	rules = list()
	for source_edge in source_node.get_child_edges(source_root):	
		if len(source_edge.tails) > max_rule_size:
			continue
		source_nt_count = len([node for node in source_edge.tails if not node.is_terminal(source_root)])
		for target_edge in target_node.get_child_edges(target_root):	
			if len(target_edge.tails) > max_rule_size:
				continue
			target_nt_count = len([node for node in target_edge.tails if not node.is_terminal(target_root)])
			if source_nt_count != target_nt_count:
				continue

			# Work out which of the source_node's children correspond to which of the target_node's
			# Note: Each node should have at most 1 alignment
			target_terminals = set([node for node in target_edge.tails if node.is_terminal(target_root)])
			target_nonterminals = set([node for node in target_edge.tails if not node.is_terminal(target_root)])
			child_alignments = []
			for node in source_edge.tails:	
				possible_alignments = target_terminals if node.is_terminal(source_root) else target_nonterminals 
				child_alignments.append(list(possible_alignments & s2t_node_alignments[node]))
				assert len(child_alignments[-1]) <= 1

			# Build up lists of the parts that make up the source and target RHS
			# index is the number of NT-NT pairs that have been added to the rule so far
			# the rule_part_maps give, for each NT, what opposite-side NT it corresponds to along with its index
			index = 1
			s2t_rule_part_map = {}
			t2s_rule_part_map = {}
			unused_target_children = set(target_edge.tails)
			has_unaligned_nt = False
			for i, s in enumerate(source_edge.tails):
				if len(child_alignments[i]) == 1:
					t = child_alignments[i][0]
					unused_target_children.remove(t)
					if not s.is_terminal(source_root):
						s2t_rule_part_map[s] = (t, index)
						t2s_rule_part_map[t] = (s, index)
						index += 1
				elif not source_edge.tails[i].is_terminal(source_root):
						has_unaligned_nt = True
						break

			has_unaligned_nt |= False in [node.is_terminal(target_root) for node in unused_target_children]
			if has_unaligned_nt:
				continue

			# At this point, all information defining the rule has been calculated.
			# All that remains is to output it, which requires turning various bits
			# of information into string form.

			# Turn the source and target RHS's into strings
			source_rhs = []
			for node in source_edge.tails:
				if node.is_terminal(source_root):
					source_rhs.append(node.label)
				else:
					target, index = s2t_rule_part_map[node]
					source_rhs.append('[%s::%s,%d]' % (node.label, target.label, index))

			target_rhs = []
			for node in target_edge.tails:
				if node.is_terminal(target_root):
					target_rhs.append(node.label)
				else:
					source, index = t2s_rule_part_map[node]
					target_rhs.append('[%s::%s,%d]' % (source.label, node.label, index))

			# Work out the rule type
			is_abstract_source = True
			is_abstract_target = True
			is_phrase = True
			for node in source_edge.tails:
				if node.is_terminal(source_root):
					is_abstract_source = False
				else:
					is_phrase = False

			for node in target_edge.tails:
				if node.is_terminal(target_root):
					is_abstract_target = False

			node_types = []
			s = 'V' if '-' in source_node.label else 'O'
			t = 'V' if '-' in target_node.label else 'O'
			node_types.append(s + t)
			for source in source_edge.tails:
				if source.is_terminal(source_root):
					continue
				target, _ = s2t_rule_part_map[source]	
				s = 'V' if '-' in source.label else 'O'
				t = 'V' if '-' in target.label else 'O'
				node_types.append(s + t)

			# Calculate the node-to-node and word-to-word alignments within this rule
			alignments = []
			for i, source_part in enumerate(source_edge.tails):
				for j, target_part in enumerate(target_edge.tails):
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
			weight = source_root.weights[source_edge] * target_root.weights[target_edge]

			#debug_str = ' ||| '.join([' '.join([node.label for node in source_children]), ' '.join([node.label for node in target_edge.tails]), str(source_root.weights[source_edge]), str(target_root.weights[target_edge])])
			source_composed_tree = ''
			if source_edge.is_composed:
				source_composed_tree = build_mini_hypergraph(source_edge.composed_edges).to_tree_string()
			target_composed_tree = ''
			if target_edge.is_composed:
				target_composed_hg = build_mini_hypergraph(target_edge.composed_edges)
				label_map = {node: str(index) for (node, (_, index)) in t2s_rule_part_map.iteritems()}	
				target_composed_tree = target_composed_hg.to_tree_string(label_map)
			#debug_str = ' '.join([str(source_edge.is_composed), str(target_edge.is_composed)]) + ' ||| ' + (source_composed_tree) + ' ||| ' + target_composed_tree 

			parts = [rule_type, lhs, ' '.join(source_rhs), ' '.join(target_rhs), ' '.join('%d-%d' % link for link in alignments), ' '.join(node_types)]
			if not args.suppress_counts:
				parts.append('1.0 ' + str(weight))
			#parts.append(debug_str)
			rule = ' ||| '.join(parts)
			rules.append(rule)
	for rule in rules:
		print rule.encode('utf-8')
	sys.stdout.flush()

def find_best_minimal_alignment(node, target_nodes, taken_target_nodes, target_generations):
	target_nodes = [target_node for target_node in target_nodes if target_node not in taken_target_nodes and target_node.is_terminal_flag == node.is_terminal_flag]
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
		# Add virtual nodes and edges to the tree structures
		source_tree.add_virtual_nodes(args.virtual_size, False)
		target_tree.add_virtual_nodes(args.virtual_size, False)

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
			minimize_alignments(source_tree, target_tree, s2t_node_alignments, t2s_node_alignments)
			minimize_alignments(target_tree, source_tree, t2s_node_alignments, s2t_node_alignments)

		# Add composed edges to the tree structures
		if args.minimal_rules:
			s2t_aligned_nodes = set(node for node, alignments in s2t_node_alignments.iteritems() if len(alignments) > 0)
			t2s_aligned_nodes = set(node for node, alignments in t2s_node_alignments.iteritems() if len(alignments) > 0)
			source_tree.add_minimal_composed_edges(args.max_rule_size, s2t_aligned_nodes)
			target_tree.add_minimal_composed_edges(args.max_rule_size, t2s_aligned_nodes)
		else:
			source_tree.add_composed_edges(args.max_rule_size)
			target_tree.add_composed_edges(args.max_rule_size)


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
	parser.add_argument('--suppress_counts', action='store_true', help='Emulate old rule extractor by hiding the count field(s)')
        parser.add_argument('--virtual_size', '-v', type=int, default=1, help='Maximum number of components in a virtual node')
        parser.add_argument('--minimal_rules', '-m', action='store_true', help='Only extract minimal rules')
        parser.add_argument('--max_rule_size', '-s', type=int, default=5, help='Maximum number of parts (terminal or non-terminal) in the RHS of a rule')
	group = parser.add_mutually_exclusive_group()
        group.add_argument('--s2t', action='store_true', help='String-to-tree mode. Target side file should contain (tokenized) sentences instead of trees.')
        group.add_argument('--t2s', action='store_true', help='Tree-to-string mode. Source side file should contain (tokenized) sentences instead of trees.')
        parser.add_argument('--debug', '-d', action='store_true', help='Debug mode')
        args = parser.parse_args()

	source_tree_file = open(args.source_trees)
	target_tree_file = open(args.target_trees)
	alignment_file = open(args.alignments)

	read_source_file = read_tree_file if not args.s2t else read_string_file
	read_target_file = read_tree_file if not args.t2s else read_string_file

	sentence_number = 1
	for source_tree, target_tree, alignment in izip(read_source_file(source_tree_file), read_target_file(target_tree_file), read_alignment_file(alignment_file)):
		print 'Sentence', sentence_number	
		# Can happen if Berkeley gives borked trees	
		if source_tree == None or target_tree == None:
			pass
		else:
			try:
				handle_sentence(source_tree, target_tree, alignment)
			except Exception as e:
				if args.debug:
					raise
		sys.stdout.flush()
		sentence_number += 1	
