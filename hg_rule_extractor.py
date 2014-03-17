#coding: utf8
import sys
import heapq
from tree import TreeNode, NonTerminalNode, TerminalNode
from collections import defaultdict
from helpers import computeSpans, Alignment, getSpanMask
from hypergraph import Hypergraph
from itertools import izip

# Reads the next line from a file stream representing input trees.
# TODO: Allow this to read in Hypergraphs or Berkeley k-best lists
def read_tree_file(stream):
	line = stream.readline()
	if line:
		line = line.decode('utf-8').strip()
		# Sometimes Berkeley Parser likes to given broken trees
		if line == '(())':
			yield None

		# Otherwise, turn the string into a tree, and then into a hypergraph
		tree = TreeNode.from_string(line)
		computeSpans(tree)
		hg = Hypergraph.from_tree(tree)
		yield hg

# Reads the next line from a file stream representing alignments.
# TODO: Allow this to add probabilities to alignment links.
def read_alignment_file(stream):
	line = stream.readline()
	if  line:
		line = line.decode('utf-8').strip()
		alignment = Alignment.from_string(line)
		yield alignment

# Determines whether source_node and target_node are node-aligned.
def are_aligned(source_node, target_node, source_hg, target_hg, s2t_word_alignments, t2s_word_alignments):
	source_terminals = set(source_node.find_terminals(source_hg))
	target_terminals = set(target_node.find_terminals(target_hg))

	has_alignments = True in [len(s2t_word_alignments[terminal]) > 0 for terminal in source_terminals]
	if not has_alignments:
		return False

	for source_terminal in source_terminals:
		for target_terminal in s2t_word_alignments[source_terminal]:
			if target_terminal not in target_terminals:
				return False

	for target_terminal in target_terminals:
		for source_terminal in t2s_word_alignments[target_terminal]:
			if source_terminal not in source_terminals:
				return False

	return True

# Extracts all available rules from the node pair source_node and target_node
def extract_rules(source_node, target_node, s2t_node_alignments, t2s_node_alignments, source_root, target_root, max_rule_size):
	rules = set()
	for source_children in source_node.get_child_sets(source_root):
		for target_children in target_node.get_child_sets(target_root):
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
					for terminal in source_children[i].find_terminals(source_root):
						source_parts.append(terminal)
			for node in unused_target_children:
				for terminal in node.find_terminals(target_root):
					target_parts.append(terminal)

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
			is_abstract = True
			is_phrase = True
			for node in source_parts:
				if node.is_terminal(source_root):
					is_abstract = False
				else:
					is_phrase = False

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

			# output !
			rule_type = 'A' if is_abstract else 'P' if is_phrase else 'G'
			lhs = '[%s::%s]' % (source_node.label, target_node.label)
			alignments = []
			rule = ' ||| '.join([rule_type, lhs, ' '.join(source_rhs), ' '.join(target_rhs), ' '.join('%d-%d' % link for link in alignments), ' '.join(node_types)])
			rules.add(rule)
	for rule in rules:
		print rule.encode('utf-8')
	sys.stdout.flush()

# Computes a dictionary whos keys are nodes in a hypergraph, and whose keys
# are integers, representing distance from the root of the hypergraph.
def compute_generations(node, root, current_generation=0, generations={}):
	if node not in generations or current_generation < generations[node]:
		generations[node] = current_generation
	for children in node.get_child_sets(root):
		for child in children:
				compute_generations(child, root, current_generation + 1, generations)
	return generations

# Removes non-minimal edges from a hypergraph. We define an edge as
# non-minimal if it skips (i.e. composes over) a node that is
# node-aligned to something in the opposite tree.
def minimize_alignments(s2t, t2s, source_root, target_root):
	def node_key(node):
		type_match = 1 if source_node.is_terminal(source_root) == node.is_terminal(target_root) else 0
		span_size = node.span.end - node.span.start
		generation = generations[node]
		return (-type_match, span_size, -generation)
	generations = compute_generations(target_root.start, target_root)
	
	for source_node, target_nodes in s2t.items():
		min_span_size = min([target_node.span.end - target_node.span.start for target_node in target_nodes])
		best_alignment = heapq.nsmallest(1, target_nodes, key=node_key)[0]
		for target_node in target_nodes.copy():
			if target_node != best_alignment:
				s2t[source_node].remove(target_node)
				t2s[target_node].remove(source_node)
				if len(t2s[target_node]) == 0:
					del t2s[target_node]

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
			for t_node in target_tree.nodes:
					if are_aligned(s_node, t_node, source_tree, target_tree, s2t_word_alignments, t2s_word_alignments):
						s2t_node_alignments[s_node].add(t_node)
						t2s_node_alignments[t_node].add(s_node)

		if args.minimal_rules:
			minimize_alignments(s2t_node_alignments, t2s_node_alignments, source_tree, target_tree)
			minimize_alignments(t2s_node_alignments, s2t_node_alignments, target_tree, source_tree)

		# Finally extract rules
		for source_node, target_nodes in s2t_node_alignments.iteritems():
			for target_node in target_nodes:
				if not source_node.is_terminal(source_tree) and not target_node.is_terminal(target_tree):
					extract_rules(source_node, target_node, s2t_node_alignments.copy(), t2s_node_alignments, source_tree, target_tree, args.max_rule_size)

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('source_trees')
        parser.add_argument('target_trees')
        parser.add_argument('alignments')
        parser.add_argument('--virtual_size', '-v', type=int, default=1, help='Maximum number of components in a virtual node')
        parser.add_argument('--minimal_rules', '-m', action='store_true', help='Only extract minimal rules')
        parser.add_argument('--max_rule_size', '-s', type=int, default=5, help='Maximum number of parts (terminal or non-terminal) in the RHS of a rule')
        args = parser.parse_args()

	source_tree_file = open(args.source_trees)
	target_tree_file = open(args.target_trees)
	alignment_file = open(args.alignments)

	sentence_number = 1
	for source_tree, target_tree, alignment in izip(read_tree_file(source_tree_file), read_tree_file(target_tree_file), read_alignment_file(alignment_file)):
		print 'Sentence', sentence_number
		# Can happen if Berkeley gives borked trees
		if source_tree == None or target_tree == None:
			continue
		handle_sentence(source_tree, target_tree, alignment)
		sys.stdout.flush()
		sentence_number += 1
