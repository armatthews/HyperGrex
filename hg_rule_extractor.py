#coding: utf8
import sys
import math
import heapq
from rule_formatters import RuleFormatter, GrexRuleFormatter, CdecT2SRuleFormatter
from tree import TreeNode, NonTerminalNode, TerminalNode
from collections import defaultdict
from helpers import computeSpans, Alignment, compute_generations, Rule, Span
from hypergraph import Hypergraph, NodeWithSpan, Edge
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
			yield None
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
def are_aligned(source_span, target_span, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments):
	source_node_terminals = source_terminals[source_span.start : source_span.end]
	target_node_terminals = target_terminals[target_span.start : target_span.end]

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
# Each rule will come from a pair of edges, one with head at source_node and
# the other with head at target_node.
# These two edges must have the same number of non-terminal children.
def extract_rules(source_node, target_node, s2t_node_alignments, t2s_node_alignments, source_root, target_root, max_rule_size, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments):
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

			weight = source_root.weights[source_edge] * target_root.weights[target_edge]
			# At this point, all information defining the rule has been calculated.
			# All that remains is to output it, which requires turning various bits
			# of information into string form.

			# Calculate the node-to-node and word-to-word alignments within this rule
			alignments = []
			for i, source_part in enumerate(source_edge.tails):	
				for j, target_part in enumerate(target_edge.tails):
					if target_part in s2t_node_alignments[source_part] or source_node in t2s_node_alignments[target_part] or \
					   target_part in s2t_word_alignments[source_part] or source_node in t2s_word_alignments[target_part]:
						alignments.append((i, j))

			yield Rule(source_edge, target_edge, s2t_rule_part_map, t2s_rule_part_map, alignments, weight)

def find_best_minimal_alignment(node, target_nodes, taken_target_nodes, target_generations):
	# only allow terminals to align to terminals, and NTs to align to NTs
	target_nodes = set([target_node for target_node in target_nodes
					if target_node.is_terminal_flag == node.is_terminal_flag])
	if len(target_nodes) == 0:
		return None

	# A minimal alignment will have the smallest span of all aligned nodes
	min_span_size = min([target_node.span.end - target_node.span.start for target_node in target_nodes])
	target_nodes = set([target_node for target_node in target_nodes
					if target_node.span.end - target_node.span.start == min_span_size])

	# Whatever node we're aligned to must not have been previously aligned to something else
	target_nodes = target_nodes - taken_target_nodes
	if len(target_nodes) == 0:
		return None

	# If there is more than one, such as in a unary chain
	# return the aligned node lowest in the tree, i.e. with max generation
	return max(target_nodes, key=lambda target_node: target_generations[target_node])

# Removes non-minimal node alignments from a hypergraph.
# If a node is aligned to multiple opposite nodes,
# its minimal alignment is the smallest one span-wise.
# If that's a tie, the minimal alignment is the one lowest
# in the tree, i.e. the one with the highest generation.
def minimize_alignments(source_root, target_root, s2t, t2s):
	if args.t2s:
		target_generations = {}
		for node in target_root.nodes:
			target_generations[node] = 2 if node.is_terminal_flag else 1
		target_generations[target_root.start] = 0
	else:
		target_generations = compute_generations(target_root)

	taken_target_nodes = set()
	minimal_s2t = defaultdict(set)
	minimal_t2s = defaultdict(set)
	for source_node in source_root.topsort():
		target_nodes = s2t[source_node]
		target_node = find_best_minimal_alignment(source_node, target_nodes, taken_target_nodes, target_generations)
		if target_node is not None:
			taken_target_nodes.add(target_node)
			minimal_s2t[source_node].add(target_node)
			minimal_t2s[target_node].add(source_node)

	# Ensure the root nodes aren't aligned to anything other than each other
	for node in minimal_s2t.iterkeys():
		if target_tree.start in minimal_s2t[node]:
			minimal_s2t[node].remove(target_root.start)

	for node in minimal_t2s.iterkeys():
		if source_tree.start in minimal_t2s[node]:
			minimal_t2s[node].remove(source_root.start)

	minimal_s2t[source_tree.start] = set([target_root.start])
	minimal_t2s[target_tree.start] = set([source_root.start])


	for k, v in minimal_s2t.iteritems():
		assert len(v) <= 1
	for k, v in minimal_t2s.iteritems():
		assert len(v) <= 1

	return minimal_s2t, minimal_t2s
	

def build_node_alignment_maps(source_tree, target_tree, are_aligned, minimal_only=False):
	s2t_node_alignments = defaultdict(set)
	t2s_node_alignments = defaultdict(set)
	for s_node in source_tree.nodes:
		s2t_node_alignments[s_node] = set()
		for t_node in target_tree.nodes:
				if are_aligned(s_node.span, t_node.span):
					s2t_node_alignments[s_node].add(t_node)
					t2s_node_alignments[t_node].add(s_node)	

	# The roots of the two trees are always node-aligned, even when there are no alignment links
	s2t_node_alignments[source_tree.start].add(target_tree.start)
	t2s_node_alignments[target_tree.start].add(source_tree.start)

	if minimal_only:
		s2t_node_alignments, t2s_node_alignments = minimize_alignments(source_tree, target_tree, s2t_node_alignments, t2s_node_alignments)

	return s2t_node_alignments, t2s_node_alignments


def find_aligned_spans(target_tree, source_nodes, are_aligned):
	aligned_spans = set()
	target_len = max(node.span.end for node in target_tree.nodes)
	for span_size in range(1, target_len + 1):
		for span_start in range(target_len - span_size + 1):
			for node in source_nodes:
				if are_aligned(node.span, Span(span_start, span_start + span_size)):
					aligned_spans.add(Span(span_start, span_start + span_size))
					break
	return aligned_spans

def find_child_sets(nodes, span, max_nt_count):
	for node in nodes:
		if node.span.start != span.start:
			continue

		if max_nt_count == 0 and not node.is_terminal_flag:
			continue

		if node.span.end == span.end:
			if max_nt_count >= 1 or node.is_terminal_flag:
				yield [node]
		else:
			if max_nt_count > 0:
				max_nt_siblings = max_nt_count - (0 if node.is_terminal_flag else 1)
			else:
				max_nt_siblings = 0

			for sibling_set in find_child_sets(nodes, Span(node.span.end, span.end), max_nt_siblings):
				assert len([1 for child in sibling_set if not child.is_terminal_flag]) <= max_nt_siblings
				child_set = [node] + sibling_set
				nt_count = len([1 for child in child_set if not child.is_terminal_flag])
				assert nt_count <= max_nt_count
				yield child_set

def add_t2s_virtual_nodes(target_tree, source_tree, are_aligned):
	aligned_spans = find_aligned_spans(target_tree, source_tree.nodes, are_aligned)
	for aligned_span in sorted(aligned_spans, key=lambda span: span.end - span.start):
		# Note: virtual_node should == target_tree.root when aligned_span covers the whole sentence
		virtual_node = NodeWithSpan('X', aligned_span, False, True)
		target_tree.nodes.add(virtual_node)

def add_t2s_virtual_edges(target_tree, source_tree, are_aligned):
	aligned_spans = find_aligned_spans(target_tree, source_tree.nodes, are_aligned)
	for aligned_span in sorted(aligned_spans, key=lambda span: span.end - span.start):
		# Note: virtual_node should == target_tree.root when aligned_span covers the whole sentence
		virtual_node = NodeWithSpan('X', aligned_span, False, True)
		valid_nodes = set([node for node in target_tree.nodes if node.is_terminal_flag or node.span in aligned_spans])
		for child_set in find_child_sets(valid_nodes, aligned_span, args.max_rule_size + 1000):
			if virtual_node in child_set:
				continue
			virtual_edge = Edge(virtual_node, tuple(child_set))
			target_tree.add(virtual_edge)

def build_word_alignment_maps(source_terminals, target_terminals, alignment):
	s2t_word_alignments = defaultdict(list)
	t2s_word_alignments = defaultdict(list)
	for s, t in alignment.links:
		s_node = source_terminals[s]
		t_node = target_terminals[t]
		s2t_word_alignments[s_node].append(t_node)
		t2s_word_alignments[t_node].append(s_node)
	return s2t_word_alignments, t2s_word_alignments

# Takes two hypergraphs representing source and target trees, as well as a word
# alignment, and finds all rules extractable there from.
def handle_sentence(source_tree, target_tree, alignment):
		# Identify the terminal nodes in both trees
		source_terminals = sorted([node for node in source_tree.nodes if node.is_terminal_flag], key=lambda node: node.span.start)
		target_terminals = sorted([node for node in target_tree.nodes if node.is_terminal_flag], key=lambda node: node.span.start)

		# Build word alignment maps
		s2t_word_alignments, t2s_word_alignments = build_word_alignment_maps(source_terminals, target_terminals, alignment)

		# Define this little helper function
		spans_are_aligned = lambda source_span, target_span: are_aligned(source_span, target_span, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments)
	
		source_tree.add_virtual_nodes_only(args.virtual_size, False)
		if not args.t2s:
			target_tree.add_virtual_nodes_only(args.virtual_size, False)
		else:	
			add_t2s_virtual_nodes(target_tree, source_tree, spans_are_aligned)
			#target_tree.add_virtual_nodes_only(1000, True, lambda nodes: 'X')

		# Build node alignments maps
		s2t_node_alignments, t2s_node_alignments = build_node_alignment_maps(source_tree, target_tree, spans_are_aligned, args.minimal_rules)

		# Add virtual nodes and edges to the tree structures
		source_tree.add_virtual_nodes(args.virtual_size, False)
		if not args.t2s:
			target_tree.add_virtual_nodes(args.virtual_size, False)
		else:
			aligned_spans = set([node.span for node in t2s_node_alignments.keys()])
			in_aligned_spans = lambda source_span, target_span: target_span in aligned_spans
			add_t2s_virtual_edges(target_tree, source_tree, in_aligned_spans)

		# Add composed edges to the tree structures
		if args.minimal_rules:
			s2t_aligned_nodes = set(node for node, alignments in s2t_node_alignments.iteritems() if len(alignments) > 0)
			t2s_aligned_nodes = set(node for node, alignments in t2s_node_alignments.iteritems() if len(alignments) > 0)
			source_tree.add_minimal_composed_edges(args.max_rule_size, s2t_aligned_nodes)
			if not args.t2s:
				target_tree.add_minimal_composed_edges(args.max_rule_size, t2s_aligned_nodes)
		else:
			source_tree.add_composed_edges(args.max_rule_size)
			if not args.t2s:
				target_tree.add_composed_edges(args.max_rule_size)

		# Finally extract rules
		formatter = GrexRuleFormatter() if not args.t2s else CdecT2SRuleFormatter()
		for source_node, target_nodes in s2t_node_alignments.copy().iteritems():
			for target_node in target_nodes:
				if not source_node.is_terminal(source_tree) and not target_node.is_terminal(target_tree):
					for rule in extract_rules(source_node, target_node, s2t_node_alignments, t2s_node_alignments, source_tree, target_tree, args.max_rule_size, source_terminals, target_terminals, s2t_word_alignments, t2s_word_alignments):
						print formatter.format_rule(rule).encode('utf-8')
		sys.stdout.flush()

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
