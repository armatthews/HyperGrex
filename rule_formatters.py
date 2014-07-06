import sys
from hypergraph import Hypergraph
from helpers import Rule

class RuleFormatter:
	pass

class GrexRuleFormatter(RuleFormatter):
	def __init__(self, suppress_counts=False):
		self.suppress_counts = suppress_counts

	def format_rule(self, rule):
		# Unpack the rule structure
		source_edge, target_edge, s2t_rule_part_map, t2s_rule_part_map, alignments, weight = rule

		# Turn the source and target LHS into a string
		lhs = '[%s::%s]' % (source_edge.head.label, target_edge.head.label)

		# Turn the source and target RHS's into strings
		source_rhs = []
		for node in source_edge.tails:
			if node.is_terminal_flag:
				source_rhs.append(node.label)
			else:
				target, index = s2t_rule_part_map[node]
				source_rhs.append('[%s::%s,%d]' % (node.label, target.label, index))

		target_rhs = []
		for node in target_edge.tails:
			if node.is_terminal_flag:
				target_rhs.append(node.label)
			else:
				source, index = t2s_rule_part_map[node]
				target_rhs.append('[%s::%s,%d]' % (source.label, node.label, index))

		# Work out the rule type
		is_abstract_source = True
		is_abstract_target = True
		is_phrase = True
		for node in source_edge.tails:
			if node.is_terminal_flag:
				is_abstract_source = False
			else:
				is_phrase = False

		for node in target_edge.tails:
			if node.is_terminal_flag:
				is_abstract_target = False

		node_types = []
		s = 'V' if source_edge.head.is_virtual else 'O'
		t = 'V' if target_edge.head.is_virtual else 'O'
		node_types.append(s + t)
		for source in source_edge.tails:
			if source.is_terminal_flag:
				continue
			target, _ = s2t_rule_part_map[source]
			s = 'V' if source.is_virtual else 'O'
			t = 'V' if target.is_virtual else 'O'
			node_types.append(s + t)

		# Figure out what type of rule this is.
		# P=phrase pair, G=partially lexicallized, A=fully abstract
		# TODO: break 'A' into three categories?
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

		parts = [rule_type, lhs, ' '.join(source_rhs), ' '.join(target_rhs), ' '.join('%d-%d' % link for link in alignments), ' '.join(node_types)]
		if not self.suppress_counts:
			feats = [1, weight]
			parts.append(' '.join(str(feat) for feat in feats))
		return ' ||| '.join(parts)

class CdecT2SRuleFormatter:
	def __init__(self):
		pass

	@staticmethod
	def build_mini_hypergraph(edges):	
		hg = Hypergraph(edges[0].head)
		edges = list(edges[:])
		while len(edges) > 0:
			edge = edges.pop()
			if len(edge.composed_edges) == 0:
				hg.add(edge)
			else:
				edges += edge.composed_edges
		return hg

	@staticmethod
	def find_terminals(hg, node=None):
		if node == None:
			node = hg.start
		incoming_edges = hg.head_index[node]
		assert len(incoming_edges) <= 1
		if len(incoming_edges) == 0:
			return [node]
		else:
			terminals = []
			for edge in incoming_edges:
				for tail in edge.tails:
					terminals += CdecT2SRuleFormatter.find_terminals(hg, tail)
			return terminals

	def format_rule(self, rule):
		# Unpack the rule structure
		source_edge, target_edge, s2t_rule_part_map, t2s_rule_part_map, alignments, weight = rule

		source_structure = source_edge.composed_edges if source_edge.is_composed else [source_edge]
		source_hg = CdecT2SRuleFormatter.build_mini_hypergraph(source_structure)
		source_tree = source_hg.to_tree_string()

		target_structure = target_edge.composed_edges if target_edge.is_composed else [target_edge]
		target_hg = CdecT2SRuleFormatter.build_mini_hypergraph(target_structure)	
		label_map = {node: '[%d]' % index for (node, (_, index)) in t2s_rule_part_map.iteritems()}
		target_side = ' '.join(label_map[node] if node in label_map else node.label for node in CdecT2SRuleFormatter.find_terminals(target_hg))
		
		feats = {'count': weight, 'sent_count': 1}
		feats = ' '.join('%s=%s' % (name, str(value)) for name, value in feats.iteritems())
		return ' ||| '.join([source_tree, target_side, ' '.join('%d-%d' % link for link in alignments), feats])

class CdecT2TRuleFormatter(CdecT2SRuleFormatter):
	def format_rule(self, rule):
		# Unpack the rule structure
		source_edge, target_edge, s2t_rule_part_map, t2s_rule_part_map, alignments, weight = rule

		source_structure = source_edge.composed_edges if source_edge.is_composed else [source_edge]
		source_hg = CdecT2SRuleFormatter.build_mini_hypergraph(source_structure)
		source_tree = source_hg.to_tree_string()

		target_label = '[%s]' % target_edge.head.label
		target_structure = target_edge.composed_edges if target_edge.is_composed else [target_edge]
		target_hg = CdecT2SRuleFormatter.build_mini_hypergraph(target_structure)
		label_map = {}
		for node, (_, index) in t2s_rule_part_map.iteritems():
			label_map[node] = '[%s,%d]' % (node.label, index)
		target_side = ' '.join(label_map[node] if node in label_map else node.label for node in CdecT2SRuleFormatter.find_terminals(target_hg))
		
		feats = {'count': weight, 'sent_count': 1}
		feats = ' '.join('%s=%s' % (name, str(value)) for name, value in feats.iteritems())
		return ' ||| '.join([target_label, source_tree, target_side, ' '.join('%d-%d' % link for link in alignments), feats])

