#coding: utf8
import sys
import json
from tree import TreeNode, NonTerminalNode, TerminalNode
from collections import defaultdict
from helpers import computeSpans, Span, enumerate_subsets, compute_generations

class Edge:
	def __init__(self, head, tails, weight=1.0):
		self.head = head
		self.tails = tails
		self.weight = weight

	def __eq__(self, other):
		return isinstance(other, Edge) and self.head == other.head and self.tails == other.tails

	def __hash__(self):
		return hash((self.head, self.tails))

	def __repr__(self):
		return "Edge(%s,%s,%s)" % (repr(self.head), repr(self.tails), repr(self.weight))

class NodeWithSpan:
	def __init__(self, label, span, generation, is_terminal=False, is_virtual=False):
		self.label = label
		self.span = span
		self.generation = generation
		self.is_terminal_flag = is_terminal
		self.is_virtual = is_virtual

	def __str__(self):
		return ('%s_%d [%d,%d)' % (self.label, self.generation, self.span.start, self.span.end)).encode('utf-8')

	def __repr__(self):
		return 'NodeWithSpan(%s,%s,%s)' % (repr(self.label), repr(self.span), repr(self.generation))

	def __eq__(self, other):
		return isinstance(other, NodeWithSpan) and self.label == other.label and self.span == other.span and self.generation == other.generation and self.is_virtual == other.is_virtual

	def __hash__(self):
		return hash((self.label, self.span, self.generation, self.is_virtual))

	def is_terminal(self, hg):
		r = (self in hg.nodes) and len(hg.head_index[self]) == 0
		assert self.is_terminal_flag == r
		return r

	def get_child_sets(self, hg):
		child_sets = hg.head_index[self]
        	for child_set in child_sets:
                	children = child_set.tails
			yield children

	def find_terminals(self, hg):
		terminals = []
		if self.is_terminal(hg):
			terminals.append(self)
		else:
			for children in self.get_child_sets(hg):
				for child in children:
					terminals += child.find_terminals(hg)
				# The terminals covered by this node should be the same no matter
				# which child set we examine, so we need only examine one.
				break
		return terminals

class Hypergraph:
	def __init__(self, start):
		self.start = start
		self.nodes = set([start])
		self.edges = set()
		self.head_index = defaultdict(set)
		self.tail_index = defaultdict(set)

	def add(self, e):
		self.nodes.add(e.head)
		self.nodes.update(e.tails)
		if e not in self.edges:
			self.edges.add(e)
		else:
			if e.weight != None and e.weight != 0.0:
				for edge in self.edges:
					if hash(edge) == hash(e):
						edge.weight += e.weight
						break
		self.head_index[e.head].add(e)
		for tail in e.tails:
			self.tail_index[tail].add(e)

	def combine(self, other):
		self.nodes.update(other.nodes)
		for edge in other.edges:
			self.add(edge)
		self.head_index.update(other.head_index)
		self.tail_index.update(other.tail_index)

	@staticmethod
	def from_tree(root, generation=0):
		if isinstance(root, NonTerminalNode):
			hg = Hypergraph(NodeWithSpan(root.label, root.span, generation, False))
			child_nodes = []
			for child in root.children:
				child_hypergraph = Hypergraph.from_tree(child, generation + 1)
				child_nodes.append(child_hypergraph.start)
				hg.combine(child_hypergraph)
			hg.add(Edge(hg.start, tuple(child_nodes)))
		else:
			hg = Hypergraph(NodeWithSpan(root.word, root.span, generation, True))
		return hg

	@staticmethod
	def from_json(s):
		o = json.load(s)
		nodes = [Node(node['label']) for node in o['nodes']]
		hg = Hypergraph(start=nodes[o['root']])
		for edge in o['edges']:
			hg.add(Edge(nodes[edge['head']], tuple(nodes[tail] for tail in edge['tails'])))
		return hg

	def to_json(self):
		o = {}
		nodes = []
		nodeindex = {}
		for ni, node in enumerate(self.nodes):
			nodes.append({'label': str(node)})
			nodeindex[node] = ni
		edges = []
		for edge in self.edges:
			edges.append({'head': nodeindex[edge.head], \
				      'tails': [nodeindex[tail] for tail in edge.tails]})
		return json.dumps({'root': nodeindex[self.start], \
				   'nodes': nodes,
				   'edges': edges})

	def add_virtual_nodes(self, max_size, allow_unary_vns):
		virtual_edges = []
		for edge in self.edges:
			for size in range(2, max_size + 1):
				for start in range(len(edge.tails) - size + 1):	
					if start == 0 and start + size == len(edge.tails) and not allow_unary_vns:
						continue
					covered_children = edge.tails[start : start + size]
					virtual_label = '-'.join(child.label for child in covered_children)
					virtual_span = Span(covered_children[0].span.start, covered_children[-1].span.end)
					virtual_node = NodeWithSpan(virtual_label, virtual_span, True)
					virtual_tails = edge.tails[:start] + (virtual_node,) + edge.tails[start + size:]
					virtual_edge1 = Edge(edge.head, virtual_tails)
					virtual_edge2 = Edge(virtual_node, covered_children)
					self.nodes.add(virtual_node)
					virtual_edges.append(virtual_edge1)
					virtual_edges.append(virtual_edge2)
		for edge in virtual_edges:	
			self.add(edge)

	def compose_edge(self, edge, max_size, already_composed=set()):
		if edge in already_composed:
			return
		for tail in edge.tails:
			for child_edge in self.head_index[tail].copy():
				self.compose_edge(child_edge, max_size, already_composed)

		tail_choices = []
		for tail in edge.tails:
			choices = []
			choices.append((tail,))
			for child_edge in self.head_index[tail]:
				if len(child_edge.tails) <= max_size - len(edge.tails) + 1:
					choices.append(child_edge.tails)
			tail_choices.append(choices)

		if len(tail_choices) > max_size:
			return

		for chosen_tails in enumerate_subsets(tail_choices):
			new_tails = []
			for tail in chosen_tails:
				assert len(tail) >= 1
				new_tails += tail

			if len(new_tails) <= max_size:
				new_edge = Edge(edge.head, tuple(new_tails))
				if new_edge not in self.edges:
					self.add(new_edge)
		already_composed.add(edge)
		
	def add_composed_edges(self, max_size):	
		composed_edges = set()
		for edge in self.head_index[self.start].copy():
			self.compose_edge(edge, max_size)

	def covering_sets(self, i, j, nodes_by_start):
		r = set()
		for node in nodes_by_start[i]:
			if node.span.end > j:
				return r
			if node.span.end == j:
				r.add((node,))
			else:
				for tail in self.covering_sets(node.span.end, j, nodes_by_start):
					r.add(tuple([node] + list(tail)))
		return r

	def add_composed_edges2(self, max_size):
		sentence_length = len(self.start.find_terminals(self))
		nodes_by_span = defaultdict(set)
		nodes_by_start = defaultdict(set)
		generation = compute_generations(self.start, self)
		for node in self.nodes:
			nodes_by_span[node.span].add(node)
			nodes_by_start[node.span.start].add(node)

		for i in range(0, sentence_length):
			for j in range(i + 1, sentence_length + 1):
				if len(nodes_by_span[Span(i, j)]) == 0:
					continue
				for covering_set in self.covering_sets(i, j, nodes_by_start):
					for head in nodes_by_span[Span(i, j)]:
						valid = True
						for tail in covering_set:
							if generation[tail] <= generation[head]:
								valid = False
								break
						if valid:
							new_edge = Edge(head, tuple(covering_set))
							if new_edge not in self.edges:
								self.add(new_edge)

if __name__ == "__main__":
	s = u'(S (NP (DT le) (JJ petit) (NN garçon)) (VP (VB a) (VBN marché) (PP (P à) (NP (DT l\') (NN école)))) (. .))'
	tree = TreeNode.from_string(s)
	computeSpans(tree)
	hg = Hypergraph.from_tree(tree)
	print 'HG has %d nodes and %d edges' % (len(hg.nodes), len(hg.edges))
	#print hg.to_json()
	hg.add_virtual_nodes(2, False)
	print 'HG has %d nodes and %d edges' % (len(hg.nodes), len(hg.edges))
	#print hg.to_json()
	hg.add_composed_edges(5)
	print 'HG has %d nodes and %d edges' % (len(hg.nodes), len(hg.edges))
	#print hg.to_json()

