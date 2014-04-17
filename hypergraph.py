#coding: utf8
import sys
import json
import Queue
from tree import TreeNode, NonTerminalNode, TerminalNode
from collections import defaultdict
from helpers import computeSpans, Span, enumerate_subsets, compute_generations

class Edge:
	def __init__(self, head, tails, is_composed=False):
		self.head = head
		self.tails = tails
		self.is_composed = is_composed
		self.composed_edges = tuple()

	def __eq__(self, other):
		return isinstance(other, Edge) and self.head == other.head and self.tails == other.tails and self.is_composed == other.is_composed and self.composed_edges == other.composed_edges

	def __ne__(self, other):
		return not (self == other)

	def __hash__(self):
		return hash((self.head, self.tails, self.is_composed, self.composed_edges))

	def __repr__(self):
		return "%sEdge(%s,%s,%s)" % ('' if not self.is_composed else 'Composed', repr(self.head), repr(self.tails), repr(self.composed_edges))

class NodeWithSpan:
	def __init__(self, label, span, is_terminal=False, is_virtual=False):
		self.label = label
		self.span = span
		self.is_terminal_flag = is_terminal
		self.is_virtual = is_virtual

	def __str__(self):
		if not self.is_terminal_flag:
			return ('%s [%d,%d)' % (self.label, self.span.start, self.span.end)).encode('utf-8')
		else:
			return ('%s' % (self.label)).encode('utf-8')

	def __repr__(self):
		return 'NodeWithSpan(%s,%s)' % (repr(self.label), repr(self.span))

	def __eq__(self, other):
		return isinstance(other, NodeWithSpan) and self.label == other.label and self.span == other.span and self.is_terminal_flag == other.is_terminal_flag and self.is_virtual == other.is_virtual

	def __ne__(self, other):
		return not (self == other)

	def __hash__(self):
		return hash((self.label, self.span, self.is_virtual, self.is_terminal_flag))

	def is_terminal(self, hg):
		return self.is_terminal_flag
		r = (self in hg.nodes) and len(hg.head_index[self]) == 0
		assert self.is_terminal_flag == r
		return r

	def get_child_edges(self, hg):
		return hg.head_index[self]

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
		self.weights = defaultdict(float)
		self.head_index = defaultdict(set)
		self.tail_index = defaultdict(set)

	def __str__(self):
		return ' '.join(str(term) for term in self.start.find_terminals(self))

	def to_tree_string_helper(self, node, alternate_labels):
		children = []
		assert(len(self.head_index[node]) <= 1)
		for edge in self.head_index[node]:
			for tail in edge.tails:
				children.append(tail)

		label = node.label if node not in alternate_labels else alternate_labels[node]
		if len(children) > 0:
			return '(%s %s)' % (label, ' '.join([self.to_tree_string_helper(child, alternate_labels) for child in children]))
		else:
			if node.is_terminal_flag:
				return label
			else:
				return '[%s]' % label

	def to_tree_string(self, alternate_labels={}):
		for node in self.nodes:
			if len(self.head_index[node]) >= 2:
				raise Exception('Invalid call to to_tree_str() on a hypergraph not representable by a tree!')
		return self.to_tree_string_helper(self.start, alternate_labels)

	def add(self, e, weight=1.0):
		assert weight >= 0.0 and weight <= 1.0
		self.nodes.add(e.head)
		self.nodes.update(e.tails)
		self.edges.add(e)
		self.head_index[e.head].add(e)
		for tail in e.tails:
			self.tail_index[tail].add(e)
		self.weights[e] += weight
		assert self.weights[e] >= 0.0 and self.weights[e] <= 1.0001
		if self.weights[e] > 1.0:
			self.weights[e] = 1.0

	def sanity_check(self):
		for edge in self.edges:
			assert edge in self.weights
			assert edge.head in self.nodes
			assert edge in self.head_index[edge.head]
			for tail in edge.tails:
				assert tail in self.nodes
				assert edge in self.tail_index[tail]

	def combine(self, other):
		self.nodes.update(other.nodes)
		for edge in other.edges:
			self.add(edge, other.weights[edge])
		#self.head_index.update(other.head_index)
		#self.tail_index.update(other.tail_index)

	@staticmethod
	def from_surface_string(string):
		words = string.strip().split()
		root_node = NodeWithSpan('ROOT', Span(0, len(words) + 1))
		x_nodes = tuple(NodeWithSpan('X', Span(i, i + 1)) for i, word in enumerate(words))
		child_nodes = tuple(NodeWithSpan(word, Span(i, i + 1), True) for i, word in enumerate(words))

		hg = Hypergraph(root_node)
		hg.nodes.update(x_nodes)
		hg.nodes.update(child_nodes)

		edge = Edge(root_node, x_nodes)
		hg.add(edge)
		for x_node, child_node in zip(x_nodes, child_nodes):
			edge = Edge(x_node, (child_node,))
			hg.add(edge)
		return hg	

	@staticmethod
	def from_tree(root, weight=1.0):
		if isinstance(root, NonTerminalNode):
			hg = Hypergraph(NodeWithSpan(root.label, root.span, False))
			child_nodes = []
			for child in root.children:
				child_hypergraph = Hypergraph.from_tree(child, weight)
				child_nodes.append(child_hypergraph.start)
				hg.combine(child_hypergraph)
			hg.add(Edge(hg.start, tuple(child_nodes)), weight)
		else:
			hg = Hypergraph(NodeWithSpan(root.word, root.span, True))
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

	def add_virtual_nodes_helper(self, edge, max_size, allow_unary_vns, label_generator):
		virtual_edges = []
		for size in range(2, max_size + 1):
			for start in range(len(edge.tails) - size + 1):	
				if start == 0 and start + size == len(edge.tails) and not allow_unary_vns:
					continue	
				covered_children = edge.tails[start : start + size]
				virtual_label = label_generator(covered_children)
				virtual_span = Span(covered_children[0].span.start, covered_children[-1].span.end)
				virtual_node = NodeWithSpan(virtual_label, virtual_span, False, True)
				#self.nodes.add(virtual_node)
				virtual_tails = edge.tails[:start] + (virtual_node,) + edge.tails[start + size:]
				virtual_edge1 = Edge(edge.head, virtual_tails)	
				virtual_edge2 = Edge(virtual_node, covered_children)
				virtual_edges.append((virtual_edge1, self.weights[edge]))
				virtual_edges.append((virtual_edge2, self.weights[edge]))
		return virtual_edges

	def add_virtual_nodes(self, max_size, allow_unary_vns, allow_recursive_vns=False, label_generator=lambda nodes:'-'.join(node.label for node in nodes)):
		done = set()
		while len(self.edges - done) > 0:
			virtual_edges = []
			for edge in self.edges - done:
				if not edge.head.is_virtual:
					virtual_edges += self.add_virtual_nodes_helper(edge, max_size, allow_unary_vns, label_generator)
				done.add(edge)

			for edge, weight in virtual_edges:
				if edge not in self.edges:
					self.add(edge, weight)

			if not allow_recursive_vns:
				break

	def compose_edge(self, edge, max_size, minimal_only=False, aligned_nodes=[]):
		tail_choices = []
		for tail in edge.tails:
			choices = []
			choices.append(((tail,), 1.0, None))

			# Composing through virtual nodes leads to us counting a ton of stuff twice.
			# Note that composing through virtual nodes will never add any extra edges either.
			if not tail.is_virtual and (not minimal_only or tail not in aligned_nodes):
				# In minimal mode, we don't compose through aligned nodes
				if not minimal_only or not tail in aligned_nodes:
					for child_edge in tail.get_child_edges(self):
						# In minimal mode, each tail must be an aligned node.
						is_valid = True
						if minimal_only:
							for t in child_edge.tails:
								if t not in aligned_nodes and not t.is_terminal_flag:
									is_valid = False
									break
						if not is_valid:
							continue
				
						if len(child_edge.tails) <= max_size - len(edge.tails) + 1:
							choices.append((child_edge.tails, self.weights[child_edge], child_edge))
			tail_choices.append(choices)

		if len(tail_choices) > max_size:
			return

		for chosen_child_edges in enumerate_subsets(tail_choices):
			new_tails = []
			new_weight = self.weights[edge]
			composed_edges = [edge]
			for tail, weight, internal_edge in chosen_child_edges:
				assert len(tail) >= 1
				new_tails += tail
				new_weight *= weight
				if internal_edge is not None:
					composed_edges.append(internal_edge)

			if len(new_tails) <= max_size:
				new_edge = Edge(edge.head, tuple(new_tails), True)
				assert len(composed_edges) > 0
				new_edge.composed_edges = tuple(composed_edges)
				if edge.tails != new_edge.tails:					
					self.add(new_edge, new_weight)	
		
	def add_composed_edges(self, max_size):
		for node in self.topsort():
			for edge in self.head_index[node].copy():
				self.compose_edge(edge, max_size)


	def add_minimal_composed_edges(self, max_size, aligned_nodes):
		for node in self.topsort():
			for edge in self.head_index[node].copy():
				self.compose_edge(edge, max_size, True, aligned_nodes)

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

	def topsort(self):
		sorted_nodes = []
		terminals = Queue.Queue()
		removed_edges = set()

		for node in self.nodes:
			if node.is_terminal(self):
				terminals.put(node)

		while not terminals.empty():
			n = terminals.get()
			sorted_nodes.append(n)
			for edge in self.tail_index[n] - removed_edges:
				m = edge.head

				needs_added = True
				for tail in edge.tails:
					if tail not in sorted_nodes:
						needs_added = False
						break

				if needs_added == True:
					removed_edges.add(edge)
				if len(self.head_index[m] - removed_edges) == 0:
					terminals.put(m)

		if len([edge for edge in self.edges if edge not in removed_edges]) > 0:
			raise Exception('Invalid attempt to topsort a hypergraph with cycles')
		else:
			assert len(self.nodes) == len(sorted_nodes)
			return sorted_nodes


if __name__ == "__main__":
	s = u'(S (NP (DT le) (JJ petit) (NN garçon)) (VP (VB a) (VBN marché) (PP (P à) (NP (DT l\') (NN école)))) (. .))'
	tree = TreeNode.from_string(s)
	computeSpans(tree)
	hg = Hypergraph.from_tree(tree)
	hg.sanity_check()
	print 'HG has %d nodes and %d edges' % (len(hg.nodes), len(hg.edges))
	print hg.to_tree_string()
	for node in hg.topsort():
		print node
	#print hg.to_json()
	hg.add_virtual_nodes(2, False)
	hg.sanity_check()
	print 'HG has %d nodes and %d edges' % (len(hg.nodes), len(hg.edges))
	#print hg.to_json()
	hg.add_composed_edges(5)
	hg.sanity_check()
	print 'HG has %d nodes and %d edges' % (len(hg.nodes), len(hg.edges))
	#print hg.to_json()	
