import sys
from tree import TreeNode, NonTerminalNode, TerminalNode
from itertools import izip
from collections import defaultdict, namedtuple

def Enum(**enums):
	return type('Enum', (), enums)

NodeAlignmentType = Enum(SRC_GROWN=1, TGT_GROWN=2, T2T=4, T2S=8, S2T=16, T2TS=32, TS2T=64, TS2TS=128, T2P=256, P2T=512)
Span = namedtuple('Span', 'start, end')
NodeAlignment = namedtuple('NodeAlignment', 'source_span, target_span, from_source, types')

class Alignment:
	def __init__(self, links):
		self.links = links

	@staticmethod
	def from_string(s):
		links = []
		for pair in [link_text.strip() for link_text in s.strip().split(' ') if len(link_text) > 0]:
			f, e = pair.split('-')
			f = int(f)
			e = int(e)
			links.append((f,e))
		return Alignment(links)

	def __str__(self):
		return ' '.join('%d-%d' % link for link in self.links)

	def from_source_word(self, w):
		return [j for (i, j) in self.links if i == w]

	def from_target_word(self, w):
		return [i for (i, j) in self.links if j == w]

class ProbabilisticAlignment:
	def __init__(self):
		self.matrix = defaultdict(float) # stores p(s, t)
		self.source_sums = defaultdict(float) # stores sum over t of p(t|s)
		self.target_sums = defaultdict(float) # stores sum over s of p(s|t)

	@staticmethod
	def from_string(s):
		alignment = ProbabilisticAlignment()
		for link_text in [link_text for link_text in s.strip().split(' ') if len(link_text) > 0]:
			if '/' in link_text:
				prob = link_text[link_text.find('/') + 1:]
				prob = float(prob)
				link_text = link_text[:link_text.find('/')]
			else:
				prob = 1.0
			i, j = link_text.split('-')
			i, j = int(i), int(j)
			assert (i, j) not in alignment.matrix
			alignment.matrix[i, j] = prob
			alignment.source_sums[i] += prob
			alignment.target_sums[j] += prob
		return alignment

	def prob(self, i, j):
		assert i != None or j != None

		if j == None:
			return 1.0 - self.source_sums[i]
		elif i == None:
			return 1.0 - self.target_sums[j]
		else:
			return self.matrix[i, j]

	def __str__(self):
		return ' '.join(['%d-%d/%f' % (i, j, p) for (i, j), p in self.matrix.iteritems()])
			

def computeSpans(node, start=0):
	if isinstance(node, NonTerminalNode):
		end = start
		for child in node.children:
			computeSpans(child, end)
			end = child.span.end
		node.span = Span(start, end)
	else:
		node.span = Span(start, start + 1)

# items is a list of lists.
# this function will enumerate all the ways of choosing one item
# from each of the inner lists.
def enumerate_subsets(items):
	for item in items:
		if len(item) == 0:
			return

	indices = [0 for item in items]
	Done = False
	while not Done:
		Done = True
		yield [item[i] for (item, i) in zip(items, indices)]
		for i in range(len(items)):
			if indices[i] < len(items[i]) - 1:
				indices[i] += 1
				Done = False
				break
			else:
				indices[i] = 0

# Computes a dictionary whos keys are nodes in a hypergraph, and whose keys
# are integers, representing distance from the root of the hypergraph.
def compute_generations(root):
	generations = {}
	for node in reversed(root.topsort()):
		if len(root.tail_index[node]) == 0:
			generations[node] = 0
		else:
			generations[node] = min(generations[head] for head in [edge.head for edge in root.tail_index[node]]) + 1
	return generations
