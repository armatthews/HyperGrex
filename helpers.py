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

def computeSpans(node, start=0):
	if isinstance(node, NonTerminalNode):
		end = start
		for child in node.children:
			computeSpans(child, end)
			end = child.span.end
		node.span = Span(start, end)
	else:
		node.span = Span(start, start + 1)

# Precondition: both entire parse trees must be annotated with spans.
def computeProjectedCoverageVector(node, isSource, alignment):
	node.coverage = 0
	# For a terminal node, the vector comes directly from alignments:
	if isinstance(node, TerminalNode):
		if isSource:
			links = alignment.from_source_word(node.span.start)
		else:
			links = alignment.from_target_word(node.span.start)
		for link in links:
			node.coverage |= (1 << link)
	else:
		for child in node.children:
			computeProjectedCoverageVector(child, isSource, alignment)
			node.coverage |= child.coverage

def computeProjectedComplementVector(node, parent=None):
	if parent == None:
		node.complement = 0
	else:
		node.complement = parent.complement
		siblings = [child for child in parent.children if child != node]
		for sibling in siblings:
			node.complement |= sibling.coverage

	if isinstance(node, NonTerminalNode):
		for child in node.children:
			computeProjectedComplementVector(child, node)

def getCoverageSpanMask(coverage):
        if coverage == 0:
                return 0

        coverage_lsb = coverage & ~(coverage - 1)
        while coverage != 0:
                coverage_msb = coverage
                coverage = coverage & (coverage - 1)

        coverage_span = (2 * coverage_msb - 1) & ~(coverage_lsb - 1)
        return coverage_span

def getCoverageSpan(coverage):
        if coverage == 0:
                return Span(0, 0)

        start = 0
        while (coverage & 1) == 0:
                start += 1
                coverage >>= 1

        end = start
        while coverage != 0:
                end += 1
                coverage >>= 1

        return Span(start, end)

def getSpanMask(span):
	r = 0
	for i in range(span.start, span.end):
		r |= (1 << i)
	return r

def buildSpanToNodeMap(tree):
	span2node = defaultdict(set)
	span2node[tree.span].add(tree)
	if isinstance(tree, NonTerminalNode):
		for child in tree.children:
			child_map = buildSpanToNodeMap(child)
			for key, value in child_map.iteritems():
				span2node[key] |= value
	return span2node

def getNodesWithinSpan(node, span):
	matches = []
	if node.span.start >= span.start and node.span.end <= span.end:
		matches.append(node)
	elif node.span.end <= span.start or node.span.start >= span.end:
		pass
	else:
		if isinstance(node, NonTerminalNode):
			for child in node.children:
				matches += getNodesWithinSpan(child, span)
	return matches

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

