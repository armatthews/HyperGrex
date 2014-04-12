class TreeNode:
	@staticmethod
	def from_string(s):
		s = s.strip()
		if len(s) == 0:
			raise Exception('Invalid attempt to create a TreeNode from an empty string')
		elif s[0] == '(':
			return NonTerminalNode.from_string(s)
		else:
			return TerminalNode.from_string(s)

class NonTerminalNode:
	def __init__(self, label, children):
		self.label = label
		self.children = children
		self.parent = None

	@staticmethod
	def from_string(s):
		s = s.strip()
		assert s[0] == '('
		assert s[-1] == ')'
		s = s[1:-1]

		label, rest = s.split(' ', 1)
		rest = rest.strip()

		children = []
		while len(rest) > 0:
			if rest[0] == '(':
				open_parens = 0
				for i, c in enumerate(rest):
					if c == '(':
						open_parens += 1
					elif c == ')':
						open_parens -= 1
					if open_parens == 0:
						break
				child_text = rest[:i + 1]
				rest = rest[i + 1:].strip()
				child = NonTerminalNode.from_string(child_text)
				children.append(child)
			else:
				if ' ' in rest:
					word, rest = rest.split(' ', 1)
					rest = rest.strip()
				else:
					word = rest
					rest = ''
				child = TerminalNode.from_string(word)
				children.append(child)

		node = NonTerminalNode(label, children)
		for child in children:
			child.parent = node
		return node

	def find_terminals(self):
		terminals = []
		for child in self.children:
			if isinstance(child, TerminalNode):
				terminals.append(child)
			else:
				terminals += child.find_terminals()
		return terminals
				

	def __str__(self):
		return '(' + self.label.encode('utf-8') + ' ' + ' '.join(str(child) for child in self.children) + ')'

	def __repr__(self):
		return str(self)

class TerminalNode:
	def __init__(self, word):
		self.word = word

	@staticmethod
	def from_string(s):
		s = s.strip()		
		assert '(' not in s
		assert ')' not in s
		assert ' ' not in s
		return TerminalNode(s)

	def __str__(self):
		return self.word.encode('utf-8')

	def __repr__(self):
		return str(self)

if __name__ == "__main__":
	print TerminalNode.from_string('boy')
	print NonTerminalNode.from_string('(NN boy)')
	print NonTerminalNode.from_string('(NP (DT the) (NN boy))').find_terminals()
	print NonTerminalNode.from_string('(NP the (NN boy))')
	print NonTerminalNode.from_string('(NP (DT the) boy)')
