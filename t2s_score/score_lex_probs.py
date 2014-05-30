import re
import sys
import argparse
import math
from collections import defaultdict

def read_lex_file(stream, source_first):
	probs = defaultdict(float)
	are_log_probs = True
	lc = 0
	for line in stream:
		lc += 1
		src, tgt, prob = line.strip().split()[:3]
		if src == '<eps>': src = 'NULL'
		if tgt == '<eps>': tgt = 'NULL'
		prob = float(prob)
		if lc == 1 and prob > 0.0:
			are_log_probs = False
		if are_log_probs:
			prob = math.exp(prob)
		if not source_first:
			src, tgt = tgt, src
		probs[src, tgt] = prob
	return probs

parser = argparse.ArgumentParser()
parser.add_argument('sgt_lex_file')
parser.add_argument('tgs_lex_file')
parser.add_argument('--sgt_name', '-n', type=str, default='lex_sgt', help='Name of the lexical source given target feature (e2f)')
parser.add_argument('--tgs_name', '-m', type=str, default='lex_tgs', help='Name of the lexical target given source feature (f2e)')
parser.add_argument('--epsilon', '-e', type=float, default=0.000000001, help='Minimal probability value to give to unseen alignments')
args = parser.parse_args()

sgt_probs = read_lex_file(open(args.sgt_lex_file), True)
tgs_probs = read_lex_file(open(args.tgs_lex_file), False)

for line in sys.stdin:
	src, tgt, align, feats = [part.strip() for part in line.split('|||')]
	src_terminals = [token.split(')', 1)[0] for token in src.split() if token[0] != '(']
	tgt_terminals = tgt.split()

	sgt_probs_by_index = [0.0 for _ in src_terminals]
	tgs_probs_by_index = [0.0 for _ in tgt_terminals]
	src_counts = [0 for _ in src_terminals] + [0]
	tgt_counts = [0 for _ in tgt_terminals] + [0]
	used_src = [False for _ in src_terminals] + [False]
	used_tgt = [False for _ in tgt_terminals] + [False]
	src_is_nt = [re.match(r'\[.+\]', w) is not None for w in src_terminals] + [False]
	tgt_is_nt = [re.match(r'\[.+\]', w) is not None for w in tgt_terminals] + [False]

	for link in align.split():
		link, count = link.split('/')
		count = float(count)
		i, j = [int(i) for i in link.split('-')]
		src_word = src_terminals[i]
		tgt_word = tgt_terminals[j]

		if src_is_nt[i] != tgt_is_nt[j]:
			print >>sys.stderr, 'Mismatched arities detected in rule:'
			print >>sys.stderr, line.strip()
			print >>sys.stderr, 'Consider escaping brackets in your input to something like -LSB- and -LSB-'
		assert src_is_nt[i] == tgt_is_nt[j]
		if not src_is_nt[i]:
			if (src_word, tgt_word) not in sgt_probs or (src_word, tgt_word) not in tgs_probs:
				print >>sys.stderr, 'WARNING: Unseen word pair:', src_word, 'and', tgt_word, 'found in line "%s"' % line.strip()
			#assert ((src_word, tgt_word) in sgt_probs)
			#assert ((src_word, tgt_word) in tgs_probs)
			sgt_probs_by_index[i] += sgt_probs[src_word, tgt_word] * count
			tgs_probs_by_index[j] += tgs_probs[src_word, tgt_word] * count
	
			src_counts[i] += count
			tgt_counts[j] += count
			used_src[i] = True
			used_tgt[j] = True

	# Scale total probability for each word by how often it was aligned
	for i, v in enumerate(sgt_probs_by_index):
		if src_counts[i] != 0:
			sgt_probs_by_index[i] = sgt_probs_by_index[i] / src_counts[i]

	for j, v in enumerate(tgs_probs_by_index):
		if tgt_counts[j] != 0:
			tgs_probs_by_index[j] = tgs_probs_by_index[j] / tgt_counts[j]

	# Add in probability for each unused terminal aligning to NULL
	sgt_probs_by_index.append(0.0)
	tgs_probs_by_index.append(0.0)
	for i, src_word in enumerate(src_terminals):
		if not used_src[i] and not src_is_nt[i]:
			sgt_probs_by_index[i] += sgt_probs[src_word, 'NULL']
			tgs_probs_by_index[-1] += tgs_probs[src_word, 'NULL']
			used_src[i] = True
			used_tgt[-1] = True
	for j, tgt_word in enumerate(tgt_terminals):
		if not used_tgt[j] and not tgt_is_nt[j]:
			sgt_probs_by_index[-1] += sgt_probs['NULL', tgt_word]
			tgs_probs_by_index[j] += tgs_probs['NULL', tgt_word]
			used_src[-1] = True
			used_tgt[j] = True

	# Compute final lex probs
	total_sgt = 1.0
	for i, v in enumerate(sgt_probs_by_index):
		if used_src[i]:
			total_sgt *= v

	total_tgs = 1.0
	for j, v in enumerate(tgs_probs_by_index):
		if used_tgt[j]:
			total_tgs *= v

	if total_sgt < args.epsilon:
		print >>sys.stderr, 'WARNING: %s underflow. Using epsilon value of %g' % (args.sgt_name, args.epsilon)
		total_sgt = args.epsilon

	if total_tgs < args.epsilon:
		print >>sys.stderr, 'WARNING: %s underflow. Using epsilon value of %g' % (args.tgs_name, args.epsilon)
		total_tgs = args.epsilon

	print ' ||| '.join([src, tgt, feats]), '%s=%g' % (args.sgt_name, total_sgt), '%s=%g' % (args.tgs_name, total_tgs)
