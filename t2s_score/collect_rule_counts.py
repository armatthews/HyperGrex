import sys
import argparse
import collections
import exceptions
from collections import defaultdict

parser = argparse.ArgumentParser()
args = parser.parse_args()

def num(s):
	try:
		return int(s)
	except exceptions.ValueError:
		return float(s)

# Accumulator for multiple copies of the same rule with different word aligns:
current_rule = None
current_scores = defaultdict(float)
align_counts = collections.Counter()
warned = False

def output_current_rule():
	align_string = ' '.join('%d-%d/%d' % (i, j, count) for ((i, j), count) in sorted(align_counts.iteritems()))
	lhs, rhs = current_rule
	print ' ||| '.join((lhs, rhs, align_string, ' '.join('%s=%s' % (k, v) for k, v in current_scores.iteritems()))).encode('utf-8')

def parse_scores(s):
	scores = defaultdict(float)
	for kvp in s.split():
		kvp = kvp.strip()
		if not kvp:
			continue
		k, v = kvp.split('=')
		v = float(v)
		assert k not in scores or scores[k] == v
		scores[k] = v
	return scores

# Read rule instances from standard in, one per line:
# NOTE: This assumes rules are sorted by fields 1-4!
for line in sys.stdin:
	line = line.decode('utf-8').strip()
	parts = [part.strip() for part in line.split('|||')]

	# Break rule line into fields:
	lhs, rhs, aligns, scores = parts
	key = (lhs, rhs)	
	scores = parse_scores(scores)
	aligns = aligns.split()

	# If different rule than previous one, write out old rule with counts:
	if key != current_rule:
		if current_rule is not None:
			output_current_rule()
		# Reset accumulators:
		current_rule = key
		current_scores = defaultdict(float)
		align_counts = collections.Counter()

	for k, v in scores.iteritems():
		if not warned and k not in ['count', 'sent_count']:
			print >>sys.stderr, 'Unrecognized features in grammar. This script will sum the value of all features, per type'
			warned = True
		current_scores[k] += v

	# Add this rule's word alignment link counts to the accumulators:
	for link in aligns:
		i, j = (num(i) for i in link.split('-'))
		align_counts[(i, j)] += scores['count']

# At end, write out final rule still in accumulators:
if current_rule is not None:
	output_current_rule()
