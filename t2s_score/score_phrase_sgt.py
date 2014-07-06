import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser('Takes rules on stdin and adds p(fields | conditions) as a feature. Assumes the input is sorted by the conditions already.')
parser.add_argument('fields', help='which fields to evaluate the probability of, in the same form as a key to sort -k. e.g. 1,1')
parser.add_argument('conditions', help='Which fields to condition on, in the same form as a key to sort -k. e.g. 1,1')
parser.add_argument('--feat_name', '-n', type=str, help='Name of the new feature')
parser.add_argument('--count_feat', type=str, default='count', help='Name of \'count\' feature')
args = parser.parse_args()

def field_getter(line, i, j):
	parts = [part.strip() for part in line.split('|||') if len(part.strip()) > 0]
	assert i >= 1 and i <= len(parts)
	assert j >= 1 and j <= len(parts)
	return tuple(parts[i - 1 : j])

def get_features(line):
	feat_str = line.split('|||')[-1].strip()
	feats = {}
	for kvp in feat_str.split():
		k, v = kvp.split('=')
		feats[k] = float(v)
	return feats

def output_rules(queue, count, name=None):
	for line in queue:
		feats = get_features(line)
		value = float(feats[args.count_feat]) / count
		print line, '%s=%f' % (name, value) if name is not None else '%f' % value

fields = map(int, args.fields.split(','))
conditions = map(int, args.conditions.split(','))

prev_key = None
count = 0
queue = set()

for line in sys.stdin:
	line = line.strip()
	key = field_getter(line, *conditions)
	value = field_getter(line, *fields)
	feats = get_features(line)

	if key != prev_key:
		if prev_key != None:
			output_rules(queue, count, args.feat_name)
		prev_key = key
		count = 0
		queue = set()
	queue.add(line)
	if args.count_feat not in feats:
		print >>sys.stderr, 'No feature named \"%s\" was found in grammar on line:\n%s' % (args.count_feat, line)
		sys.exit(1)
	count += float(feats[args.count_feat])

if prev_key != None:
	output_rules(queue, count, args.feat_name)
