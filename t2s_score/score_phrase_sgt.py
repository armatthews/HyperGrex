import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser('Takes rules on stdin and adds p(fields | conditions) as a feature. Assumes the input is sorted by the conditions already.')
parser.add_argument('fields', help='which fields to evaluate the probability of, in the same form as a key to sort -k. e.g. 1,1')
parser.add_argument('conditions', help='Which fields to condition on, in the same form as a key to sort -k. e.g. 1,1')
parser.add_argument('--feat_name', '-n', type=str, help='Name of the new feature')
args = parser.parse_args()

def field_getter(line, i, j):
	parts = [part.strip() for part in line.split('|||') if len(part.strip()) > 0]
	assert i >= 1 and i <= len(parts)
	assert j >= 1 and j <= len(parts)
	return tuple(parts[i - 1 : j])

def output_rules(queue, count, name=None):
	for fields, lines in queue.iteritems():
		for line in lines:
			value = 1.0 * len(lines) / count
			print line, '%s=%f' % (name, value) if name is not None else '%f' % value

fields = map(int, args.fields.split(','))
conditions = map(int, args.conditions.split(','))

prev_key = None
count = 0
queue = defaultdict(list)

for line in sys.stdin:
	line = line.strip()
	key = field_getter(line, *conditions)
	value = field_getter(line, *fields)

	if key != prev_key:
		if prev_key != None:
			output_rules(queue, count, args.feat_name)
		prev_key = key
		count = 0
		queue = defaultdict(list)
	queue[value].append(line)
	count += 1

if prev_key != None:
	output_rules(queue, count, args.feat_name)
