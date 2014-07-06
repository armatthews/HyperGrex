import sys
import argparse
from math import log, exp

parser = argparse.ArgumentParser()
parser.add_argument('--rarity', '-r', action='store_true', help='Convert count feature into rarity')
parser.add_argument('--count-name', '-c', type=str, default='count', help='Name of count feature')
parser.add_argument('--rarity-name', '-n', type=str, default='rarity', help='Name of rarity feature, added by this script')
args = parser.parse_args()

for line in sys.stdin:
	parts = [part.strip() for part in line.strip().split('|||')]
	feat_str = parts[-1]

	feats = {}
	for kvp in feat_str.split():
		k, v = kvp.split('=')
		feats[k] = float(v)

	new_feats = {}
	for k, v in feats.iteritems():
		if k.endswith(args.count_name) and args.rarity:
			new_feats[k[:-len(args.count_name)] + args.rarity_name] = (exp(1.0 / v) - 1) / (exp(1.0) - 1)
		elif not k.endswith('?'):
			new_feats["log_" + k] = log(v)
		else:
			new_feats[k] = v

	feat_str = ' '.join('%s=%g' % (k, v) for k,v in new_feats.iteritems())
	parts[-1] = feat_str
	print ' ||| '.join(parts)
