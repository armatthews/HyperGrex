#!/bin/bash
set -e

dir=$(dirname $0)
script=$(basename $0)

if [[ $# != 3 ]]; then
  echo "Usage: $0 <extracted-rules.txt> <sgt-params.txt> <tgs-params.txt>" >&2;
  exit 1;
fi
sgt=$2
tgs=$3

tmp1=$(mktemp /tmp/$script.XXXXXXX)
tmp2=$(mktemp /tmp/$script.XXXXXXX)
tmp3=$(mktemp /tmp/$script.XXXXXXX)
tmp4=$(mktemp /tmp/$script.XXXXXXX)

echo "Collecting rule counts ..." >&2
cat $1 | grep -v ^Sentence | perl -e 'while(<>){ s/ \|\|\| /\t/g; print;}' | LC_ALL=C sort -k 1,4 | perl -e 'while(<>){ s/\t/ ||| /g; print;}' | python $dir/collect_rule_counts.py > $tmp1
echo "Computing relative frequencies target | source ..." >&2
cat $tmp1 | perl -e 'while(<>){ s/ \|\|\| /\t/g; print;}' | LC_ALL=C sort -k 1,1 | perl -e 'while(<>){ s/\t/ ||| /g; print;}' | python $dir/score_phrase_sgt.py 2,2 1,1 -n phrase_tgs > $tmp2
echo "Computing relative frequencies source | target ..." >&2
cat $tmp2 | perl -e 'while(<>){ s/ \|\|\| /\t/g; print;}' | LC_ALL=C sort -k 2,2 | perl -e 'while(<>){ s/\t/ ||| /g; print;}' | python $dir/score_phrase_sgt.py 1,1 2,2 -n phrase_sgt > $tmp3
echo "Computing lexical translation probabilities ..." >&2
cat $tmp3 | python $dir/score_lex_probs.py $sgt $tgs -n lex_sgt -m lex_tgs > $tmp4
echo "Computing rarity features ..." >&2
cat $tmp4 | python $dir/log_feats.py -r
rm -f $tmp1 $tmp2 $tmp3 $tmp4
echo "SUCCESS." >&2


