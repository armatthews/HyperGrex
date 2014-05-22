dir=$(dirname $0)
tmp1=$(mktemp)
tmp2=$(mktemp)
tmp3=$(mktemp)
tmp4=$(mktemp)

if [[ $# != 3 ]]; then
  echo "Usage: bash $0 <grammar> <e2f> <f2e>" >&2;
  exit 1;
fi

e2f=$2
f2e=$3

pv $1 | sed 's/ ||| /\t/g' | LC_ALL=C sort -k 1,4 | sed 's/\t/ ||| /g' | python $dir/collect_rule_counts.py > $tmp1
pv $tmp1 | sed 's/ ||| /\t/g' | LC_ALL=C sort -k 1,1 | sed 's/\t/ ||| /g' | python $dir/score_phrase_sgt.py 2,2 1,1 -n phrase_tgs > $tmp2
pv $tmp2 | sed 's/ ||| /\t/g' | LC_ALL=C sort -k 2,2 | sed 's/\t/ ||| /g' | python $dir/score_phrase_sgt.py 1,1 2,2 -n phrase_sgt > $tmp3
pv $tmp3 | python $dir/score_lex_probs.py $e2f $f2e -n lex_sgt -m lex_tgs > $tmp4
pv $tmp4 | python $dir/log_feats.py -r
