HyperGrex
=========

A hypergraph-based syntactic translation grammar extractor for use with [cdec](http://www.cdec-decoder.org/) and similar translation systems.

## Supported functionality

 * extract various kinds (tree-to-string, tree-to-tree, string-to-tree) of tree transduction rules from aligned parallel corpora with parses on one or both sides
 * tranduction rules can be *minimal* or *composed*, with limits on the size and complexity of the rules
 * extract rules from parse forests or *k*-kbest lists of parses
 * score extracted rules with a variety of standard features

## Example tree-to-string extraction

    python hg_rule_extractor.py test_data/test.fr test_data/test.en test_data/test.al --t2s -m -s 1000 > rules.t2s

The options are:
 * `test_data/test.fr` is the source side of the bitext, parsed, one tree per line
 * `test_data/test.en` is the target side of the bitext, one sentence per line (not parsed)
 * `--t2s` indicates that xRs rules should be extracted
 * `-m` indicates that *minimal* (non-composed) rules should be extracted
 * `-s 1000` indicates that rules my have up to 1000 symbols in them (effectively, this disables any size-based filtering)

## For further information

 * For information on tree-to-string (xRs) translation rules, see
   * [this paper by GHKM](http://www.aclweb.org/anthology/N/N04/N04-1035.pdf), and
   * [this paper by Liang Huang](http://www.cis.upenn.edu/~lhuang3/amta06-sdtedl.pdf)
   * [this follow up paper by Galley et al.](http://www.cs.columbia.edu/nlp/papers/2006/galley_al_06.pdf).
   * [this dissertation by Jon May](http://www.isi.edu/~jonmay/pubs/thesis_ss.pdf)
 * For more information on the supported tree-to-tree formalism, see
   * this paper

This software is a rewrite of the [Grex grammar extractor](https://github.com/ghanneman/grex)

