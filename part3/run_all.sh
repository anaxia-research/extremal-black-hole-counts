#!/bin/sh
# Runs every check of Part III and writes the logs.
cd "$(dirname "$0")/src"
for s in taylor_check cong stest lemma_check nscompare consts cardy apriori envelope unif lc v2check tables3 figdata; do
  echo "== $s"; python3 $s.py > ../logs/$s.log 2>&1 || echo "FAILED $s"
done
