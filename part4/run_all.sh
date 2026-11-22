#!/bin/sh
# Runs every check of Part IV and writes the logs.
cd "$(dirname "$0")/src"
run() { name=$1; shift; echo "== $name"; python3 "$@" > ../logs/$name.log 2>&1 || echo "FAILED $name"; }
run constraints constraints.py
run zero4 zero4.py 80
run sseries sseries.py 20000
run transf transf.py
run analytic analytic.py
run identities identities.py
run polarlemma polarlemma.py
run cuspcheck cuspcheck.py
run exist exist.py
run polar polar.py 1 70
run tables4 tables4.py
run genus_rows genus_rows.py
run integral integral4.py
run signs signs.py
