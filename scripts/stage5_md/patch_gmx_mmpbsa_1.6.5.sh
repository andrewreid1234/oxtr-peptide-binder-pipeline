#!/bin/bash
# Fix the LIMITATIONS.md O8 blocker in gmx_MMPBSA 1.6.5.
#
# utils.py:list2range() returns a bare '' when input_list is empty, while every
# other return path returns {'num': [...], 'string': [...]}. Callers do
# list2range(x)['string'], so an empty residue-classification list raises
# "TypeError: string indices must be integers".
#
# Idempotent. Run after any reinstall of the env.
set -eu
ENV=${1:-/scratch/drewdog/gmx_mmpbsa/env}
F="$ENV/lib/python3.10/site-packages/GMXMMPBSA/utils.py"
[ -f "$F" ] || { echo "not found: $F"; exit 1; }
if grep -q "return {'num': \[\], 'string': \[\]}" "$F"; then
  echo "already patched"; exit 0
fi
cp -n "$F" "$F.orig-1.6.5"
python3 - "$F" <<'PY'
import sys
p=sys.argv[1]; s=open(p).read()
old="    ranges = []\n    ranges_str = []\n    if not input_list:\n        return ''\n"
new=("    ranges = []\n    ranges_str = []\n    if not input_list:\n"
     "        # upstream returns a bare '' here; every other path returns a dict\n"
     "        # and callers do list2range(x)['string']. See LIMITATIONS.md O8.\n"
     "        return {'num': [], 'string': []}\n")
assert old in s, "bug signature not found -- version may differ from 1.6.5"
open(p,'w').write(s.replace(old,new,1))
PY
echo "patched $F (backup at $F.orig-1.6.5)"
