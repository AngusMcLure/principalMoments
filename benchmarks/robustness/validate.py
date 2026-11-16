"""Compare reference runs and demonstrate loss of information in double input.
Usage: python benchmarks/robustness/validate.py results_dir reference160_dir
"""
import csv
import sys
from pathlib import Path
import mpmath as mp
mp.mp.dps = 160
left, right = map(Path, sys.argv[1:3])
a = list(csv.DictReader((left/'reference.csv').open()))
b = list(csv.DictReader((right/'reference.csv').open()))
assert len(a) == len(b)
functional, node_weight = mp.mpf(0), mp.mpf(0)
same_double = True
for x, y in zip(a, b):
    assert (x['id'], x['side']) == (y['id'], y['side'])
    functional = max(functional, abs(mp.mpf(x['value'])-mp.mpf(y['value'])))
    same_double &= float(x['value']) == float(y['value'])
    for field in ['nodes', 'weights']:
        xx, yy = x[field].split(';'), y[field].split(';')
        assert len(xx) == len(yy)
        node_weight = max(node_weight, max(abs(mp.mpf(i)-mp.mpf(j))
                                          for i,j in zip(xx,yy)))
assert functional < mp.mpf('1e-60') and node_weight < mp.mpf('1e-60')
assert same_double
text = (f'Recomputed all {len(a)} reference rows at 160 decimal digits.\n'
        f'Maximum value difference: {mp.nstr(functional,20)}\n'
        f'Maximum node/weight difference: {mp.nstr(node_weight,20)}\n'
        f'All rounded-to-double reference functional values identical: {same_double}\n')
(left/'reference-validation.txt').write_text(text)
print(text)
print('Input collision, distinct distributions with identical rounded moments:')
for k in [100,101,201]:
    first = mp.mpf(1)/3*(mp.mpf(1)/3)**k + mp.mpf(2)/3*(mp.mpf(2)/3)**k
    second = mp.mpf(2)/3*(mp.mpf(2)/3)**k
    assert float(first) == float(second)
    print(k, float(first), float(second), 'exact difference', mp.nstr(first-second,12))
print('Exact difference in mean/prevalence:', mp.nstr(mp.mpf(1)/9,20))
