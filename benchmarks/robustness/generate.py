"""Generate feasible exact fixtures and 100-decimal references; requires mpmath.
Run from the repository root. References use full-precision input moments, while
R receives their nearest double values. This deliberately measures input and
solver error together; it cannot certify sharp bounds for noisy moments.
"""
import csv, math, os
from pathlib import Path
import mpmath as mp
mp.mp.dps = int(os.environ.get("ROBUSTNESS_DPS", "100"))
out = Path(os.environ.get('ROBUSTNESS_RESULTS', 'benchmarks/robustness/results'))
out.mkdir(parents=True, exist_ok=True)
meta, inputs, reference = [], [], []

def gaussian(t, n):
    h = mp.matrix([[t[i+j] for j in range(n)] for i in range(n)])
    b = mp.matrix([[t[i+j+1] for j in range(n)] for i in range(n)])
    inv = mp.cholesky(h).T**-1
    j = inv.T*b*inv
    eig, vec = mp.eigsy((j+j.T)/2)
    return [(eig[i], t[0]*vec[0,i]**2) for i in range(n)]

def principal(m):
    v = len(m)-1
    if v%2:
        n=(v+1)//2
        low=gaussian(m,n)
        g=gaussian([m[k+1]-m[k+2] for k in range(2*(n-1))], n-1) if n>1 else []
        up=[(x,w/(x*(1-x))) for x,w in g]
        right=m[1]-sum(x*w for x,w in up)
        up=[(mp.mpf(0),1-sum(w for x,w in up)-right)]+up+[(mp.mpf(1),right)]
    else:
        n=v//2
        g=gaussian(m[1:],n)
        low=[(x,w/x) for x,w in g]
        low=[(mp.mpf(0),1-sum(w for x,w in low))]+low
        g=gaussian([m[k]-m[k+1] for k in range(2*n)],n)
        up=[(x,w/(1-x)) for x,w in g]
        up += [(mp.mpf(1),1-sum(w for x,w in up))]
    return low,up

def record(family, label, api, orders, m, s=25, pools=0, units=0, refs=None):
    cid=len(meta)+1
    meta.append(dict(id=cid,family=family,label=label,api=api,d=len(orders),s=s,
                     pools=pools,units=units,max_order=max(orders),
                     spread=float(max(orders)/min(orders))))
    for k,z in zip(orders,m):
        inputs.append(dict(id=cid,order=str(k),moment=mp.nstr(z,105)))
    for side,rule in (refs or {}).items():
        # q represents positive-pool probability for the fixed-pool solver.
        # x represents negative-individual probability for the mixed solver.
        f=(lambda x: -mp.expm1(mp.log1p(-x)/s)) if api=='fixed' else (lambda x: 1-x)
        value=sum(w*f(x) for x,w in rule)
        reference.append(dict(id=cid,side=side,value=mp.nstr(value,100),
                              nodes=';'.join(mp.nstr(x,100) for x,w in rule),
                              weights=';'.join(mp.nstr(w,100) for x,w in rule)))
    return cid

# Consecutive-moment sweeps, including a prevalence distribution concentrated
# near zero and realistic pooled moments from Beta mixing for individual P.
for label, kind, args in [
    ('Q_uniform','qbeta',('1','1')),
    ('Q_beta_2_3','qbeta',('2','3')),
    ('Q_U_shaped','qbeta',('.2','.2')),
    ('Q_rare','qbeta',('.01','10')),
    ('P_beta_p01_rho01_s25','pbeta',('.09','8.91',25)),
    ('P_beta_p001_rho01_s50','pbeta',('.009','8.991',50))]:
    a,b=map(mp.mpf,args[:2]); s=args[2] if kind=='pbeta' else 25
    def moment(k):
        if kind=='qbeta': return mp.rf(a,k)/mp.rf(a+b,k)
        return sum((-1)**j*mp.binomial(k,j)*mp.rf(b,s*j)/mp.rf(a+b,s*j)
                   for j in range(k+1))
    moments=[mp.mpf(1)]+[moment(k) for k in range(1,101)]
    for v in list(range(1,25))+[28,32,40,60,100]:
        rules=principal(moments[:v+1]) if v<=32 else None
        record('fixed',label,'fixed',list(range(1,v+1)),moments[1:v+1],s,v,s*v,
               dict(zip(['lower','upper'],rules)) if rules else None)
    print('generated',label,flush=True)

# No independent maximum pool size: isolate size from moment order with a
# fixed distribution of Q and use a known three-moment rule as the oracle.
for v in [1,2,3]:
    m=[mp.mpf(1)]+[mp.mpf(1)/(k+1) for k in range(1,v+1)]
    rules=principal(m)
    for s in [1,25,1000,10**6,10**12]:
        record('pool_size','Q_uniform','fixed',list(range(1,v+1)),m[1:],s,v,s*v,
               dict(zip(['lower','upper'],rules)))

# Actual mixed designs: attainable exponents are subset sums with multiplicity.
def exponents(sizes, counts):
    a={0}
    for s,v in zip(sizes,counts): a={x+s*r for x in a for r in range(v+1)}
    return sorted(a-{0})

designs=[]
for sizes in [[2,3],[10,40],[25,26]]:
    for n in [1,2,3,4]: designs.append((sizes,[n,n]))
for sizes in [[2,101],[100,101],[2,3,7],[2,3,7,15],[2,4,8,16,32]]:
    designs.append((sizes,[1]*len(sizes)))
for sizes,counts in designs:
    orders=exponents(sizes,counts); d=len(orders); r=(d+1)//2
    for side in ['lower','upper']:
        nr=r-(side=='upper') if d%2 else r
        ends=([0,1] if side=='upper' else []) if d%2 else [int(side=='upper')]
        # Known principal measure, with rational nodes and positive weights.
        x=[mp.mpf(i)/(nr+1) for i in range(1,nr+1)]+list(map(mp.mpf,ends))
        w=[mp.mpf(i+1) for i in range(len(x))]; w=[z/sum(w) for z in w]
        m=[sum(wi*xi**k for xi,wi in zip(x,w)) for k in orders]
        rule=sorted(zip(x,w))
        record('mixed_known',f'{sizes}x{counts}_{side}','mixed',orders,m,
               pools=sum(counts),units=sum(s*v for s,v in zip(sizes,counts)),
               refs={side:rule})

# Same dimension/exponent ratios, hugely different individual totals.
for factor in [1,10,100,1000,10**6]:
    orders=[factor*k for k in [2,3,5]]
    x=[mp.mpf(1)/4,mp.mpf(3)/4]; w=[mp.mpf(1)/3,mp.mpf(2)/3]
    x=[z**(mp.mpf(1)/factor) for z in x]
    m=[sum(wi*xi**k for xi,wi in zip(x,w)) for k in orders]
    record('scaled_orders',f'scale_{factor}','mixed',orders,m,pools=2,units=5*factor,
           refs={'lower':list(zip(x,w))})

# Actual theory distribution: individual prevalence Beta(.09,8.91).
for sizes,counts in [([10,40],[1,1]),([10,40],[2,2]),([25,26],[1,1]),
                     ([25,26],[2,2]),([2,3],[3,3])]:
    orders=exponents(sizes,counts)
    m=[mp.rf(mp.mpf('8.91'),k)/mp.rf(9,k) for k in orders]
    record('mixed_beta',f'{sizes}x{counts}','mixed',orders,m,pools=sum(counts),
           units=sum(s*v for s,v in zip(sizes,counts)))

for name,rows in [('cases',meta),('inputs',inputs),('reference',reference)]:
    with (out/(name+'.csv')).open('w') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
print('cases',len(meta),flush=True)
