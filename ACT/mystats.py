
import math

def mean(l):
    if len(l):
        return float(sum(l))/len(l)
    else:
        return None

def std(l):
    s=mean(l)
    ss=0.0
    for e in l:
        ss+=(e-s)**2
    return math.sqrt(ss/len(l))

def percentile(l, pct):
    if not l:
        return None
    l=sorted(l)
    if pct==100: return l[-1]
    n=len(l)-1
    if n==0: return l[-1]
    fidx=n*pct/100.0
    idx=int(fidx)
    return (fidx-idx)*(l[idx+1]-l[idx])+l[idx]
