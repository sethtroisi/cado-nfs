load alpha4FFS.sage
load real.roots.stats.sage

def alpha_infinity(f,alim,blim,R):
    F=R.base_ring()
    p=F.cardinality()
    t,x=R.gens()
    A.<t>=F['t']
    S.<x>=A['x']
    f=S(f)
    cancelation=real_roots_stats(f,floor(blim/2),alim,blim,R)
    return -sum(i*cancelation[i] for i in range(len(cancelation)))

def complete_alpha(f,B,alim,blim,R):
    F=R.base_ring()
    p=F.cardinality()
    t,x=R.gens()
    A.<t>=F['t']
    S.<x>=A['x']
    f=S(f)
    return alpha(f,B)+alpha_infinity(f,alim,blim,R)


"""
EXAMPLES
R.<t,x>=GF(2)['t,x']
complete_alpha(x^3+t^6+1,10,32,30,R)
complete_alpha(x^3+t^6+1,10,30,30,R)
complete_alpha(x^3+t^6+1,10,30,32,R)
"""
