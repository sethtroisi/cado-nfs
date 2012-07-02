# sage: load init_rings.sage
# sage: f=S(x^5-x+t)
# sage: l,makefblist=prepare_small_primes_stats(f,10)
# Preparing root property for primes up to degree 10.
# Names l and makefblist are COMPULSORY because they are used as global
# variables by beta.
# sage: beta(f,15,10,9)
# Approximation of the probability that
# for polynomials a,b of degree in [0,10] respectively [0..9]
# the norm of a-b*theta is 15-smooth.
# sage: beta(f,25,28,28)
# Variables l and makefblist depend on f only, they are reused for different
# parameters.


load real.roots.stats.sage
load makefb.sage

def prob_val(r,pp,max_prec=None,At=GF(2)['t,x']):
    p=At.characteristic()
    if r == []:
        return []
    n=len(r)
    max_val=ceil(max_prec/pp.degree())
    if max_val == None:
        max_val=max([e[3] for e in r])
    prob=[QQ(1)]+[QQ(0) for i in range(max_val)]
    # prob[val] = prob to have the exact valuation val 
    for val in range(max_val+1):
        for i in range(n):
            if r[i][4] == val:
                e=r[i]
                pk=e[0]
                k=ZZ(pk.degree()/pp.degree())
                vali=e[3]
                prob[vali]+=1/p^(k*pp.degree())
                prob[val]-=1/p^(k*pp.degree())
    return prob



def prepare_small_primes_stats(f,MAX_DEGREE,At=GF(2)['t']):
    p=At.characteristic()
    l=Primes(At,p^(MAX_DEGREE+1)) 
    makefblist=[]
    for i in range(len(l)):
        pp=l[i]
        r=all_roots(S(f),pp,20)
        prob=prob_val(r,pp,20,At)
        makefblist.append(prob)
    return l,makefblist

# always use with l,makefblist=prepare_small_primes_stats
def beta(f,b,A,B,R=GF(2)['t,x'],lenl=0):
    # lenl = number of small primes
    f=R(f)
    t,x=R.gens()
    p=R.characteristic()
    At.<t>=GF(p)['t']
    S.<x>=At['x']
    result=0
    d=S(f).degree()
    fc=S(f).coeffs()
    top_degree=max([ZZ(i*A)+(d-i)*B+At(fc[i]).degree() for i in range(d+1)])
    deg_prob=[0]*(top_degree+1)
    new_deg_prob=[0]*(top_degree+1)
    for dega in range(A+1):
        for degb in range(B+1):
            cancelation=real_roots_stats(f,top_degree,dega,degb,R,dega,degb)
            local_max_degree=max([ZZ(j*dega)+(d-j)*degb+At(fc[j]).degree()  for j in range(d+1)])
            for i in range(min(local_max_degree+1,len(cancelation))):
                deg_prob[top_degree-local_max_degree+i]+=p^(dega+degb)/p^(A+B+2)*cancelation[i]
    deg_prob[top_degree]=1-sum([deg_prob[j] for j in range(top_degree+1)]) 
    for i in range(lenl):
        pp=At(l[i])
        probi=makefblist[i]
        new_deg_prob=[el for el in deg_prob]
        for val in range(1,len(probi)):
            for j in range(top_degree):
                new_deg_prob[min(j+val*pp.degree(),top_degree)]+=deg_prob[j]*probi[val]
                new_deg_prob[j]-=deg_prob[j]*probi[val]
        deg_prob=[el for el in new_deg_prob]
    assert abs(sum([deg_prob[j] for j in range(top_degree+1)]) - 1) < 0.01
    return sum([dickman_rho((top_degree-j)/b)*deg_prob[j] for j  in    range(top_degree+1)] ) 
    #  dickman_rho to be replaced with psi(x,y,z) = probability that a
    #  polynomial of degree x has all irreducible divisors of degree in [z,y] 
