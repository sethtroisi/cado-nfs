"""
 "estimate_<name>" makes probabilistic estimations of <name>. 
"""

def estimate_alpha_l(f, l, nt):
    """
    Should compute the same thing as alpha_l, but experimentally
    """
    s=0
    n=0
    x=f.parent().gen()
    A=f.base_ring()
    F=homogenize(f)
    dl=l.degree()
    for i in range(nt):
        a=A.random_element(ceil(log(nt,2)))
        b=A.random_element(ceil(log(nt,2)))
        c=A.random_element(f(0).degree()+f.degree()*ceil(log(nt,2)))
        if gcd(a,b) == 1:
            ds=valuation(c,l)-valuation(A(F(a,b)),l)
            s+=ZZ(ds)
            n+=1
    return float(s*dl/n) 


def estimate_alpha(f,B,nt):
    """
    same as alpha(f,B), but experimentally. Only for debugging.
    """
    A=f.base_ring()
    q=A.base_ring().cardinality()
    b=floor(log(B,q))
    result=0
    for d in range(1,b+1):
        for l in A.polynomials(of_degree=d):
            if l.is_monic() and l.is_irreducible():
                result+=estimate_alpha_l(f,l,nt)
    return result

"""
alpha(S(f),4)
estimate(S(f),4,1000)
"""

