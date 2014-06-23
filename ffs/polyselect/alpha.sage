""" 
 See estimate_alpha.sage for a probabilistic approach 
 and approx_alpha for a factor base approach.
"""


load init_rings.sage 
load l_roots.sage

def average_valuation_affine(f,l,depth=0,prec=10):
    """
    returns the average l-valuation of the polynomial f. Works recursively.
    """
    S= f.parent();A=f.base_ring();x = S.gen();
    F=A.base_ring();q=F.cardinality();t=A.gen()
    l=A(l)
    nl=q^l.degree()
    if not A(l).is_irreducible():
        print("Error: l="+str(l)+" must be irreducible")
        return 0
    fcont=gcd(f.coefficients())
    v = valuation (fcont, l)
    fv = S(f/l^v)
    df=fv.derivative()
    result=v
    for r in roots_mod_l(fv,A(l)):
        assert fv(A(r)) % l == 0
        if df(A(r)) % l != 0:
            result += 1/(nl-1)
        else:
            f2 = fv(r + l*x)
            if depth < prec:
                result+=average_valuation_affine(S(f2), A(l),depth+1)/nl
    return result


def average_valuation_homogeneous_coprime(f,l):
    """
    returns the average l-valuation of F(a,b) for coprime a,b . Projective
    roots are counted as well.
    """
    S= f.parent();A=f.base_ring();x = S.gen();
    F=A.base_ring();q=F.cardinality();t=A.gen()
    l=A(l)
    nl=q^l.degree()
    aff=average_valuation_affine(S(f), A(l))
    proj=average_valuation_affine(S((S(f).reverse())(l*x)), A(l))
    return aff*nl/(nl+1)+proj/(nl+1)
    # Lemma: homogeneous_val = Prob(b=0)*proj + Prob(b!=0)*affine

def alpha_l(f,l):
    """
    computes the contribution of l in alpha(f)
    """
    A=f.base_ring();q=A.base_ring().cardinality()
    l=A(l)
    nl=q^l.degree()
    return l.degree()*(1/(nl-1)-average_valuation_homogeneous_coprime(f,l))

def alpha(f,B,B0=1):
    """
    computes alpha(f), by summing contribution of l such that B0+1<=Norm(l)<= B
    """
    A=f.base_ring();q=A.base_ring().cardinality()
    b=floor(log(B,q))
    b0=ceil(log(B0,q))
    result=0
    for d in range(b0+1,b+1):
        for l in A.polynomials(of_degree=d):
            if l.is_monic() and l.is_irreducible():
                result+=alpha_l(f,A(l)).n()
    return result

"""
 A.<t>=GF(2)['t']
 S.<x>=A['x']
 # or p=2; load init_rings.sage
 f=x^6+(t+t^7)*x+(t^2+t+1)
 alpha(S(f),1000)
 p=7
 load init_rings.sage
 f=x^6+3*t
 alpha(S(f),1000)
"""
