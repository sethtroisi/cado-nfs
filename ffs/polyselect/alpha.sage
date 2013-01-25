load init_rings.sage 
load tools.sage  

"""
 The main function is alpha(f,B); put B=p^k for alpha_l with l of degree  <=k.
 See estimate_alpha.sage for a probabilistic approach 
 and approx_alpha for a factor base approach.
"""

def average_valuation_affine(f,l,depth=0,MAX_PREC=10):
    """
    returns the average l-valuation of the polynomial f. Works recursively.
    """
    if len(l.factor()) != 1:
        print("Error: l="+str(l)+" must be irreducible")
        return 0
    fcont=gcd(f.coefficients())
    v = valuation (fcont, l)
    S= f.parent()
    A=f.base_ring()
    x = S.gen()
    F=A.base_ring()
    q=F.cardinality()
    K.<wK>=GF(q^l.degree())
    Ky.<y>=K['y']
    fv = S(f/l^v)
    Q=fv.derivative()
    for r in fmodp_roots(fv.coeffs(),l,K,Ky):
        assert fv(r) % l == 0
        if Q(r) % l != 0:
            v += QQ(1/(Norm(l)-1))
        else:
            f2 = fv(r + l*x)
            if depth < MAX_PREC:
                v += QQ(average_valuation_affine(f2, l,depth+1)/Norm(l))
    return v


def average_valuation_homogeneous_coprime(f,l):
    """
    returns the average l-valuation of F(a,b) for coprime a,b . Projective
    roots are counted as well.
    """
    S= f.parent()
    x = S.gen()
    affine_average=average_valuation_affine(f, l)
    proj_average=QQ(average_valuation_affine((f.reverse())(l*x), l))
    return (affine_average*Norm(l)/(Norm(l)+1)+proj_average/(Norm(l)+1))
    # Lemma: homogenous_val = Prob(b=0)*proj + Prob(b!=0)*affine

def alpha_l(f,l):
    """
    computes the contribution of l in alpha(f)
    """
    return float((1/(Norm(l)-1)-average_valuation_homogeneous_coprime(f,l))*l.degree())

def alpha(f,B,B0=1):
    """
    computes alpha(f), by summing contribution of l such that B0+1<=Norm(l)<= B
    """
    A=f.base_ring()
    q=A.base_ring().cardinality()
    b=floor(log(B,q))
    b0=ceil(log(B0,q))
    result=0
    for d in range(b0+1,b+1):
        for l in A.polynomials(of_degree=d):
            if l.is_monic() and l.is_irreducible():
                result+=alpha_l(f,l)
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
