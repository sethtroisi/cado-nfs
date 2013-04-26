"""
The main function is approx_alpha(f,B). It is a good approximation of alpha
using the factor base.
"""

load makefb.sage


"""
format of makefb 
ex [(t^5, 7,4, t^4+t^2+t)]
stays for val(f(t^4+t^2+t),t)=7 and val(f(t^2+t),t)=4
"""
def approx_alpha_l(f,l,prec=20):
    A=f.base_ring();t=A.gen();q=A.base_ring().cardinality()
    l=A(l)
    dl=l.degree()
    nl=q^dl
    result=0
    lr=all_roots(S(f),A(l),prec)
    ram,simple=makefb_format(S(f),A(l),lr)
    result=0
    for r_ in ram:
        k=valuation(r_[0],l)
        v1=r_[1]
        result+=v1*1/nl^k
    for r in simple:
        result+=1/(nl-1) 
    return dl*(1/(nl-1)-nl/(nl+1)*result) 

def approx_alpha(f,B,prec=10,B0=1):
    """
    computes the alpha(f), use B=p^k for polys l with deg <=k. 
    """
    A=f.base_ring()
    b=floor(log(B,A.base_ring().cardinality()))
    b0=ceil(log(B0,A.base_ring().cardinality()))
    result=0
    for d in range(b0+1,b+1):
        for l in A.polynomials(of_degree=d):
            if l.is_monic() and l.is_irreducible():
                result+=approx_alpha_l(f,l)
    return result.n()

"""
 #  EXAMPLE
 load init_rings.sage
 f=S(x^6+(t^2+t)*x^3+(t^5 + t + 1)*x+t^6 + t^4 + t^3 + t)
 sage: time approx_alpha(f,2000)                     
 sage: time alpha(f,2000)       
"""
