load alpha.sage
load alpha_infty.sage 

def sigma(f,s,e=0):
    """
    s= skewness ,i.e. deg(a)-deg(b)
    e= max degree of pairs a,b
    """
    d=f.degree()
    fc=f.coeffs()
    return d*(e-s/2)+max([fc[i].degree()+i*s for i in range(d+1)])

def completed_alpha(f,s,B=1000):
    """
    s= skewness
    """
    return alpha(f,B)+alpha_infty(f,s)

def epsilon(f,s,e=0,B=1000):
    """
    s= skewness ,i.e. deg(a)-deg(b)
    e= max degree of pairs a,b
    """
    return alpha(f,B)+alpha_infty(f,s)+sigma(f,s,e)

"""
load init_rings.sage
f=x^5+t^40*x+t
sigma(S(f),10)
sigma(S(f),0)
epsilon(S(f),10,7)
"""
