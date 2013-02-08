load init_rings.sage
load tools.sage
load real_roots_stats.sage


# The main function is alpha(f,B) which computes the root property of f with
# respect to irreducible polynomials of position less than B,e.g. B=2^5 for
# polynomials of degree below 5 if characteristic=2
#The "estimate_<name_of_function>" makes probabilistic estimates of
#<name_of_function> 

# Polynomials must be in (GF(p)['t'])['x']
# A.<t>=GF(2)['t']
# S.<x>=A['x']
# alpha(S((x-t)*(x-t^2)),100)



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
    #Q = bar(fv,p).derivative()
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
    returns the average l-valuation of f(a,b) for a,b coprime. Projective
    roots are counted as well.
    """
    S= f.parent()
    x = S.gen()
    affine_average=average_valuation_affine(f, l)
    proj_average=QQ(average_valuation_affine((f.reverse())(l*x), l))
    return (affine_average*Norm(l)/(Norm(l)+1)+proj_average/(Norm(l)+1))
    # hom = Prob(b=0)*proj + Prob(b!=0)*affine

def alpha_p(f,l):
    """
    Computes the contribution at l of the alpha value of f
    """
    return float((1/(Norm(l)-1)-average_valuation_homogeneous_coprime(f,l))*l.degree())

def estimate_alpha_p(f, p, nt):
    """
    Should compute the same thing as alpha_p, but experimentally
    """
    s=0
    n=0
    x=f.parent().gen()
    A=f.base_ring()
    F=homogenize(f)
    l=p.degree()
    for i in range(nt):
        a=A.random_element(ceil(log(nt,2)))
        b=A.random_element(ceil(log(nt,2)))
        c=A.random_element(f(0).degree()+f.degree()*ceil(log(nt,2)))
        if gcd(a,b) == 1:
            ds=valuation(c,p)-valuation(A(F(a,b)), p)
            s+=ZZ(ds)
            n+=1
    #sys.stdout.write("%f\r" % float(s*l/n))
    return float(s*l/n) 


def alpha(f,B,B0=1):
    """
    Computes the alpha value of f, for irreducible polynomials l such that
    B0+1<=Norm(l)<= B
    """
    A=f.base_ring()
    q=A.base_ring().cardinality()
    b=ceil(log(B,q))
    b0=ceil(log(B0,q))
    result=0
    for d in range(b0+1,b+1):
        for l in A.polynomials(of_degree=d):
            if l.is_monic() and l.is_irreducible():
                result+=alpha_p(f,l)
    return result

def completed_alpha(f,B,skewness):
    return alpha(f,B)+alpha_infty(f,skewness)


def estimate_alpha(f,B,nt):
    """
    Same as alpha(f,B), but experimentally.  Only for debugging.
    """
    A=f.base_ring()
    return sum([estimate_alpha_p(f, p, nt) for p in Primes(A,B+1)])



def estimate_average_valuation_homogeneous_coprime(f, p, nt,x0=0, x1=0,y0=0,y1=0):
    """
    this computes the same thing as average_valuation_homogeneous_coprime(f,p),
    but does so experimentally.  Only for debugging.
    """
    if x1==0:
        x1=nt
        if y1==0:
            y1=nt
    s=0
    n=0
    x=f.parent().gen()
    A=f.base_ring()
    F=homogenize(f)
    for a in element_list(A,x0+x1,x0):
        for b in element_list(A,y0+y1,y0):
            if gcd(a,b) == 1:
                #print a,b
                s+=valuation(A(F(a,b)), p)
                n+=1
	return float(s/n) if n > 0 else Infinity






def estimate_random_valuation_p( p, nt):
	"""
	Should compute the same thing as alpha_p, but experimentally
	"""
	s=0
	n=0
	A=p.parent()
	list=element_list(A,nt^2,0)
	for i in range(nt):
		ia=randrange(0,nt^2)
		ib=randrange(0,nt^2)
		ic=randrange(0,nt^2)
		a,b,c=list[ia],list[ib],list[ic]
		if gcd(a,b) == 1:
			s+=valuation(c,p)
			n+=1
	sys.stdout.write("%f\r" % float(s/n))
	return float(s/n) if n > 0 else Infinity


