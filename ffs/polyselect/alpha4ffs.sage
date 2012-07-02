load tools.sage

# The main function is alpha(f,B) which computes the root property of f with
# respect to irreducible polynomials of position less than B,e.g. B=2^5 for
# polynomials of degree below 5 if characteristic=2
#The "estimate_<name_of_function>" makes probabilistic estimates of
#<name_of_function> 

# Polynomials must be in (GF(p)['t'])['x']
# A.<t>=GF(2)['t']
# S.<x>=A['x']
# alpha(S((x-t)*(x-t^2)),100)



def number_of_roots(f,p):
    """
    Counts the roots of f mod p, without multiplicities. Projective roots
    are also counted (without multiplicities).
    """
    fp= bar(f,p)
    if fp != 0:
        s=len([r[1] for r in fp.roots()])
        if (f.degree()>fp.degree()):
            s+=1
    else:
        # f reducing to zero mod p is a degenerate case. Not clear what
        # we should return...
        print "Warning, counting roots of zero polynomial\n"
        s=f.degree()
    return s


def average_valuation_affine(f,p,depth=0):
    """
    returns the average p-valuation of the polynomial f. Works recursively.
    """
    MAX_PRECISION=10
    fcont=gcd(f.coefficients())
    v = valuation (fcont, p)
    ZP= f.parent()
    A=ZP.base_ring()
    x = ZP.gen()
    fv = ZP(f/p^v)
    Q = bar(fv,p).derivative()
    for r in bar(fv,p).roots():
        assert fv(A(r[0].polynomial())) % p ==0
        if Q(r[0]) != 0:
            v += QQ(1/(Norm(p)-1))
        else:
            f2 = fv(A(r[0].polynomial())+p*x)
            if depth<MAX_PRECISION:
                v += QQ(average_valuation_affine(f2, p,depth+1)/Norm(p))
    return v


def average_valuation_homogeneous_coprime(f,p):
    """
    returns the average p-valuation of f(a,b) for a,b coprime. Projective
    roots are counted as well.
    """
    ZP= f.parent()
    x = ZP.gen()
    affine_average=average_valuation_affine(f, p)
    proj_average=QQ(average_valuation_affine((f.reverse())(p*x), p))
    if proj_average ==0:
        return affine_average
    else:
        return (affine_average*(Norm(p)-1)/Norm(p)+proj_average/(Norm(p)))


def alpha_p(f,p):
    """
    Computes the contribution at p of the alpha value of f
    """
    return float((1/(Norm(p)-1)-average_valuation_homogeneous_coprime(f,p))*p.degree())

def Estimate_alpha_p(f,p,nt):
    """
    Same as alpha_p(f,p), but experimentally.  Only for debugging.
    """
    return float((1/(Norm(p)-1)-estimate_average_valuation_homogeneous_coprime(f,p,nt))*p.degree())


def alpha(f,B):
    """
    Computes the alpha value of f, up to prime bound B
    """
    A=f.base_ring()
    return sum([alpha_p(f, p) for p in Primes(A,B+1)])

def estimate_alpha(f,B,nt):
    """
    Same as alpha(f,B), but experimentally.  Only for debugging.
    """
    A=f.base_ring()
    return sum([Estimate_alpha_p(f, p, nt) for p in Primes(A,B+1)])



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
	list=element_list(A,nt^2,0)
	for i in range(nt):
		ia=randrange(0,nt)
		ib=randrange(0,nt)
		ic=randrange(0,nt^2)
		a=list[ia]
                b=list[ib]
                c=list[ic]
		if gcd(a,b) == 1:
			ds=valuation(c,p)-valuation(A(F(a,b)), p)
			s+=ZZ(ds)
			n+=1
	#sys.stdout.write("%f\r" % float(s*l/n))
	return float(s*l/n) 



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


