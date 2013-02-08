load ../makefb/tools.sage

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


def average_valuation_affine(f,p):
    """
    returns the average p-valuation of the polynomial f. Works recursively.
    """
    fcont=gcd(f.coefficients())
    v = valuation (fcont, p)
    ZP= f.parent()
    A=ZP.base_ring()
    x = ZP.gen()
    fv = ZP(f/p^v)
    Q = bar(fv,p).derivative()
    for r in bar(fv,p).roots():
        if Q(r[0]) != 0:
            v += QQ(Norm(p)/(Norm(p)-1))
        else:
            f2 = fv(A(r[0].polynomial())+p*x)
            v += QQ(average_valuation_affine(f2, p)/Norm(p))
    return v


def average_valuation_homogeneous_coprime(f,p):
    """
    returns the average p-valuation of f(a,b) for a,b coprime. Projective
    roots are counted as well.
    """
    r = average_valuation_affine(f, p) * Norm(p)
    ZP= f.parent()
    x = ZP.gen()
    # Infinity expands as 0 when swapped. So we already know the
    # first-order term. For this reason, we we may drop one p, since it's
    # already counting one order ahead.
    r += QQ(average_valuation_affine((f.reverse())(p*x), p))
    r /= QQ(Norm(p)+1)
    return r


def alpha_p(f,disc,p):
    """
    Computes the contribution at p of the alpha value of f
    """
    return float((1/(Norm(p)-1)-average_valuation_homogeneous_coprime(f,p))*p.degree())


def alpha(f,B):
    """
    Computes the alpha value of f, up to prime bound B
    """
    A=f.base_ring()
    disc = f.discriminant()
    return sum([alpha_p(f, disc, p) for p in Primes(A,B+1)])



def estimate_average_valuation_homogeneous_coprime(f, p, x0, x1, y0, y1):
	"""
	this computes the same thing as average_valuation_homogeneous_coprime(f,p),
	but does so experimentally.  Only for debugging.
	"""
	s=0
	n=0
	x=f.parent().gen()
	A=f.base_ring()
	F=homogenize(f)
	for a in element_list(A,x0+x1,x0):
		for b in element_list(A,y0+y1,y0):
			if gcd(a,b) == 1:
				s+=valuation(F(a,b), p)
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
		ia=randrange(0,nt^2)
		ib=randrange(0,nt^2)
		ic=randrange(0,nt^2)
		a,b,c=list[ia],list[ib],list[ic]
		if gcd(a,b) == 1:
			ds=valuation(c,p)-valuation(F(a,b), p)
			s+=ZZ(ds)
			n+=1
	sys.stdout.write("%f\r" % float(s*l/n))
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
