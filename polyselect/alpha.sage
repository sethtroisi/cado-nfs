
# Code for computing alpha, and some variations around it.
# The main function exported here is alpha(f,p)

# We also provide splits alpha_affine and alpha_projective.
# alpha_projective, roughly said, depends only on the top coefficients
# (except for very special cases where touching a low degree coefficients
# might change the value).

def number_of_roots(f,p):
    """
    Counts the roots of f mod p, with multiplicities. Projective roots
    are also counted (with multiplicities).
    """
    fp= GF(p)['x'](f)
    if fp != 0:
        s=sum([r[1] for r in fp.roots()])
        s+=f.degree()-fp.degree()
    else:
        s=f.degree()
    return s

def average_valuation_affine(f,p):
    """
    returns the average p-valuation of the polynomial f. Works recursively.
    """
    v = valuation (f.content(), p)
    ZP= f.parent()
    x = ZP.gen()
    fv = ZP(f/p^v)
    Q = GF(p)['x'](fv.derivative())
    for r in fv.roots(GF(p)):
        if Q(r[0]) != 0:
            v += 1/(p-1)
        else:
            f2 = fv(Integers()(r[0])+p*x)
            v += average_valuation_affine(f2, p)/p
    return v

def average_valuation_homogeneous_coprime(f,p):
    """
    returns the average p-valuation of f(a,b) for a,b coprime. Projective
    roots are counted as well.
    """
    # modulo p^n, roots touch one class having p^n*(1-1/p) representatives,
    # amongst the p^2n*(1-1/p^2) coprime representative pairs for
    # P^1(Z/p^n). So the contribution is, for p^n, p/(p+1) * p^-n
    r = average_valuation_affine(f, p) * p
    ZP= f.parent()
    x = ZP.gen()
    # Infinity expands as 0 when swapped. So we already know the
    # first-order term. For this reason, we we may drop one p, since it's
    # already counting one order ahead.
    r += average_valuation_affine((f.reverse())(p*x), p)
    r /= p+1
    return r

def alpha_p(f,disc,p):
    """
    Computes the contribution at p of the alpha value of f
    """
    if disc % p == 0:
        s = float((1/(p-1)-average_valuation_homogeneous_coprime(f,p))*log(p))
    else:
        # Then we know that the average valuation is
        # number_of_roots/(p-1), multiplied by p/(p+1) so as to take into
        # account the non-coprime pairs.
        s = float((1-number_of_roots(f,p)*p/(p+1))*log(p)/(p-1))
    return s

def alpha_p_affine(f,disc,p):
    """
    Computes the contribution at p of the alpha value of f (affine part
    only)
    """
    if disc % p == 0:
        e = average_valuation_affine(f, p) * p / (p+1)
        s = float((1/(p-1)-e)*log(p))
    else:
        s = float((1-len(f.roots(GF(p)))*p/(p+1))*log(p)/(p-1))
    return s

def alpha_p_projective(f,disc,p):
    """
    Computes the contribution at p of the alpha value of f (projective part
    only). Always <= 0.
    """
    x=f.parent().gen()
    q=f.reverse()(p*x)
    if disc % p == 0:
        e = average_valuation_affine(q,p)/(p+1)
        s = float((-e)*log(p))
    else:
        fp= GF(p)['x'](f)
        assert(fp != 0)
        s = float((-(f.degree()-fp.degree())*p/(p+1))*log(p)/(p-1))
    return s

def alpha(f,B):
    """
    Computes the alpha value of f, up to prime bound B
    """
    disc = f.discriminant()
    return sum([alpha_p(f, disc, p) for p in prime_range(2,B)])

def alpha_affine(f,B):
    """
    Computes the affine part of the alpha value of f, up to prime bound B
    """
    disc = f.discriminant()
    return sum([alpha_p_affine(f, disc, p) for p in prime_range(2,B)])

def alpha_projective(f,B):
    """
    Computes the affine part of the alpha value of f, up to prime bound B.
    Always <= 0.
    """
    disc = f.discriminant()
    return sum([alpha_p_projective(f, disc, p) for p in prime_range(2,B)])

# -*-*- debug -*-*-
def estimate_average_valuation_homogeneous_coprime(f, p, x0, x1, y0, y1):
    """
    this computes the same thing as average_valuation_homogeneous_coprime(f,p),
    but does so experimentally.  Only for debugging.
    """
    s=0
    n=0
    x=f.parent().gen()
    for a in range(x0,x0+x1):
        for b in range(y0,y0+y1):
            if gcd(a,b) == 1:
                s+=valuation(f.resultant(a*x-b), p)
                n+=1
    return float(s/n) if n > 0 else Infinity

def estimate_alpha_p(f, p, nt):
    """
    Should compute the same thing as alpha_p, but experimentally
    """
    s=0
    n=0
    x=f.parent().gen()
    l=log(p)
    for i in range(nt):
        a=randrange(0,nt^2)
        b=randrange(0,nt^2)
        c=randrange(0,nt^2)
        if gcd(a,b) == 1:
             s+=valuation(c,p)-valuation(f.resultant(a*x-b), p)
             n+=1
             sys.stdout.write("%f\r" % float(s*l/n))
    return float(s*l/n) if n > 0 else Infinity

# auxiliary
def alpha_simplistic(f,B):
    """
    This ignores the multiple roots altogether.
    """
    return sum([float((1-number_of_roots(f,p)*p/(p+1))*log(p)/(p-1)) for p in prime_range(2,B)])
