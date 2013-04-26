load l_roots.sage

def deg(r):
    """
    extends degree to rational fractions
    """
    if r == 0:
        return -Infinity
    return r.numerator().degree()-r.denominator().degree()

def gap(f,r,m):
    """
    r=a_(-deg(r)) t^degr(r)+...+a_m 1/t^m with a_m possibly 0;
    computes the gap of (r,m)
    """
    F=f.base_ring().base_ring();A.<t>=F['t']
    K=FractionField(A);S.<x>=K['x']
    f=S(f)
    d=f.degree()
    cf=f.coeffs()
    r=K(r)
    result=max(deg(cf[i]*r^i) for i in [0..d])-deg(R(f)(t,r))
    for a_m1 in F:
        r1=r+a_m1/t^(m+1)
        r1=K(r1)
        gap_r1=max(deg(cf[i]*r1^i) for i in [0..d])-deg(R(f)(t,r1))
        if gap_r1 < result:
            result =gap_r1
    return result


def Laurent_roots(f,prec=5):
    F=f.base_ring().base_ring();q=F.cardinality()
    A.<t>=F['t'];K=FractionField(A);S.<x>=K['x']
    R_.<t,x>=F['t,x']; R=FractionField(R_)
    f=S(f)
    d=f.degree()
    cf=f.coeffs()
    s=ceil(max((deg(cf[i])-deg(cf[d]))/(d-i) for i in [0..d-1]))
    f_bar=R(f)(1/t,x)
    z=R(f_bar)(t,x*t^(-s)).denominator().degree()
    f_tilde=t^z*R(f_bar)(t,x*t^(-s))
    f_tilde=S(R_(f_tilde))
    v=valuation(f_tilde.coeffs()[0],K(t))
    rr=affine_roots(S(f_tilde),A(t),prec) 
    # for a proven code, replace prec by prec+v 
    result=[]
    for r_ in rr:
        r_tilde=A(r_[0])
        k_tilde=r_[1]
        r=r_tilde(1/t)*t^s
        m=k_tilde-s-1
        r=K(R(r).numerator())/K(R(r).denominator())
        gap_=gap(S(f),r,m)
        if gap_ in [1..prec]:
            result.append([r,m,gap_])
    return result

def alpha_infty(f,skew,prec=5):
    rr=Laurent_roots(f,prec)
    q=f.base_ring().base_ring().cardinality()
    result=0
    for r_ in rr:
        r=r_[0]
        m=r_[1]
        k=r_[2]
        Nr=deg(r)+m
        result-=1/(q^2-1)*q^-(Nr+abs(skew-deg(r)))
    return result

