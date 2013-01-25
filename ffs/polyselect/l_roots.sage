# The main function is affine_roots. It is used in ... .



# input:  irreducible i in GF(q)[t] and a poly f in GF(q)[t][x]
# output: r in GF(q)[t] such that f(r)=0 mod l and deg(r)<=l. 
def roots_mod_l(f,l):
    """
    ! roots() not implemented for extension fields
    """
    A=f.base_ring()
    if not A(l).is_irreducible():
        print "l must be irreducible"
        return []
    Fq=A.base_ring(); R.<t,x>=Fq['t,x']; q=Fq.cardinality() 
    dl=l.degree()
    Kl.<w>=GF(q^l.degree())
    wl=A(l).change_ring(Kl).roots(multiplicities=false)[0]
    U.<xx>=Kl['x']
    fl=R(f)(wl,xx)
    if fl == 0:
        result=Kl.list()
    elif U(fl).degree() == 0:
        return []
    else:
        result=U(fl).roots(multiplicities=false)
    if A(l).degree() > 1:
        M=(matrix([vector(Kl(wl^i)) for i in range(dl)])^(-1)).transpose()
        return [A(list(M*vector(Kl(r)))) for r in result]
    else:
        return result

def inverse_mod(a,pk):
    g,_,u=xgcd(pk,a)
    assert (not g.is_zero())
    return u

# input : a simple root r of precision prec 
# output: a root r_new of precision 2*prec
def hensel_lift(f,l,r,prec,df=None):
    A=r.parent();t=A.gen();F=A.base_ring(); S.<x>=A['x'];
    if df==None:
        df=S(f).derivative()
    assert S(f)(r) % l^prec == 0
    assert df(r) % l != 0
    r_new=r-f(r)*A(df(r)).inverse_mod(A(l)^(2*prec))
    r_new= r_new % l^(2*prec)
    assert (r-r_new) % l^prec == 0
    assert f(r_new) % l^(2*prec) ==0
    return r_new

# input : a root r in GF(q)[t] s.t. f(r) = 0 mod l^prec
# output: a list "final" of all the roots r_new of precision prec+1
# which extend r.
def naive_lift(f,l,r,prec):
    A=r.parent();t=A.gen();F=A.base_ring();R.<t,x>=F['t,x'];S.<x>=A['x'];
    f=S(f)
    assert f(r) % l^prec == 0
    v=prec
    S=f.parent()
    f_new=S(f(r+ l^v*x)/l^v)
    final=[]
    for e in roots_mod_l(f_new,A(l)):
        r_new=r+l^v * e
        assert f(r_new) % l^(v+1) == 0
        final.append(r_new)
    return final

# input : f in GF(q)[t][x], irreducible l in GF(q)[t], prec positive integers 
# output: the list of [r,k] such that f(r)=0 mod l^k, for k<=prec. 
def affine_roots(f,l,prec):
    A=l.parent();t=A.gen();F=A.base_ring();R.<t,x>=F['t,x'];S.<x>=A['x'];
    f=S(f)
    df=f.derivative()
    ll=roots_mod_l(S(f),A(l))
    if prec == 1:
        return [[e,1] for e in ll]
    final=[]
    for r in ll:
        if A(df(r)) % A(l) != 0:
            pr=1
            while pr<prec:
                r=hensel_lift(S(f),A(l),A(r),pr)
                pr*=2
            assert f(r) % l^prec == 0
            final+=[[r % l^i,i] for i in [1..prec]]
        else:
            pr=1
            lr=[r]
            final.append([r,1])
            while pr<prec:
                lr_new=[]
                for rr in lr:
                    lr_new+=naive_lift(S(f),A(l),A(rr),pr)
                pr+=1
                final+=[[e,pr] for e in lr]
                lr=lr_new
    return final

# output: affine roots i.e. F(r,1) % l^k = 0 and 
#         projective roots i.e. F(1,r) % l^k = 0 with r % l = 0 
def all_roots(f,l,prec):
    A=l.parent();t=A.gen();F=A.base_ring();R.<t,x>=F['t,x'];S.<x>=A['x'];
    f=S(f) 
    aff=affine_roots(f,l,prec)
    cf=f.coeffs()
    cf.reverse()
    f_proj=S(cf) # f_proj(x)=F(t,1,x), f=F(t,x,1)
    if A(S(f_proj)(0)) % A(l) == 0:
        proj=[(A(0),1)]+[(r_[0]*l,r_[1]+1) for r_ in affine_roots(S(S(f_proj)(x*l)/l),A(l),prec-1)]
    else:
        proj=[]
    return aff,proj

""" EXAMPLE
load init_rings.sage
affine_roots(S(x^2-t^3+t),A(t^4+t+1),10)
affine_roots(S(x-t^3+t),A(t^2+t+1),10)
all_roots(S((t^3+t+1)*x^2+t^2*x+(t+1)),A(t^3+t+1),3)
"""
