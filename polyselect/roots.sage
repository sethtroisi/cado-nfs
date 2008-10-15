# An example:
# all_roots(((x-3)*(x-19)^2+2^16)*(x-1234567)+2^40,2)
# studies the exponent of 2 in F(a,b), where F(x,y) = y^deg(f)*f(x/y)
# which returns (among others):
# 2^52 : (2576298330889987:1), delta=62-61
# This result means that for a,b fixed mod 2^52, we have to 
# add 62-61 = 1 * log(2) to the concerned cells, which are those
# in the lattice a=2576298330889987*i+j, b=2^52*i.

# For roots at infinity, the 2nd part of the output is (1:x), which means
# the lattice is a=i+x*j, b=2^52*j.

# The delta value is not always 1, and it is given by a difference wrt
# smaller exponents.

# Remark: to avoid rounding errors, instead of adding delta*log(2) in the
# sieving, it is preferable to add round(hi*log(2)) - round(lo*log(2)),
# where delta = hi - lo.

def lift_root_unramified(f,df,r,p,kmax):
    assert f(r) % p == 0
    assert df(r) % p != 0
    k = 1
    while k < kmax:
        r = Integers()(Integers(p^(2*k))(r - f(r)/df(r)))
        k = k * 2;
    return r

# The roots are of the form phi(x), and phi is now to constrain roots mod
# p^m ; the content of original_f(phi(x)) is known to be equal to k0
def all_roots_affine(f,p,kmax,k0,m,phi):
    res=[]
    assert valuation (f.content(), p) == 0
    if (k0 >= kmax):
        return res
    Z=Integers();
    ZP=f.parent()
    x=ZP.gen()
    K = GF(p)
    KP=K['x']
    Q = f.derivative()
    for r in f.roots(K):
        r=Z(r[0])
        if Q(r) % p != 0:
            rr = lift_root_unramified(f,Q,r,p,kmax-k0)
            for l in range(1,kmax-k0):
                res.append((m+l,phi(rr) % p^(m+l),k0+l-1,k0+l))
        else:
            ff=f(r+p*x)
            v = valuation(ff.content(), p)
            res.append((m+1,phi(r) % p^(m+1), k0, k0 + v))
            nphi=phi(r+p*x)
            nm=m+1
            res.extend(all_roots_affine(ZP(ff/p^v),p,kmax,k0+v,m+1,nphi))
    return res

def all_roots(f,p):
    Z=Integers();
    ZP=f.parent()
    x=ZP.gen()
    kmax=ceil(64*log(2)/log(p))
    aff=all_roots_affine(f,p,kmax,0,0,x)
    final=[]
    # affine
    for r in aff:
        final.append((r[0],r[1],1,r[3],r[2]))
        # print "%d^%d : (%d:%d), delta=%d-%d" % (p,r[0],r[1],1,r[3],r[2])
    # projective
    # That's a special precaution, as the all_roots_affine code assumes
    # we've got no content on input. Note that this value v is used later
    # on.
    fh = ZP((p*x)^f.degree()*f(1/(p*x)))
    v = valuation (fh.content(), p)
    if v > 0:
        final.append((1,1,0,v,0))
        fh=ZP(fh/p^v)
    proj=all_roots_affine(fh,p,kmax-1,0,0,x)
    for r in proj:
        final.append((1+r[0],1,p*r[1],v+r[3],v+r[2]))
        # print "%d^%d : (%d:%d), delta=%d-%d" % (p,1+r[0],1,p*r[1],v+r[3],v+r[2])
    # Now print and check:
    ZP2=PolynomialRing(Z,['X','Y'])
    X,Y=ZP2.gens()
    hom_f=ZP2(f(X/Y)*Y^(f.degree()))
    for r in final:
        rr = tuple([p]+list(r))
        print "%d^%d : (%r:%r), delta=%d-%d" % rr
        assert valuation(hom_f(rr[2]+p^rr[1]*X,rr[3]+p^rr[1]*Y).content(),p) == rr[4]
