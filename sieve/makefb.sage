#!/usr/bin/env sage

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

def all_roots(f,p,maxbits):
    Z=Integers();
    ZP=f.parent()
    x=ZP.gen()
    kmax=ceil(maxbits*log(2)/log(p))
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
    final2=[]
    for r in final:
        rr = tuple([p] + list(r))
        final2.append((p^r[0],r[1],r[2],r[3],r[4]))
        assert valuation(hom_f(rr[2]+p^rr[1]*X,rr[3]+p^rr[1]*Y).content(),p) == rr[4]
    return final2

def rewrite_roots(lr):
    if lr == []:
        return lr
    rr = []
    for r in lr:
        if r[2] == 1:
            roo = r[1]
        else:   # projective root, encoded as p^k + 1/r.
            roo = r[0] + ((r[2]/r[1]) % r[0])
        rr.append((r[0], r[3], r[4], roo))
    rr.sort()
    # merge lines with identical first three entries
    ss = []
    old0 = 1
    old1 = 0
    old2 = 0
    lr = []
    for r in rr:
        if r[0] == old0 and r[1] == old1 and r[2] == old2:
            lr.append(r[3])
        else:
            if lr != []:
                ss.append((old0, old1, old2, lr))
            old0 = r[0]
            old1 = r[1]
            old2 = r[2]
            lr = [r[3]]
    ss.append((old0, old1, old2, lr))
    return ss

x=PolynomialRing(Integers(),['x']).gen()


# c158
f=69042960*x^5 -150777388929552*x^4 +129632089360584232586*x^3 +263540073072166885859228735*x^2 -48947696216829079528262929524278*x -32290899357812757163668740209282119536

#f=1008593880*x^5 - 47389790327*x^4 - 84256212127259029352*x^3 + 3474222647711706240332297*x^2 + 76764659243128790828718944401*x + 62435925692971697863740890240
#alim=40000000
alim=100
powerlim=12    # in bits




print "# f={0}".format(f)
print "# alim={0}".format(alim)
print "# powerlim={0}".format(powerlim)
p = 2
while p <= alim:
    xx = all_roots(f, p, powerlim)
    xx = rewrite_roots(xx)
    for r in xx:
        if (r[1] != 1) or (r[2] != 0):
            stri = "{0}:{1},{2}: ".format(r[0],r[1],r[2]) + r[3][0].str()
        else:
            stri = "{0}: ".format(r[0]) + r[3][0].str()
        for i in range(1,len(r[3])):
            stri = stri + ",{0}".format(r[3][i])
        print stri
    p = p.next_prime()
