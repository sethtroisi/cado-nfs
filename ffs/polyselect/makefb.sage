#!/usr/bin/env sage
load tools.sage

def inverse_mod(a,pk):
    g,_,u=xgcd(pk,a)
    assert (not g.is_zero())
    return u

def lift_root_unramified(f,df,r,p,kmax):
    assert f(r) % p == 0
    assert df(r) % p != 0
    k = 1
    A=f.base_ring()
    while k < kmax:
        pass
        r = A((r - f(r)*inverse_mod(df(r),p^(2*k))+p^(2*k))) % p^(2*k) 
        k = k * 2;
    return r


# The roots are of the form phi(x), and phi is now to constrain roots mod
# p^m ; the content of original_f(phi(x)) is known to be equal to k0
def all_roots_affine(f,p,kmax,k0,m,phi):
    if k0 >=kmax:
        return []
    assert valuation (gcd(f.coefficients()), p) == 0
    A=f.base_ring()
    ZP=f.parent()
    x=ZP.gen()
    df = f.derivative()
    res=[]
    for r in bar(f,p).roots():
        r=A(r[0].polynomial())
        if df(r) % p != 0:
            rr = lift_root_unramified(f,df,r,p,kmax-k0)
            for l in range(1,kmax-k0):
                res.append((m+l,phi(rr) % p^(m+l),k0+l-1,k0+l))
        else:
            ff=f(r+p*x)
            v = valuation(gcd(ff.coefficients()), p)
            res.append((m+1,phi(r) % p^(m+1), k0, k0 + v))
            nphi=phi(r+p*x)
            nm=m+1
            res.extend(all_roots_affine(ZP(ff/p^v),p,kmax,k0+v,m+1,nphi))
    return res


def all_roots(f,p,powerlim):
    ZP=f.parent()
    x=ZP.gen()
    A=ZP.base_ring();
    kmax=floor(powerlim/p.degree())
    #print f,p
    aff=all_roots_affine(f,A(p),kmax,0,0,x)
    final=[]
    # affine
    for r in aff:
        final.append((r[0],r[1],1,r[3],r[2]))
        # print "%d^%d : (%d:%d), delta=%d-%d" % (p,r[0],r[1],1,r[3],r[2])
    # projective
    # That's a special precaution, as the all_roots_affine code assumes
    # we've got no content on input. Note that this value v is used later
    # on.
    fh=(f.reverse())(p*x)
    v = valuation (gcd(fh.coefficients()), p)
    if v > 0:
        final.append((1,1,0,v,0))
        fh=ZP(fh/p^v)
    proj=all_roots_affine(fh,p,kmax-1,0,0,x)
    for r in proj:
        final.append((1+r[0],1,p*r[1],v+r[3],v+r[2]))
        # print "%d^%d : (%d:%d), delta=%d-%d" % (p,1+r[0],1,p*r[1],v+r[3],v+r[2])
    # Now print and check:
    ZP2=PolynomialRing(A,['X','Y'])
    X,Y=ZP2.gens()
    cl=f.coeffs()
    hom_f=ZP2(0)
    for i in range(f.degree()+1):
        hom_f+=cl[i]*X^i*Y^(f.degree()-i)
    final2=[]
    for r in final:
        rr = tuple([p] + list(r))
        final2.append((p^r[0],r[1],r[2],r[3],r[4]))
        vf=valuation(gcd(hom_f(rr[2]+p^rr[1]*X,rr[3]+p^rr[1]*Y).coefficients()),p)
        if vf != rr[4]:
            pass
        assert valuation(gcd(hom_f(rr[2]+p^rr[1]*X,rr[3]+p^rr[1]*Y).coefficients()),p) == rr[4]
    return final2

def rewrite_roots(lr):
    if lr == []:
        return lr
    rr = []
    for r in lr:
        if r[2] == 1:
            roo = r[1]
        else:   # projective root, encoded as p^k + 1/r.
            roo = r[0] + ((r[2]*inverse_mod(r[1],r[0])) % r[0])
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

#dlim= max degree on the sieving domain
#powerlim=maximal degree of the irreducibles powers we consider

def makefb(f,dlim,powerlim,filename="",typo="cado"): 
	if filename == "":
		gd=sys.stdout
	else:
		gd=open(filename,"w")
	gd.write("# f={0}\n".format(f))
	gd.write("# dlim={0}\n".format(dlim))
	gd.write("# powerlim={0}\n".format(powerlim))
	A=f.base_ring()
	if A!=ZZ and A!=QQ:
		F=A.base_ring()
		q=len(F.list())
		dummy=2^ceil(log(q)/log(2))
		t=A.gen()
		Zt.<t0>=ZZ['t']
	def hexify(ri):
		if typo=="cado":
			return hex(Zt(ri)(dummy))
		else:
			return format(ri)
		F=A.base_ring()
		q=F.cardinality()
	for p in Primes(A,q^min(dlim,powerlim)):
		#print 'ok'
		xx = all_roots(f, p, powerlim)
		xx = rewrite_roots(xx)
		for r in xx:
			if (r[1] != 1) or (r[2] != 0):
				stri =hexify(r[0])+":"+format(r[1])+","+format(r[2])+": "+hexify(r[3][0])
			else:
				stri = hexify(r[0])+": "+hexify(r[3][0])
			for i in range(1,len(r[3])):
				stri = stri + ","+hexify(r[3][i])
			gd.write(stri+"\n")


