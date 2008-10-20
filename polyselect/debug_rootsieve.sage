
attach rotation_bound.sage
attach kleinjung.sage
attach higher_order_roots.sage

Z=Integers()
ZP.<x>=PolynomialRing(Z)
n=110547499294735934596573920716806495615146307355718579549435266277149267619863094991687986227655049181879310236052609641971560327957264468408615218306269241
ad=1008593880; p=1628876881135933; m=161423410553519491417914681351; E=44.88
f,g=lemme_21(n,5,ad,p,m)
g1=g(x+22977)
f1=f(x+22977)+(3*x+18223)*g1

# Now alpha(f1,2000) is 0.18.
# Within rotation by (jx+k), |j|<=4, |k|<=2^16, the best alpha is reached
# for:
# sage: alpha(f1+(x-26071)*g1,2000)
# -4.6975881658522551
#
# the C version catches this in 82 seconds on nougatine.

f=f1
g=g1


ZP2=PolynomialRing(Z,['X','Y'])
X,Y=ZP2.gens()

def hom(f):
    return ZP2(f(X/Y)*Y^(f.degree()))

ZP3.<t,u,v>=PolynomialRing(Integers(),3)
h=f(t)+(u*t+v)*g(t)


def get_reference(sbound,p):
    complete=matrix(RR,sbound)
    a = flog(p)/(p-1)
    x=ZP.gen()
    t0=cputime()
    for k in range(sbound):
        for l in range(sbound):
            a=flog(p)/(p-1)
            r=f+(k*x+l)*g
            s=alpha_p_simplistic(r,p)-a
            if (valuation(r[r.degree()],p)>=1):
                s -= -p/(p+1)*a
            c=alpha_p_affine_nodisc(r,p)-a
            # simple[k,l]=s
            complete[k,l]=c
    return complete



def testp(rdict,p,reference):
    rotation_clear(rdict)
    rotation_handle_p(rdict,p)
    allcoeffs=reduce((lambda x,y: x+y),[list(x) for x in
         list(matrix(sbound,sbound,rdict['sarr'])-reference)],[])
    mi,ma=min(allcoeffs), max(allcoeffs);
    res=max(ma,-mi)
    print "Maximal inaccuracy: %e (large means bug!)" % res
    return res < 0.1

def manyp(rdict,plim):
    rotation_clear(rdict)
    t0=cputime()
    for p in prime_range(plim):
        t1=cputime()
        hits=rotation_handle_p(rdict,p)
        print "%r: %.2f seconds, %r hits" % (p,cputime()-t1,hits)
    print "total: %.2f seconds" % (cputime()-t0)



####################################
# To try:

# Fix a size for the sieving area. It's taken square here, but the
# program accepts also a rectangle. We look into polynomials h(x,u,v) for
# u and v integers within the interval [0,sbound-1]
sbound=40
p=43

print "First computing reference scores with the naive method"
refp = get_reference(sbound, p)
print "Took %.2f seconds" % (cputime()-t0)
print "Now with root sieve"
t0=cputime()
rdict=rotation_init(f,g,0,sbound,0,sbound)
testp(rdict,p,complete)
print "Took %.2f seconds" % (cputime()-t0)


K=GF(p)
KP=K['t']
KP3=K['t','u','v']
KF3=KP3.fraction_field()
ZF3=ZP3.fraction_field()

# More extensive testing:
# ref2=get_reference(sbound,2)
# ref3=get_reference(sbound,3)
# ref5=get_reference(sbound,5)
# ref7=get_reference(sbound,7)

# rdict=rotation_init(f,g,0,sbound,0,sbound)
# t0=cputime();testp(rdict,2,ref2);testp(rdict,3,ref3);testp(rdict,5,ref5);cputime()-t0


# zview(mround(matrix(sbound,sbound,sarr)-complete),1,0,p)
# zview(mround(matrix(sbound,sbound,sarr)-complete),2,1,p)


