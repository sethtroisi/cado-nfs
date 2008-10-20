
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


sbound=32
simple=matrix(RR,sbound)
complete=matrix(RR,sbound)
dz=matrix(ZZ,sbound)
p=2
a = flog(p)/(p-1)

x=ZP.gen()
print "First computing reference scores with the naive method"
t0=cputime()
for k in range(sbound):
    for l in range(sbound):
        a=flog(p)/(p-1)
        r=f+(k*x+l)*g
        s=alpha_p_simplistic(r,p)-a
        if (valuation(r[r.degree()],p)>=1):
            s -= -p/(p+1)*a
        c=alpha_p_affine_nodisc(r,p)-a
        simple[k,l]=s
        complete[k,l]=c

print "Took %.2f seconds" % (cputime()-t0)


print "Now with root sieve"
t0=cputime()
rdict=rotation_init(f,g,0,sbound,0,sbound)
rotation_clear(rdict)
rotation_handle_p(rdict,p)
allcoeffs=reduce((lambda x,y: x+y),[list(x) for x in
     list(matrix(sbound,sbound,rdict['sarr'])-complete)],[])
mi,ma=min(allcoeffs), max(allcoeffs);
print "Took %.2f seconds" % (cputime()-t0)
print "Maximal inaccuracy: %e (large means bug!)" % max(ma,-mi)

# zview(mround(matrix(sbound,sbound,sarr)-complete),1,0,p)
# zview(mround(matrix(sbound,sbound,sarr)-complete),2,1,p)


