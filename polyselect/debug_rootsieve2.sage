attach rotation_bound.sage
attach kleinjung.sage
attach higher_order_roots.sage

Z=Integers()
ZP.<x>=PolynomialRing(Z)

f = 626160*x^5 + 199078173920*x^4 + 46653788280990838*x^3 - 89368397408157908492*x^2 - 4261015116515121510492288*x + 1385334566237042810410340937
g = 2583241664039*x-5806276217823462199974022
sboundu=2
sboundv=2000

f2=f
g2=g
# Output from polyselect: Rotate by x+1036: alpha improved from -0.83 to -4.06

rdict=rotation_init(f2,g2,0,sboundu,0,sboundv)
pmax=2000
rotation_clear(rdict)
tt0=cputime()
for p in prime_range(1,pmax):
    tt=cputime()
    hits,best,u,v=rotation_handle_p(rdict,p)
    check=alpha(f2+(u*x+v)*g2,p)
    z=log(10^-20+abs(check-best))
    gooddigits=floor(-z/log(10))
    print "Done %d [%.2f, tot %.2f]." \
            " Best alpha=%.4f, for %d,%d (%d dd ok)" % \
            (p, cputime()-tt, cputime()-tt0, best,u,v,gooddigits)
