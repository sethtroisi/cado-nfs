attach "alpha.sage"
attach "rotation_bound.sage"

def lemme_21(n,d,ad,p,m):
    """Implements Kleinjung's lemma 2.1"""
    coeffs=[ad]
    deltas=[]
    rs=[]
    r=n
    a=ad
    Z=Integers()
    mm=[1]
    immp=[Integers(p)(1)]
    imp=1/Integers(p)(m)
    for i in [1..d]:
        mm.append(mm[-1]*m)
        immp.append(immp[-1]*imp)
    for i in reversed([0..d-1]):
        r=Z((r-a*mm[i+1])/p)
        ai_mod_p=Z(r*immp[i])
        k=(r/mm[i]-ai_mod_p)/p
        kr=round(k)
        delta=p*(kr-k)
        a=p*kr+ai_mod_p
        rs.append(r)
        coeffs.append(a)
        deltas.append(delta)
    coeffs.reverse()
    rs.reverse()
    deltas.reverse()
    f=Integers()['x'](coeffs)
    return (f,deltas,rs)

# For example:

# n=1230186684530117755130494958384962720772853569595334792197322452151726400507263657518745202199786469389956474942774063845925192557326303453731548268507917026122142913461670429214311602221240479274737794080665351419597459856902143413
# f=ZZ['x']([-277565266791543881995216199713801103343120,-18185779352088594356726018862434803054,6525437261935989397109667371894785,-46477854471727854271772677450,-5006815697800138351796828,1276509360768321888,265482057982680])
# g=ZZ['x']([-1291187456580021223163547791574810881,34661003550492501851445829])
# f0,ds,rs=lemme_21(n,6,f[6],g[1],-g[0])
# (f-f0)/g

def rotation (n,d,p,m):
   """Computes optimal rotation for degree-d polynomial for integer n with
      linear polynomial p*x-m"""
   B = 2000
   ad = n/m^d % p
   f = lemme_21(n,d,ad,p,m)[0]
   g = parent(f)(p*x-m)
   s = skew_l2norm_tk (f)
   lognorm = flog (l2norm_tk (f, s))
   alp = alpha (f, B)
   E = lognorm + alp
   print "original polynomial: lognorm=", lognorm, " alpha=", alp, " E=", E
   w = 2
   r = lognorm_plus_alpha_rot_scons_linear (f, g, l2norm_tk, s, flog10(w))
   oldr = 2 * r
   while r < oldr:
      oldr = r
      w = 2 * w
      r = lognorm_plus_alpha_rot_scons_linear (f, g, l2norm_tk, s, flog10(w))
   # the constant polynomial varies up to V = sqrt(w*s), and the linear
   # polynomial up to U = sqrt(w/s)
   U = ceil(sqrt(w/s))
   V  = ceil(sqrt(w*s))
   Emin = E
   umin = 0
   vmin = 0
   print "rotation bounds:", U, V
   sys.stdout.flush()
   for u in range(-U,U+1):
      print "u=", u
      for v in range(-V,V+1):
         frot = f + parent(f)(u*x+v)*g
         lognorm = flog (best_l2norm_tk (frot))
         alp = alpha (frot, B)
         E = lognorm + alp
         if E < Emin:
            Emin = E
            umin = u
            vmin = v
            print "u=", umin, "v=", vmin, "lognorm=", lognorm, "alpha=", alp, "E=", Emin
            sys.stdout.flush()
   return p, m, umin*x+vmin, Emin


