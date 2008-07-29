attach "alpha.sage"
attach "rotation_bound.sage"

# for reading kleinjung.res
import re

def lemme_21(n,d,ad,p,m):
    """Implements Kleinjung's lemma 2.1. ad must satisfy ad*m^d==n mod p.
    ad may also be equal to zero, in which case it is recovered
    automatically"""
    if ad==0:
        ad=ZZ(Integers(p)(n/m^d))
    assert(ad*Integers(p)(m)^d-n==0)
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
    g=Integers()['x']([-m,p])
    #return (f,deltas,rs)
    return (f,g)

# For example:

# n=1230186684530117755130494958384962720772853569595334792197322452151726400507263657518745202199786469389956474942774063845925192557326303453731548268507917026122142913461670429214311602221240479274737794080665351419597459856902143413
# f=ZZ['x']([-277565266791543881995216199713801103343120,-18185779352088594356726018862434803054,6525437261935989397109667371894785,-46477854471727854271772677450,-5006815697800138351796828,1276509360768321888,265482057982680])
# g=ZZ['x']([-1291187456580021223163547791574810881,34661003550492501851445829])
# f0,_=lemme_21(n,6,f[6],g[1],-g[0])
# (f-f0)/g

# Do for instance:
# tuples=import_kleinjung_dot_res(n,d,"/net/tiramisu/cado1/cado/Examples/35_143+_c136/kleinjung.res", RR('1.0e+19'))
def import_kleinjung_dot_res(n,d,filename,normbound):
    """This function reads and parses the kleinjung.res file. Only
    reports below the specified normbound are processed."""
    tuples=[]
    fi=open(filename,"r")
    for s in fi.readlines():
        ma=re.match("^p=(\d+) m=(\d+) norm=([\d\.]+e\+\d+)$",s)
        # No match is something normal
        if ma == None: continue
        tu=ma.groups()
        norm=RR(tu[2])
        if norm >= normbound: continue
        p=ZZ(tu[0])
        m=ZZ(tu[1])
        s0= "p=%r m=%r supnorm=%.2e"%(p,m,norm)
        print s0
        f,g=lemme_21(n,d,0,p,m)
        #tuples.append((p,m,norm))
        tuples.append((f,g))
    # tuples.sort()
    fi.close()
    return tuples

def optimize_tuple(f,g):
    """Optimizes the polynomial pair f,g. Return the best logmu+alpha"""
    p=g[1]
    m=-g[0]
    s=skew_l2norm_tk(f)
    a0=alpha_projective(f,2000)
    a1=alpha_affine(f,2000)
    toto=lognorm_plus_alpha_rot_scons_linear
    c=[(a0+toto(f,g,l2norm_tk,s,i),i) for i in [0..14]]
    cm=min(c)
    s0= "p=%r m=%r"%(p,m)
    s1= "%.2f -> %.2f [10^%d]" %(flog(best_l2norm_tk(f))+a1,cm[0],cm[1])
    print s0 + " " + s1
    return cm

def optimize_tuples(l):
    """Given a list as returned by import_kleinjung_dot_res, optimizes
    all one by one"""
    cms=[(optimize_tuple(l[i][0],l[i][1]),i) for i in range(len(l))]
    mi=min(cms)
    i=mi[1]
    f=l[i][0]
    g=l[i][1]
    y=mi[0][0]
    w=mi[0][1]
    print "min: [#%d] p=%r m=%r %.2f [10^%d]" % (i,p,m,y,w)

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

