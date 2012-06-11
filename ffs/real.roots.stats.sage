# usage: real_roots_stats(x^2+t^4+t+1,12,10,10,R)
# 12 is the maximal number of conditions we impose on a pair of polynomials
# (a,b) such that b^deg(f)*f(a/b) has some decrease in degree
# next two parameters are maximal degrees of a and b
# R=GF(2)['t,x']
# OUTPUT: a vector "cancelation" such that cancelation[k]=probabilitie that
# the cancelation (gap) is exactly k



load real_roots.sage


def all_roots_at_t_in_makefb_style(f,prec,R):
    F=R.base_ring()
    t,x=R.gens()
    At.<t>=F['t']
    S.<x>=At['x']
    f=S(f)
    df=f.derivative()
    l=roots_mod_t(S(f))
    truncroots=[]
    # truncroots= list of pairs 
    #(truncated root, father, precision, gap=0, proba=0)
    # where father(r)=index of r mod t^deg(r)
    for r in l:
        if At(df(r)) % t != 0:
            truncroots.append([r,None,1,0,0])
            rr=r
            pr=1
            while pr<prec:
                rr=At(hensel_lift(S(f),rr,pr,R))
                pr*=2
            for pri in [2..prec]:
                for j in range(len(truncroots)):
                    if (truncroots[j][0] == rr.truncate(pri-1)) and truncroots[j][2] == pri-1:
                        fatherind=j
                        break
                truncroots.append([rr.truncate(pri),fatherind,pri,0,0])
        else:
            pr=1
            lr=[[r,None,1,0,0]]
            rind=len(truncroots)
            while pr<prec:
                lr_new=[]
                for eind in range(len(lr)):
                    e=lr[eind]
                    if e[2] == pr:
                        rr=e[0]
                        lr_new+=[[rrr,eind+rind,pr+1,0,0] for rrr in naive_lift(f,rr,pr,R)]
                pr+=1
                lr+=lr_new
            truncroots+=lr
    return truncroots

def truncated_root_to_statistics(r,prec,A,B,At):
    # r is a truncated root which guarantees precision prec
    # presented as r=(r1,r2), root=r1(t)+r2(1/t)
    # At(resp. B) is the maximal degree of a (resp. b)
    r1,r2=r[0],r[1]
    assert prec>=1
    if r1 != 0:
		ubdb=min(B,A-r1.degree())
		lbdb=prec-r1.degree()
		# (u)lbdb = (upper) lower bound degree of b
		return sum([2^(i+(i+1+r1.degree()-prec)) for i in [lbdb..ubdb]])/2^(A+B+2) 
		# a*t^(prec-(1+r1.degree()))= b*(r1*t^(prec-(1+r1.degree()))+ rev(r2,prec-(1+r1.degree()))) +c
		# with deg(b)<=B0
		# deg(c)<=B0 and (c = -b*(r1*t^prec+rev(r2)) mod t^(prec-(1+r1.degree())
    else:
		if r2 == 0:
			v2=prec-1
		else:
			v2=valuation(At(r2),At(t))
		# v2= deg(b)-deg(a)
		ubdb=min(B, A+v2)
		lbdb=prec-1
		return sum([2^(i+(i-v2- (prec-1-v2))) for i in [lbdb..ubdb]])/2^(A+B+2)

def r_tilde_to_gap(r_tilde,prec,f_bar,s,R):
    F=R.base_ring()
    p=F.cardinality()
    t,x=R.gens()
    At.<t>=F['t']
    K=FractionField(At)
    S.<x>=K['x']
    r_tilde0=r_tilde+t^prec
    r_bar=K(r_tilde*At(t)^s)
    r_bar0=K(r_tilde0*At(t)^s)
    if r_bar != 0:
        value=K(S(f_bar)(r_bar))
    else:
        value=S(f_bar).constant_coefficient()
    if r_bar0 != 0:
        value0=K(S(f_bar)(r_bar0))
    else:
        value0=S(f_bar).constant_coefficient()
    # sage has a bug: cannot evaluate polynomials of S at 0
    fbcoeffs=[ K(e) for e in S(f_bar).coeffs()]
    gap=min(valuation(value,At(t)),valuation(value0,At(t)))-min([valuation(fbcoeffs[0],At(t))]+[valuation(fbcoeffs[i],At(t))+i*valuation(r_bar,At(t))
        for i in [1..S(f_bar).degree()]])
    return gap


def root_of_reversed_poly_to_real_root(s,r_tilde,R):
    t,x=R.gens()
    F=R.base_ring() 
    At.<t>=F['t']
    # s,m are such that f_tilde=t^m*f_bar(x*t^s)
    # parameter
    r_bar=r_tilde*t^s
    r=r_bar(1/t)
    num=r.numerator()
    den=r.denominator()
    r1=num.quo_rem(den)[0]
    aux=(num % den)/den
    r2= At(aux(1/t))
    return [r1,r2]

def tree_to_cancelation(truncroots,N):
    cancelation=[0]*(N+1)
    for gap in range(1,N+1):
        for i in range(len(truncroots)):
            e=truncroots[i]
            if e[3] != gap or (e[1] != None and truncroots[e[1]][3] == gap):
                continue
            goodsons=[]
            for j in range(len(truncroots)):
                if truncroots[j][1] == i and truncroots[j][3]> gap:
                    goodsons.append(j)
            cancelation[gap]+=e[4]-sum([truncroots[j][4] for j in goodsons])
    cancelation[0]=1-sum(cancelation)
    return cancelation

def real_roots_stats(f,max_prec,A,B,R):
    max_prec=min(max_prec,A,B)
    F=R.base_ring()
    p=F.cardinality()
    t,x=R.gens()
    At.<t>=F['t']
    K=FractionField(At)
    S.<x>=K['x']
    f=S(f)
    f_bar=S(R(f)(1/t,x))
    f_tilde,s,m=rotate(f_bar,R)
    truncroots=all_roots_at_t_in_makefb_style(f_tilde,max_prec,R)
    N=S(f).degree()*max(A,B)+max([At(e).degree() for e in S(f).coefficients()])+1
    for i in range(len(truncroots)):
        r_tilde=truncroots[i][0]
        prec=truncroots[i][2]
        r=root_of_reversed_poly_to_real_root(s,r_tilde,R)
        truncroots[i][3]=r_tilde_to_gap(r_tilde,prec,f_bar,s,R)
        truncroots[i][4]=truncated_root_to_statistics(r,prec,A,B,At).numerical_approx()
    cancelation=tree_to_cancelation(truncroots,N)
    return cancelation

