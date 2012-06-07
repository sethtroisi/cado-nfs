load real_roots.sage



def all_roots_at_t_in_makefb_style(f,prec,R):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    S.<x>=A['x']
    f=S(f)
    df=f.derivative()
    l=roots_mod_t(S(f))
    unramified=[]
    ramified=[]
    for r in l:
        if A(df(r)) % t != 0:
            unramified.append(r)
        else:
            pr=1
            lr=[[r,None,1]]
            # lr= list of pairs (truncated root, father, precision)
            # where father(r)=r mod t^deg(r)
            while pr<prec:
                lr_new=[]
                for e in lr:
                    if e[2] == pr:
                        rr=e[0]
                        lr_new+=[[rrr,rr,pr+1] for rrr in naive_lift(f,rr,pr,R)]
                pr+=1
                lr+=lr_new
            ramified+=lr
    return unramified,ramified

def delete_fathers(ramified):
    # an element of ramified is of type
    # (r,father(r),precesion)
    # if r has a unique son, then we replace all the triples
    # (r',son(r),precision') by (r',r,precision')
    clean_ramified=[]
    for e1 in ramified:
        flag=true
        for e2 in ramified:
            if e2[1] == e1[0]:
                flag=false
                break
        if flag:
            clean_ramified.append(e1)
    return clean_ramified

def truncated_root_to_statistics(r,prec,AA,B):
    # r is a truncated root which guarantees precision prec
    # presented as r=(r1,r2), root=r1(t)+r2(1/t)
    # A (resp. B) is the maximal degree of a (resp. b)
    # we approximate the proportion of pairs (a,b) with  
    r1,r2=r[0],r[1]
    if prec<=0:
        return 0
    if r1 != 0:
        B0=min(B,AA-r1.degree())
    else:
        B0=min(AA,B-(prec-r2.degree()))
    # a*t^prec= b*(r1*t^prec+ rev(r2)) +c
    # with deg(b)<=B0
    # deg(c)<=B0 and (c = -b*(r1*t^prec+rev(r2)) mod t^prec
    if B0<prec:
        return 0
    else:
        return 2^(B0+(B0-prec))/2^(AA+B)

def root_of_reversed_poly_to_real_root(s,r_tilde,A):
    t=A.gen()
    # s is such that f_tilde=t^m*f_bar(x*t^s) for some m
    r_bar=r_tilde*t^s
    r=r_bar(1/t)
    num=r.numerator()
    den=r.denominator()
    r1=num.quo_rem(den)[0]
    aux=(num % den)/den
    r2= A(aux(1/t))
    return [r1,r2]

def real_roots_stats(f,max_prec,AA,B,R):
    F=R.base_ring()
    p=F.cardinality()
    t,x=R.gens()
    A.<t>=F['t']
    S.<x>=A['x']
    f=S(f)
    f_bar=R(f)(1/t,x)
    f_tilde,s=rotate(f_bar,R)
    u,ram=all_roots_at_t_in_makefb_style(f_tilde,max_prec+abs(s),R)
    N=f.degree()*max(AA,B)
    cancelation0=[0]*N
    for e in ram:
        r_tilde=e[0]
        prec=e[2]-abs(s)
        r=root_of_reversed_poly_to_real_root(s,r_tilde,A)
        cancelation0[prec]+= truncated_root_to_statistics(r,prec,AA,B)
    for prec in range(1,max_prec+1):
        r=root_of_reversed_poly_to_real_root(s,A(1+t^prec),A)
        # all the unramified roots behave the same,
        # take 1+t^prec as a surrogate
        cancelation0[prec]+=len(u)* truncated_root_to_statistics(r,prec,AA,B)
        # cancelation 0 = proba(at least a cancelation of prec)
    cancelation=[0]+[cancelation0[i]-cancelation0[i+1] for i in  range(1,N-1)]
    cancelation[0]=1-sum(cancelation)
    return cancelation


