load alpha.sage
load alpha_infty.sage

"""
input: lr the list of Laurent roots, each root given as [r,m,gap_]
       prec an integer
       q = base field cardinality
       s= skewness
output: c such that c[i]=Prob(i=max f_i a^i b^d-i - deg F(a,b))
        on a set of pairs in GF(q)[t] such that 
        << deg a- deg b =s >>.
        Exceptionally, c[prec] =Prob ( gap >=prec).
"""
def cancelations(lr,prec,q,s):
    c=[1]+[0]*prec
    for r_ in lr:
        r=r_[0]
        m=r_[1]
        Nr=m+deg(r)
        gap_=r_[2]
        if deg(r) == s:
            prob=1/(q+1)*q^(-Nr)
            c[gap_]+=prob
            c[gap_-1]-=prob
    return c


"""
 Murphi's E.
 input: f,g pair of polynomials in GF(q)[t][x]
        e= sieve paramete i.e. the pairs considered by the sieve have
        deg(a)+deg(b)<=2e
        s= skewness, i.e. max(deg(a))-max(deg(b))=s
        alphaf=alpha(f)
        beta= smoothness bound (in degree)
"""
def E(f,g,e,s,beta,alphaf):
    A=f.base_ring(); t=A.gen(); F=A.base_ring()
    q=F.cardinality(); S.<x>=A['x']
    f=S(f)
    d=f.degree()
    fc=f.coeffs()
    fcd=[ZZ(A(pp).degree()) for pp in fc]
    dg=S(g).coeffs()[0].degree()
    result=0
    drho=[]
    lr=Laurent_roots(f,prec)
    db0=ceil(e-s/2)
    da0=ceil(e+s/2)
    for j in range(d*(da0+10)+max(fcd)):
        drho.append(dickman_rho(max(0,j+alphaf)/beta))
    for skew in range(-db0,da0+1):
        c=cancelations(lr,d*e+max(fcd)+20,q,skew)
        for da in range(max(skew,0),min(db0-skew,da0)+1):
            db=da-skew
            sd=max([ZZ(j*da)+(d-j)*db+fcd[j]  for j in range(d+1)])
            # sd= symbolic degree
            for i in range(max(0,sd-len(drho)+1),min(sd,len(c))):
                old=result
                result+=((q-1)^2*q^(da+db)*c[i])*drho[sd-i]*dickman_rho((dg+db+1/(q-1))/beta)
                #print da,db,old-result,i
    return result 


""" EXAMPLE
p=2
load init_rings.sage
f=t^12 + t^10 + t^8 + t^2*x^5 + t*x^5 + x^6 + t^5 + x^5 + t^3 + t^2*x + t*x+t
g=t^104 + t^14 + t^13 + t^11 + t^10 + t^8 + t^7 + t^5 + t^4 + t^3 + t + x+1
f=S(f)
g=S(g)
alphaf=alpha(f,2000)
for s in [-1..7]:
   10^-6*E(f,g,27,s,25,alphaf)
# lpb0=22, ( I=J=12, deg(sq)=30 implies deg(a),deg(b)=27
"""

