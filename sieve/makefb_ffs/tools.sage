def Primes(A,alim):
	if A==ZZ:
		return primes(alim)
	else:
		res=[]
		for a in element_list(A,alim):
			if (not a.is_zero()) and (a.is_irreducible()) and (a.is_monic()):
				res.append(a)
	return res
    

def bar(f,p):
    A=f.base_ring()
    if A==ZZ:
        res=f.change_ring(GF(p))
    else:
        t=A.gen()
        q=A.base_ring().characteristic()
        kp.<w>=GF(q^p.degree(),name='w',modulus=p)
        if w==GF(q)(1):
            w=p.roots(GF(q))[0][0] # convention of sage is w=1 when the residual field is GF(q)
        newX.<x>=kp['x']
        res=0
        i=0
        for c in f.coeffs():
            res+=c.change_ring(kp)(t=w)*x^i
            i+=1
        res=newX(res)
    return res


def log_norm(p):
    A=p.parent()
    if A==ZZ:
        res=log(p)
    else:
        F=A.base_ring()
	q=len(F)
        res=p.degree()*log(q)
    return res

def element_list(A,alim):
	if A==ZZ:
		return range(alim)
	else:
		q=len(A.base_ring())
		F=A.base_ring()
		t=A.gen()
		l=F.list()
		res=[]
		for i in range(alim):
			iPol=A(0)			
			k=0
			dig=ZZ(i).digits(q)
			for j in dig:
				iPol+=l[j]*t^k
				k+=1
			res+=[iPol]
		return res			

def cast(d,A,q): # GF(q) to ZZ for small q
	l=element_list(A,q)
	i=0
	for j in l:
		if j==d:
			return i
		i+=1
		if i>q:
			j='error'
	print 'cast error'
	return 0	
	
def position(a):
	if (a+0).parent()==ZZ:#int+0=integer
		return a
	else:
		A=(a+0).parent()
		q=ZZ(len(A.base_ring()))
		res=0
		l=a.coeffs()
		k=0
		for d in l:
			res+=cast(d,A,q)*q^k
			k+=1 		
		return res

def parse_line(l,A=ZZ):
	l=l.split('\n')[0] # eliminating '\n'	
	words=l.split(':')
	q=A(words[0])
	if len(words)==3:
		n=[ZZ(i) for i in words[1].split(',')]
		r=[A(ri) for ri in words[2].split(',')] 
	else:
		n=[1,0]
		r=[A(ri) for ri in words[1].split(',')]
	return q,n,r 	


