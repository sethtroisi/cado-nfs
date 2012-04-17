load ../makefb/tools.sage
def sieve(alim,A,filename,thr=3,q=2):
    domain=[[ZZ(0) for j in range(alim)] for i in range(alim)]
    fd=open(filename,"r").readlines()
    p=0
    for line in fd:
        if line[0]=='#':
            continue
        pk,n,r=parse_line(line,A)
        if p==0 or pk % p !=0:
            p=pk
            if q^(p.degree())>alim:
                line="error"
                continue
        a0_b0=[]
        for ri in r:
	    	if ri.degree()<pk.degree():
		    	for b0 in element_list(A,q^pk.degree()): # deg(b0)<deg(pk)
    				a0_b0.append([ri*b0 % pk, b0])
	    	else:
	    		for a0 in element_list(A,q^pk.degree()):
	    			a0_b0.append([a0,ri*a0 % pk]) #(ri-pk)%pk
    	for a0b0 in a0_b0:
    		a0,b0=a0b0[0],a0b0[1]
    		for i in element_list(A,floor(alim/q^pk.degree())):			
    			a=a0+i*pk				
    			for j in element_list(A,floor(alim/q^pk.degree())):
    				b=b0+j*pk
    				if q^(a.degree()+1)<alim and  q^(b.degree()+1)<alim:
    					domain[position(a)][position(b)]+=ZZ(n[0]-n[1])*p.degree()
    res=[]
    l=GF(q).list()
    for i in range(alim):
        for j in range(alim):
            if domain[i][j]>thr:
                res.append([int2pol(q,ZZ(i).digits(q),l,A),int2pol(q,ZZ(j).digits(q),l,A)])
    return res
