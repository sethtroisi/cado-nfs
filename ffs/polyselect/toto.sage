for a in element_list(A,x0+x1,x0):
	for b in element_list(A,y0+y1,y0):
		if gcd(a,b) == 1:
			v=valuation(f.resultant(a-b*x), p)
			if v>0:
				print a,b
			s+=v 
			n+=1
