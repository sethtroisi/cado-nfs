# finds the best LLL-reduced form for skewness 'skew' and translation 'k'
def size_optimization_aux(f,g,skew,k=0):
   d = f.degree()
   max_rot = d-2
   m = matrix (ZZ, d+1, max_rot+2)
   x = f.variables()[0]
   ft = f(x=x+k)
   gt = g(x=x+k)
   for i in [0..d]:
      m[i,0] = ft[i]*skew^i
   for j in [0..max_rot]:
      m[j,j+1] = gt[0]*skew^j
      m[j+1,j+1] = gt[1]*skew^(j+1)
   m = m.transpose().LLL().transpose()
   minnorm = infinity
   for j in range(max_rot+2):
      v = m.column(j)
      if v[d] == 0:
         continue
      norm = sum(v[i]^2 for i in [0..d])
      if norm < minnorm:
         minnorm = norm 
         bestf = sum((v[i]//skew^i)*x^i for i in [0..d])
   return bestf, gt, sqrt(1.0*minnorm)

# once a good LLL-reduced form is found, tries to minimize the degree d-2
# coefficient using translation, and then iterates until a minimum is found
def size_optimization_aux2(f,g,skew,k=0):
   x = f.variables()[0]
   f = f(x=x+k)
   g = g(x=x+k)
   d = f.degree()
   l = size_optimization_aux(f,g,skew)
   var('t')
   while True:
      best = l
      fopt = l[0]
      eq = fopt[d-2] + (d-1)*t*fopt[d-1] + d*(d-1)/2*t^2*fopt[d]
      for (x,_) in SR(eq).roots(ring=RR):
         for kk in [floor(x),ceil(x)]:
            ll = size_optimization_aux(f,g,skew,kk)
            if ll[2] < l[2]:
               best = ll
      if best == l:
         break
      l = best
   return l
