def listB1 (n=30):
  B1 = 105.
  L = []
  for _ in xrange(n):
    L.append(B1.floor())
    B1 += sqrt(B1)
  return L

print listB1()
