load sieve.sage

q=2
F=GF(q)
A.<t>=F['t']
S.<x>=A['x']

print "# f=x^5 + x^2*t + x^2 + t^2"
f=x^5 + x^2*t + x^2 + t^2
print "making fb"
makefb(f,10,3,"F2_607","human")
print "sieving on 100 x 100 pairs with a trheshold 5"
l=sieve(100,A,"F2_607",5)
print "number of pairs found:",len(l)
print "\n# f=(t^2 + t + 1)*x^5 + t^2*x^4 + t^2*x^3 + t^2*x + t^2 + 1"
f=(t^2 + t + 1)*x^5 + t^2*x^4 + t^2*x^3 + t^2*x + t^2 + 1
print "making fb"
makefb(f,10,3,"paul_poly","human")
print "sieving on 100 x 100 pairs with a trheshold 5"
l=sieve(100,A,"paul_poly",5)
print "number of pairs found:",len(l)
