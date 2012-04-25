q=3
F=GF(q)
A.<t>=F['t']
S.<x>=A['x']
if q==2:
    F = x^5 + x^4 + x^3 + x^2 + x + t^2
    G = x*t^25 + x*t^24 + x*t^23 + x*t^21 + x*t^17 + x*t^16 + x*t^15 + x*t^13 + x*t^10 + x*t^9 + x*t^7 + x*t^5 + x*t^3 + x + t^26 + t^25 + t^24 + t^23 + t^21 + t^20 + t^17 + t^14 + t^13 + t^11 + t^9 + t^6 + t^4 + t^3 + t^2;
else:
    F = x^5 + x*t + x + t^2
    G = x*t^25 + x*t^23 + 2*x*t^22 + 2*x*t^21 + 2*x*t^20 + x*t^18 + 2*x*t^17 + x*t^16 + 2*x*t^14 + 2*x*t^13 + 2*x*t^11 + x*t^10 + 2*x*t^9 + 2*x*t^8 + x*t^7 + 2*x*t^6 + x*t^5 + x*t^4 + 2*x*t^2 + x*t + 2*x + t^26 + 2*t^25 + t^23 + t^22 + 2*t^20 + t^19 + t^18 + 2*t^17 + t^16 + t^12 + t^11 + 2*t^10 + t^8 + t^7 + 2*t^6 + t^4 + 2*t^3 + t + 2


if q==2: 
    dlim = 6
elif q==3:
    dlim= 4
else: 
    print "Error: p must be 2 or 3\n"

powerlim=1 #no powers
load makefb.sage
makefb(F,dlim,powerlim,"Aroots","cado")
makefb(G,dlim,powerlim,"Rroots","cado")
