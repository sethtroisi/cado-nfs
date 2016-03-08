q=2
F=GF(q)
A.<t>=F['t']
S.<x>=A['x']
f=x^3+t+1
p=t+1
dlim=4
powerlim=7
load makefb.sage
makefb(f,dlim,powerlim,"","human")
