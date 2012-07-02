try:
    F=GF(p)
except NameError:
    F=GF(2)
A.<t>=F['t']
K.<t>=FractionField(A)
S.<x>=A['x']
R.<t,x>=F['t,x']
