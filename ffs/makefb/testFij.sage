load "readparam.sage"

def skewgauss(q, rho, S):
    Fpt = q.parent()
    t = Fpt.gen()
    a0 = q
    b0 = Fpt(0)
    a1 = rho
    b1 = t^S
    while True:
        qq = a0 // a1
        if qq.degree()+b1.degree() > a0.degree():
            break
        a0 -= qq*a1
        b0 -= qq*b1
        qq = a1 // a0
        if qq.degree()+b0.degree() > a1.degree():
            break
        a1 -= qq*a0
        b1 -= qq*b0
        if a0.degree() <= b0.degree():
            break
    return a0, b0 // t^S, a1, b1 // t^S

# this will define, among others, pol0
readparam("../param.2.607")
t = q0.parent().gen()
x = pol0.parent().gen()

q = t^25 + t^3 + 1
rho = t^24 + t^23 + t^22 + t^21 + t^19 + t^16 + t^15 + t^13 + t^11 + t^10 + t^7 + t^3

a0, b0, a1, b1 = skewgauss(q, rho, 0)

a = x*a0+a1
b = x*b0+b1

Fij = pol0.parent()(0)
C = pol0.coeffs()
d = len(C)-1
for k in range(0, len(C)):
    Fij += C[k]*a^k*b^(d-k)
