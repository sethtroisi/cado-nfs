
"""
writes f in GF(q)[t][x] in hexadecimal form
"""
def hexify(f):
    A=f.base_ring()
    q=A.base_ring().cardinality()
    dummy=2^ceil(log(q,2))
    fc=S(f).coeffs()
    result=''
    aux=0
    if f.degree() < 1:
        return hex(ZZ(f))
    l=fc[0].coeffs()
    for j in range(len(l)):
        aux+=ZZ(l[j])*dummy^j
    result=result+hex(aux)
    for i in range(1,len(fc)):
        aux=0
        l=fc[i].coeffs()
        for j in range(len(l)):
            aux+=ZZ(l[j])*dummy^j
        result=result+','+hex(aux)
    return result

"""
associates an element of GF(q)[t] to a string
"""
def sigma_hex(hexpol,q=2):
    A.<t>=GF(q)['t']
    dummy=2^ceil(log(q,2))
    l=A.base_ring().list()
    return A(ZZ(int(hexpol,16)).digits(dummy))

""" EXAMPLES
p=49
load init_rings.sage
hexify(t^2+4*t+3)
sigma_hex('21',2)

"""
