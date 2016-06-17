R.<t> = GF(2)[]
S.<x> = R[]
F = x^5 + x^2*t + x^2 + t^2
G = x*t^121 + x*t^120 + x*t^117 + x*t^116 + x*t^115 + x*t^114 + x*t^113 + x*t^112 + x*t^108 + x*t^106 + x*t^105 + x*t^103 + x*t^102 + x*t^101 + x*t^99 + x*t^90 + x*t^89 + x*t^88 + x*t^87 + x*t^85 + x*t^80 + x*t^79 + x*t^78 + x*t^77 + x*t^76 + x*t^74 + x*t^70 + x*t^66 + x*t^65 + x*t^62 + x*t^58 + x*t^56 + x*t^52 + x*t^51 + x*t^46 + x*t^44 + x*t^42 + x*t^41 + x*t^40 + x*t^38 + x*t^37 + x*t^36 + x*t^35 + x*t^33 + x*t^30 + x*t^27 + x*t^26 + x*t^21 + x*t^18 + x*t^17 + x*t^16 + x*t^15 + x*t^14 + x*t^11 + x*t^10 + x*t^7 + x*t^6 + x*t^5 + x + t^122 + t^120 + t^117 + t^116 + t^115 + t^113 + t^110 + t^105 + t^104 + t^102 + t^101 + t^97 + t^95 + t^94 + t^93 + t^90 + t^82 + t^80 + t^79 + t^77 + t^76 + t^75 + t^74 + t^73 + t^72 + t^70 + t^68 + t^67 + t^64 + t^58 + t^54 + t^53 + t^52 + t^51 + t^49 + t^45 + t^44 + t^43 + t^42 + t^41 + t^40 + t^39 + t^38 + t^37 + t^35 + t^31 + t^30 + t^28 + t^21 + t^19 + t^18 + t^15 + t^13 + t^11 + t^10 + t^9 + t^8 + t^7 + t^6

def int2pol(q):
   l = ZZ(q).bits()
   return sum([l[i]*t^i for i in range(len(l))])

def roots(q):
   Q = int2pol(q)
   U.<tbar> = GF(2^Q.degree(),name='tbar',modulus=Q)
   V.<x> = U[]
   FF = V(F)
   l = FF.roots()
   for r in l:
      z = r[0].polynomial().coeffs()
      rho = sum([ZZ(z[i])*2^i for i in range(len(z))])
      print hex(ZZ(q)), hex(rho)
      os.system("./sieve2_607 %s %s" % (hex(ZZ(q)), hex(rho)))

def fb(q0,q1):
   for q in range(q0,q1):
      Q = int2pol(q)
      if Q.is_irreducible():
         roots(q)
