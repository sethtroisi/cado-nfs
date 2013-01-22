load alpha4ffs.sage
q=2
F=GF(q)
A.<t>=F['t']
S.<x>=A['x']
print "# A polynomial f has a good root property if alpha(f) is negative and small. For instance if the leading or the constant coefficients have many small divisors, f has a good root property."
f=x^5+x+t+1
print "# f=",f
g=sigma("1b92c17dec4c4cf4f5ab9c1e86f")*x+sigma("d0e134790925d9e08")
print "# g=sigma(1b92c17dec4c4cf4f5ab9c1e86f)*x+sigma(d0e134790925d9e08)"
print "\nf and g is the pair used by Joux and Lercier in 2002 for GF(2^531)"
print "alpha(f,100)=",alpha(f,100)
print "alpha(g,100)=",alpha(g,100)
print "alpha_p(f,t+1)=",alpha_p(f,t+1)
print "estimate_alpha_p(f,t+1,120)=",estimate_alpha_p(f,t+1,120)
