load alpha.sage

# implements the formula on pages 86-87 of Murphy's thesis
# s is the skewness
def MurphyE(f,g,s):
    K = 1000
    # Bf = 4*10^5; Bg = 10^5; area = 1e13
    Bf = 1e7; Bg = 5e6; area = 1e16 # values used by pol51opt.c
    df = f.degree()
    dg = g.degree()
    alpha_f = alpha(f,2000)
    alpha_g = alpha(g,2000) # pol51opt.c uses alpha=0 for the linear polynomial
    E = 0
    sx = sqrt(area*s)
    sy = sqrt(area/s)
    for i in range(K):
       theta_i = float(pi/K*(i+1/2))
       xi = cos(theta_i)*sx
       yi = sin(theta_i)*sy
       fi = f(x=xi/yi)*yi^df
       gi = g(x=xi/yi)*yi^dg
       ui = (log(abs(fi))+alpha_f)/log(Bf)
       vi = (log(abs(gi))+alpha_g)/log(Bg)
       E += dickman_rho(ui)*dickman_rho(vi)
    return E/K

# example: RSA-768 polynomials
#skewness 44204.72 norm 1.35e+28 alpha -7.30 Murphy_E 3.79e-09
R.<x> = PolynomialRing(ZZ)
f = 265482057982680*x^6+1276509360768321888*x^5-5006815697800138351796828*x^4-46477854471727854271772677450*x^3+6525437261935989397109667371894785*x^2-18185779352088594356726018862434803054*x-277565266791543881995216199713801103343120
g=34661003550492501851445829*x-1291187456580021223163547791574810881
s=44204.72

f1=314460*x^5-6781777312402*x^4-1897893500688827450*x^3+18803371566755928198084581*x^2+2993935114144585946202328612634*x-6740134391082766311453285355967829910
g1=80795283786995723*x-258705022998682066594976199123
s1=1779785.90
# gives MurphyE(f1,g1,s1) = 2.21357103514175e-12


