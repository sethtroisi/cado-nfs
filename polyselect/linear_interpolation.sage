
# This file works with a convex (or concave) function initially sampled
# at integer values. It aims at providing a good piecewise linear
# approximation.

def quick_eval(t,x):
    """ Given a table whose values are some f(x) at integer points
    between 0 and len(t)-1, returns f(x) at an arbitratry real number x,
    assuming piecewise linear interpolation.
    The value is returned as a vector.
    """
    f=floor(x)
    m=len(t)-2
    if f < m:
        v=[t[f][i]+(x-f)*(t[f+1][i]-t[f][i]) for i in range(len(t[f]))]
    else:
        f=m-1
        v=[t[f][i]+(x-f)*(t[f][i]-t[f-1][i]) for i in range(len(t[f]))]
    return vector(v)

def max_trichotomy(f,a,b,epsilon):
    """
    Computes the maximum of f within the interval [a,b]. Assumes that the
    derivative of f stays negative. Returns the maxpoint index, and the
    corresponding value. Result is accurate to precision epsilon.
    """
    ea=f(a); eb=f(b)
    while b-a > epsilon:
        third=(b-a)/3
        c=a+third; ec=f(c)
        d=c+third; ed=f(d)
        if ec > ed:
            b=d; eb=ed
        else:
            a=c; ea=ec
    c=(a+b)/2
    return (float(c),float(f(c)))

def min_trichotomy(f,a,b,epsilon):
    c=max_trichotomy(lambda x: -f(x), a,b,epsilon)
    return (c[0],-c[1])

def worst_approximation(table,a,b):
    """ Given a table whose values are some f(x) at integer points
    between 0 and len(t)-1, returns a value such that the linear
    interpolation between points a and b gives the largest error.
    This assumes that the worst approximation is reached just once.
    """
    a0=a; la=1; wa=quick_eval(table,a)
    b0=b; lb=0; wb=quick_eval(table,b)
    # a == la * a0 + (1-la) * b0
    # b == lb * a0 + (1-lb) * b0
    while b-a > 0.01:
        third=(b-a)/3
        c=a+third; lc=(la+la+lb)/3; vc=quick_eval(table,c)
        d=c+third; ld=(la+lb+lb)/3; vd=quick_eval(table,d)
        # linear interpolations:
        wc=(wa+wa+wb)/3; ec=abs(wc-vc)
        wd=(wa+wb+wb)/3; ed=abs(wd-vd)
        if ec > ed:
            b=d; lb=ld; wb=wd
        else:
            a=c; la=lc; wa=wc
    c=(a+b)/2; vc=quick_eval(table,c); wc=(wa+wb)/2
    return (float(c),float(abs(wc-vc)))

def worst_approximation_piecewise(table,lc):
    """ Returns the worst approximation from the piecewise segments given
    by lc"""
    w=[worst_approximation(table,lc[i],lc[i+1]) for i in range(len(lc)-1)]
    return max([(z[1],z[0]) for z in w])

def cutoffs_inner(table,a,b,e):
    """ Internal function. """
    c=worst_approximation(table,a,b)
    if (c[1] < e):
        v=[]
    else:
        v=cutoffs_inner(table,a,c[0],e) + [c[0]] + cutoffs_inner(table,c[0],b,e)
    return v

def find_interpolation_points(t,e):
    """
    Given a table whose values are some f(x) at integer points between 0
    and len(t)-1, returns the list of interpolation points necessary for
    obtaining a piecewise linear approximation of accuracy e.
    """
    return [0] + cutoffs_inner(t,0,len(t)-1,e) + [len(t)-1]

def reduced_interpolation_table(table,epsilon):
    """
    Gives a reduced interpolation table, with integer interpolation
    points.
    """
    lc=find_interpolation_points(table,epsilon)
    lc0=[floor(round(z)) for z in lc]
    return [(z,quick_eval(table,z)) for z in lc0]
