# The characteristic in global variable.
# Give a stupid name so that there is no collision with other libs.
__p_charac_in_ffs = 2

def str2fppol(s):
    p = __p_charac_in_ffs
    i = Integers()(s, 16)
    if p == 2:
        dig = i.digits(p)
    else:
        assert (p == 3)
        dig = i.digits(4)
    return GF(p)['t'](dig)

def str2ffspol(s):
    p = __p_charac_in_ffs
    Fpt.<t> = GF(p)['t']
    Fptx.<x> = Fpt['x']
    ss = s.split(',')
    R = Fptx(0)
    for i in range(0,len(ss)):
        R += str2fppol(ss[i])*x^i
    return R

def readparam(paramfile):
    intparam = [ "I", "J", "sqside", "firstside", "S",
            "fbb0", "lpb0", "thresh0", "powerlim0",
            "fbb1", "lpb1", "thresh1", "powerlim1",
            ]
    fppolparam = [ "q", "rho", "q0", "q1" ]
    ffspolparam = [ "pol0", "pol1" ]
    f = open(paramfile, "r")
    lines = f.readlines()
    f.close()
    for l in lines:
        if l[0] == '#' or l.rstrip() == "":
            continue;
        w = l.rstrip().split('=')
        if len(w) != 2:
            print "Error while parsing " + l
        else:
            # For known parameters, we convert them appropriately. 
            # The others are set as strings.
            if w[0] in intparam:
                globals()[w[0]] = Integers()(w[1])
            elif w[0] in fppolparam:
                globals()[w[0]] = str2fppol(w[1])
            elif w[0] in ffspolparam:
                globals()[w[0]] = str2ffspol(w[1])
            else:
                globals()[w[0]] = w[1]

