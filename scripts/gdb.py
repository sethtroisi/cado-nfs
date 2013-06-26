# This script allows printing MPZ and MPFR types in a readable way in gdb.
# Written by E. Thomé, with minor changes by A. Kruppa
# Use in gdb by saying:
#   source path/to/this/script/gdb.py
# Then all mpz_t etc. should magically be printed nicely.

# Reference on python scripting:
# http://sourceware.org/gdb/current/onlinedocs/gdb/Python.html

# Prefer gdb 7.2 for this, because of the following bug (which bites
# you if you have mpz_t's in structs):
# http://comments.gmane.org/gmane.comp.gdb.patches/59238
# Otherwise gdb 7.1 will work.

import gdb
from gmpy import mpz,mpf

__all__=['update_display_limit', 'make_mp_printer_objects', 'hook_mp_printers', 'unhook_mp_printers' ]

try:
    gdb.pretty_printers.remove(make_mp_printer_objects)
except NameError:
    pass

# arrange so that at most 40 characters of each field are printed.
display_limit=40

def gmpy_from_mpn(d, nlimbs, wordsize):
    res=mpz(0)
    word=pow(mpz(2),wordsize)
    for i in range(0,nlimbs):
        a=mpz(int(d[nlimbs-1-i]))
        # python ints are signed, beware.
        if a < 0: a+=word
        res=res*word+a
    return res

def truncate_output(r):
    # not completely satisfactory
    if display_limit == 0:
        return r
    l=len(r)
    if l >= display_limit:
        r=r[0:display_limit-15] +"...<%d>..."%(l-(display_limit-10)) +r[l-5:l]
    return r
    
def update_display_limit(l):
    global display_limit
    display_limit=l

class mpfr_printer:
    def __init__(self, val):
        self.val = val
    def mantissa_exponent(self):
        X = self.val
        sign = int(X['_mpfr_sign'])
        exp =  int(X['_mpfr_exp'])
        prec = int(X['_mpfr_prec'])
        d =    X['_mpfr_d']
        wordsize=64
        nlimbs=(prec+wordsize-1)/wordsize
        # print "(via %s) Have prec %d, %d limbs\n" % (self.vt, prec,nlimbs)
        # try:
        mantissa=gmpy_from_mpn(d, nlimbs, wordsize)
        mantissa*=sign
        e=exp-nlimbs*wordsize
        return mantissa, e
    def to_string(self):
        X = self.val
        mantissa,e = self.mantissa_exponent()
        prec = int(X['_mpfr_prec'])
        exp =  int(X['_mpfr_exp'])
        wordsize=64
        special=-pow(2,wordsize-1)
        sign = int(X['_mpfr_sign'])
        if exp == special+2:
            return "NaN"
        if exp == special+1:
            if sign < 0:
                return "-0"
            else:
                return "0"
        if exp == special+3:
            if sign < 0:
                return "-inf"
            else:
                return "+inf"
        res=mpf(mantissa, prec)
        if e>0: res*=pow(mpz(2),e)
        else: res/=pow(mpz(2),-e)
        return truncate_output(str(res))

class mpz_printer:
    def __init__(self, val):
        self.val = val
    def to_string(self):
        X = self.val
        # There's apparently a bug in array member of structs. Their
        # starting address seems to be constantly equal to the starting
        # address of the struct itself...
        # print "X at %s\n" % X.address
        alloc = int(X['_mp_alloc'])
        size = int(X['_mp_size'])
        d =    X['_mp_d']
        wordsize=64
        nlimbs=size
        if size<0: nlimbs=-int(size)
        # try:
        mantissa=gmpy_from_mpn(d, nlimbs, wordsize)
        # except RuntimeError:
        # # it's not necessarily a good idea to do this.
        # return "<error>"
        if size<0: mantissa=-mantissa
        return "mpz_t(alloc=%d, size=%d, d=%x) = %s" % \
          (alloc, size, long(d), truncate_output(str(mantissa)))

class mpc_printer:
    def __init__(self, val):
        # print val['re'].address
        self.re = mpfr_printer(val['re'].dereference())
        self.im = mpfr_printer(val['im'].dereference())
    def to_string(self):
        return self.re.to_string() + "+i*" + self.im.to_string()

class mpfrx_printer:
    def __init__(self, val):
        self.val = val
    def to_string(self):
        X = self.val
        n=int(X['deg'])+1
        coeffs=[mpfr_printer(X['coeff'][i]).to_string() for i in range(n)]
        res="(%s %s)" % (n-1, " ".join(coeffs))
        return res

def make_mp_printer_objects(val):
    try:
        # beware. If we hook on mpfr_t to remove the extra braces, then
        # some bugs pop up -- apparently dereferencing this pointer does
        # not work as it should in arrays of mpfr_t's. Haven't checked,
        # but I assume it's similar for other types.
        t=str(val.type)
        # print "[request for %s]\n" % t
        if (t == 'mpz_ptr' or t == 'mpz_srcptr'):
            return mpz_printer(val.dereference())
        if (t == 'mpz_t' or t == '__mpz_struct [1]'):
            return mpz_printer(val.dereference())
        if (t == '__mpz_struct'):
            return mpz_printer(val)

        if (t == 'mpfr_ptr' or t == 'mpfr_srcptr'):
            return mpfr_printer(val.dereference())
        if (t == 'mpfr_t' or t == '__mpfr_struct [1]'):
            return mpfr_printer(val.dereference())
        if (t == '__mpfr_struct'):
            return mpfr_printer(val)

        if (t == 'mpc_ptr' or t == 'mpc_srcptr'):
            return mpc_printer(val.dereference())
        if (t == 'mpc_t' or t == '__mpc_struct [1]'):
            return mpc_printer(val.dereference())
        if (t == '__mpc_struct'):
            return mpc_printer(val)

        if (t == 'mpfrx_ptr' or t == 'mpfrx_srcptr'):
            return mpfrx_printer(val.dereference())
        if (t == 'mpfrx_t' or t == '__mpfrx_struct [1]'):
            return mpfrx_printer(val.dereference())
        if (t == '__mpfrx_struct'):
            return mpfrx_printer(val)

    except RuntimeError:
    # constructors may abandon building if the object looks too complicated.
        return None
    return None

def remove_all_printers():
    while len(gdb.pretty_printers):
        gdb.pretty_printers.pop(0)

def hook_mp_printers():
    gdb.pretty_printers.append(make_mp_printer_objects)

def unhook_mp_printers():
    gdb.pretty_printers.remove(make_mp_printer_objects)

# this is just for easily sourcing this again and again while debugging.
remove_all_printers()
hook_mp_printers()

