To construct a factor base, type the following in sage:

  load "makefb.sage"
 
  makefb(f,dlim,powerlim,filename="",typo="cado")

where the parameters are:
  f        = the FFS polynomial
  dlim     = max degree of the factor base
  powerlim = max degree of the irreducible powers which occur in makefb
  filename =[optional] name of the output file, put "" or no parameter for
                       standard output
  typo     =[optional] output format: "human" for human readable,
                       "cado" [by default] is like cado but in hexadecimal.


Example:
  
  q=2
  A.<t>=GF(2)['t']
  x=A['x'].gen()
  f= x^3+t+1
  dlim=13
  powerlim=12

  makefb(f,dlim,powerlim, typo="human")


---------------------------------------------------------------------------

We recall here the CADO format for the factor base: 
(cado-nfs/sieve/README.fb_format)


Factor base file format:
------------------------

An entry is of the form:

q:n1,n2: r1,r2,r3

In the (frequent) case where n1,n2=1,0 this can be abridged with:

q: r1,r2,r3

Here, q is a irreducible or a irreducible power, ri are the corresponding
roots and the contribution that must be subtracted at these positions is
(n1-n2)*degree(p) (assuming smaller powers of this irreducible have
alredy been taken care of).  By position, we mean (a,b) such that a -
b*ri = 0 mod q.

The roots ri must be sorted in lexicographical order.  If a root ri is
greater or equal to q, it means that this is a projective root:
subtracting q gives a root for the reciprocal polynomial (or
equivalently, (1:(ri-q)) is the projective root).

It is allowed to have several lines with the same q, but there must be
only one line for a given (q,n1,n2) triple. 

It is not assumed that the entries are sorted in increasing degree of q.
For small prime and prime powers, this is clearly not the case, since we
might want to have all powers of the same prime at the same place in the
file. However, it might be that for large primes las.c assumes some
sorting. TODO: check this!

Lines starting with a # are comments.
For freerels.c, it is currently mandatory to have one line in the header
of the form:
# DEGREE: 5
(yes, this is a bit strange to rely on a comment...)


Rem:
----
For reference, the old factor base format was based on lines of the form:
p: r1,r2,r3
with no powers and no indication of the multiple of log p to be
subtracted (and no projective roots).
