rootfind(f, n) = {
local(r, l, p, q, e, logn, i, powers, powlen, nextpow);
logn = log (n);

/* powers[i] contains the smallest q^(i+1) > p, q prime. nextpow contains
   the minimum value in powers[] */

print ("# Roots for polynomial ", f);

powlen = floor (logn / log(2));
powers = vector (powlen, i, 2^(i+1));
nextpow = vecmin (powers[1]); /* = 4, but lets do it the generic way */

forprime(p = 2, n, 
  while (p > nextpow, /* Between the last p and this one, 
                         there was a (or another) prime power */
    i = 1;
    while (powers[i] != nextpow, i++);
    q = round (sqrtn (nextpow, i + 1));
    r = List (lift (polrootspadic (f, q, i+1)));

    l = length (r);
    if (l > 0, 
      print1 (nextpow, ": ");
      for (i = 1, l-1, 
        print1 (r[i], ",")
      );
      print (r[l], " # Prime power ", q, "^", i+1);
    );

    q = nextprime (q + 1);
    powers[i] = q^(i + 1);
    nextpow = vecmin (powers);
  );

  r = List (lift (polrootsmod (f, p)));
  l = length (r);
  if (l > 0, 
    print1 (p,": ");
    for (i = 1, l-1, 
      print1 (r[i],",")
    );
    print (r[l]);
  )
)
}
