rootfind(f, n) = {
local(r, l, p, q, e, logn, i, powers, powlen, nextpow);
logn = log (n);

/* powers[i] contains the smallest q^(i+1) > p, q prime. nextpow contains
   the minimum value in powers[] */

print ("# Roots for polynomial ", f);
print ("# DEGREE: ", poldegree(f));

powlen = floor (logn / log(2));
powers = vector (powlen, i, 2^(i+1));
nextpow = vecmin (powers[1]); /* = 4, but lets do it the generic way */

forprime(p = 2, n, 
  while (0 && p > nextpow, /* Between the last p and this one, 
                         there was a (or another) prime power */
    /* Roots mod p^k with polrootspadic simply does not work the way I
       thought it would. Currently disabled. */
    i = 1;
    while (powers[i] != nextpow, i++);
    q = round (sqrtn (nextpow, i + 1));
    r = lift (polrootspadic (f, q, i+1)~);
    l = length (List(r));

    /* If p divides leading coefficient, we may get strange extra roots
       which are rational. Discard them */
    while (l > 0 && type(r[1]) == "t_FRAC",
      /* print ("# Discarding ", r[1], " from ", r); */
      r = vecextract (r, "2..");
      l--;
    );

    if (l > 0, 
      print1 (nextpow, ": ");
      for (i = 1, l-1, 
        print1 (r[i] % (q^(i+1)), ",")
      );
      print (r[l] % (q^(i+1)), " # Prime power ", q, "^", i+1);
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
