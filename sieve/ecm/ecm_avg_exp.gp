/****************************************************************
 *                      General functions                       *
 ****************************************************************/


/* Print factors of n as powers of primes, e.g., n=360 prints 8 9 5 */
printfac(n)={
  local(i, f); 
  f = Vec(factorint(n)); 
  for(i = 1, length(f[1]),
    print1 (f[1][i]^f[2][i], " ")
  ); 
  print();
}


/* Divide vectors elementwise */
vecdiv (a,b) = {
  local(l);
  l = length(a);
  if (length(b) != l, error("Vectors differ in length"));
  vector(l, i, a[i] / b[i])
}


/* Computes lcm(1, 2, 3, ..., B1) */
big_lcm (B1) = {
  local (p, e);
  e = 1;
  forprime (p = 2, B1, 
    e *= p^floor(log (B1) / log (p))
  );
  return (e);
}


/* Enumerates primes p in [pmin, pmax] where p==r (mod m) holds 
   and p+c is B1,B2-smooth. The stage 2 part of the order must be 
   a prime greater than B1.
   Returns [s1, s2] where s1 is the number p+c values that are B1-smooth,
   and s2 is the number of p+c values that are B1,B2 smooth but not B1-smooth.
   If v==1, each prime where p+c is B1,B2-smooth is printed. */
   
list_smooth (B1, B2, pmin, pmax, c, r, m, v) =
{
  local (p, E, q, s1, s2);

  E = big_lcm (B1);
  p = nextprime (pmin);
  s1 = 0;
  s2 = 0;
  if (r < 0, r = m - r);

  while (p <= pmax,
    if (m == 0 || p % m == r,
      q = (p + c) / gcd (E, p + c);
      if (q == 1, 
        s1++;
	if (v == 1, print (p))
      ,
        if (B1 < q && q <= B2 && isprime (q), 
	  s2++;
	  if (v == 1, print (p))
	);
      );
    );
    p = nextprime (p + 1);
  );
  return ([s1, s2]);
}


/****************************************************************
 *                      Functions for P+1                       *
 ****************************************************************/

/* Compute value of i-th Chebyshev polynomial defined by
   V_i(x+1/x) = x^i + 1/x^i */
Chebyshev(i,X) = {
  local(n); 
  if (i == 0, return (2)); 
  if (i == 1, return(X)); 
  if(i % 2 == 0, return (Chebyshev (i/2, X)^2-2));
  if(i % 3 == 0, n = Chebyshev (i/3, X); return (n^3 - 3 * n));
  if(i % 5 == 0, n = Chebyshev (i/5, X); return (n^5 - 5 * n^3 + 5*n));
  if(i % 7 == 0, n = Chebyshev (i/7, X); return (n^7 - 7*n^5 + 14*n^3 - 7*n));
  return (Chebyshev ((i+1)/2, X) * Chebyshev ((i-1)/2, X) - X)
}


/* Determine the order of X in F_p[sqrt(X^2-4)]
   It divides p-1 or p+1, depending on whether X^2-4 is a QR mod p or not */
V_order (X, o) = {
  local(m, d, f, po, i, j, k, p, Xp);
  if (X == 0, return 4);
  m = lift(X) + lift(-X);
  if (o == 0, o=m-kronecker (lift(X)^2-4, m));
  /* if (Chebyshev(o,X) != 2, error("V_order: group order should be ", o, " but isn't")); */
  f = factorint (o)~;
  po = o;
  for (i = 1, length(f), 
    p = f[1,i];
    k = f[2,i];
    po /= p^k; /* All p's removed */
    j = 0;
    Xp = Chebyshev(po, X);
    while (Xp != 2, po *= p; j++; if (j == k, break); Xp = Chebyshev(p, Xp));
  );
  return (po);
}


list_pp1_order (pmin, pmax, res, mod, verbose) =
{
  local (p, o, d);

  p = nextprime (pmin);
  n = 0;

  while (p <= pmax,
    if (mod != 0 && p % mod != res,
      p = nextprime (p + 1);
      next;
    );
    /* Compute the order of alpha where alpha+1/alpha = 2/7 */
    /* o = p - kronecker(lift(Mod(2/7, p))^2-4, p); */
    /* p == 1 (mod 6): o = p-1,  p == 5 (mod 6): p = p+1 */
    o = p + (p % 6 - 3) / 2;
    po = V_order (Mod(2/7, p), o);
    if (o % po != 0, error ("Order of element does not divide order of group"));
    if (verbose != 0, print(p, " ", o, " ", po, " ", o/po));
    p = nextprime (p + 1);
  );
}


/****************************************************************
 *                      Functions for ECM                       *
 ****************************************************************/

/* Generate an elliptic curve E, using the parameterization for
   ratinal torsion 12 from Montgomery's thesis.
   Also generates a starting point.
   Both are returned as an array [E, P] */

ecmtorsion12(n,p,verbose) = {
  local (E, x, y, t, t2, a, A, B, u, v, c, pnt, nmod2);

  /* If the n parameter is a negative rational, use that as the u value */
  if (n < 0,
    u = Mod(-n,p);
  , 
    /* We want a rational u so that v^2 = u^3-12u is a rational square. 
       u=-2, v=4 satisfies this and together with the torsion point (0,0)
       generates all rational points on the curve (which has rank 1).
       Adding the torsion point (0,0) does not seem to affect the
       ECM curve we get out in the end, so we simply ignore it.
       So we get all other (u,v) pairs as multiples of the point (-2,4) 
       on the curve v^2=u^3-12*u. */
    
    if (p != 0,
      c = ellinit([Mod(0,p),Mod(0,p),Mod(0,p),Mod(-12,p),Mod(0,p)]);
      pnt = ellpow (c, [Mod(-2,p), Mod(4,p)], n);
    ,
      c = ellinit([0,0,0,-12,0]);
      pnt = ellpow (c, [-2, 4], n);
    );
    if (verbose, print ("Resulting point on v^2=u^3-12*u : ", pnt));
    u = pnt[1];
    v = pnt[2];
  );

  if (u == 0, 
    print ("u = ", u);
    return(0);
  );

  /* We want t^2 = (u^2-12)/(4u) = v^2/(4*u^2) (so t=v/(2*u))
     and a=(t^2-1)/(t^2+3). */

  t = v/u/2;
  t2 = t^2;
  /* Avoid a = 1/0 */
  if (t2 == -3, return(0));
  a = (t2 - 1)/(t2 + 3);
  /* Avoid a=0 */
  if (a == 0, 
    error ("a == ", a);
  );

  A = (-3*a^4 - 6*a^2 + 1)/(4*a^3);
  if (A == 2 || A == -2, return(0));

  B = (a^2 - 1)^2/(4*a^3);
  if (B == 0, return(0));

  x = (3*a^2+1)/(4*a);
  /* y = sqrt(3*a^2+1)/(4*a), 
     
     or y = sqrt(4*(t^4+3) / (t^2+3)^2) / (4*a)

          = 2 * sqrt(t^4 + 3) / (t^2+3) / (4*a)
 
        with t^2 = (u^2-12)/(4*u) we have
        t^4 + 3 = ((u^2+12)/(4*u))^2
  */
  y = 2 * (u^2+12)/(4*u)  / (t^2+3) / (4*a);

  if (B*y^2 != x^3 + A*x^2 + x, error ("Point is wrong"));

  /* Now we have the curve By^2 = x^3 + Ax^2 + x.
     Transform into a curve of form 
     y^2 = x^3 + a*x^2 + b*x
     by multiplying the equation by B^3 to receive
     B^4 * y^2 = B^3*x^3 + A*B^3*x^2 + B^3*x
     and change of variables (x,y) -> (x'*B, y'*B^2) to get the curve
     y'^2 = x'^3 + A*B*x'^2 + B^2*x' */

  if (verbose, 
    print("u = ", u, ", t^2 = ", t2, ", a = ", a, ", A = ", A, ", B = ", B, 
          ", (x,y) = (", x, ",", y, ") on By^2 = x^3 + Ax^2 + x");
  );

  E = ellinit ([0, B*A, 0, B^2, 0]);
  x *= B;
  y *= B^2;
  
  if (y^2 != x^3 + A*B*x^2 + B^2*x, error ("New point is wrong"));

  if (!ellisoncurve(E,[x,y]), error ("Point is not on curve"));

  return ([E, [x, y]]);
}


ecmtorsion16(n,p,verbose) = {
  local (E, A, B, x, y);

  if (n != 1, error ("Only one curve with rational torsion 16 implemented so far"));

  if (p == 0, A = 1, A = Mod (1, p));
  x = A;
  y = A;
  A *= 54721/14400;
  x *= 8/15;
  B = x^3 + A*x^2 + x;
  
  if (B*y^2 != x^3 + A*x^2 + x, error ("Point is wrong"));
  
  /* Now we have the curve By^2 = x^3 + Ax^2 + x.
     Transform into a curve of form 
     y^2 = x^3 + a*x^2 + b*x
     by multiplying the equation by B^3 to receive
     B^4 * y^2 = B^3*x^3 + A*B^3*x^2 + B^3*x
     and change of variables (x,y) -> (x'*B, y'*B^2) to get the curve
     y'^2 = x'^3 + A*B*x'^2 + B^2*x' */

  if (verbose, 
    print("u = ", u, ", t^2 = ", t2, ", a = ", a, ", A = ", A, ", B = ", B, 
          ", (x,y) = (", x, ",", y, ") on By^2 = x^3 + Ax^2 + x");
  );

  E = ellinit ([0, B*A, 0, B^2, 0]);
  x *= B;
  y *= B^2;
  
  if (y^2 != x^3 + A*B*x^2 + B^2*x, error ("New point is wrong"));

  if (!ellisoncurve(E,[x,y]), error ("Point is not on curve"));

  return ([E, [x, y]]);
}


/* Returns a curve and a point on the curve, generated by the 
   Brent-Suyama parameterization with order divisible by 12 */

ecmsigma(s, p)={
  local (t, u, v, x, w, a, b, A, E, X);
  
  if (p != 0,
    t = Mod(s, p);
  ,
    t = s;
  );
  v = 4*t;
  u = t^2 - 5;
  a = (v-u)^3*(3*u+v);
  b = 4*u^3*v;
  if (gcd (b, p) != 1, return(0));
  A = a/b-2;
/*
  if (poldegree(gcd(3*z^2 + 2*A*z + 1, z^3 + A*z^2 + z)) > 0, 
    print("Skipped singular curve sigma = ",s,", A = ",A);
    return(0);
  );
*/
  if (A == 2 || A == -2, 
/*    print("Skipped singular curve, sigma = ", s,", A = ", A); */
    return(0);
  );

  x = u^3/v^3;
  w = x^3 + A*x^2 + x;
  if (w == 0, return(0));
  /* print ("sigma = ", s, ", p = ",p,", A = ", A, ", w = ", w); */

  /* Now we use the curve  E : wy^2 = x^3 + A*x^2 + x 
     which has the point (x, 1).
     This curve is equivalent to 
     Y^2 = X^3 + w*A*X^2 + w^2 * X
     with X=w*x, Y=w^2*y = w^2 (because we made y=1)
     The trick is that if w is a quadratic residue, then
     E is isomorphic to y^2 = x^3 + A*x^2 + x
     whereas if w is a non-residue, then E is the twist curve instead.
  */

  /* X = w*x; */
  E = ellinit ([0, w*A, 0, w^2, 0]); 

  return([E, [x * w, w^2]]);
}


/* Prints order of curve, of starting point, and their quotient
   parameterization = 0 is Brent-Suyama,
   parameterization = 1 is Montgomery torsion 12 */
ecm_printorder (parameterization, parameter, p) =
{
  local (E, P, o, po);
  if (parameterization != 0 && parameterization != 1 
      && parameterization != 2, error ("Invalid parameterization"));
  if (parameterization == 0, EP = ecmsigma (parameter, p));
  if (parameterization == 1, EP = ecmtorsion12 (parameter, p));
  if (parameterization == 2, EP = ecmtorsion16 (parameter, p));
  if (EP == 0, print (p, " 0 0 0"); return);
  E = EP[1];
  P = EP[2];
  o = p + 1 - ellap (E);
  po = ellorder (E, P, o);
  print (p, " ", o, " ", po, " ", o / po);
}


/* Count those prime in [pmin, pmax] where the starting point
   on the elliptic curve with Brent-Suyama (parameterization=0)
   or Montgomery torsion 12 (parameterization=1)
   parameterization has order which is B1,B2-smooth. 
   If v=1, print each such prime */

list_smooth_ellorder (B1, B2, pmin, pmax, parameterization, parameter, d, v) =
{
  local(E, C, P, o, p, q, s1, s2, s3);

  if (parameterization != 0 && parameterization != 1,
    error ("Invalid parameterization"));

  E = big_lcm (B1);
  p = nextprime (pmin);
  s1 = 0;
  s2 = 0;
  while (p <= pmax,
    if (parameterization == 0,
      EP = ecmsigma (parameter, p);
    ,
      EP = ecmtorsion12 (parameter, p);
    );
    if (EP == 0, next);
    E = EP[1];
    P = EP[2];
    o = ellorder (E, P);
    q = o / gcd (E, o);
      if (q == 1, 
        s1++;
        if (v == 1, print (p));
      );
      if (B1 < q && q <= B2 && isprime(q),
        s2++;
        if (v == 1, print (p));
      );
      
      /* These are the points precomputed for stage 2. If any of them has
         z-coordinate 0, the factor is found during normalization */
      if (q > 1 && d > 0 && d % q == 0,
        s3++;
        /* print ("q|d: p = ", p, ", q = ", q); */
        if (v == 1, print (p));
      );
      if (1 < q && q <= B1 && q < d / 2 && gcd(q, d) == 1,
        s3++;
        /* print ("p = ", p, ", q = ", q); */
        if (v == 1, print (p));
      );
    p = nextprime (p + 1);
  );
  return ([s1, s2, s3]);
}



/*
   For each prime p in [pmin, pmax] and p == r (mod m),
   and each curve modulo p with parameter in the set s 
   (Brent or Montgomery parameterization with order divisible by 12),
   compute the number of times that valuation(q, o) = i
   for q=2,3,5,7,11, 0 <= i < 15, where o is the order
   of the curve 
*/

ecm_dist_exp(parameterization,s,pmin,pmax,r,m) = {
  local(n, p2, p3, p5, p7, p11, EP, singular, len, i);
  n = 0; singular = 0; len = 15; svec = Vec(s);
  p2 = listcreate(len); for(i=1,len,listput(p2, 0));
  p3 = listcreate(len); for(i=1,len,listput(p3, 0));
  p5 = listcreate(len); for(i=1,len,listput(p5, 0));
  p7 = listcreate(len); for(i=1,len,listput(p7, 0));
  p11 = listcreate(len); for(i=1,len,listput(p11, 0));
  for (i = 1, length(svec),  
    forprime (p = pmin, pmax,
      if (m == 0 || p % m == r,
        if (parameterization == 0,
          EP = ecmsigma (svec[i], p);
        ,
          EP = ecmtorsion12 (svec[i], p);
        );
        if (EP != 0,
          n++;
          o = p + 1 - ellap(EP[1]); 
          p2[min (valuation(o,2)+1, len)]++; 
          p3[min (valuation(o,3)+1, len)]++; 
          p5[min (valuation(o,5)+1, len)]++; 
          p7[min (valuation(o,7)+1, len)]++;
          p11[min (valuation(o,11)+1, len)]++;
        , singular++);
      );
    );
  );
/*  printp("2: ",precision(1.*p2/n,9),", 3: ",precision(1.*p3/n,9), \
           ", 5: ",precision(1.*p5/n,9),", 7: ",precision(1.*p7/n,9)); */
  return([Vec(p2),Vec(p3),Vec(p5),Vec(p7),Vec(p11), Vec([n,singular])]);
}


ecm_avg_exp(d) = {
  local(n, s, i, j, avg, l);
  l = length(d);
  n = d[l][1]; /* The total number of good curves we tried */
  avg = listcreate(l - 1);
  for(i = 1, l - 1,
    s = 0;
    for (j = 1, length(d[i]),
      s += (j - 1) * d[i][j];
    );
    listput (avg, s / n);
  );
  return(Vec(avg));
}

/* Computes the order of the curve, the order of the starting point and
   the index of the starting point for curves with Montgomery torsion 12 
n=0;v=vector(20);forprime(p=100000, 200000, EP=ecmtorsion12(2, p); E=EP[1]; P=EP[2]; o=p + 1 - ellap(E); po = ellorder(E,P,o); n++; v[valuation(o/po,2)+1]++; print(p, ": ", o, " / ", po, " = ", o/po));
*/
