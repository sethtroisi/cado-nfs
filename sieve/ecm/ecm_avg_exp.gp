ecmtorsion12(n,p) = {
  local (t2, a, A, B);

  /* If the n parameter is a negative rational, use that as the u value */
  if (n < 0,
    u = Mod(-n,p);
  ,
    /* We want a rational u so that u^3-12u is a rational square v^2. 
       u=4, v=4 satisfies this, and we get more such pairs as multiples
       of this point on the curve v^2=u^3-12*u. */
    c = ellinit([Mod(0,p),Mod(0,p),Mod(0,p),Mod(-12,p),Mod(0,p)]);
    u = ellpow (c, [Mod(4,p),Mod(4,p)], n)[1];
  );

  if (u == Mod(0, p), 
    print ("u = ", u);
    return(0);
  );
  t2 = (u^2 - Mod(12,p))/(4*u);
  if (t2 == Mod (-3, p) || t2 == Mod(1,p) || t2 == Mod(-1,p),
     print("t2 = ", t2);
     return(0);
  );
  a = (t2 - Mod(1,p))/(t2 + Mod(3,p));
  if (a == Mod (0, p), return(0));
  A = (-3*a^4 - 6*a^2 + Mod(1,p))/(4*a^3);
  if (A == Mod(2, p) || A == Mod(-2, p), 
    print ("A = ", A);
    return(0);
  );
  B = (a^2 - Mod(1,p))^2/(4*a^3);
  if (B == Mod(0, p), 
    print ("B = ", B);
    return(0);
  );
/* print("u = ", u, ", t^2 = ", t2, ", a = ", a, ", A = ", A, ", B = ", B); */
  E = ellinit ([0, B*A, 0, B^2, 0]);
}

ecmsigma(s, p)={
  local (t, u, v, x, w, a, b, A, E, X);
  
  t = Mod(s, p);
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
  if (A == Mod(2, p) || A == Mod(-2, p), 
/*    print("Skipped singular curve, sigma = ", s,", A = ", A); */
    return(0);
  );

 x = u^3/v^3;
  w = x^3 + A*x^2 + x;
  if (w == 0, return(0));
  /* print ("sigma = ", s, ", p = ",p,", A = ", A, ", w = ", w); */

  /* Now we use the curve  E : wy^2 = x^3 + A*x^2 + x.
     This curve is equivalent to 
     Y^2 = X^3 + w*A*X^2 + w^2 * X
     with X=w*x, Y=w^2*y 
     The trick is that if w is a quadratic residue, then
     E is isomorphic to y^2 = x^3 + A*x^2 + x
     whereas if w is a non-residue, then E is the twist curve instead.
  */

  /* X = w*x; */
  E = ellinit ([0, w*A, 0, w^2, 0]); 

  return(E);
}

primeexp(p, n) = {
  local (r, i);
  r = n;
  i = 0;
  while (r%p == 0, i++; r/=p);

  return(i);
}

ecm_avg_exp(s,pmin,pmax,r,m) = {
  local(n, p2, p3, p5, p7, c, singular);
  n = 0; singular = 0;
  p2 = listcreate(10); for(i=1,10,listput(p2, 0));
  p3 = listcreate(10); for(i=1,10,listput(p3, 0));
  p5 = listcreate(10); for(i=1,10,listput(p5, 0));
  p7 = listcreate(10); for(i=1,10,listput(p7, 0));
  forprime (p = pmin, pmax,
    if (p % m == r,
      c = ecmsigma (s, p);
      if (c != 0,
        n++;
        o = ellsea(c,p); 
        p2[min (primeexp(2,o)+1, 10)]++; 
        p3[min (primeexp(3,o)+1, 10)]++; 
        p5[min (primeexp(5,o)+1, 10)]++; 
        p7[min (primeexp(7,o)+1, 10)]++;
      , singular++);
    );
  );
/*  printp("2: ",precision(1.*p2/n,9),", 3: ",precision(1.*p3/n,9), \
           ", 5: ",precision(1.*p5/n,9),", 7: ",precision(1.*p7/n,9)); */
  return([Vec(p2),Vec(p3),Vec(p5),Vec(p7),Vec([n,singular])]);
}

ecm_avg_exp_t12(s,pmin,pmax,r,m) = {
  local(n, p2, p3, p5, p7, c, singular);
  n = 0; singular = 0;
  p2 = listcreate(10); for(i=1,10,listput(p2, 0));
  p3 = listcreate(10); for(i=1,10,listput(p3, 0));
  p5 = listcreate(10); for(i=1,10,listput(p5, 0));
  p7 = listcreate(10); for(i=1,10,listput(p7, 0));
  forprime (p = pmin, pmax,
    if (p % m == r,
      c = ecmtorsion12 (s, p);
      if (c != 0,
        n++;
        o = ellsea(c,p); 
        p2[min (primeexp(2,o)+1, 10)]++; 
        p3[min (primeexp(3,o)+1, 10)]++; 
        p5[min (primeexp(5,o)+1, 10)]++; 
        p7[min (primeexp(7,o)+1, 10)]++;
      , singular++);
    );
  );
/*  printp("2: ",(1.*p2/n,9),", 3: ",precision(1.*p3/n,9), \
           ", 5: ",precision(1.*p5/n,9),", 7: ",precision(1.*p7/n,9)); */
  return([Vec(p2),Vec(p3),Vec(p5),Vec(p7),Vec([n,singular])]);
}
