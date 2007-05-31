f:=191024405174092*x^4-21845598385191*x^3-116985015868315*x^2-3606923667439*x+
167656188134287:

get_alpha := proc(f, B) local s, p, lc, e, q;
   s := 0;
   p := 2;
   lc := lcoeff(f,x);
   while p <= B do
      if lc mod p <> 0 then
         q := 0;
         for e in Roots(f) mod p do
            q:=q+e[2];
         od;
         s := s + evalf((1-q*p/(p+1))*log(p)/(p-1))
      fi;
      p := nextprime(p);
   od;
   s;
end:
