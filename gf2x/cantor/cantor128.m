/* Magma code that checks the cantor_main output: */

// w := 64;

load "/tmp/toto";
PP := PolynomialRing(GF(2));
F := PP!Intseq(Seqint(f, 2^w), 2);
G := PP!Intseq(Seqint(g, 2^w), 2);
FG := PP!Intseq(Seqint(fg, 2^w), 2);

FG eq F*G;
if FG ne F*G then
    exit;
end if;


// w := 64;
load "/tmp/toto";
PP := PolynomialRing(GF(2));
x := PP.1;
F128<z> := ext<GF(2) | x^128 + x^7 + x^2 + x + 1>;
P128<x> := PolynomialRing(F128);

fi:=f;
gi:=g;
fgi:=fg;

function convertF128(F128, x)
  return F128!PP!Intseq(x, 2);
end function;

function convertP128(F128, f)
  return Polynomial([ convertF128(F128, z): z in Intseq(Seqint(f, 2^w), 2^64)]);
end function;

function reglue(ffi)
    return &+[ PP!Eltseq(Coefficient(ffi,i))*x^(64*i) : i in [0..Degree(ffi)]];
end function;

ffi := convertP128(F128,fi);
ggi := convertP128(F128,gi);

reglue(ffi) eq F;
reglue(ggi) eq G;
reglue(ffi*ggi) eq F*G;

h := convertP128(F128,fgi);
reglue(h) eq F*G;  




beta_0 := F128!1;
beta_1 := F128!Intseq((7451958335100816136 +2^64*2979905974001933049 ), 2);
beta_2 := F128!Intseq((18379359142562139938+2^64*11678753174964950217), 2);
beta_3 := F128!Intseq((6032672185962850376 +2^64*8744256601846606146 ), 2);
beta_4 := F128!Intseq((14608789766813931076+2^64*645900494672607458  ), 2);
beta_5 := F128!Intseq((5568895992531017964 +2^64*5316835906050072984 ), 2);
beta_6 := F128!Intseq((11619261390503595532+2^64*377604546732988956  ), 2);
beta_7 := F128!Intseq((6679137017075335448 +2^64*5571574281931689094 ), 2);
beta_8 := F128!Intseq(102777404866455936327711046477332611352, 2);
beta:=[beta_0,beta_1,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7];

function next_beta(b)
    x:=PolynomialRing(Parent(b)).1;
    return Roots(x^2+x+b)[1][1];
end function;

// can be used to print omega values.
function F128_to_int(b)
    Z:=IntegerRing();
    return Seqint(ChangeUniverse(Eltseq(b),Z),2);
end function;

s1:=x^2+x;
s2:=Evaluate(s1,s1);
s3:=Evaluate(s2,s1);
s4:=Evaluate(s2,s2);
s5:=Evaluate(s1,s4);
s6:=Evaluate(s3,s3);
s7:=Evaluate(s3,s4);
spoly:=[s1,s2,s3,s4,s5,s6,s7];

// throw in some more, while we're at it.
Append(~beta,next_beta(beta[#beta]));
Append(~spoly,Evaluate(s1,spoly[#spoly]));
Append(~beta,next_beta(beta[#beta]));
Append(~spoly,Evaluate(s1,spoly[#spoly]));

&and [Evaluate(spoly[i],beta[i+1]) eq 1:i in [1..#spoly]];             
&and [Evaluate(s1,beta[i+1]) eq beta[i]:i in [1..#spoly]];


function one_omega(i)
    n:=Intseq(i,2);
    if #n eq 0 then return F128!0; end if;
    return &+ [beta[j]*n[j]:j in [1..#n]];
end function;

&and[spoly[j] eq PP!&*[x-one_omega(i):i in [0..2^j-1]]:j in [1..#spoly]];


///////////////////////////////////////////////////////////////////////////


function split(fi)
    return [Intseq(Seqint(ChangeUniverse(Eltseq(Coefficient(fi,k)),Integers()),2),2^w) : k in [0..Degree(fi)] ];
end function;

///////////////////////////////////////////////////////////////////////////

procedure quotremSi(i, pf, ii)
  printf "QuotremSi(%o, %o, %o)\n", i, pf, ii;
end procedure;

procedure mulSi(i, pf, ii)
  printf "mulSi(%o, %o, %o)\n", i, pf, ii;
end procedure;



/*
procedure reduceModTrunc(k, length)
  pf := 0;
  len := length;

  // go down, computing quotrems
  printf "Going down...\n";
  ii := 2^k;
  for i := k-1 to 0 by -1 do
    if (len ge (2^i)) then
      quotremSi(i, pf, ii);
      ii -:= 2^i;
      pf +:= 2^i;
      len -:= 2^i;
    end if;
  end for;

  printf "Going up...\n";
  len := length;
  i := 0;
  while IsEven(len) do
    len div:= 2;
    i +:= 1;
  end while;
  pf := length;
  ii := 2^i;

  i +:= 1;
  len div:= 2;
  while i le k-1 do
    if ( (len mod 2) eq 1) then
      mulSi(i, pf - ii, ii);
      ii +:= 2^i;
    end if;
    len div:= 2;
    i +:= 1;
  end while;
end procedure;
*/

ind_S:=[
    [1],                        // S1
    [1],                        // S2
    [4, 2, 1],                  // S3
    [1],                        // S4
    [16, 2, 1],                 // S5
    [16, 4, 1],                 // S6
    [64, 32, 16, 8, 4, 2, 1],   // S7
    [1],                        // S8
    [256, 2, 1],                // S9
    [256, 4, 1],                // S10
    [1024, 512, 256, 8, 4, 2, 1],       // S11
    [256, 16, 1],               // S12
    [4096, 512, 256, 32, 16, 2, 1],     // S13
    [4096, 1024, 256, 64, 16, 4, 1],    // S14
    [16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1],     // S15
    [1],                        // S16
    [65536, 2, 1],              // S17
    [65536, 4, 1],              // S18
    [262144, 131072, 65536, 8, 4, 2, 1],        // S19
    [65536, 16, 1],             // S20
    [1048576, 131072, 65536, 32, 16, 2, 1],     // S21
    [1048576, 262144, 65536, 64, 16, 4, 1],     // S22
    [4194304, 2097152, 1048576, 524288, 262144, 131072, 65536, 128, 64, 32, 16, 8, 4, 2, 1],    // S23
    [65536, 256, 1],            // S24
    [16777216, 131072, 65536, 512, 256, 2, 1],  // S25
    [16777216, 262144, 65536, 1024, 256, 4, 1], // S26
    [67108864, 33554432, 16777216, 524288, 262144, 131072, 65536, 2048, 1024, 512, 256, 8, 4, 2, 1],    // S27
    [16777216, 1048576, 65536, 4096, 256, 16, 1],       // S28
    [268435456, 33554432, 16777216, 2097152, 1048576, 131072, 65536, 8192, 4096, 512, 256, 32, 16, 2, 1],       // S29
    [268435456, 67108864, 16777216, 4194304, 1048576, 262144, 65536, 16384, 4096, 1024, 256, 64, 16, 4, 1]     // S30
];

// reduce mod Sk + omega and Sk + omega + 1
procedure reduceSi_inplace(k, ~f, omega)
    K:=2^k;
    assert #f eq 2*K;
    a:=P128!f;

    // put rest mod Sk+omega in low part of f and quo in high part.
    for i in [2*K-1..K by -1] do
        fi:=f[1+i];
        for j in ind_S[k] do
            f[1+j+i-K] +:= fi;
        end for;
        f[1+i-K] +:= fi * omega;
    end for;

    print P128!(f[K+1..2*K]) eq a div (spoly[k]+omega);
    print P128!(f[1..K]) eq a mod (spoly[k]+omega);

    for i in [0..K-1] do
        f[1+i+K] +:= f[1+i];
    end for;

    print P128!(f[K+1..2*K]) eq a mod (spoly[k]+omega+1);

end procedure;

// reduce mod Sk + omega and Sk + omega + 1
function reduceSi(k,a,omega)
    K:=2^k;
    f:=[F128!0:i in [1..2*K]];
    for i in [0..Degree(a)] do f[1+i]:=Coefficient(a,i); end for;
    for i in [2*K-1..K by -1] do
        fi:=f[1+i];
        if k gt 0 then
            for j in ind_S[k] do
                f[1+j+i-K] +:= fi;
            end for;
        end if;
        f[1+i-K] +:= fi * omega;
    end for;
    for i in [0..K-1] do
        f[1+i+K] +:= f[1+i];
    end for;
    h:=P128!(f[K+1..2*K]);
    if k gt 0 then assert h eq a mod (spoly[k]+omega+1); end if;
    l:=P128!(f[1..K]);
    if k gt 0 then assert l eq a mod (spoly[k]+omega); end if;
    return [l,h];
end function;

function multieval(k,f,i_omega)
    // our input is defined modulo s_k + omega_{i_omega}
    assert k ge 1;
    // the end case is where we have output defined modulo s_1+b, for some
    // omega value b. We have to compute f mod x+c and x+c+1, where c is the
    // root of x^2+x+b, hence of index 2*i_omega
    i_omega *:= 2;
    children := reduceSi(k-1, f, one_omega(i_omega));
    if k eq 1 then
        return [Coefficient(children[1],0),Coefficient(children[2],0)];
    else
        // then we recurse, since the children are not leaves.
        return
            multieval(k-1,children[1],i_omega)
            cat multieval(k-1,children[2],i_omega+1);
    end if;
end function;

// f has 2^k coefficients. Returns chunks of 2^kappa coefficients,
// corresponding to reductions mod x^(2^kappa)-x-omega, for all the
// 2^(k-kappa) relevant values of omega.
//
// The ``truncation'' by n actually computes only the first ceiling(n/2^kappa)
// terms of the sequence defined above. n0 is initially zero, but recursive
// calls update it to indicate the number of terms computed on the left side
// of the tree. If n==2^k, it's a full computation.
function reduce_top_cantor(k,kappa,f,i_omega,n,n0)
    assert Degree(f)+1 le 2^k;
    if n0 ge n then
        return [Parent(f)|];
    elif k eq kappa then
        return [f];
    else
        i_omega *:= 2;
        children := reduceSi(k-1, f, one_omega(i_omega));
        return
            reduce_top_cantor(k-1,kappa,children[1],i_omega,n,n0)
            cat
            reduce_top_cantor(k-1,kappa,children[2],i_omega+1,n,n0+2^(k-1));
    end if;
end function;

// must hold.
print multieval(7,ffi,0) eq [Evaluate(ffi,one_omega(i)):i in [0..127]]; 

function expand(a,t,t0)
    // write a, having at most 2^t terms, in base x^(2^t0)+x
    // a has at most 2^t coefficients. Split in pieces of 2^t0 coefficients.
    if t eq t0 then return [a]; end if;
    x:=Parent(a).1;
    K:=2^(t-1);
    K0:=2^(t-1-t0);
    // split a modulo x^(2^(t-1))-x^(2^(t-1-t0))
    f:=[CoefficientRing(Parent(a))|0:i in [1..2*K]];
    for i in [0..Degree(a)] do f[1+i]:=Coefficient(a,i); end for;
    // compute quotient.
    for i in [2*K-1..K by -1] do
        f[1+i+K0-K] +:= f[1+i];
    end for;
    s:=x^K-x^K0;
    h:=P128!(f[K+1..2*K]);
    l:=P128!(f[1..K]);
    assert h eq a div s;
    assert l eq a mod s;
    return expand(l,t-1,t0) cat expand(h,t-1,t0);
end function;


function expand_trunc(a,t,t0,n)
    K:=2^(t-1);
    if t eq t0 then return [a]; end if;
    if n le K then
        return expand_trunc(a,t-1,t0,n);
    end if;
    x:=Parent(a).1;
    K0:=2^(t-1-t0);
    f:=[CoefficientRing(Parent(a))|0:i in [1..n]];
    for i in [0..Degree(a)] do f[1+i]:=Coefficient(a,i); end for;
    for i in [n-1..K by -1] do
        f[1+i+K0-K] +:= f[1+i];
    end for;
    h:=P128!(f[K+1..n]);
    l:=P128!(f[1..K]);
    return expand_trunc(l,t-1,t0,K) cat expand_trunc(h,t-1,t0,n-K);
end function;


function gm_trick(two_t,f,j)
    // Compute f on omega_{j + i}, i in [0..2^two_t-1]
    assert IsPowerOf(two_t,2);
    eta:=2^two_t;
    // f has eta coefficients.
    if eta eq 2 then
        omega := one_omega(2*j);
        f0:=Coefficient(f,0);
        f1:=Coefficient(f,1);
        z := f1 * omega + f0;
        return [ z, z + f1 ];
    end if;
    t:=two_t div 2;
    tau:=2^t;
    // we'll first write f in base x^(2^t)+x ; that makes 2^t terms.
    g:=expand(f,two_t,t);
    T:=[0..tau-1];
    h:=[Parent(f)![Coefficient(g[1+i],l):i in T]:l in T];
    evals_h:=[gm_trick(t, h[1+l], j): l in [0..tau-1]];
    f_mod_xtau_x_omega_phi:=[Parent(f)![evals_h[1+i][1+o]:i in T]:o in T];
    // f_mod_xtau_x_omega_phi eq [f mod (x^tau-x-one_omega(i)):i in T];
    return &cat[gm_trick(t,f_mod_xtau_x_omega_phi[1+i],j*tau+i):i in
  T];
    // gm_trick(two_t,f,j) eq [Evaluate(f,one_omega(i)):i in [0..2^(two_t)-1]];
end function;

/*
function gm_trick_trunc(two_t,f,j,n)
    // Compute f on omega_{j * eta + i}, i in [0..n-1]
    // n must be at most 2^two_t
    assert IsPowerOf(two_t,2);
    eta:=2^two_t;
    assert n le eta;
    // f has eta coefficients.
    if n eq 0 then return [ CoefficientRing(f) | ]; end if;
    if eta eq 2 then
        omega := one_omega(2*j);
        f0:=Coefficient(f,0);
        f1:=Coefficient(f,1);
        z := f1 * omega + f0;
        if n eq 1 then return [z]; end if;
        return [ z, z + f1 ];
    end if;
    t:=two_t div 2;
    tau:=2^t;
    // we'll first write f in base x^(2^t)+x ; that makes 2^t terms.
    g:=expand_trunc(f,two_t,t,n);
    T:=[0..tau-1];
    h:=[Parent(f)![Coefficient(g[1+i],l):i in [0..#g-1]]:l in T];
    // h_lambda has Ceiling(n/2^t) coeffs at most. More exactly, it's
    // n div 2^t + (lambda gt (n mod 2^t) select 1 else 0);
    // &and [1+Degree(h[lambda+1]) eq n div 2^t + (lambda lt (n mod 2^t)
    // select 1 else 0) : lambda in T];
    ncoeffs_h:=func<l|n div 2^t + (l lt (n mod 2^t) select 1 else 0)>;
    // here we know that the evaluation is never trivial, but we forcibly
    // compute only the values we're interested in.
    evals_h:=[gm_trick_trunc(t, h[1+l], j, ncoeffs_h(l)): l in T];
    
    // and then here we fail. The f's make no sense if we don't have all the
    // values.
    f_mod_xtau_x_omega_phi:=[Parent(f)![evals_h[1+i][1+o]:i in T]:o in T];
    // f_mod_xtau_x_omega_phi eq [f mod (x^tau-x-one_omega(i)):i in T];
    return &cat[gm_trick_trunc(t,f_mod_xtau_x_omega_phi[1+i],j*tau+i):i in
  T];
    // gm_trick_trunc(two_t,f,j) eq [Evaluate(f,one_omega(i)):i in [0..2^(two_t)-1]];
end function;
*/

// children:=reduce_top_cantor(7,4,ffi,0);f:=children[1][1];j:=0;two_t:=4;

// n:=5000;f:=Parent(f)![Random(F128):i in [0..n-1]];j:=0;two_t:=16;


// f has n coeffs, which is at most 2^k.
function multieval_toplevel(f)
    n:=1+Degree(f);
    k:=Ilog2(n)+1;
    assert 1 + Degree(f) le 2^k;
    r:=Ilog2(k);
    kappa:=2^r;
    children:=reduce_top_cantor(k,kappa,f,0,2^k,0);
    return &cat [ gm_trick(kappa,children[i],i):i in [1..#children]];
end function;

// f:=ffi;k:=7;multieval_toplevel(k,f) eq [Evaluate(f,one_omega(i)):i in [0..2^k-1]];
