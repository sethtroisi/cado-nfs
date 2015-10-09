/* This check code is vastly incomplete, and full of old clutter */
/* commit bf225ce125bc21c5ae5f4c3b2f8acc652e248f46 reportedly had more working
 * stuff, but the code has evolved somewhat since then.
 */

v2i:=func<r|Seqint(ChangeUniverse(Eltseq(r),Integers()),2)>;
m2i:=func<m|[v2i(r):r in Rows(m)]>;
i2m:=func<var|Matrix(GF(2),[Intseq(x,2,64):x in var])>;
load "mn.m";
load "t.m"; M:=Matrix(GF(2),Matrix(var));
nr:=Nrows(M);
nc:=Ncols(M);

nr_orig:=nr;
nc_orig:=nc;

assert nc le nr;
nc:=nr;

x:=Matrix(GF(2),nr,nc,[]);InsertBlock(~x,M,1,1);M:=x;
load "placemats.m";

load "rw.m"; rw:=var;
load "cw.m"; cw:=var;

load "b.m"; 

Cr:=Matrix(GF(2),nr,nrpad,[]); InsertBlock(~Cr,IdentityMatrix(GF(2),nr),1,1);
Cc:=Matrix(GF(2),nc,ncpad,[]); InsertBlock(~Cc,IdentityMatrix(GF(2),nc),1,1);

assert IsOne(Cr*Transpose(Cr));
assert IsOne(Cc*Transpose(Cc));

assert ncpad eq nv*ncp;
assert nrpad eq nh*nrp;
nz:=nrp div nv;
assert nz eq ncp div nh;

pr:=func<x|(qr*nh+qq)*nz+r+1 where qq,qr is Quotrem(q, nv) where q,r is Quotrem(x-1,nz)>;
prinv:=func<x|(qr*nv+qq)*nz+r+1 where qq,qr is Quotrem(q, nh) where q,r is Quotrem(x-1,nz)>;
prp:=SymmetricGroup(nrpad)![pr(x):x in [1..nrpad]];
Pr:=PermutationMatrix(GF(2),[pr(x):x in [1..nrpad]]);
assert Pr eq PermutationMatrix(GF(2),prp);



sc:=SymmetricGroup(ncpad)!colperm;
if assigned rowperm then
    sr:=SymmetricGroup(nrpad)!rowperm;
else
    sr:=SymmetricGroup(nrpad)!colperm;
end if;
sr_inv:=sr^-1;
sc_inv:=sc^-1;
Sr:=PermutationMatrix(GF(2),sr);
Sc:=PermutationMatrix(GF(2),sc);
if assigned rowperm then
    Sr:=Pr^-1*Sr;
    sr:=prp^-1*sr;
end if;

Mx:=M;

Mx[2] eq (Sr*Mx)[2^(sr^-1)];
Transpose(Mx)[2] eq (Sc*Transpose(Mx))[2^(sc^-1)];

S:=Sr;
P:=Pr;

// de-correlation permutation:
q:=SymmetricGroup(ncpad)![preshuf(x):x in [1..ncpad]];
Q:=PermutationMatrix(GF(2),q);
q0:=SymmetricGroup(nc)![preshuf(x):x in [1..nc]];
Q0:=PermutationMatrix(GF(2),q0);
Transpose(Pr*Sr)*Mt*Sc eq Mx*Q;
// t means twisted.

/* Note that Mt being equal to Pr times Sr*M*Sc^-1 is a
 * direct consequence of the fact that balancing_workhorse.c permutes the
 * rows which are being sent.
 */
if Mt eq P*S*Mx*Q*S^-1 then
    print "Shuffled product detected";
elif  Mt eq S*Mx*Q*S^-1 then
    print "non-shuffled product detected";
    assert false;
    // this has been deprecated.
    Pr:=Identity(Parent(Pr));
else
    assert false;
end if;
assert Transpose(Pr*Sr)*Mt*Sc eq Mx*Q;
assert Mt eq Pr*Sr*Mx*Q*Sc^-1;


colweights:=func<M|[#[i:i in [1..Nrows(M)]|M[i,j] ne 0]: j in [1..Ncols(M)]]>;
rowweights:=func<M|[#[j:j in [1..Ncols(M)]|M[i,j] ne 0]: i in [1..Nrows(M)]]>;
// colsums:=func<M|[&+[Integers()!M[i,j]:i in [1..Nrows(M)]]: j in [1..Ncols(M)]]>;
// rowsums:=func<M|[&+[Integers()!M[i,j]:j in [1..Ncols(M)]]: i in [1..Nrows(M)]]>;

assert colweights(Mt*Sc) eq colweights(Mx*Q);
assert rowweights(Transpose(Pr*Sr)*Mt) eq rowweights(M);



if nullspace eq "left" then
    nV:=nr;
elif nullspace eq "right" then
    nV:=nc;
else
    assert false;
end if;


load "x.m";
C0:=Matrix(GF(2),nV,64,[]);
Xt:=Matrix(GF(2),nV,n,[]);
for i in [1..64] do for u in var[i] do C0[u][i] +:=1; end for; end for;
for i in [1..n] do for u in var[i] do Xt[u][i] +:=1; end for; end for;


function vblock(nr,var, nc)
    assert #var le nr;
    R:=Matrix(GF(2),nr,nc,[]);
    nb:=nc div 64;
    for i in [0..nr-1] do
        R[i+1]:=Vector(GF(2),Intseq(u,2,nc) where u is Seqint(var[1+i*nb..(i+1)*nb],2^64));
    end for;
  return R;
end function;


load "Y.0.m";      Y0:= vblock(nV, var, n);
load "V0-64.0.m";  V0:= vblock(nV, var, 64);
load "V0-64.1.m";  V1:= vblock(nV, var, 64);
assert Transpose(V0)*Cr*Mx*Q eq Transpose(V1)*Cr;

// load "V0-64.2.m";  V2:= vblock(nV, var, 64);
// load "V0-64.3.m";  V3:= vblock(nV, var, 64);
// load "V0-64.4.m";  V4:= vblock(nV, var, 64);
// load "V0-64.5.m";  V5:= vblock(nV, var, 64);
// load "V0-64.6.m";  V6:= vblock(nV, var, 64);
// load "V0-64.7.m";  V7:= vblock(nV, var, 64);
// load "V0-64.8.m";  V8:= vblock(nV, var, 64);
// load "V0-64.9.m";  V9:= vblock(nV, var, 64);
load "V0-64.10.m"; V10:=vblock(nV, var, 64);
// load "V0-64.20.m"; V20:=vblock(nV, var, 64);
//

z:=Transpose(V0)*Cr;
for i in [1..10] do z*:=Mx*Q; end for;
assert z eq Transpose(V10)*Cr;


// load "V64-128.0.m";  V0b:= vblock(nV, var, 64);
// load "V64-128.10.m"; V10b:=vblock(nV, var, 64);
// load "V64-128.20.m"; V20b:=vblock(nV, var, 64);


// load "V0-64.12.m"; V12:=vblock(nc, var, n);
// load "V0-64.13.m"; V13:=vblock(nc, var, n);
// load "V0-64.14.m"; V14:=vblock(nc, var, n);
// load "V0-64.15.m"; V15:=vblock(nc, var, n);
// load "V0-64.16.m"; V16:=vblock(nc, var, n);
// load "V0-64.17.m"; V17:=vblock(nc, var, n);
// load "V0-64.18.m"; V18:=vblock(nc, var, n);
load "C.0.m";    C0:=vblock(nV, var, 64);
load "C.10.m";   C10:=vblock(nV, var, 64);

load "H1.m";    H1:=vblock(nr, var, 64);
var:=[(p div j+q*j) mod 2^64:j in [1..nc]]
    where p is PreviousPrime(2^64)
    where q is PreviousPrime(2^63);
H0:=vblock(nc, var, 64);
assert Cr*M*Q*Transpose(Cc)*H0 eq H1;

C0 eq Xt or C0 eq Submatrix(Xt,1,1,Nrows(Xt),64);
C1 eq Cr*Mx*Q*Transpose(Cr)*Xt or C1 eq Submatrix(Cr*Mx*Q*Transpose(Cr)*Xt,1,1,Nrows(Xt),64);

load "A0-64.0-18.m";
A:=vblock(18*64,var,64);
Transpose(Xt)*V0 eq Submatrix(A,1,1,64,64);

// v0z:=[Seqint(ChangeUniverse((Eltseq(V0[i])),Integers()),2):i in [1..nr]];
// v0bz:=[Seqint(ChangeUniverse((Eltseq(V0b[i])),Integers()),2):i in [1..nr]];

// it's somewhat misleading sometimes.
// [Index(Rows(Mt*V0),r):r in Rows(V1)];

Z:=Integers();

/*
Sort(rowweights(Matrix(Z,Transpose(Mt)*V0))) eq 
Sort(rowweights(Matrix(Z,V1)));

([w[pr(i)]: i in [1..#w]] where w is rowweights(Matrix(Z,Transpose(Mt)*V0)))\
 eq rowweights(Matrix(Z,V1));

rowweights(Matrix(Z,Transpose(Mt)*V0)) eq
[w[prinv(i)]: i in [1..#w]] where w is rowweights(Matrix(Z,V1));

rowweights(Matrix(Z,Pr*Transpose(Mt)*V0)) eq rowweights(Matrix(Z,V1));;
*/

TMt:=xtr(Pr^-1)*xtr(Mt);

TM:=xtr(M*Q0);

// equivalent to TMt:=Transpose(Mt*Pr^-1) for nullspace=left
// For nullspace=left:
// TMt eq Transpose(Mt*Pr^-1);
// TMt eq Transpose(Pr*Sr*Cr*M*Transpose(Cc)*Sc^-1*Pr^-1);
// TMt eq Transpose(Pr*Sr*Mx*Sc^-1*Pr^-1);
// TMt eq Transpose(Pr*Sc)^-1*Transpose(Mx)*Transpose(Pr*Sr);
// assert Sr eq Sc
// Transpose(Mx)*(Pr*Sr)^-1*V0 eq (Pr*Sr)^-1*V1;
//
// For nullspace=right
// TMt eq Pr^-1*Mt
// TMt eq Sr*Cr*M*Transpose(Cc)*Sc^-1
// TMt eq Sr*Mx*Sc^-1
// assert Sr eq Sc
// Mx*Sr^-1*V0 eq Sr^-1*V1

conj_if_left:=func<x|xtr(xtr(x)*Pr^-1)*xtr(Pr)>;
// for nullspace == left this is Pr*x*Pr^-1

// it's normal if one of the two is false. It would help chasing bugs, in
// fact !
TMt eq conj_if_left(Sc*Cc*TM*Transpose(Cr)*Sr^-1);

// Mtt eq Pr*Sc*Cr*M*Transpose(Cc)*Sr^-1*Pr^-1;

TM*V0 eq V1; assert TM*V0 eq V1;
TM*V1 eq V2; assert TM*V1 eq V2;
TM*V2 eq V3; assert TM*V2 eq V3;
TM*V3 eq V4; assert TM*V3 eq V4;
TM^10*V0 eq V10; assert TM^10*V0 eq V10;
TM^10*V0b eq V10b; assert TM^10*V0b eq V10b;

// (Pr^-1*Mt)^10*C0 eq C10;
//
// This one only valid for nullspace==left
// Pr*(Pr^-1*Mt)^10*Pr^-1*C0 eq C10;
//
print "check vectors:";
C0 eq Submatrix(Xt,1,1,Nrows(Xt),64); 
C10 eq Submatrix(Transpose(TM)^10*Xt,1,1,Nrows(Xt),64);

Transpose(TM)^10*C0 eq C10;


Transpose(V0)*C10 eq Transpose(V10)*C0;
Transpose(C10)*V0 eq Transpose(C0)*V10;
Transpose(C10)*V0b eq Transpose(C0)*V10b;

print "Done checking krylov stuff";

load "A0-128.0-30.m";   // suited for matrix t100b

stride:=m*n div 64;
for i in [0..30-1] do
    vblock(m,var[i*stride+1..i*stride+stride],n) eq
    Transpose(Xt)*TM^i*Y0;
end for;

K:=GF(2);
b:=m+n;
N:=nr;
// KR<x>:=PowerSeriesRing(K);
KR<x>:=PolynomialAlgebra(K);
KRmn:=RMatrixSpace(KR,m,n);
KRbb:=RMatrixSpace(KR,b,b);
KRmb:=RMatrixSpace(KR,m,b);
KRnb:=RMatrixSpace(KR,n,b);
KRn1:=RMatrixSpace(KR,n,1);
KRnn:=RMatrixSpace(KR,n,n);
KRmn:=RMatrixSpace(KR,m,n);
KRNm:=RMatrixSpace(KR,N,m);
KRNn:=RMatrixSpace(KR,N,n);
KRN:=RMatrixSpace(KR,N,N);
KNn:=KMatrixSpace(K,N,n);

function mycoeff(M,k)
        s:=KMatrixSpace(K,Nrows(M), Ncols(M));
        return s![Coefficient(P,k):P in Eltseq(M)];
end function;

// beware -- we're using the sequence starting at 1.
// A_sequence:=&+[x^(i-1)*KRmn!(Transpose(Xt)*TM^i*V0):i in [1..100]];
A_sequence:=KRmn!0;
z:=Y0;
for i in [1..100] do
    A_sequence +:= x^(i-1)*KRmn!(Transpose(Xt)*z);
    z:=TM*z;
end for;
// assert A_sequence eq &+[x^(i-1)*KRmn!(Transpose(Xt)*TM^(i-1)*V0):i in [1..100]];


load "F.m";

degf:=#var div stride - 1;

F:=&+[x^(degf-i)*KRnn!vblock(n,var[i*stride+1..(i+1)*stride],n):i in [0..degf]];
Fr:=&+[x^i*      KRnn!vblock(n,var[i*stride+1..(i+1)*stride],n):i in [0..degf]];

printf "FYE: %o\n", IsZero(mycoeff(A_sequence*Transpose(F),degf));
for i in [degf+1..degf+5] cat [97..99] do
    // degf or degf+1 ? I'm seeing non-null at degf... It's pretty harmless.
    IsZero(mycoeff(A_sequence*Transpose(F),i));
end for;

not IsZero(mycoeff(A_sequence*Transpose(F),degf-1));

excess:=10;

V_sequence:=&+[x^i*KRNn!(TM^i*Y0):i in [0..degf+excess]];

assert IsZero(Transpose(Xt)*mycoeff(V_sequence*Transpose(F),degf+1));

// Transpose(Xt)*(KRN!TM)^4*(KRN!TM*V_sequence) eq mycoeff(A_sequence,3);

IsZero(mycoeff(Transpose(KRNm!Xt)*(&+[x^(i-1)*KRNn!(TM^(i-1)*Y0):i in [1..10\
0]])*Transpose(F),degf+1));


for i in [0..excess-1] do
    IsZero(Transpose(Xt)*TM*mycoeff(V_sequence*Transpose(F),degf + i));
    IsZero(TM*mycoeff(V_sequence*Transpose(F),degf + i));
end for;

S:=&+[TM^i*Y0*mycoeff(Transpose(Fr),i):i in [0..degf]];
IsZero(TM*S);

// S:=&+[TM^i*V0*mycoeff(Transpose(F),i):i in [0..2]];

load "S0-64.10.m"; S10:=vblock(nc, var, n);
load "S0-64.20.m"; S20:=vblock(nc, var, n);
load "S64-128.10.m"; S10b:=vblock(nc, var, n);
load "S64-128.20.m"; S20b:=vblock(nc, var, n);
load "K.0.m";    K0:=vblock(nc, var, 64);

S eq S10+S20;

// K0 is basically S in column echelon form.
Rowspace(Transpose(S)) eq Rowspace(Transpose(K0));
IsZero(TM*K0);



load "W.m";    W:=vblock(nr, var, 64);    


// For our simple test cases, K0==W. In more subtle situations, things
// might differ a bit.

print "FYI: ", K0 eq W;
IsZero(TM*W);
