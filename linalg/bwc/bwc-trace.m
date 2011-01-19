
m:=64;
n:=64;

load "/tmp/bwc/t.m"; M:=Matrix(GF(2),Matrix(var));
nr:=Nrows(M);
nc:=Ncols(M);

assert nc le nr;
nc:=nr;
x:=Matrix(GF(2),nr,nc,[]);InsertBlock(~x,M,1,1);M:=x;
load "/tmp/bwc/placemats.m";

load "/tmp/bwc/b.m"; 

//L:=[]; for i in [0..nh-1] do
//t:=nr div nh + (i lt nr mod nh select 1 else 0); 
//s:=i*nrp;
//L cat:=[s+1..s+t];
//end for;Lr:=L;
//for j in [1..#Lr] do Cr[Lr[j]][j]:=1; end for;
Cr:=Matrix(GF(2),nh*nrp,nr,[]);
InsertBlock(~Cr,IdentityMatrix(GF(2),nr),1,1);
//L:=[]; for j in [0..nv-1] do
//t:=nc  div nv + (j lt nc mod nv select 1 else 0);
//s:=j*ncp;
//L cat:=[s+1..s+t];
//end for;Lc:=L;
// for j in [1..#Lc] do Cc[Lc[j]][j]:=1; end for;
Cc:=Matrix(GF(2),nv*ncp,nc,[]);
InsertBlock(~Cc,IdentityMatrix(GF(2),nc),1,1);

assert IsOne(Transpose(Cr)*Cr);
assert IsOne(Transpose(Cc)*Cc);

ncx := nv*ncp;
nrx := nh*nrp;

sc_inv:=SymmetricGroup(nv*ncp)!var;
sr_inv:=SymmetricGroup(nh*nrp)!var;
sr:=sr_inv^-1;
sc:=sc_inv^-1;
Sr:=PermutationMatrix(GF(2),sr_inv)/*^-1*/;     // magma thinks backwards
Sc:=PermutationMatrix(GF(2),sc_inv)/*^-1*/;     // magma thinks backwards

Mx:=Cr*M*Transpose(Cc);

Mx[2] eq (Sr*Mx)[Eltseq(sr)[2]];
Transpose(Mx)[2] eq (Sc*Transpose(Mx))[Eltseq(sc)[2]];

nz:=nrp div nv;
assert nz eq ncp div nh;
pr:=func<x|(qr*nh+qq)*nz+r+1 where qq,qr is Quotrem(q, nv) where q,r is Quotrem(x-1,nz)>;
prinv:=func<x|(qr*nv+qq)*nz+r+1 where qq,qr is Quotrem(q, nh) where q,r is Quotrem(x-1,nz)>;
Pr:=PermutationMatrix(GF(2),[pr(x):x in [1..nrx]]);


// t means twisted.

/* Note that Mt being equal to Pr times Sr*Cr*M*Transpose(Cc)*Sc^-1 is a
 * direct consequence of the fact that balancing_workhorse.c permutes the
 * rows which are being sent.
 */
assert Transpose(Pr*Sr*Cr)*Mt*Sc*Cc eq M;
assert Mt eq Pr*Sr*Cr*M*Transpose(Cc)*Sc^-1;

colweights:=func<M|[#[i:i in [1..Nrows(M)]|M[i,j] ne 0]: j in [1..Ncols(M)]]>;
rowweights:=func<M|[#[j:j in [1..Ncols(M)]|M[i,j] ne 0]: i in [1..Nrows(M)]]>;
// colsums:=func<M|[&+[Integers()!M[i,j]:i in [1..Nrows(M)]]: j in [1..Ncols(M)]]>;
// rowsums:=func<M|[&+[Integers()!M[i,j]:j in [1..Ncols(M)]]: i in [1..Nrows(M)]]>;

assert colweights(Mt*Sc*Cc) eq colweights(M);
assert rowweights(Transpose(Pr*Sr*Cr)*Mt) eq rowweights(M);





load "/tmp/bwc/x.m";
C0:=Matrix(GF(2),nrx,64,[]);
Xt:=Matrix(GF(2),nrx,n,[]);
for i in [1..64] do for u in var[i] do C0[u][i] +:=1; end for; end for;
for i in [1..n] do for u in var[i] do Xt[u][i] +:=1; end for; end for;

function vblock(nr,var, nc)
    R:=Matrix(GF(2),nr,nc,[]);
    nb:=nc div 64;
    for i in [0..nr-1] do
        R[i+1]:=Vector(GF(2),Intseq(u,2,nc) where u is Seqint(var[1+i*nb..(i+1)*nb],2^64));
    end for;
  return R;
end function;


load "/tmp/bwc/Y.0.m";     Y0:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.0.m"; V0:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.1.m"; V1:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.2.m"; V2:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.3.m"; V3:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.4.m"; V4:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.5.m"; V5:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.6.m"; V6:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.7.m"; V7:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.8.m"; V8:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.9.m"; V9:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.10.m"; V10:=vblock(ncx, var, n);
load "/tmp/bwc/V0-64.20.m"; V20:=vblock(ncx, var, n);
// load "/tmp/bwc/V0-64.12.m"; V12:=vblock(ncx, var, n);
// load "/tmp/bwc/V0-64.13.m"; V13:=vblock(ncx, var, n);
// load "/tmp/bwc/V0-64.14.m"; V14:=vblock(ncx, var, n);
// load "/tmp/bwc/V0-64.15.m"; V15:=vblock(ncx, var, n);
// load "/tmp/bwc/V0-64.16.m"; V16:=vblock(ncx, var, n);
// load "/tmp/bwc/V0-64.17.m"; V17:=vblock(ncx, var, n);
// load "/tmp/bwc/V0-64.18.m"; V18:=vblock(ncx, var, n);
load "/tmp/bwc/C.10.m";    C10:=vblock(nrx, var, 64);

v0z:=[Seqint(ChangeUniverse((Eltseq(V0[i])),Integers()),2):i in [1..nr]];

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
TMt eq conj_if_left(Sc*Cc*xtr(M)*Transpose(Cr)*Sr^-1);

// Mtt eq Pr*Sc*Cr*M*Transpose(Cc)*Sr^-1*Pr^-1;

TMt*V0 eq V1; assert TMt*V0 eq V1;
TMt*V1 eq V2; assert TMt*V1 eq V2;
TMt*V2 eq V3; assert TMt*V2 eq V3;
TMt*V3 eq V4; assert TMt*V3 eq V4;

TMt^10*V0 eq V10; assert TMt^10*V0 eq V10;

// (Pr^-1*Mt)^10*C0 eq C10;
//
// This one only valid for nullspace==left
// Pr*(Pr^-1*Mt)^10*Pr^-1*C0 eq C10;
Transpose(TMt)^10*C0 eq C10;


Transpose(V0)*C10 eq Transpose(V10)*C0;
Transpose(C10)*V0 eq Transpose(C0)*V10;

// D10:=Pr*(Pr^-1*Mt)^10*Pr^-1*C0;
// Transpose(V0)*D10 eq Transpose(V10)*C0;
// Transpose(D10)*V0 eq Transpose(C0)*V10;

print "Done checking krylov stuff";


load "/tmp/bwc/A0-64.0-20.m";

stride:=m*n div 64;
for i in [0..20-1] do
    vblock(m,var[i*stride+1..i*stride+stride],n) eq Transpose(Xt)*TMt^i*V0;
end for;

K:=GF(2);
b:=m+n;
N:=nrx;
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
// A_sequence:=&+[x^(i-1)*KRmn!(Transpose(Xt)*TMt^i*V0):i in [1..100]];
A_sequence:=KRmn!0;
z:=V0;
for i in [1..100] do
    A_sequence +:= x^(i-1)*KRmn!(Transpose(Xt)*z);
    z:=TMt*z;
end for;
// assert A_sequence eq &+[x^(i-1)*KRmn!(Transpose(Xt)*TMt^(i-1)*V0):i in [1..100]];


load "/tmp/bwc/F0-64.m";

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

V_sequence:=&+[x^i*KRNn!(TMt^i*V0):i in [0..degf+excess]];

assert IsZero(Transpose(Xt)*mycoeff(V_sequence*Transpose(F),degf+1));

// Transpose(Xt)*(KRN!TMt)^4*(KRN!TMt*V_sequence) eq mycoeff(A_sequence,3);

IsZero(mycoeff(Transpose(KRNm!Xt)*(&+[x^(i-1)*KRNn!(TMt^(i-1)*V0):i in [1..10\
0]])*Transpose(F),degf+1));


for i in [0..excess-1] do
    IsZero(Transpose(Xt)*TMt*mycoeff(V_sequence*Transpose(F),degf + i));
    IsZero(TMt*mycoeff(V_sequence*Transpose(F),degf + i));
end for;

S:=&+[TMt^i*V0*mycoeff(Transpose(Fr),i):i in [0..degf]];
IsZero(TMt*S);

// S:=&+[TMt^i*V0*mycoeff(Transpose(F),i):i in [0..2]];

load "/tmp/bwc/S0-64.10.m"; S10:=vblock(ncx, var, n);

S eq S10;

b:=Basis(sub<VectorSpace(K,N)|Rows(Transpose(S))>);

IsZero(b[1] * Transpose(TMt));

load "/tmp/bwc/W.m";    W:=vblock(nrx, var, 64);    
load "/tmp/bwc/Wu.m";    Wu:=vblock(nr, var, 64);

print "The following checks are for nullspace=left only";
IsZero(Transpose(W)*Mt);
IsZero(Transpose(W)*Pr*Sr*Cr*M*Transpose(Cc)*Sc^-1);
IsZero(Transpose(W)*Pr*Sr*Cr*M);
Matrix(((Pr*Sr)^-1*W)[1..Nrows(M)]) eq Wu;
IsZero(Transpose(Wu)*M);


