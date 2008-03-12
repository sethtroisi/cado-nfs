#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include "struct.h"
#include "basic.h"
#include "random.h"
#include "alloc.h"
#include "mat_ops.h"
#include "echelon.h"
#include "ReadWrite.h"
#include "timing.h"

#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
#define	WBITS	(CHAR_BIT * sizeof(unsigned long))


/* FIXME !!! Get rid of jumbo functions */

/* function passed to sub-processes */
static void MyXORfunction(unsigned long *invec, unsigned long *inoutvec,
			  int *len, MPI_Datatype * dtype)
{
    unsigned long i;
    for (i = 0; i < *len; ++i) {
	inoutvec[i] ^= invec[i];
    }
};

/*

The Lanczos iterative process 


Remark: Needs MPI installed and running

*/



unsigned long *LanczosIterations(unsigned long *a, unsigned long *Y,
				 unsigned long m, unsigned long n,
				 unsigned long Block,
				 unsigned long *Resultado)
{

    float SumTimeSD = 0;
    float SumTimeTSD = 0;
    float SumTimeOper = 0;
    float SumTimeCom = 0;
    float SumTimeComBcast = 0;
    float SumTimeComRed = 0;
    unsigned long i;
    unsigned long SizeABlock;

    int size, p;
    unsigned long src, dest, tag, t0, t1;
    dest = 0;
    tag = 123;

    t0 = microseconds();

    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks



    if (p == size - 1) {
	SizeABlock = m / size + (m % size);
    } else {
	SizeABlock = m / size;
    }



// To all Processes




    MPI_Op newop;
    MPI_Op_create((MPI_User_function *) MyXORfunction, 0, &newop);

    unsigned long N = Block;
    unsigned long *Final;
    Final = malloc(sizeof(unsigned long));
    Final[0] = 0;

    unsigned long SizeRowV;

    SizeRowV = iceildiv(N, WBITS);


    unsigned long *ZnN;
    ZnN = Allocmn(n, N);

    memset(ZnN, 0, n * SizeRowV * sizeof(unsigned long));

    unsigned long *X;
    X = Allocmn(n, N);


    unsigned long **V;
    V = Allocmn3(n, N);




    unsigned long *ATAC;
    unsigned long *ATAC_Dist;
    unsigned long *ResmN_Dist;


    ATAC = Allocmn(n, N);
    ResmN_Dist = Allocmn(m, N);
    ATAC_Dist = Allocmn(n, N);


    SMultDmatrixBit(m, n, N, a, Y, ResmN_Dist, p * (m / size), SizeABlock);
    TSMultDmatrixBit(SizeABlock, n, N, a, ResmN_Dist, ATAC_Dist,
		     p * (m / size), SizeABlock);

    MPI_Reduce(ATAC_Dist, ATAC, iceildiv(Block, WBITS) * n, MPI_UNSIGNED_LONG,
	       newop, 0, MPI_COMM_WORLD);

    MPI_Bcast(ATAC, n * iceildiv(Block, WBITS), MPI_UNSIGNED_LONG, 0,
	      MPI_COMM_WORLD);

    free(ResmN_Dist);
    free(ATAC_Dist);


    for (i = 0; i < n * iceildiv(Block, WBITS); ++i) {
	V[0][i] = ATAC[i];
	V[1][i] = ZnN[i];
	V[2][i] = ZnN[i];
	X[i] = ZnN[i];
    }



    //initial conditions data prodution



    unsigned long *IN;

    IN = Allocmn(N, N);

    memset(IN, 0, N * SizeRowV * sizeof(unsigned long));

    for (i = 0; i < N; ++i) {
	IN[i * SizeRowV + i / WBITS] = (1UL << i % WBITS);
    }


    unsigned long *ZN;
    ZN = Allocmn(N, N);

    memset(ZN, 0, N * SizeRowV * sizeof(unsigned long));

    unsigned long **Winv;
    Winv = Allocmn3(N, N);
    unsigned long **S;
    S = Allocmn3(N, N);

    for (i = 0; i < N * SizeRowV; ++i) {
	Winv[0][i] = ZN[i];
	Winv[1][i] = ZN[i];
	Winv[2][i] = ZN[i];
	S[0][i] = IN[i];
	S[1][i] = IN[i];
	S[2][i] = IN[i];
    }


    unsigned long *SizeS;
    SizeS = malloc(3 * sizeof(unsigned long));
    SizeS[0] = N;
    SizeS[1] = N;
    SizeS[2] = N;

    unsigned long **Set;
    Set = malloc(3 * sizeof(unsigned long *));
    Set[0] = malloc((N + 1) * sizeof(unsigned long));
    Set[1] = malloc((N + 1) * sizeof(unsigned long));
    Set[2] = malloc((N + 1) * sizeof(unsigned long));

    for (i = 0; i < N; ++i) {
	Set[0][i] = i;
	Set[1][i] = i;
	Set[2][i] = i;
    }




// To all Processes




    unsigned long *AVi_Dist;


    ResmN_Dist = Allocmn(m, N);
    AVi_Dist = Allocmn(n, N);


    unsigned long *ResNn;
    ResNn = Allocmn(N, n);

    unsigned long *AVi;
    AVi = Allocmn(n, N);
    unsigned long *ViAVi;
    ViAVi = Allocmn(N, N);
    unsigned long *ViA2Vi;
    ViA2Vi = Allocmn(N, N);


    unsigned long *UseSet;
    UseSet = malloc((N + 1) * sizeof(unsigned long));

    UseSet[0] = N;
    for (i = 1; i <= N; ++i) {
	UseSet[i] = i - 1;
    }



    unsigned long *Set0;
    Set0 = malloc(sizeof(unsigned long));

    unsigned long *nW_inv;
    nW_inv = Allocmn(N, N);




    unsigned long *ntaken;

    ntaken = malloc((N + 1) * sizeof(unsigned long));



    unsigned long *nSi;
    nSi = Allocmn(N, N);



    SMultDmatrixBit(m, n, N, a, V[0], ResmN_Dist, p * (m / size), SizeABlock);
    TSMultDmatrixBit(SizeABlock, n, N, a, ResmN_Dist, AVi_Dist,
		     p * (m / size), SizeABlock);


    MPI_Reduce(AVi_Dist, AVi, iceildiv(Block, WBITS) * n, MPI_UNSIGNED_LONG,
	       newop, 0, MPI_COMM_WORLD);


    if (p == 0) {

/*    
    TransposeBit(n, N, V[0], ResNn);
    DMultBit( N, n, N, ResNn, AVi, ViAVi);
    TransposeBit(n, N, AVi, ResNn);
    DMultBit( N,n, N, ResNn, AVi, ViA2Vi); 
 */

	TVUBit(n, N, V[0], AVi, ViAVi);
	TVUBit(n, N, AVi, AVi, ViA2Vi);

	Construction_S_WinvBit(N, ViAVi, UseSet, nSi, Set0, nW_inv, ntaken);

	Refresh_ArrayBit(N, N, Winv, nW_inv);
	Refresh_ArrayBit(N, N, S, nSi);


	SizeS[2] = SizeS[1];
	SizeS[1] = SizeS[0];
	SizeS[0] = ntaken[0];

	for (i = 0; i < N + 1; ++i) {
	    Set[2][i] = Set[1][i];
	    Set[1][i] = Set[0][i];
	    Set[0][i] = ntaken[i];
	}

// Construct X0



	unsigned long *ViWiinv;
	unsigned long *ViTV0;
	unsigned long *Xi;
	unsigned long *ResNnX;

	ViWiinv = Allocmn(n, N);
	ViTV0 = Allocmn(N, N);
	Xi = Allocmn(n, N);
	ResNnX = Allocmn(N, n);

//    TransposeBit(n, N, V[0], ResNnX);    
//    DMultBit( N, n, N, ResNnX, ATAC ,ViTV0);

	TVUBit(n, N, V[0], ATAC, ViTV0);

//    DMultBit( n, N, N, V[0], Winv[0] ,ViWiinv);

	VUBit(n, N, V[0], Winv[0], ViWiinv);

//    DMultBit( N, n, N, ResNnX, ATAC ,ViTV0);

//    DMultBit( n, N, N, ViWiinv, ViTV0 ,Xi);

	VUBit(n, N, ViWiinv, ViTV0, Xi);

	for (i = 0; i < iceildiv(Block, WBITS) * n; ++i) {
	    X[i] ^= Xi[i];
	}

	free(ViWiinv);
	free(ViTV0);
	free(Xi);
	free(ResNnX);


// end construct X0


	printf("Step= %d      SizeWi=  %lu     Number of processes= %d\n", 0,
	       SizeS[0], size);


    }

// To all processes

    unsigned long Step = 0;
    unsigned long *AVi1;
    AVi1 = Allocmn(n, N);
    unsigned long *Vi1AVi1;
    Vi1AVi1 = Allocmn(N, N);
    unsigned long *Vi1A2Vi1;
    Vi1A2Vi1 = Allocmn(N, N);

    unsigned long *AVi1_Dist;
    AVi1_Dist = Allocmn(n, N);



    SMultDmatrixBit(m, n, N, a, V[1], ResmN_Dist, p * (m / size), SizeABlock);
    TSMultDmatrixBit(SizeABlock, n, N, a, ResmN_Dist, AVi1_Dist,
		     p * (m / size), SizeABlock);

    MPI_Reduce(AVi1_Dist, AVi1, iceildiv(Block, WBITS) * n, MPI_UNSIGNED_LONG,
	       newop, 0, MPI_COMM_WORLD);


    if (p == 0) {

/*
    TransposeBit(n, N, V[1], ResNn);
    DMultBit( N, n, N, ResNn, AVi1, Vi1AVi1);
    TransposeBit(n, N, AVi1, ResNn);
    DMultBit( N,n, N, ResNn, AVi1, Vi1A2Vi1);
*/
	TVUBit(n, N, V[1], AVi1, Vi1AVi1);
	TVUBit(n, N, AVi1, AVi1, Vi1A2Vi1);

    }


//Desalloc 

    free(AVi1);
    free(AVi1_Dist);







    //unsigned long Step = 1;
    unsigned long SumRank = SizeS[0];


//Alloc



    unsigned long *V1i;
    V1i = Allocmn(n, N);
    unsigned long *D1i;
    D1i = Allocmn(N, N);
    unsigned long *ResNN;
    ResNN = Allocmn(N, N);
    unsigned long *ResNN2;
    ResNN2 = Allocmn(N, N);
    unsigned long *E1i;
    E1i = Allocmn(N, N);
    unsigned long *F1i;
    F1i = Allocmn(N, N);
    unsigned long *ResNN3;
    ResNN3 = Allocmn(N, N);
    unsigned long *ResnN;
    ResnN = Allocmn(n, N);


    t1 = microseconds();

//the iteration



    while (((TestZero(N, N, ViAVi) == 0) || (p != dest)) && (Final[0] == 0)) {

	unsigned long T00;

	T00 = microseconds();

	if (p == 0)
	{

	    for (i = 1; i <= Set0[0]; ++i)
		UseSet[i] = Set[0][i];
	    UseSet[0] = Set0[0];

// for V1i

	    for (i = 0; i < iceildiv(N, WBITS) * n; ++i)
		V1i[i] = AVi[i];
// end for V1i


//        D1i:=IN-Winv`i*(ViA2Vi*S`i*Transpose(S`i)+ViAVi);

	    unsigned long *ResVN;
	    ResVN = Allocmn(Set0[0], N);
	    unsigned long *ResNV;
	    ResNV = Allocmn(N, Set0[0]);

	    SelectLinesBit(N, S[0], SizeS[0], ResVN);
	    TransposeBit(SizeS[0], N, ResVN, ResNV);
	    DMultBit(N, SizeS[0], N, ResNV, ResVN, ResNN);


// for V1i

//      DMultBit(n, N, N, V1i, ResNN, ResnN);
	    VUBit(n, N, V1i, ResNN, ResnN);

	    for (i = 0; i < iceildiv(N, WBITS) * n; ++i)
		V1i[i] = ResnN[i];

// end for V1i



	    DMultBit(N, N, N, ViA2Vi, ResNN, ResNN2);
	    DMatrixSumBit(N, N, ResNN2, ViAVi, ResNN);
	    DMultBit(N, N, N, Winv[0], ResNN, ResNN2);
	    DMatrixSumBit(N, N, IN, ResNN2, D1i);



//      E1i:=-Winv`i1*ViAVi*S`i*Transpose(S`i);

	    SelectLinesBit(N, S[0], SizeS[0], ResVN);
	    TransposeBit(SizeS[0], N, ResVN, ResNV);
	    DMultBit(N, SizeS[0], N, ResNV, ResVN, ResNN);
	    DMultBit(N, N, N, ViAVi, ResNN, ResNN2);
	    DMultBit(N, N, N, Winv[1], ResNN2, E1i);



// F1i:=-Winv`i2*(IN-Vi1AVi1*Winv`i1)*(Vi1A2Vi1*S`i1*Transpose(S`i1)+Vi1AVi1)*S`i*Transpose(S`i);


	    unsigned long *ResVN2;
	    ResVN2 = Allocmn(SizeS[1], N);
	    unsigned long *ResNV2;
	    ResNV2 = Allocmn(N, SizeS[1]);



	    SelectLinesBit(N, S[1], SizeS[1], ResVN2);
	    TransposeBit(SizeS[1], N, ResVN2, ResNV2);
	    DMultBit(N, SizeS[1], N, ResNV2, ResVN2, ResNN);
	    DMultBit(N, N, N, Vi1A2Vi1, ResNN, ResNN2);
	    DMatrixSumBit(N, N, ResNN2, Vi1AVi1, ResNN);
	    DMultBit(N, SizeS[0], N, ResNV, ResVN, ResNN2);
	    DMultBit(N, N, N, ResNN, ResNN2, ResNN3);
	    DMultBit(N, N, N, Vi1AVi1, Winv[1], ResNN);
	    DMatrixSumBit(N, N, IN, ResNN, ResNN2);
	    DMultBit(N, N, N, ResNN2, ResNN3, ResNN);
	    DMultBit(N, N, N, Winv[2], ResNN, F1i);


// V1i:=AVi*S`i*Transpose(S`i)+V`i*D1i+V`i1*E1i+V`i2*F1i;


//      DMultBit(n, N, N, V[0], D1i, ResnN);
	    VUBit(n, N, V[0], D1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);
//      DMultBit(n, N, N, V[1], E1i, ResnN);

	    VUBit(n, N, V[1], E1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);
//      DMultBit(n, N, N, V[2], F1i, ResnN);

	    VUBit(n, N, V[2], F1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);



	    Refresh_ArrayBit(n, N, V, V1i);




//Desalloc   


	    free(ResVN);
	    free(ResVN2);
	    free(ResNV2);
	    free(ResNV);



	    for (i = 0; i < iceildiv(N, WBITS) * N; ++i) {
		Vi1AVi1[i] = ViAVi[i];
		Vi1A2Vi1[i] = ViA2Vi[i];
	    }


	}

// To All Processes

	unsigned long T1, T2, T3, T4, T0, T5;

	T0 = microseconds();


	MPI_Bcast(V1i, n * iceildiv(Block, WBITS), MPI_UNSIGNED_LONG, 0,
		  MPI_COMM_WORLD);


	T1 = microseconds();

	SMultDmatrixBit(m, n, N, a, V1i, ResmN_Dist, p * (m / size),
			SizeABlock);

	T2 = microseconds();
	SumTimeSD += (T2 - T1);

	TSMultDmatrixBit(SizeABlock, n, N, a, ResmN_Dist, AVi_Dist,
			 p * (m / size), SizeABlock);

	T3 = microseconds();
	SumTimeTSD += (T3 - T2);

	//sleep(2);


	if (p != 0) {
	    Step += 1;
	}


	MPI_Reduce(AVi_Dist, AVi, iceildiv(Block, WBITS) * n,
		   MPI_UNSIGNED_LONG, newop, 0, MPI_COMM_WORLD);


	//sleep(2);

	T4 = microseconds();
	if (p == 0) {
	    SumTimeOper += (T3 - T1) + (T0 - T00);
	    SumTimeCom += (T4 - T0 - (T3 - T1));
	    SumTimeComBcast += T1 - T0;
	    SumTimeComRed += T4 - T3;
	}


	if (p == 0) {



	    unsigned long *ResNn;
	    ResNn = Allocmn(N, n);
/*
    TransposeBit(n, N, V1i, ResNn);

    unsigned long tNnN0=microseconds();

    DMultBit( N, n, N, ResNn, AVi, ViAVi);

    unsigned long tNnN1=microseconds();
*/

	    TVUBit(n, N, V1i, AVi, ViAVi);

	    TVUBit(n, N, AVi, AVi, ViA2Vi);
/*
    TransposeBit(n, N, AVi, ResNn);

    DMultBit( N,n, N, ResNn, AVi, ViA2Vi);
*/

	    Construction_S_WinvBit(N, ViAVi, UseSet, nSi, Set0, nW_inv,
				   ntaken);

	    Refresh_ArrayBit(N, N, Winv, nW_inv);
	    Refresh_ArrayBit(N, N, S, nSi);

	    SizeS[2] = SizeS[1];
	    SizeS[1] = SizeS[0];
	    SizeS[0] = ntaken[0];

	    SumRank += SizeS[0];

	    Refresh_ArrayBit(N + 1, 1, Set, ntaken);

// Construct X0

	    unsigned long *ViWiinv;
	    unsigned long *ViTV0;
	    unsigned long *Xi;
	    unsigned long *ResNnX;

	    ViWiinv = Allocmn(n, N);
	    ViTV0 = Allocmn(N, N);
	    Xi = Allocmn(n, N);
	    ResNnX = Allocmn(N, n);

/*
 
    TransposeBit(n, N, V[0], ResNnX);
    DMultBit( N, n, N, ResNnX, ATAC ,ViTV0);
*/

	    unsigned long tnNN0 = microseconds();

	    TVUBit(n, N, V[0], ATAC, ViTV0);

	    unsigned long tnNN1 = microseconds();

//    DMultBit( n, N, N, V[0], Winv[0] ,ViWiinv);

	    VUBit(n, N, V[0], Winv[0], ViWiinv);



	    //DMultBit( N, n, N, ResNnX, ATAC ,ViTV0);



//    DMultBit( n, N, N, ViWiinv, ViTV0 ,Xi);

	    VUBit(n, N, ViWiinv, ViTV0, Xi);


	    for (i = 0; i < iceildiv(Block, WBITS) * n; ++i) {
		X[i] ^= Xi[i];
	    }

	    free(ViWiinv);
	    free(ViTV0);
	    free(Xi);
	    free(ResNnX);

// end construct X0

	    T5 = microseconds();
	    SumTimeOper += (T5 - T4);

	    float SD = T2 - T1;
	    float TSD = T3 - T2;
	    float TIteration = T5 - T00;
	    float BCalc = T0 - T00;
	    float ACalc1 = T5 - T4;
	    float ACalc2 = T4 - T3;
	    //float tNnN=tnNN1-tNnN0;
	    float tnNN = tnNN1 - tnNN0;

	    printf
		("Step= %lu      SizeWi=  %lu    SumRank= %lu    Machine= %d   TimeSD=%f s  TimeTSD=%f s  TITime=%f s   SAlgebra1=%f s  Salgebra2=%f   TnNN=%f\n",
		 Step, SizeS[0], SumRank, p, SD / 1000000, TSD / 1000000,
		 TIteration / 1000000, BCalc / 1000000, ACalc1 / 1000000,
		 tnNN / 1000000);

	    //sleep(5);

	    Step += 1;


	    free(ResNn);
	}
//Stop criteria for all processes

	if (p == 0) {
	    if (SizeS[0] == 0) {
		Final[0] = 1;
	    }
	    for (src = 1; src < size; ++src) {
		MPI_Send(Final, 1, MPI_UNSIGNED_LONG, src, tag,
			 MPI_COMM_WORLD);
	    }
	} else {
	    MPI_Recv(Final, 1, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD,
		     &status);
	}




    }



    if (p == 0) {


	char *f;
	f = malloc(10 * sizeof(unsigned long));

	sprintf(f, "ControlFile_%lu_%lu_%lu.txt", m, n, Block);

	FILE *File;
	File = fopen(f, "a");

	printf("ACompT = %f ms AComT = %f ms ",
	       SumTimeOper / (Step - 1) / 1000,
	       SumTimeCom / (Step - 1) / 1000);
	fprintf(File, "ACompT = %f ms AComT = %f ms ",
		SumTimeOper / (Step - 1) / 1000,
		SumTimeCom / (Step - 1) / 1000);

	printf("ABcastT = %f ms  ARedT= %f ms  ",
	       SumTimeComBcast / (Step - 1) / 1000,
	       SumTimeComRed / (Step - 1) / 1000);
	fprintf(File, "ABcastT = %f ms  ARedT= %f ms  ",
		SumTimeComBcast / (Step - 1) / 1000,
		SumTimeComRed / (Step - 1) / 1000);
	printf("TPTime = %lu ms  ", (t1 - t0) / 1000);

	fprintf(File, "TPTime = %lu ms  ", (t1 - t0) / 1000);

	fclose(File);



	for (i = 0; i < iceildiv(Block, WBITS) * n; ++i) {
	    Resultado[i] = X[i] ^ Y[i];
	}

    } else {

    }


//Desalloc 

    free(AVi_Dist);
    free(ResmN_Dist);
    //free(C);
    //free(Y);
    free(ATAC);

    //free(a);
    free(Final);
    free(ResNn);
    free(ResNN);
    free(ResNN2);
    free(ResnN);
    free(ResNN3);
    free(V1i);
    free(AVi);
    free(ViAVi);
    free(ViA2Vi);
    free(D1i);
    free(E1i);
    free(F1i);

    free(IN);
    free(ZN);
    free(ZnN);


    free(nW_inv);
    free(ntaken);
    free(nSi);

    free(SizeS);
    free(UseSet);
    free(Set0);

    free(Vi1AVi1);
    free(Vi1A2Vi1);

    FreeAllocmn3(V);
    FreeAllocmn3(Winv);
    FreeAllocmn3(Set);
    FreeAllocmn3(S);
    // free(f);
    // free(f1);

    return 0;
}



/*

"Small" Linear algebra to get a 64 block of vectors in the Kernel of MTM

*/






unsigned long *KernelSparse(unsigned long *a, unsigned long *R,
			    unsigned long m, unsigned long n,
			    unsigned long Block, unsigned long *Ker,
			    unsigned long *Ind)
{


    unsigned long SizeABlock;
    int size, p;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

    if (p == size - 1) {
	SizeABlock = m / size + (m % size);
    } else {
	SizeABlock = m / size;
    }



// To all Processes



    MPI_Op newop;
    MPI_Op_create((MPI_User_function *) MyXORfunction, 0, &newop);

    unsigned long N = Block;
    unsigned long SizeRowV;

    SizeRowV = iceildiv(N, WBITS);

    unsigned long *ATAR;
    unsigned long *ATAR_Dist;
    unsigned long *ResmN_Dist;


// To do Transpose(a)*a*R


    ATAR = Allocmn(n, N);
    ResmN_Dist = Allocmn(m, N);
    ATAR_Dist = Allocmn(n, N);



    SMultDmatrixBit(m, n, N, a, R, ResmN_Dist, p * (m / size), SizeABlock);
    TSMultDmatrixBit(SizeABlock, n, N, a, ResmN_Dist, ATAR_Dist,
		     p * (m / size), SizeABlock);
    MPI_Reduce(ATAR_Dist, ATAR, iceildiv(Block, WBITS) * n, MPI_UNSIGNED_LONG,
	       newop, 0, MPI_COMM_WORLD);

    free(ResmN_Dist);
    free(ATAR_Dist);


    if (p == 0) {

	unsigned long *ResVn, *ResNN, *E, *Aout, *ResnN, *ResNn, *ResNn2,
	    *ResnV, *ListLines, *LengListLines;
	Aout = Allocmn(Block, n);
	ResNn = Allocmn(Block, n);
	ResNn2 = Allocmn(Block, n);

	ResnN = Allocmn(n, Block);
	E = Allocmn(Block, Block);
	ResNN = Allocmn(Block, Block);
	ListLines = malloc((N + 1) * sizeof(unsigned long));
	LengListLines = malloc(sizeof(unsigned long));

	TransposeBit(n, N, ATAR, ResNn);

	GaussElimBit(Block, n, ResNn, Aout, E, ListLines, LengListLines);

	TransposeBit(n, N, R, ResNn);

	DMultBit(N, N, n, E, ResNn, ResNn2);

	ResVn = Allocmn(Block - LengListLines[0], n);

	SelectLinesListBit(n, ResNn2, ListLines, Block - LengListLines[0],
			   ResVn);

	ResnV = Allocmn(n, Block - LengListLines[0]);

	TransposeBit(Block - LengListLines[0], n, ResVn, Ker);

	Ind[0] = n;
	Ind[1] = Block - LengListLines[0];

	free(ResVn);
	free(ResNN);
	free(E);
	free(Aout);
	free(ResnN);
	free(ResNn);
	free(ResNn2);
	free(ResnV);
	free(ListLines);
	free(LengListLines);
    }



    free(ATAR);



    return 0;
}


/*

Complete Lanczos

*/




struct DenseMatrix Lanczos(struct SparseMatrix M, unsigned long Block)
{

    int size, p;



    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks

    struct DenseMatrix Kernel;
    Kernel.Data = Allocmn(M.Ncols, Block);

    unsigned long *Y;
    Y = Allocmn(M.Ncols + 1, Block);

    unsigned long T0, T1;

    if (p == 0) {
	RandomDMatrixBitTest(M.Ncols, Block, Y);
    }

    MPI_Bcast(Y, M.Ncols * iceildiv(Block, WBITS), MPI_UNSIGNED_LONG, 0,
	      MPI_COMM_WORLD);


// The Lanczos iteration

    unsigned long *Result, *Index;
    Result = Allocmn(M.Ncols, Block);
    Index = malloc(2 * sizeof(unsigned long));

    T0 = microseconds();

    LanczosIterations(M.Data, Y, M.Nrows, M.Ncols, Block, Result);

    T1 = microseconds();


// Small linear algebra to compute the  elements of the Kernel 


    KernelSparse(M.Data, Result, M.Nrows, M.Ncols, Block, Kernel.Data, Index);
    Kernel.Nrows = Index[0];
    Kernel.Ncols = Index[1];

    //displayMatrix(Kernel.Data,Index[0],Index[1],'d');
    //displaySMatrixNew(M.Data,M.Nrows,M.Ncols,'M');

    free(Y);
    free(Result);


    return Kernel;

}


#if 0
unsigned long Lanczos(struct SparseMatrix M, unsigned long Block,
		      struct DenseMatrix Kernel)
{

    int size, p;



    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks


    unsigned long *Y;
    Y = Allocmn(M.Ncols + 1, Block);

    unsigned long T0, T1;

    if (p == 0) {
	RandomDMatrixBitTest(M.Ncols, Block, Y);
    }

    MPI_Bcast(Y, M.Ncols * iceildiv(Block, WBITS), MPI_UNSIGNED_LONG, 0,
	      MPI_COMM_WORLD);


// The Lanczos iteration

    unsigned long *Result, *Index;
    Result = Allocmn(M.Ncols, Block);
    Index = malloc(2 * sizeof(unsigned long));

    T0 = microseconds();

    LanczosIterations(M.Data, Y, M.Nrows, M.Ncols, Block, Result);

    T1 = microseconds();


// Small linear algebra to compute the  elements of the Kernel 


    KernelSparse(M.Data, Result, M.Nrows, M.Ncols, Block, Kernel.Data, Index);
    Kernel.Nrows = Index[0];
    Kernel.Ncols = Index[1];

    //displayMatrix(Kernel.Data,Index[0],Index[1],'d');
    //displaySMatrixNew(M.Data,M.Nrows,M.Ncols,'M');

    free(Y);
    free(Result);


    return 0;

}
#endif
