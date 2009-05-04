#include "mpi_select.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
//#include "struct.h"
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




void ToConstructX(DenseMatrix X, DenseMatrix Y, DenseMatrix V, DenseMatrix W)
{
        DenseMatrix ViWiinv,ViTV0,Xi,ResNnX;
        unsigned long i;
        unsigned long n=X->Nrows;
        unsigned long N=X->Ncols;
	ViWiinv->Data = Allocmn(n, N);ViWiinv->Nrows=n;ViWiinv->Ncols=N;
	ViTV0->Data = Allocmn(N, N);ViTV0->Nrows=N;ViTV0->Ncols=N;
	Xi->Data = Allocmn(n, N);Xi->Nrows=n;Xi->Ncols=N;
	ResNnX->Data = Allocmn(N, n);ResNnX->Nrows=N;ResNnX->Ncols=n;

	TVUBit_new(V,Y, ViTV0);
	VUBit_new(V,W, ViWiinv);
	VUBit_new(ViWiinv, ViTV0, Xi);
	for (i = 0; i < iceildiv(N, WBITS) * n; ++i) {
	    X->Data[i] ^= Xi->Data[i];
	}

	free(ViWiinv->Data);
	free(ViTV0->Data);
	free(Xi->Data);
	free(ResNnX->Data);

}



void ToConstructD(DenseMatrix D1i, DenseMatrix S, DenseMatrix W, DenseMatrix V, DenseMatrix ViAVi, DenseMatrix ViA2Vi, DenseMatrix IN,unsigned long s)
{
            unsigned long N=V->Ncols;
            unsigned long n=V->Nrows;
            unsigned long i;
	    DenseMatrix ResVN,ResNV,ResNN,ResNN2,ResnN;
            ResnN->Data = Allocmn(n, N);ResnN->Nrows=n;ResnN->Ncols=N;
            ResNN->Data = Allocmn(N, N);ResNN->Nrows=N;ResNN->Ncols=N;
            ResNN2->Data = Allocmn(N, N);ResNN2->Nrows=N;ResNN2->Ncols=N;
	    ResVN->Data = Allocmn(s, N);ResVN->Nrows=s;ResVN->Ncols=N;
	    ResNV->Data = Allocmn(N, s);ResNV->Ncols=s;ResVN->Nrows=N;
	    SelectLinesBit_new(S,s,ResVN);
	    TransposeBit_new(ResVN, ResNV);
	    DMultBit_new(ResNV, ResVN, ResNN);

// for V1i
	    VUBit_new(V, ResNN, ResnN);
	    for (i = 0; i < iceildiv(N, WBITS) * n; ++i)
		V->Data[i] = ResnN->Data[i];
// end for V1i

	    DMultBit_new(ViA2Vi, ResNN, ResNN2);
	    DMatrixSumBit_new(ResNN2, ViAVi, ResNN);
	    DMultBit_new(W, ResNN, ResNN2);
	    DMatrixSumBit_new(IN, ResNN2, D1i);
free(ResVN->Data);
free(ResNV->Data);
free(ResNN->Data);
free(ResNN2->Data);
free(ResnN->Data);
}



void ToConstructE(DenseMatrix E1i, DenseMatrix S, DenseMatrix W, DenseMatrix ViAVi, DenseMatrix ResNN, unsigned long s)
{
	    unsigned long N=ViAVi->Ncols;
            DenseMatrix ResVN,ResNV,ResNN2;
	    ResVN->Data = Allocmn(s, N);ResVN->Nrows=s;ResVN->Ncols=N;
	    ResNV->Data = Allocmn(N, s);ResNV->Ncols=s;ResVN->Nrows=N;
            ResNN2->Data = Allocmn(N, N);ResNN2->Nrows=N;ResNN2->Ncols=N;
	    SelectLinesBit_new(S, s, ResVN);
	    TransposeBit_new(ResVN, ResNV);
	    DMultBit_new(ResNV, ResVN, ResNN);
	    DMultBit_new(ViAVi, ResNN, ResNN2);
	    DMultBit_new(W, ResNN2, E1i);
free(ResVN->Data);
free(ResNV->Data);
free(ResNN2->Data);
}






void ToConstructF(DenseMatrix F1i, DenseMatrix S, DenseMatrix W[3], DenseMatrix Vi1AVi1, DenseMatrix Vi1A2Vi1,DenseMatrix IN, DenseMatrix Res, unsigned long *Set0, unsigned long *SizeS)
{
            unsigned long N=Vi1AVi1->Ncols;
            DenseMatrix ResVN,ResNV,ResVN2,ResNV2,ResNN,ResNN2,ResNN3;
	    ResVN->Data = Allocmn(Set0[0], N);ResVN->Nrows=Set0[0];ResVN->Ncols=N;
	    ResNV->Data = Allocmn(N, Set0[0]);ResNV->Ncols=Set0[0];ResVN->Nrows=N;         
	    ResVN2->Data = Allocmn(SizeS[1], N);ResVN2->Nrows=SizeS[1];ResVN2->Ncols=N;
	    ResNV2->Data = Allocmn(N, SizeS[1]);ResNV2->Ncols=SizeS[1];ResNV2->Nrows=N;
            ResNN->Data = Allocmn(N, N);ResNN->Nrows=N;ResNN->Ncols=N;
            ResNN2->Data = Allocmn(N, N);ResNN2->Nrows=N;ResNN2->Ncols=N;
            ResNN3->Data = Allocmn(N, N);ResNN3->Nrows=N;ResNN3->Ncols=N;
	    SelectLinesBit_new(S, SizeS[1], ResVN2);
	    TransposeBit_new(ResVN2, ResNV2);
	    DMultBit_new(ResNV2, ResVN2, ResNN);
	    DMultBit_new(Vi1A2Vi1, ResNN, ResNN2);
	    DMatrixSumBit_new(ResNN2, Vi1AVi1, ResNN);
	    DMultBit_new(ResNN, Res, ResNN3);
	    DMultBit_new(Vi1AVi1, W[1], ResNN);
	    DMatrixSumBit_new(IN,ResNN, ResNN2);
	    DMultBit_new(ResNN2, ResNN3, ResNN);
	    DMultBit_new(W[2], ResNN, F1i);

free(ResVN->Data);
free(ResNV->Data);
free(ResVN2->Data);
free(ResNV2->Data);


}



void ToConstructV(DenseMatrix V1i, DenseMatrix V[3], DenseMatrix D1i, DenseMatrix E1i, DenseMatrix F1i)
{
            unsigned long N=V1i->Ncols;
            unsigned long n=V1i->Nrows;
            DenseMatrix ResnN;ResnN->Data = Allocmn(n, N);ResnN->Nrows=n;ResnN->Ncols=N;
	    VUBit_new(V[0], D1i, ResnN);
	    DMatrixSumBit_new(V1i, ResnN, V1i);
	    VUBit_new(V[1], E1i, ResnN);
	    DMatrixSumBit_new(V1i, ResnN, V1i);
	    VUBit_new(V[2], F1i, ResnN);
	    DMatrixSumBit_new(V1i, ResnN, V1i);
free(ResnN->Data);
}





/*

The Lanczos iterative process 


*/



#if 1
void LanczosIterations_new(DenseMatrix Result, SparseMatrix M, DenseMatrix Y)
{

//    float SumTimeSD = 0;
//    float SumTimeTSD = 0;
    float SumTimeOper = 0;
    float SumTimeCom = 0;
//    float SumTimeComBcast = 0;
 //   float SumTimeComRed = 0;
    unsigned long i,Block;
    unsigned long SizeABlock;

    int size, p;
    unsigned long src, dest, tag, t0, t1,m,n;
    dest = 0;
    tag = 123;
    m = M->Nrows;
    n = M->Ncols;
    Block = Y->Ncols;
    t0 = microseconds();

    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks

// Displaying some information of whats happening

if (p==0){
	printf("\n");
	printf("Rank - It will increase until reaches Rank of Transposed sparse times sparse  \n");
	printf("STSVTime - Time for Transpose sparse times sparse matrix times block  \n");
	printf("TSLA - Time for small linear algebra   \n");
	printf("TI - Time per iteration\n");
	printf("ETT - Estimated total time to the end of process\n");
	printf("\n");
	printf("Detailed information will be displayed the first 50 steps and when completed 5 percent of total work\n");
}

    if (p == size - 1) {
	SizeABlock = m / size + (m % size);
    } else {
	SizeABlock = m / size;
    }

// To all Processes

    MPI_Op newop;
    MPI_Op_create((MPI_User_function *) MyXORfunction, 0, &newop);

    unsigned long N = Block;
    unsigned long *Final= malloc(sizeof(unsigned long));
    Final[0] = 0;
    unsigned long SizeRowV= iceildiv(N, WBITS);

    DenseMatrix ZnN;
    ZnN->Data = Allocmn(n, N);ZnN->Nrows=n;ZnN->Ncols=N;

    memset(ZnN->Data, 0, n * SizeRowV * sizeof(unsigned long));

    DenseMatrix X;
    X->Data = Allocmn(n, N);X->Nrows=n;X->Ncols=N;

    DenseMatrix V[3];

    for (i=0; i<3; i++){V[i]->Data=Allocmn(n,N);V[i]->Nrows=n;V[i]->Ncols=N;}

    DenseMatrix ATAC;
    ATAC->Data = Allocmn(n, N);

    STSMatrix_Vector(ATAC,M,Y);

    for (i = 0; i < n * iceildiv(Block, WBITS); ++i) {
	V[0]->Data[i] = ATAC->Data[i];
	V[1]->Data[i] = ZnN->Data[i];
	V[2]->Data[i] = ZnN->Data[i];
	X->Data[i] = ZnN->Data[i];
    }

    //initial conditions data prodution

    DenseMatrix IN;
    IN->Data = Allocmn(N, N);IN->Nrows=N;IN->Ncols=N;

    memset(IN->Data, 0, N * SizeRowV * sizeof(unsigned long));

    for (i = 0; i < N; ++i) {
	IN->Data[i * SizeRowV + i / WBITS] = (1UL << i % WBITS);
    }

    DenseMatrix ZN;
    ZN->Data = Allocmn(N, N);ZN->Nrows=n;ZN->Ncols=N;

    memset(ZN->Data, 0, N * SizeRowV * sizeof(unsigned long));

    DenseMatrix Winv[3],S[3];

    for (i=0; i<3; i++){
	Winv[i]->Data=Allocmn(N,N);Winv[i]->Nrows=N;Winv[i]->Ncols=N;
	S[i]->Data=Allocmn(N,N);S[i]->Nrows=N;S[i]->Ncols=N;}

    for (i = 0; i < N * SizeRowV; ++i) {
	Winv[0]->Data[i] = ZN->Data[i];
	Winv[1]->Data[i] = ZN->Data[i];
	Winv[2]->Data[i] = ZN->Data[i];
	S[0]->Data[i] = IN->Data[i];
	S[1]->Data[i] = IN->Data[i];
	S[2]->Data[i] = IN->Data[i];
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


    DenseMatrix ResNn,AVi,ViAVi,ViA2Vi,nW_inv;
    ResNn->Data = Allocmn(N, n);ResNn->Nrows=N;ResNn->Ncols=n;
    AVi->Data = Allocmn(n, N);AVi->Nrows=n;AVi->Ncols=N;
    ViAVi->Data = Allocmn(N, N);ViAVi->Nrows=N;ViAVi->Ncols=N;
    ViA2Vi->Data = Allocmn(N, N);ViA2Vi->Nrows=N;ViA2Vi->Ncols=N;
    nW_inv->Data = Allocmn(N, N);nW_inv->Nrows=N;nW_inv->Ncols=N;


    unsigned long *UseSet;
    UseSet = malloc((N + 1) * sizeof(unsigned long));

    UseSet[0] = N;
    for (i = 1; i <= N; ++i) {
	UseSet[i] = i - 1;
    }

    unsigned long *Set0;
    Set0 = malloc(sizeof(unsigned long));

    unsigned long *ntaken;
    ntaken = malloc((N + 1) * sizeof(unsigned long));
    DenseMatrix nSi;
    nSi->Data= Allocmn(N, N);nSi->Nrows=N;nSi->Ncols=N;


    STSMatrix_Vector(AVi,M,V[0]);


    if (p == 0) {


	TVUBit_new(V[0], AVi, ViAVi);
	TVUBit_new(AVi, AVi, ViA2Vi);
	Construction_S_WinvBit_new(ViAVi, UseSet, nSi, Set0, nW_inv, ntaken);
	Refresh_ArrayBit_new(Winv, nW_inv);
	Refresh_ArrayBit_new(S, nSi);

	SizeS[2] = SizeS[1];
	SizeS[1] = SizeS[0];
	SizeS[0] = ntaken[0];

	for (i = 0; i < N + 1; ++i) {
	    Set[2][i] = Set[1][i];
	    Set[1][i] = Set[0][i];
	    Set[0][i] = ntaken[i];

	}


	ToConstructX(X,ATAC,V[0],Winv[0]);
      
	printf("Step= %d Rank= %lu  Number of processes= %d \n", 0,SizeS[0], size);

    }

// To all processes

    unsigned long Step = 1;
    DenseMatrix AVi1,Vi1AVi1,Vi1A2Vi1;

    AVi1->Data = Allocmn(n, N);AVi1->Nrows=n;AVi1->Ncols=N;
    Vi1AVi1->Data = Allocmn(N, N);Vi1AVi1->Nrows=N;Vi1AVi1->Ncols=N;
    Vi1A2Vi1->Data = Allocmn(N, N);Vi1A2Vi1->Nrows=N;Vi1A2Vi1->Ncols=N;
  
    STSMatrix_Vector(AVi1,M,V[1]);

    if (p == 0) {

	TVUBit_new(V[1], AVi1, Vi1AVi1);
	TVUBit_new(AVi1, AVi1, Vi1A2Vi1);

    }

    free(AVi1->Data);

    unsigned long SumRank = SizeS[0];

//Alloc

    DenseMatrix V1i,D1i,E1i,F1i,ResNN,ResNN2,ResNN3,ResnN;
    V1i->Data  = Allocmn(n, N);V1i->Nrows=n;V1i->Ncols=N;
    D1i->Data = Allocmn(N, N);D1i->Nrows=N;D1i->Ncols=N;
    ResNN->Data = Allocmn(N, N);ResNN->Nrows=N;ResNN->Ncols=N;
    ResNN2->Data =Allocmn(N, N);ResNN2->Nrows=N;ResNN2->Ncols=N;
    E1i->Data = Allocmn(N, N);E1i->Nrows=N;E1i->Ncols=N;
    F1i->Data = Allocmn(N, N);F1i->Nrows=N;F1i->Ncols=N;
    ResNN3->Data = Allocmn(N, N);ResNN3->Nrows=N;ResNN3->Ncols=N;
    ResnN->Data = Allocmn(n, N);ResnN->Nrows=n;ResnN->Ncols=N;

    unsigned long T00, T0;

    float TTI=0;

    t1 = microseconds();

//the iteration

    while (((TestZero_new(ViAVi) == 0) || (p != dest)) && (Final[0] == 0)) {

	T00 = microseconds();

	if (p == 0) {

	    for (i = 1; i <= Set0[0]; ++i)
		{UseSet[i] = Set[0][i];}

	    UseSet[0] = Set0[0];

// V1i:=AVi*S`i*Transpose(S`i)+V`i*D1i+V`i1*E1i+V`i2*F1i;

	    for (i = 0; i < iceildiv(N, WBITS) * n; ++i)
		{V1i->Data[i] = AVi->Data[i];}

            DenseMatrix Res;Res->Data = Allocmn(N, N);Res->Nrows=N;Res->Ncols=N;

            ToConstructD(D1i,S[0],Winv[0],V1i,ViAVi,ViA2Vi,IN,Set0[0]);
	    ToConstructE(E1i,S[0],Winv[1],ViAVi,Res,Set0[0]);
	    ToConstructF(F1i,S[1],Winv,Vi1AVi1,Vi1A2Vi1,IN,Res,Set0,SizeS);
	    ToConstructV(V1i,V,D1i,E1i,F1i);
	    Refresh_ArrayBit_new(V, V1i);

            free(Res->Data);

	    for (i = 0; i < iceildiv(N, WBITS) * N; ++i) {
		Vi1AVi1->Data[i] = ViAVi->Data[i];
		Vi1A2Vi1->Data[i] = ViA2Vi->Data[i];
	    }
	}

// To All Processes

	unsigned long  T4, T5;
	T0 = microseconds();

        STSMatrix_Vector(AVi,M,V1i);

	T4 = microseconds();

	if (p == 0) {
	    SumTimeOper +=(T0 - T00);
	    SumTimeCom += (T4 - T0 );
	}
	if (p == 0) {
	    DenseMatrix ResNn;
	    ResNn->Data = Allocmn(N, n);ResNn->Ncols=n;ResNn->Nrows=N;
	    TVUBit_new(V1i, AVi, ViAVi);
	    TVUBit_new(AVi, AVi, ViA2Vi);
	    Construction_S_WinvBit_new(ViAVi, UseSet, nSi, Set0, nW_inv,
				   ntaken);
	    Refresh_ArrayBit_new(Winv, nW_inv);
	    Refresh_ArrayBit_new(S, nSi);
	    SizeS[2] = SizeS[1];
	    SizeS[1] = SizeS[0];
	    SizeS[0] = ntaken[0];

	    SumRank += SizeS[0];
	    Refresh_ArrayBit(N + 1, 1, Set, ntaken);

	    ToConstructX(X, ATAC,V[0],Winv[0]);

	    T5 = microseconds();

	    SumTimeOper += (T5 - T4);
	    float SD = T4 - T0;
	    float TIteration = T5 - T00;
	    float BCalc = T0 - T00;
	    float ACalc1 = T5 - T4;
            TTI+=TIteration;


if ((Step<=50) || ((Step>50) &&  (Step % (m/(SumRank/(Step+1)*20)))==0)){

            float SR=SumRank;
            float ST=Step+1;
            
            char *f=malloc(10*sizeof(unsigned long));
            TimeConvert(f,(m/(SR/ST)-(Step+1))*(TTI/(Step+1)/ 1000000));
	    printf
		("Step= %lu Rank= %lu  STSVTime=%.3fs TSLA=%.3fs TI=%.3fs  ETT=%s\n",
		 Step, SumRank, SD / 1000000,
		 (BCalc / 1000000) + (ACalc1 / 1000000),
		 TIteration / 1000000,f);
            free(f);
	    }
	    Step += 1;
	    free(ResNn->Data);
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
       
        printf("\n");
	printf("Total Computation out of Big Multiplication Time = %f s \n",SumTimeOper / 1000000);
        printf("Average Computation out of Big Multiplication Time = %f ms \n",SumTimeOper / (Step-1) / 1000);
        printf("\n");
        printf("\n");
        printf("Total Big Multiplications Time = %f s \n",SumTimeCom /1000000);
        printf("Average Big Multiplications Time = %f ms \n",SumTimeCom / (Step-1) / 1000);
	printf("\n");	
	printf("Total Preparation Time before iterative process = %lu ms\n",(t1 - t0) / 1000); 
        printf("Average Iteration Time = %.3f s \n",TTI/1000000/(Step-1));
        printf("\n");

	for (i = 0; i < iceildiv(Block, WBITS) * n; ++i) {
	    Result->Data[i] = X->Data[i] ^ Y->Data[i];
	}
    }


//  Free 

    free(X->Data);
    free(ATAC->Data);
    free(Final);
    free(ResNn->Data);
    free(ResNN->Data);
    free(ResNN2->Data);
    free(ResnN->Data);
    free(ResNN3->Data);
    free(V1i->Data);
    free(AVi->Data);
    free(ViAVi->Data);
    free(ViA2Vi->Data);
    free(D1i->Data);
    free(E1i->Data);
    free(F1i->Data);
    free(IN->Data);
    free(ZN->Data);
    free(ZnN->Data);
    free(nW_inv->Data);
    free(ntaken);
    free(nSi->Data);
    free(SizeS);
    free(UseSet);
    free(Set0);
    free(Vi1AVi1->Data);
    free(Vi1A2Vi1->Data);
    FreeAllocmn3_new(V);
    FreeAllocmn3_new(Winv);
    FreeAllocmn3(Set);
    FreeAllocmn3_new(S);

}

#endif






/*
"Small" Linear algebra to get a 64 block of vectors in the Kernel of MTM
*/




void KernelSparse_new(DenseMatrix Ker, SparseMatrix M, DenseMatrix R)
{
    
    unsigned long m = M->Nrows;
    unsigned long n = M->Ncols;
    unsigned long Block=R->Ncols;
    unsigned long Sizea;
    int size, p;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

// To all Processes

    MPI_Op newop;
    MPI_Op_create((MPI_User_function *) MyXORfunction, 0, &newop);

    unsigned long N = Block;
    unsigned long *ListLines, *LengListLines;
    
    DenseMatrix ResVn2, ResNN, E, Aout, ResnN, ResNn, ResNn2,
	 ResVn3, ResnV3;

    ResnV3->Data = Allocmn(n, Block);ResnV3->Nrows=n;ResnV3->Ncols=Block;
    ResNn->Data = Allocmn(Block, n);ResNn->Nrows=Block;ResNn->Ncols=n;
    ResNn2->Data = Allocmn(Block, n);ResNn2->Nrows=Block;ResNn2->Ncols=n;
    ResnN->Data = Allocmn(n, Block);ResnN->Nrows=n;ResnN->Ncols=Block;
    ResNN->Data = Allocmn(Block, Block);ResNN->Nrows=Block;ResNN->Ncols=Block;
    ResVn3->Data = Allocmn(Block, n);ResVn3->Nrows=Block;ResVn3->Ncols=n;
    Aout->Data = Allocmn(Block, n);Aout->Nrows=Block;Aout->Ncols=n;
    E->Data = Allocmn(Block, Block);E->Nrows=Block;E->Ncols=Block;
    ListLines = malloc((N + 1) * sizeof(unsigned long));
    LengListLines = malloc(sizeof(unsigned long));
    ResVn2->Data = Allocmn(Block, n);ResVn2->Nrows=Block;ResVn2->Ncols=n;

    unsigned long *NC= Allocmn(1, 1);

    DenseMatrix ATAR;
    ATAR->Nrows=M->Ncols;
    ATAR->Ncols=Block;
    ATAR->Data = Allocmn(ATAR->Nrows,ATAR->Ncols);
    
/* To Make Transpose(A)*A*R   */
  
   STSMatrix_Vector(ATAR,M,R);


/*  small gauss in Result to obtain independent vectors  */

    if (p == 0) {
	TransposeBit_new(ATAR, ResNn);
	GaussElimBit_new(ResNn, Aout, E, ListLines, LengListLines);
	TransposeBit_new(R, ResNn);
	DMultBit_new(E, ResNn, ResNn2);
	SelectLinesListBit_new(ResNn2, ListLines, Block - LengListLines[0],
			   ResVn2);
	GaussElimBit_new(ResVn2, Aout, E, ListLines, LengListLines);
	Sizea = LengListLines[0];
	TransposeBit_new(Aout, ResnV3);
        
    }

/* To Make A*R  */

    DenseMatrix d;
    d->Data = malloc(m * sizeof(unsigned long));d->Nrows=m;d->Ncols=Block;
    unsigned long i;

    SMatrix_Vector(d,M,ResnV3);



/* Analise the resulting matrix  */

if (p == 0) {

	if (TestZero_new(d)) {
	    for (i = 0; i < n; i++) {
		Ker->Data[i] = ResnV3->Data[i];
	    };
           
	    Ker->Nrows = n;
	    Ker->Ncols = Sizea;	    
	} 
	else {
	    
	    TransposeBit_new(d, ResVn3);
	    GaussElimBit_new(ResVn3, Aout, E, ListLines, LengListLines);
	    TransposeBit_new(ResnV3, ResNn);
	    DMultBit_new(E, ResNn, ResNn2);
	    for (i = LengListLines[0]; i < Sizea; i++) {
		ListLines[i] = i;
	    };
	    SelectLinesListBit_new(ResNn2, ListLines, Block, ResVn3);
	    GaussElimBit_new(ResVn3, Aout, E,
			 ListLines, LengListLines);
	    for (i = 0; i < LengListLines[0]; i++) {
		ListLines[i] = i;
	    };
	    SelectLinesListBit_new(Aout, ListLines, Block, ResVn3);
	    TransposeBit_new(ResVn3,Ker);
	    Ker->Nrows = n;
	    Ker->Ncols = LengListLines[0];
	}
        free(ResVn2->Data);
    }
// Free

	free(ATAR->Data);
	free(d->Data);
        free(NC);
	free(ResnV3->Data);
	free(ResVn3->Data);
	free(ResNN->Data);
	free(E->Data);
	free(Aout->Data);
	free(ResnN->Data);
	free(ResNn->Data);
	free(ResNn2->Data);
	free(ListLines);
	free(LengListLines);
 
    
}






/*

Complete Lanczos

*/

/*

The function Lanczos returns in "Kernel", "Block" independent vectors in the kernel of the sparse matrix M

*/


void Lanczos_new(DenseMatrix Kernel, SparseMatrix M,
	     unsigned long Block)
{

    int size, p;

    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks

    DenseMatrix Y;
    Y->Data = Allocmn(M->Ncols + 1, Block);

    unsigned long T0, T1,T2;

    if (p == 0) {
	RandomDMatrixBitTest_new(M->Ncols, Block, Y);
    }

    MPI_Bcast(Y->Data, M->Ncols * iceildiv(Block, WBITS), MPI_UNSIGNED_LONG, 0,
	      MPI_COMM_WORLD);
    Y->Ncols=Block;
    Y->Nrows=M->Ncols;

// The Lanczos iteration

    DenseMatrix Result;
    Result->Data = Allocmn(M->Ncols, Block);Result->Nrows=M->Ncols;Result->Ncols=Block;

    T0 = microseconds();

    LanczosIterations_new(Result,M, Y);

    T1 = microseconds();

// Small linear algebra to compute the  elements of the Kernel 

    Kernel->Nrows = M->Ncols;
    Kernel->Ncols = Block;

    KernelSparse_new(Kernel,M, Result);
    
    T2 = microseconds();

if (p==0) {
	float T10=(T1-T0);
	float T21=(T2-T1);
	printf("\n");
	printf("Total time for iterated process = %f s \n",T10/1000000);
	printf("Time for post processing data = %f s \n",T21/1000000);
	printf("\n");
}
    free(Y->Data);
    free(Result->Data);

}






///////////////////////////////////////////////////////////////////////////////
//
//
//
//
//
//
//
///////////////////////////////////////////////////////////////////////////////





#if 1      //Last working function


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

// Displaying some information of whats happening

if (p==0){
   printf("\n");
   printf("Rank - It will increase until reaches Rank of Transposed sparse times sparse  \n");
   printf("TSD - Time for sparse matrix times block  \n");
   printf("TTSD - Time for transposed sparse matrix times block   \n");
   printf("TSLA - Time for small linear algebra   \n");
   printf("TI - Time per iteration\n");
   printf("ETT - Estimated total time to the end of process\n");
   printf("\n");
   printf("Detailed information will be displayed the first 50 steps and when completed 5 percent of total work\n");
}




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


	TVUBit(n, N, V[0], ATAC, ViTV0);

	VUBit(n, N, V[0], Winv[0], ViWiinv);

	VUBit(n, N, ViWiinv, ViTV0, Xi);

	for (i = 0; i < iceildiv(Block, WBITS) * n; ++i) {
	    X[i] ^= Xi[i];
	}
	free(ViWiinv);
	free(ViTV0);
	free(Xi);
	free(ResNnX);


// end construct X0

	printf("Step= %d Rank= %lu  Number of processes= %d \n", 0,SizeS[0], size);

    }
// To all processes

    unsigned long Step = 1;
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

    unsigned long T00, T0;

    float TTI=0;

    t1 = microseconds();

//the iteration



    while (((TestZero(N, N, ViAVi) == 0) || (p != dest)) && (Final[0] == 0)) {


	T00 = microseconds();

	if (p == 0) {

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



	    VUBit(n, N, V[0], D1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);

	    VUBit(n, N, V[1], E1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);

	    VUBit(n, N, V[2], F1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);



	    Refresh_ArrayBit(n, N, V, V1i);



//T0 = microseconds();

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

	unsigned long T1, T2, T3, T4, T5;

	T0 = microseconds();
#if 1

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



	if (p != 0) {
	    Step += 1;
	}


	MPI_Reduce(AVi_Dist, AVi, iceildiv(Block, WBITS) * n,
		   MPI_UNSIGNED_LONG, newop, 0, MPI_COMM_WORLD);

#endif




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



	    TVUBit(n, N, V1i, AVi, ViAVi);

	    TVUBit(n, N, AVi, AVi, ViA2Vi);




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

//          unsigned long tnNN0 = microseconds();

	    TVUBit(n, N, V[0], ATAC, ViTV0);

//          unsigned long tnNN1 = microseconds();

	    VUBit(n, N, V[0], Winv[0], ViWiinv);

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
            TTI+=TIteration;
	    //float ACalc2 = T4 - T3;
	    //float tNnN=tnNN1-tNnN0;
	    //float tnNN = tnNN1 - tnNN0;


if ((Step<=50) || ((Step>50) &&  (Step % (m/(SumRank/(Step+1)*20)))==0)){
            float SR=SumRank;
            float ST=Step+1;
            char *f=malloc(10*sizeof(unsigned long));
            TimeConvert(f,(m/(SR/ST)-(Step+1))*(TTI/(Step+1)/ 1000000));
	    printf
		("Step= %lu Rank= %lu TSD=%.3fs TTSD=%.3fs TSLA=%.3fs TI=%.3fs TBcast=%.3fs TRed=%.3fs ETT=%s\n",
		 Step, SumRank, SD / 1000000, TSD / 1000000,
		 (BCalc / 1000000) + (ACalc1 / 1000000),
		 TIteration / 1000000,SumTimeComBcast/(1000000*(Step+1)),SumTimeComRed/(1000000*(Step+1)),f);
            free(f);
	    }
            
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
       

        printf("\n");
 //       printf("Computation times\n");
        printf("\n");
	printf("Total Computation Time = %f s \n",SumTimeOper / 1000000);
        printf("Average Computation Time = %f ms \n",SumTimeOper / (Step-1) / 1000);
        printf("\n");
  //      printf("Comunication times\n");
        printf("\n");
        printf("Total Comunication Time = %f s \n",SumTimeCom /1000000);
        printf("Average Comunication Time = %f ms \n",SumTimeCom / (Step-1) / 1000);
	printf("Average Broadcasting Time = %f ms \n",SumTimeComBcast / (Step-1) / 1000);
        printf("Average Reducing Time= %f ms  \n",SumTimeComRed / (Step-1) / 1000);

printf("\n");
//printf("Preparation and iterative process times\n");
printf("\n");

	printf("Total Preparation Time before iterative process = %lu ms\n",(t1 - t0) / 1000); 
        printf("Average Iteration Time = %.3f s \n",TTI/1000000/(Step-1));
        printf("\n");
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
    free(X);
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

    // free(f1);

    return 0;
}

#endif








void KernelSparse(SparseMatrix M, unsigned long *R,
		  unsigned long Block, DenseMatrix Ker)
{
    unsigned long * a = M->Data;
    unsigned long m = M->Nrows;
    unsigned long n = M->Ncols;

    unsigned long Sizea;
    int size, p;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    MPI_Status status;


    unsigned long SizeABlock = M->slices[p]->i1 - M->slices[p]->i0;



// To all Processes



    MPI_Op newop;
    MPI_Op_create((MPI_User_function *) MyXORfunction, 0, &newop);

    unsigned long N = Block;
    unsigned long SizeRowV;

    SizeRowV = iceildiv(N, WBITS);

    
    unsigned long *ATAR_Dist;
    unsigned long *ResmN_Dist;
    unsigned long *ResVn2, *ResNN, *E, *Aout, *ResnN, *ResNn, *ResNn2,
	*ListLines, *LengListLines, *ResVn3;
    unsigned long *ResnV3;


    ResnV3 = Allocmn(n, Block);
    ResNn = Allocmn(Block, n);
    ResNn2 = Allocmn(Block, n);
    ResnN = Allocmn(n, Block);
    ResNN = Allocmn(Block, Block);
    ResVn3 = Allocmn(Block, n);
    Aout = Allocmn(Block, n);
    E = Allocmn(Block, Block);
    ListLines = malloc((N + 1) * sizeof(unsigned long));
    LengListLines = malloc(sizeof(unsigned long));

    ResVn2 = Allocmn(Block, n);


    unsigned long *NC;
    NC = Allocmn(1, 1);

// To Make Transpose(A)*A*R



    unsigned long *ATAR;
    ATAR = Allocmn(n, N);
    

    ResmN_Dist = Allocmn(m, N);
    ATAR_Dist = Allocmn(n, N);

     

    SMultDmatrixBit(m, n, N, a, R, ResmN_Dist, M->slices[p]->i0, SizeABlock);
    TSMultDmatrixBit(SizeABlock, n, N, a, ResmN_Dist, ATAR_Dist,
		     M->slices[p]->i0, SizeABlock);
    MPI_Reduce(ATAR_Dist, ATAR, iceildiv(Block, WBITS) * n, MPI_UNSIGNED_LONG,
	       newop, 0, MPI_COMM_WORLD);

    free(ResmN_Dist);
    free(ATAR_Dist);


   


    if (p == 0) {
	TransposeBit(n, N, ATAR, ResNn);
	GaussElimBit(Block, n, ResNn, Aout, E, ListLines, LengListLines);
	TransposeBit(n, N, R, ResNn);
	DMultBit(N, N, n, E, ResNn, ResNn2);


/*  small gauss in Result to obtain independent vectors  */

	
	SelectLinesListBit(n, ResNn2, ListLines, Block - LengListLines[0],
			   ResVn2);
	unsigned long Size1 = Block - LengListLines[0];
	GaussElimBit(Size1, n, ResVn2, Aout, E, ListLines, LengListLines);
	Sizea = LengListLines[0];
	TransposeBit(LengListLines[0], n, Aout, ResnV3);
        
    }

    

    //To make A*R


    unsigned long *d;
    d = malloc(m * sizeof(unsigned long));

    unsigned long *ResmN_Dist1;
   ResmN_Dist1 = Allocmn(m, Block);


    unsigned long j, i;

    if (p==0) {
           
    }

    MPI_Bcast(ResnV3, n, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);


    MPI_Bcast(NC, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);


    SMultDmatrixBit(m, n, Block, a, ResnV3, ResmN_Dist1, M->slices[p]->i0,
		    SizeABlock);


    if (p != 0) {
	MPI_Send(ResmN_Dist1, SizeABlock, MPI_UNSIGNED_LONG, 0, 1,
		 MPI_COMM_WORLD);

  };

    

    if (p == 0) {
	for (i = 0; i < SizeABlock ; i++) {
	    d[i] = ResmN_Dist1[i];
	}

	for (i = 1; i < size; i++) {

	    MPI_Recv(ResmN_Dist1,
                    M->slices[i]->i1 - M->slices[i]->i0,
                    MPI_UNSIGNED_LONG, i,
		     1, MPI_COMM_WORLD, &status);

	    for (j = M->slices[i]->i0; j < M->slices[i]->i1 ; j++) {
		d[j] = ResmN_Dist1[j - M->slices[i]->i0];
	    }
	}

      

	if (TestZero(m, Block, d)) {
	    for (i = 0; i < n; i++) {
		Ker->Data[i] = ResnV3[i];
	    };
           
	    Ker->Nrows = n;
	    Ker->Ncols = Sizea;
	    //printf("\n No second Small Lin Alg  \n");
	} else {
	    //printf("\n With second Small Lin Alg  \n");

	    TransposeBit(m, Block, d, ResVn3);

            

	    GaussElimBit(Block, m, ResVn3, Aout, E, ListLines, LengListLines);

	    TransposeBit(n, Block, ResnV3, ResNn);

	    DMultBit(Block, Block, n, E, ResNn, ResNn2);

	    for (i = LengListLines[0]; i < Sizea; i++) {
		ListLines[i] = i;
	    };


	    SelectLinesListBit(n, ResNn2, ListLines, Block, ResVn3);
	    GaussElimBit(Sizea - LengListLines[0], n, ResVn3, Aout, E,
			 ListLines, LengListLines);
	    for (i = 0; i < LengListLines[0]; i++) {
		ListLines[i] = i;
	    };
	    SelectLinesListBit(n, Aout, ListLines, Block, ResVn3);
	    TransposeBit(LengListLines[0], n, ResVn3, Ker->Data);
	    Ker->Nrows = n;
	    Ker->Ncols = LengListLines[0];
           

	}
        
        //free(d);
        free(ResVn2);
       

    }


// See where to put them here is good for one proc

free(ResmN_Dist1);
free(ATAR);
free(d);
       
        free(NC);
	free(ResnV3);
	free(ResVn3);
	free(ResNN);
	free(E);
	free(Aout);
	free(ResnN);
	free(ResNn);
	free(ResNn2);
	free(ListLines);
	free(LengListLines);
 
    
}














#if 0

void Lanczos(DenseMatrix Kernel, SparseMatrix M,
	     unsigned long Block)
{

    int size, p;



    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks

    //DenseMatrix Kernel;
    //Kernel->Data = Allocmn(M->Ncols, Block);

    DenseMatrix Y;
    Y->Data = Allocmn(M->Ncols + 1, Block);

    unsigned long T0, T1,T2;

    if (p == 0) {
	RandomDMatrixBitTest_new(M->Ncols, Block, Y);
    }

    MPI_Bcast(Y->Data, M->Ncols * iceildiv(Block, WBITS), MPI_UNSIGNED_LONG, 0,
	      MPI_COMM_WORLD);


// The Lanczos iteration

    DenseMatrix Result;
    // , *Index;
    Result->Data = Allocmn(M->Ncols, Block);
    // Index = malloc(2 * sizeof(unsigned long));

    T0 = microseconds();

   

 //   LanczosIterations(M->Data, Y, M->Nrows, M->Ncols, Block, Result);

     LanczosIterations(Result,M, Y);

    

    T1 = microseconds();


// Small linear algebra to compute the  elements of the Kernel 

    Kernel->Nrows = M->Ncols;
    Kernel->Ncols = Block;

    KernelSparse(Kernel,M, Result);


   
    // MPI_Bcast(Index, 2, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    
    T2 = microseconds();

if (p==0) {
    float T10=(T1-T0);
    float T21=(T2-T1);
printf("\n");
printf("Total time for iterated process = %f s \n",T10/1000000);
printf("Time for post processing data = %f s \n",T21/1000000);
printf("\n");
}
//    displayMatrixScreen(Kernel->Data,M->Ncols,Kernel->Ncols);
//    displaySMatrixNew(M->Data,M->Nrows,M->Ncols,'M');

    free(Y);
    free(Result);

}




#endif











#if 1
unsigned long Lanczos(SparseMatrix M, unsigned long Block,
		      DenseMatrix Kernel)
{

    int size, p;



    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks


    unsigned long *Y;
    Y = Allocmn(M->Ncols + 1, Block);

    unsigned long T0, T1;

    if (p == 0) {
	RandomDMatrixBitTest(M->Ncols, Block, Y);
    }

    MPI_Bcast(Y, M->Ncols * iceildiv(Block, WBITS), MPI_UNSIGNED_LONG, 0,
	      MPI_COMM_WORLD);


// The Lanczos iteration

    unsigned long *Result, *Index;
    Result = Allocmn(M->Ncols, Block);
    Index = malloc(2 * sizeof(unsigned long));

    T0 = microseconds();

    LanczosIterations(M->Data,Y, M->Nrows, M->Ncols, Block, Result);

    T1 = microseconds();


// Small linear algebra to compute the  elements of the Kernel 


    KernelSparse(M, Result, Block, Kernel);
    Kernel->Nrows = M->Ncols;
    Kernel->Ncols = Block;

    //displayMatrix(Kernel->Data,Index[0],Index[1],'d');
    //displaySMatrixNew(M->Data,M->Nrows,M->Ncols,'M');

    free(Y);
    free(Result);


    return 0;

}
#endif









#if  0

void *LanczosIterations(DenseMatrix Resultado, SparseMatrix M,
			DenseMatrix Y)
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





//initial conditions data prodution



    unsigned long *ZnN;
    ZnN = Allocmn(n, N);

    memset(ZnN, 0, n * SizeRowV * sizeof(unsigned long));

    unsigned long *X;
    X = Allocmn(n, N);


    unsigned long **V;
    V = Allocmn3(n, N);



/*
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

    MPI_Bcast(ATAC, n * iceildiv(Block,WBITS), MPI_UNSIGNED_LONG, 0,
	      MPI_COMM_WORLD);

    free(ResmN_Dist);
    free(ATAC_Dist);
*/

    DenseMatrix ATAC;
    ATAC->Data = Allocmn(a->Ncols, Y->Ncols);

    STSMatrix_Vector(ATAC, a, Y);


    for (i = 0; i < n * iceildiv(Block, WBITS); ++i) {
	V[0][i] = ATAC->Data[i];
	V[1][i] = ZnN[i];
	V[2][i] = ZnN[i];
	X[i] = ZnN[i];
    }


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

// End initial conditions data prodution



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

///// Aqui

    SMultDmatrixBit(m, n, N, a, V[0], ResmN_Dist, p * (m / size), SizeABlock);
    TSMultDmatrixBit(SizeABlock, n, N, a, ResmN_Dist, AVi_Dist,
		     p * (m / size), SizeABlock);


    MPI_Reduce(AVi_Dist, AVi, iceildiv(Block, WBITS) * n, MPI_UNSIGNED_LONG,
	       newop, 0, MPI_COMM_WORLD);


    if (p == 0) {


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


	TVUBit(n, N, V[0], ATAC, ViTV0);

	VUBit(n, N, V[0], Winv[0], ViWiinv);

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

    unsigned long T00, T0;

    t1 = microseconds();

//the iteration



    while (((TestZero(N, N, ViAVi) == 0) || (p != dest)) && (Final[0] == 0)) {


	T00 = microseconds();

	if (p == 0) {

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



	    VUBit(n, N, V[0], D1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);

	    VUBit(n, N, V[1], E1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);

	    VUBit(n, N, V[2], F1i, ResnN);

	    DMatrixSumBit(n, N, V1i, ResnN, V1i);



	    Refresh_ArrayBit(n, N, V, V1i);



//T0 = microseconds();

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

	unsigned long T1, T2, T3, T4, T5;

	T0 = microseconds();
#if 1

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



	if (p != 0) {
	    Step += 1;
	}


	MPI_Reduce(AVi_Dist, AVi, iceildiv(Block, WBITS) * n,
		   MPI_UNSIGNED_LONG, newop, 0, MPI_COMM_WORLD);

#endif




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

	    TVUBit(n, N, V1i, AVi, ViAVi);

	    TVUBit(n, N, AVi, AVi, ViA2Vi);

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



//          unsigned long tnNN0 = microseconds();

	    TVUBit(n, N, V[0], ATAC, ViTV0);

//          unsigned long tnNN1 = microseconds();

	    VUBit(n, N, V[0], Winv[0], ViWiinv);

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
	    //float ACalc2 = T4 - T3;
	    //float tNnN=tnNN1-tNnN0;
	    //float tnNN = tnNN1 - tnNN0;

	    printf
		("Step= %lu      SumRank= %lu    TimeSD=%f s  TimeTSD=%f s   TSLAlgebra=%f s   TITime=%f s \n",
		 Step, SumRank, SD / 1000000, TSD / 1000000,
		 (BCalc / 1000000) + (ACalc1 / 1000000),
		 TIteration / 1000000);

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




#endif


/*

    MPI_Bcast(ResnV3->Data, n, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);


    MPI_Bcast(NC, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);


    SMultDmatrixBit(m, n, Block, M->Data, ResnV3->Data, ResmN_Dist1->Data, M->slices[p]->i0,
		    SizeABlock);


    if (p != 0) {
	MPI_Send(ResmN_Dist1->Data, SizeABlock, MPI_UNSIGNED_LONG, 0, 1,
		 MPI_COMM_WORLD);

  };

    

    if (p == 0) {
	for (i = 0; i < SizeABlock ; i++) {
	    d->Data[i] = ResmN_Dist1->Data[i];
	}

	for (i = 1; i < size; i++) {

	    MPI_Recv(ResmN_Dist1->Data,
                    M->slices[i]->i1 - M->slices[i]->i0,
                    MPI_UNSIGNED_LONG, i,
		     1, MPI_COMM_WORLD, &status);

	    for (j = M->slices[i]->i0; j < M->slices[i]->i1 ; j++) {
		d->Data[j] = ResmN_Dist1->Data[j - M->slices[i]->i0];
	    }
	}

*/      

