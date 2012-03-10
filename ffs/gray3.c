// Copy-paste from:
//   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.119.1344

#include <stdio.h>
#define N 3
#define K 4
void printgray3()
{
    int n[K+1];
    int g[K+1];
    int u[K+1];
    int i, j, k;
    /* the maximum for each digit */
    /* the Gray code */
    /* +1 or −1 */
    for (i=0; i<=K; i++) {g[i]=0; u[i]=1; n[i]=N;};
    while (g[K]==0) {
        printf("(");
        for (j=K-1; j>=0; j--) printf (" %d", g[j]);
        printf (")\n");
        i=0;
        /* enumerate next Gray code */
        k=g[0]+u[0];
        while ((k>=n[i]) || (k<0)) {
            u[i]=-u[i];
            i++;
            k=g[i]+u[i];
        };
        g[i]=k;
    };
}

// Recursive version, to highlight the changing digit.

int C[4] = {0,0,0,0};
void printC() {
    printf("( %d %d %d %d )", C[3], C[2], C[1], C[0]);
}

void gray3_rec(int k, int way) {
    if (k == 0)
        return;
    gray3_rec(k-1, way);
    // C[k-1] += way;
    if (way>0) { C[k-1]++; printf ("+"); }
    else { C[k-1]--; printf ("-"); }
    printf("%d ", k-1); printC(); printf("\n");
    gray3_rec(k-1, -way);
    // C[k-1] += way;
    if (way>0) { C[k-1]++; printf ("+"); }
    else { C[k-1]--; printf ("-"); }
    printf("%d ", k-1); printC(); printf("\n");
    gray3_rec(k-1, way);
}

#ifdef GRAY3_MAIN
int main() {
    printC(); printf("\n");
    rec(4, 1);
}
#endif
