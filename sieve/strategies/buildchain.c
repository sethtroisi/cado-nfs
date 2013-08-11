#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "methods.h"

#define MAX(a,b) ((a)<=(b))?(b):(a)



float prior_checksum(histogram_t *P) {
    float p = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int b = 15; b <= METHOD_MAXBITS; ++b)
            p += P[i][b];
    return p;
}

// TODO: put an appropriate value, here.
// ATM, we assume that probability of having a factor of b bits is 1/b,
// and we divide by 4, assuming the repartition in classes mod 12 is
// uniform.
void get_prior(histogram_t *P) {
    for (int i = 0; i < 4; ++i)
        for (int b = 15; b <= METHOD_MAXBITS; ++b)
            P[i][b] = 1/((float)(4*b));
}

// Return max bit-size b for which 95% of the primes of b bit have been
// removed.
int cleaned_primes(histogram_t *P) {
    histogram_t prior[4];
    get_prior(prior);
    for (int b = 15; b <= METHOD_MAXBITS; ++b)
        for (int i = 1; i < 4; ++i)
            prior[0][b] += prior[i][b];

    int res = -1;
    float PP;
    for (int b = 15; b <= METHOD_MAXBITS; ++b) {
        PP = 0;
        for (int i = 0; i < 4; ++i)
            PP += P[i][b];
        PP /= prior[0][b];
        if (PP < 0.05)
            res = b;
    }
    return res;
}


// History of which p-1 / p+1 where already run.
// These are histogram of success proba according to class mod 12.
// They are only pointers to the corresponding histograms of existing
// methods (should not involve allocations).
typedef struct {
    float * pm1[4];
    float * pp1[4];
    float * ecm11[4];
} ppm1_history_struct;

// Evaluate the score of running a method, knowing the histogram of
// probability of having a prime of given bit size for each class mod 12.
// The formula is - Log(Proba_of_failure) / time.
// We actually compute the score for the 3 possible input sizes: 1, 2 and
// 3 words.
void method_score(float *score, histogram_t *success, float *ms,
        histogram_t *prob)
{
    float p = 1;
    for (int b = 15; b <= METHOD_MAXBITS; ++b) {  // for each bitsize
        for (int i = 0; i < 4; ++i)               // for each class
            p *= 1 - prob[i][b] * success[i][b];
    }
    p = 1-p;
    assert (p >= 0 && p <= 1);
    p = -logf(1-p);
    for (int i = 0; i < 3; ++i)
        score[i] = p/ms[i];
}


// Given:
//   - prob, a bit-histogram of probability of occurence of prime
//   - succes_prob, a bit-histogram of probability of finding such a prime
//     with a method
// Compute:
//   - the new probability of occurence of prime for each bit-size,
//     assuming the method failed.
// Result is stored in place in prob.
//
// Formula (Bayesian thm):
//     p := p*prob(failure) / ( 1 - p*prob(success))
// NB: this formula seems to be correct, because updating twice with the
// same failure probability is equivalent to updating once with the
// squared failure probability.
void update_prob(histogram_t pr, histogram_t succ_pr)
{
    // I'm always afraid of numerical instability...
    float eps = 0.00001;
    for (int i = 15; i <= METHOD_MAXBITS; ++i) {
        assert (pr[i] >= 0 && pr[i] <= 1);
        assert (succ_pr[i] >= 0 && succ_pr[i] <= 1);
        float num = pr[i]*(1-succ_pr[i]);
        float den = 1 - pr[i]*succ_pr[i];
        if (num < eps)
            pr[i] = 0;
        else if (den < eps)
            pr[i] = 1;
        else
            pr[i] = num/den;
        assert (pr[i] >= 0 && pr[i] <= 1);
    }
}


// Given a history for p-1 and p+1, the type of the method and its
// success probability, compute the new probability of success, obtained
// by subtracting the probabilities for things that have already been run
// (and failed!)
void take_ppm1_history_into_account(histogram_t *P, int type,
        histogram_t *p, ppm1_history_struct *history)
{
    // ECM does not redo p+1 or p-1. Just copy p to P.
    if (type == ECMM12) {
        for (int i = 0; i < 4; ++i)
            for (int j = 15; j <= METHOD_MAXBITS; ++j)
                P[i][j] = p[i][j];
        return;
    }
    // ECM11 is just one curve.
    if (type == ECM11) {
        for (int i = 0; i < 4; ++i)
            for (int j = 15; j <= METHOD_MAXBITS; ++j)
                P[i][j] = p[i][j] - history->ecm11[i][j];
    }
    // p-1 redo p-1. That's also easy: just a subtract.
    else if (type == PM1) {
        for (int i = 0; i < 4; ++i)
            for (int j = 15; j <= METHOD_MAXBITS; ++j)
                P[i][j] = p[i][j] - history->pm1[i][j];
    }
    // for p+1 with 2/7 as a starting point, we do
    //   p-1 for p = 1, 7  mod 12
    //   p+1 for p = 5, 11 mod 12
    else if (type == PP1_27) {
        for (int j = 15; j <= METHOD_MAXBITS; ++j) {
            P[MOD12_1][j]  = p[MOD12_1][j]  - history->pm1[MOD12_1][j];
            P[MOD12_5][j]  = p[MOD12_5][j]  - history->pp1[MOD12_5][j];
            P[MOD12_7][j]  = p[MOD12_7][j]  - history->pm1[MOD12_7][j];
            P[MOD12_11][j] = p[MOD12_11][j] - history->pp1[MOD12_11][j];
        }
    }
    // for p+1 with 6/5 as a starting point, we do
    //   p-1 for p = 1, 5  mod 12
    //   p+1 for p = 7, 11 mod 12
    else if (type == PP1_65) {
        for (int j = 15; j <= METHOD_MAXBITS; ++j) {
            P[MOD12_1][j]  = p[MOD12_1][j]  - history->pm1[MOD12_1][j];
            P[MOD12_5][j]  = p[MOD12_5][j]  - history->pm1[MOD12_5][j];
            P[MOD12_7][j]  = p[MOD12_7][j]  - history->pp1[MOD12_7][j];
            P[MOD12_11][j] = p[MOD12_11][j] - history->pp1[MOD12_11][j];
        }
    }
    else
        abort();
    // It could be that this subtraction put us in the negative part.
    for (int i = 0; i < 4; ++i)
            for (int j = 15; j <= METHOD_MAXBITS; ++j)
                P[i][j] = MAX(0, P[i][j]);
}

float * max_prob(float *p1, float *p2) {
    float s1, s2;
    s1 = 0;
    s2 = 0;
    for(int i = 15; i <= METHOD_MAXBITS; ++i) {
        s1 += p1[i];
        s2 += p2[i];
    }
    return (s1 < s2) ? p2 : p1;
}

// So, we have run the method (and it failed).
// The ppm1 history must be updated accordingly.
void update_ppm1_history(ppm1_history_struct *hist, cofac_method_ptr meth)
{
  if (meth->type == ECMM12)
    return;
  if (meth->type == ECM11) {
    for (int i = 0; i < 4; ++i)
      hist->ecm11[i] = max_prob(hist->ecm11[i], meth->success[i]);
    return;
  }
  if (meth->type == PM1) {
    for (int i = 0; i < 4; ++i)
      hist->pm1[i] = max_prob(hist->pm1[i], meth->success[i]);
    return;
  }
  if (meth->type == PP1_27) {
    hist->pm1[MOD12_1] = max_prob(hist->pm1[MOD12_1], meth->success[MOD12_1]);
    hist->pp1[MOD12_5] = max_prob(hist->pp1[MOD12_5], meth->success[MOD12_5]);
    hist->pm1[MOD12_7] = max_prob(hist->pm1[MOD12_7], meth->success[MOD12_7]);
    hist->pp1[MOD12_11]= max_prob(hist->pp1[MOD12_11],meth->success[MOD12_11]);
    return;
  }
  if (meth->type == PP1_65) {
    hist->pm1[MOD12_1] = max_prob(hist->pm1[MOD12_1], meth->success[MOD12_1]);
    hist->pm1[MOD12_5] = max_prob(hist->pm1[MOD12_5], meth->success[MOD12_5]);
    hist->pp1[MOD12_7] = max_prob(hist->pp1[MOD12_7], meth->success[MOD12_7]);
    hist->pp1[MOD12_11]= max_prob(hist->pp1[MOD12_11],meth->success[MOD12_11]);
    return;
  }
}

int
get_best_method(cofac_method_t * list_meth, int nm, histogram_t *prior,
        ppm1_history_struct *ppm1_history, int wordsize, float *sco)
{
    int best_meth = -1;
    float best_score = 0;

    // current success probability for the given method.
    histogram_t P[4];
    float sc[3];

    for (int i = 0; i < nm; ++i) {
        take_ppm1_history_into_account(P, list_meth[i]->type,
                list_meth[i]->success, ppm1_history);
        method_score(sc, P, list_meth[i]->ms, prior);
        float score = sc[wordsize-1];
    //    fprintf(stderr, "    Score of method %d: %f\n", i, score);
        if ((best_meth == -1) || (score > best_score)) {
            best_meth = i;
            best_score = score;
        }
    }
    *sco = best_score;
    return best_meth;
}

// usage: ./buildchain <method_filename> <word_size> <chain_length>

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr, "usage: ./buildchain <method_filename> " 
                "<word_size> <chain_length>\n");
        return EXIT_FAILURE;
    }
    cofac_method_t *list_meth;
    int nm;
    int wordsize = atoi(argv[2]);
    int len = atoi(argv[3]);
    int chain[len];
    nm = methods_read(&list_meth, argv[1]);

    // The current knowledge about probabilities of occurence of a prime.
    histogram_t P[4];
    get_prior(P);

    // The history of p+1 and p-1. Initialize with 0.
    ppm1_history_struct ppm1_hist;
    float zero_hist[METHOD_MAXBITS+1] = {0,}; // fill in with zeros.
    for (int i = 0; i < 4; ++i) {
        ppm1_hist.pm1[i] = zero_hist;
        ppm1_hist.pp1[i] = zero_hist;
        ppm1_hist.ecm11[i] = zero_hist;
    }

    for (int i = 0; i < len; ++i) {
        float sc;
        fprintf(stderr, "Choosing method for step %d...\n", i);
        int m = get_best_method(list_meth, nm, P, &ppm1_hist, wordsize, &sc);
        fprintf(stderr, "  The winner for step %d is: "
                "method number %d with score %.3f\n", i, m, sc);
        chain[i] = m;

        // Update the prior
        histogram_t PP[4];
        take_ppm1_history_into_account(PP, list_meth[m]->type,
                list_meth[m]->success, &ppm1_hist);
        for (int j = 0; j < 4; ++j)
            update_prob(P[j], PP[j]);

        // Update the ppm1 history
        update_ppm1_history(&ppm1_hist, list_meth[m]);

        fprintf(stderr, "  Current prior checksum = %.4f\n", prior_checksum(P));
        fprintf(stderr, "  Bit size that has been cleaned = %d\n",
                cleaned_primes(P));
    }

    for (int i = 0; i < len; ++i) {
        printf("%d: ", i);
        method_print(list_meth[chain[i]], stdout);
    }

    free(list_meth);
    return EXIT_SUCCESS;
}
