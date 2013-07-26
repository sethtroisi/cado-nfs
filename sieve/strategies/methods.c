#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"

static void complain(const char *msg)
{
  fprintf(stderr, "Error while reading method: %s\n", msg);
  abort();
}

int method_read_stream(cofac_method_ptr meth, FILE *file)
{
#define MAXLINE 1000
  char str[MAXLINE];
  char *ret;

  // Parse first line
  {
    ret = fgets(str, MAXLINE, file);
    if (ret == NULL) complain("could not read first line");

    char type[MAXLINE];
    int B1, B2, param;
    int n;
    n = sscanf(str, "METHOD=%s B1=%d B2=%d param=%d", type, &B1, &B2, &param);
    if (n != 4) complain("first line in wrong");

    if (strcmp(type, "ECM") == 0)
      meth->type = ECM;
    else if (strcmp(type, "ECMM12") == 0)
      meth->type = ECMM12;
    else if (strcmp(type, "PM1") == 0)
      meth->type = PM1;
    else if (strcmp(type, "PP1_27") == 0)
      meth->type = PP1_27;
    else if (strcmp(type, "PP1_65") == 0)
      meth->type = PP1_65;
    else
      complain("unknown method");
    meth->B1 = B1;
    meth->B2 = B2;
    meth->param = param;
  }

  // Parse benchs
  {
    ret = fgets(str, MAXLINE, file);
    if (ret == NULL) complain("could not read second line");
    char *tmp;
    tmp = strstr(str, "Bench:");
    if (tmp == NULL) complain("second line is not Bench:");

    for (int i = 1; i < 4; ++i) {
      int ii;
      float tm;
      ret = fgets(str, MAXLINE, file);
      if (ret == NULL) complain("could not read line of bench");
      int n = sscanf(str, "%d words: %f", &ii, &tm);
      if (n != 2 || ii != i) complain("problem when parsing bench");
      meth->ms[i-1] = tm;
    }
  }

  // Parse probas
  {
    ret = fgets(str, MAXLINE, file);
    if (ret == NULL) complain("could not read Proba header line");
    char *tmp;
    tmp = strstr(str, "Success Probability");
    if (tmp == NULL) complain("proba header line is wrong");
    for (int i = 15; i <= METHOD_MAXBITS; ++i) {
      int ii;
      float tm1, tm5, tm7, tm11;
      ret = fgets(str, MAXLINE, file);
      if (ret == NULL) complain("could not read line of proba");
      int n = sscanf(str, "%d: %f %f %f %f", &ii, &tm1, &tm5, &tm7, &tm11);
      if (n != 5 || ii != i) complain("problem when parsing proba line");
      meth->success[MOD12_1][i] = tm1;
      meth->success[MOD12_5][i] = tm5;
      meth->success[MOD12_7][i] = tm7;
      meth->success[MOD12_11][i] = tm11;
    }
  }

  return 1;
#undef MAXLINE
}

int methods_read(cofac_method_t **meths, const char *f)
{
  cofac_method_t * list_meth = NULL;
  int nm = 0;
  FILE *file;
  file = fopen(f, "r");

  int c;
  while ((c = getc(file)) != EOF) {
    ungetc(c, file);
    list_meth = (cofac_method_t *)realloc(list_meth,
        (nm+1)*sizeof(cofac_method_t));
    method_read_stream(list_meth[nm], file);
    nm++;
  }
  fclose(file);
  fprintf(stderr, "Successfully read %d methods.\n", nm);
  *meths = list_meth;
  return nm;
}

void method_print_full(cofac_method_srcptr meth, FILE *file)
{
  fprintf(file, "METHOD=");
  if (meth->type == ECM)
    fprintf(file, "ECM");
  else if (meth->type == ECMM12)
    fprintf(file, "ECMM12");
  else if (meth->type == PM1)
    fprintf(file, "PM1");
  else if (meth->type == PP1_27)
    fprintf(file, "PP1_27");
  else if (meth->type == PP1_65)
    fprintf(file, "PP1_65");
  else
    abort();
  fprintf(file, " B1=%d B2=%d param=%d\n  Bench:\n",
      meth->B1, meth->B2, meth->param);
  for (int i = 0; i < 3; ++i)
    fprintf(file, "    %d words: %.3f\n", i+1, meth->ms[i]);
  fprintf(file, "  Success Probability for p of given bit size for 1,5,7,11 mod 12\n");
  for (int i = 15; i <= METHOD_MAXBITS; ++i)
    fprintf(file, "    %d: %.2f %.2f %.2f %.2f\n", i,
        meth->success[MOD12_1][i],
        meth->success[MOD12_5][i],
        meth->success[MOD12_7][i],
        meth->success[MOD12_11][i]);
}

void method_print(cofac_method_srcptr meth, FILE *file)
{
  if (meth->type == ECM)
    fprintf(file, "ECM");
  else if (meth->type == ECMM12)
    fprintf(file, "ECMM12");
  else if (meth->type == PM1)
    fprintf(file, "PM1");
  else if (meth->type == PP1_27)
    fprintf(file, "PP1_27");
  else if (meth->type == PP1_65)
    fprintf(file, "PP1_65");
  else
    abort();
  fprintf(file, " %d %d %d\n", meth->B1, meth->B2, meth->param);
}

#if 0

int main(int argc, char **argv) {
  cofac_method_t * m;
  int nm;
  nm = methods_read(&m, argv[1]);

  for (int i = 0; i < nm; ++i)
    method_print(m[i], stdout);

  free(m);
}

#endif
