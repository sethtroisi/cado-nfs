#include "mat_int64.h"
#include "timing.h"

int main()
{
  /*//[[0, 2, 3],[4, 5, 6],[7, 8, 10]]*/
  /*mat_Z_t A_Z;*/
  /*mat_Z_init(A_Z, 3, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 1, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 2, 1, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 3, 1, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 4, 2, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 5, 2, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 6, 2, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 7, 3, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 8, 3, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 10, 3, 3);*/
  /*mat_Z_t C_Z;*/
  /*mat_Z_init(C_Z, 3, 3);*/
  /*mat_Z_LLL(C_Z, A_Z);*/
  /*mat_Z_clear(A_Z);*/

  /*mat_int64_t A;*/
  /*mat_int64_init(A, 3, 3);*/
  /*mat_int64_set_identity(A);*/
  /*A->coeff[1][1] = 0;*/
  /*A->coeff[1][2] = 2;*/
  /*A->coeff[1][3] = 3;*/
  /*A->coeff[2][1] = 4;*/
  /*A->coeff[2][2] = 5;*/
  /*A->coeff[2][3] = 6;*/
  /*A->coeff[3][1] = 7;*/
  /*A->coeff[3][2] = 8;*/
  /*A->coeff[3][3] = 10;*/

  /*mat_int64_t C;*/
  /*mat_int64_init(C, 3, 3);*/
  /*mat_int64_set_zero(C);*/
  /*mat_int64_LLL(C, A);*/

  /*mat_Z_to_mat_int64(A, C_Z);*/

  /*ASSERT_ALWAYS(mat_int64_equal(A, C));*/

  /*mat_int64_clear(A);*/
  /*mat_int64_clear(C);*/
  /*mat_Z_clear(C_Z);*/

  /*//[[1, 0, 0],[0, 1, 0],[0, 0, 1]]*/
  /*mat_Z_init(A_Z, 3, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 1, 1, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 1, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 1, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 2, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 1, 2, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 2, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 3, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 3, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 1, 3, 3);*/
  /*mat_Z_init(C_Z, 3, 3);*/
  /*mat_Z_LLL(C_Z, A_Z);*/
  /*mat_Z_clear(A_Z);*/

  /*mat_int64_init(A, 3, 3);*/
  /*mat_int64_set_identity(A);*/
  /*A->coeff[1][1] = 1;*/
  /*A->coeff[1][2] = 0;*/
  /*A->coeff[1][3] = 0;*/
  /*A->coeff[2][1] = 0;*/
  /*A->coeff[2][2] = 1;*/
  /*A->coeff[2][3] = 0;*/
  /*A->coeff[3][1] = 0;*/
  /*A->coeff[3][2] = 0;*/
  /*A->coeff[3][3] = 1;*/

  /*mat_int64_init(C, 3, 3);*/
  /*mat_int64_set_zero(C);*/
  /*mat_int64_LLL(C, A);*/

  /*mat_Z_to_mat_int64(A, C_Z);*/

  /*ASSERT_ALWAYS(mat_int64_equal(A, C));*/

  /*mat_int64_clear(A);*/
  /*mat_int64_clear(C);*/
  /*mat_Z_clear(C_Z);*/

  /*//[[-1, -1, -3, 7, 0],[-1, 0, -1, -1, 1],[0, 0, 2, 1, -9],[-1, -1, -1, 7, 1],[0, -2, 0, 0, 8]]*/
  /*mat_Z_init(A_Z, 5, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, -1, 1, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, -1, 1, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, -3, 1, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 7, 1, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 1, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, -1, 2, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 2, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, -1, 2, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, -1, 2, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 1, 2, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 3, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 3, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 2, 3, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 1, 3, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, -9, 3, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, -1, 4, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, -1, 4, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, -1, 4, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 7, 4, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 1, 4, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 5, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, -2, 5, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 5, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 0, 5, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 8, 5, 5);*/

  /*mat_Z_init(C_Z, 5, 5);*/
  /*mat_Z_LLL(C_Z, A_Z);*/
  /*mat_Z_clear(A_Z);*/
  /*mat_int64_init(A, 5, 5);*/
  /*mat_int64_set_identity(A);*/
  /*A->coeff[1][1] = -1;*/
  /*A->coeff[1][2] = -1;*/
  /*A->coeff[1][3] = -3;*/
  /*A->coeff[1][4] = 7;*/
  /*A->coeff[1][5] = 0;*/
  /*A->coeff[2][1] = -1;*/
  /*A->coeff[2][2] = 0;*/
  /*A->coeff[2][3] = -1;*/
  /*A->coeff[2][4] = -1;*/
  /*A->coeff[2][5] = 1;*/
  /*A->coeff[3][1] = 0;*/
  /*A->coeff[3][2] = 0;*/
  /*A->coeff[3][3] = 2;*/
  /*A->coeff[3][4] = 1;*/
  /*A->coeff[3][5] = -9;*/
  /*A->coeff[4][1] = -1;*/
  /*A->coeff[4][2] = -1;*/
  /*A->coeff[4][3] = -1;*/
  /*A->coeff[4][4] = 7;*/
  /*A->coeff[4][5] = 1;*/
  /*A->coeff[5][1] = 0;*/
  /*A->coeff[5][2] = -2;*/
  /*A->coeff[5][3] = 0;*/
  /*A->coeff[5][4] = 0;*/
  /*A->coeff[5][5] = 8;*/

  /*mat_int64_init(C, 5, 5);*/
  /*mat_int64_set_zero(C);*/
  /*mat_int64_LLL(C, A);*/
  /*mat_Z_to_mat_int64(A, C_Z);*/
  /*ASSERT_ALWAYS(mat_int64_equal(A, C));*/
  /*mat_int64_clear(A);*/
  /*mat_int64_clear(C);*/
  /*mat_Z_clear(C_Z);*/

  /*//[[42, 32, 32, 42, 48],[37, 46, 52, 35, 43],[36, 52, 61, 37, 60],[55, 55, 51, 54, 62],[55, 44, 58, 37, 51]]*/
  /*mat_Z_init(A_Z, 5, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, 42, 1, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 32, 1, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 32, 1, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 42, 1, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 48, 1, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, 37, 2, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 46, 2, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 52, 2, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 35, 2, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 43, 2, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, 36, 3, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 52, 3, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 61, 3, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 37, 3, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 60, 3, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, 55, 4, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 55, 4, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 51, 4, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 54, 4, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 62, 4, 5);*/
  /*mat_Z_set_coeff_int64(A_Z, 55, 5, 1);*/
  /*mat_Z_set_coeff_int64(A_Z, 44, 5, 2);*/
  /*mat_Z_set_coeff_int64(A_Z, 58, 5, 3);*/
  /*mat_Z_set_coeff_int64(A_Z, 37, 5, 4);*/
  /*mat_Z_set_coeff_int64(A_Z, 51, 5, 5);*/

  /*mat_Z_init(C_Z, 5, 5);*/
  /*mat_Z_LLL(C_Z, A_Z);*/
  /*printf("Output of mat_Z:\n");*/
  /*mat_Z_fprintf(stdout, C_Z);*/
  /*mat_Z_clear(A_Z);*/
  /*mat_int64_init(A, 5, 5);*/
  /*mat_int64_set_identity(A);*/
  /*A->coeff[1][1] = 42;*/
  /*A->coeff[1][2] = 32;*/
  /*A->coeff[1][3] = 32;*/
  /*A->coeff[1][4] = 42;*/
  /*A->coeff[1][5] = 48;*/
  /*A->coeff[2][1] = 37;*/
  /*A->coeff[2][2] = 46;*/
  /*A->coeff[2][3] = 52;*/
  /*A->coeff[2][4] = 35;*/
  /*A->coeff[2][5] = 43;*/
  /*A->coeff[3][1] = 36;*/
  /*A->coeff[3][2] = 52;*/
  /*A->coeff[3][3] = 61;*/
  /*A->coeff[3][4] = 37;*/
  /*A->coeff[3][5] = 60;*/
  /*A->coeff[4][1] = 55;*/
  /*A->coeff[4][2] = 55;*/
  /*A->coeff[4][3] = 51;*/
  /*A->coeff[4][4] = 54;*/
  /*A->coeff[4][5] = 62;*/
  /*A->coeff[5][1] = 55;*/
  /*A->coeff[5][2] = 44;*/
  /*A->coeff[5][3] = 58;*/
  /*A->coeff[5][4] = 37;*/
  /*A->coeff[5][5] = 51;*/

  /*mat_int64_init(C, 5, 5);*/
  /*mat_int64_set_zero(C);*/
  /*mat_int64_LLL(C, A);*/
  /*printf("Output of mat_int64:\n");*/
  /*mat_int64_fprintf(stdout, C);*/
  /*mat_int64_clear(A);*/
  /*mat_int64_clear(C);*/
  /*mat_Z_clear(C_Z);*/

  mat_Z_t A_Z, C_Z;
  mat_int64_t A, C;

  //[[52, 63, 48, 48, 57, 53, 49, 58],[46, 34, 40, 34, 53, 36, 32, 60],[42, 39, 49, 39, 34, 49, 46, 36],[39, 50, 60, 57, 39, 58, 39, 51],[32, 36, 50, 59, 43, 42, 35, 50],[42, 38, 52, 41, 37, 51, 52, 63],[58, 32, 63, 58, 47, 62, 41, 46],[41, 49, 60, 47, 42, 63, 39, 59]]
  mat_Z_init(A_Z, 8, 8);
  mat_Z_set_coeff_int64(A_Z, 52, 1, 1);
  mat_Z_set_coeff_int64(A_Z, 63, 1, 2);
  mat_Z_set_coeff_int64(A_Z, 48, 1, 3);
  mat_Z_set_coeff_int64(A_Z, 48, 1, 4);
  mat_Z_set_coeff_int64(A_Z, 57, 1, 5);
  mat_Z_set_coeff_int64(A_Z, 53, 1, 6);
  mat_Z_set_coeff_int64(A_Z, 49, 1, 7);
  mat_Z_set_coeff_int64(A_Z, 58, 1, 8);
  mat_Z_set_coeff_int64(A_Z, 46, 2, 1);
  mat_Z_set_coeff_int64(A_Z, 34, 2, 2);
  mat_Z_set_coeff_int64(A_Z, 40, 2, 3);
  mat_Z_set_coeff_int64(A_Z, 34, 2, 4);
  mat_Z_set_coeff_int64(A_Z, 53, 2, 5);
  mat_Z_set_coeff_int64(A_Z, 36, 2, 6);
  mat_Z_set_coeff_int64(A_Z, 32, 2, 7);
  mat_Z_set_coeff_int64(A_Z, 60, 2, 8);
  mat_Z_set_coeff_int64(A_Z, 42, 3, 1);
  mat_Z_set_coeff_int64(A_Z, 39, 3, 2);
  mat_Z_set_coeff_int64(A_Z, 49, 3, 3);
  mat_Z_set_coeff_int64(A_Z, 39, 3, 4);
  mat_Z_set_coeff_int64(A_Z, 34, 3, 5);
  mat_Z_set_coeff_int64(A_Z, 49, 3, 6);
  mat_Z_set_coeff_int64(A_Z, 46, 3, 7);
  mat_Z_set_coeff_int64(A_Z, 36, 3, 8);
  mat_Z_set_coeff_int64(A_Z, 39, 4, 1);
  mat_Z_set_coeff_int64(A_Z, 50, 4, 2);
  mat_Z_set_coeff_int64(A_Z, 60, 4, 3);
  mat_Z_set_coeff_int64(A_Z, 57, 4, 4);
  mat_Z_set_coeff_int64(A_Z, 39, 4, 5);
  mat_Z_set_coeff_int64(A_Z, 58, 4, 6);
  mat_Z_set_coeff_int64(A_Z, 39, 4, 7);
  mat_Z_set_coeff_int64(A_Z, 51, 4, 8);
  mat_Z_set_coeff_int64(A_Z, 32, 5, 1);
  mat_Z_set_coeff_int64(A_Z, 36, 5, 2);
  mat_Z_set_coeff_int64(A_Z, 50, 5, 3);
  mat_Z_set_coeff_int64(A_Z, 59, 5, 4);
  mat_Z_set_coeff_int64(A_Z, 43, 5, 5);
  mat_Z_set_coeff_int64(A_Z, 42, 5, 6);
  mat_Z_set_coeff_int64(A_Z, 35, 5, 7);
  mat_Z_set_coeff_int64(A_Z, 50, 5, 8);
  mat_Z_set_coeff_int64(A_Z, 42, 6, 1);
  mat_Z_set_coeff_int64(A_Z, 38, 6, 2);
  mat_Z_set_coeff_int64(A_Z, 52, 6, 3);
  mat_Z_set_coeff_int64(A_Z, 41, 6, 4);
  mat_Z_set_coeff_int64(A_Z, 37, 6, 5);
  mat_Z_set_coeff_int64(A_Z, 51, 6, 6);
  mat_Z_set_coeff_int64(A_Z, 52, 6, 7);
  mat_Z_set_coeff_int64(A_Z, 63, 6, 8);
  mat_Z_set_coeff_int64(A_Z, 58, 7, 1);
  mat_Z_set_coeff_int64(A_Z, 32, 7, 2);
  mat_Z_set_coeff_int64(A_Z, 63, 7, 3);
  mat_Z_set_coeff_int64(A_Z, 58, 7, 4);
  mat_Z_set_coeff_int64(A_Z, 47, 7, 5);
  mat_Z_set_coeff_int64(A_Z, 62, 7, 6);
  mat_Z_set_coeff_int64(A_Z, 41, 7, 7);
  mat_Z_set_coeff_int64(A_Z, 46, 7, 8);
  mat_Z_set_coeff_int64(A_Z, 41, 8, 1);
  mat_Z_set_coeff_int64(A_Z, 49, 8, 2);
  mat_Z_set_coeff_int64(A_Z, 60, 8, 3);
  mat_Z_set_coeff_int64(A_Z, 47, 8, 4);
  mat_Z_set_coeff_int64(A_Z, 42, 8, 5);
  mat_Z_set_coeff_int64(A_Z, 63, 8, 6);
  mat_Z_set_coeff_int64(A_Z, 39, 8, 7);
  mat_Z_set_coeff_int64(A_Z, 59, 8, 8);

  mat_Z_init(C_Z, 8, 8);
  mat_Z_LLL(C_Z, A_Z);
  printf("Output of mat_Z:\n");
  mat_Z_fprintf(stdout, C_Z);
  mat_Z_clear(A_Z);
  mat_int64_init(A, 8, 8);
  mat_int64_set_identity(A);
  A->coeff[1][1] = 52;
  A->coeff[1][2] = 63;
  A->coeff[1][3] = 48;
  A->coeff[1][4] = 48;
  A->coeff[1][5] = 57;
  A->coeff[1][6] = 53;
  A->coeff[1][7] = 49;
  A->coeff[1][8] = 58;
  A->coeff[2][1] = 46;
  A->coeff[2][2] = 34;
  A->coeff[2][3] = 40;
  A->coeff[2][4] = 34;
  A->coeff[2][5] = 53;
  A->coeff[2][6] = 36;
  A->coeff[2][7] = 32;
  A->coeff[2][8] = 60;
  A->coeff[3][1] = 42;
  A->coeff[3][2] = 39;
  A->coeff[3][3] = 49;
  A->coeff[3][4] = 39;
  A->coeff[3][5] = 34;
  A->coeff[3][6] = 49;
  A->coeff[3][7] = 46;
  A->coeff[3][8] = 36;
  A->coeff[4][1] = 39;
  A->coeff[4][2] = 50;
  A->coeff[4][3] = 60;
  A->coeff[4][4] = 57;
  A->coeff[4][5] = 39;
  A->coeff[4][6] = 58;
  A->coeff[4][7] = 39;
  A->coeff[4][8] = 51;
  A->coeff[5][1] = 32;
  A->coeff[5][2] = 36;
  A->coeff[5][3] = 50;
  A->coeff[5][4] = 59;
  A->coeff[5][5] = 43;
  A->coeff[5][6] = 42;
  A->coeff[5][7] = 35;
  A->coeff[5][8] = 50;
  A->coeff[6][1] = 42;
  A->coeff[6][2] = 38;
  A->coeff[6][3] = 52;
  A->coeff[6][4] = 41;
  A->coeff[6][5] = 37;
  A->coeff[6][6] = 51;
  A->coeff[6][7] = 52;
  A->coeff[6][8] = 63;
  A->coeff[7][1] = 58;
  A->coeff[7][2] = 32;
  A->coeff[7][3] = 63;
  A->coeff[7][4] = 58;
  A->coeff[7][5] = 47;
  A->coeff[7][6] = 62;
  A->coeff[7][7] = 41;
  A->coeff[7][8] = 46;
  A->coeff[8][1] = 41;
  A->coeff[8][2] = 49;
  A->coeff[8][3] = 60;
  A->coeff[8][4] = 47;
  A->coeff[8][5] = 42;
  A->coeff[8][6] = 63;
  A->coeff[8][7] = 39;
  A->coeff[8][8] = 59;

  printf("New?\n");
  getchar();

  mat_int64_init(C, 8, 8);
  mat_int64_set_zero(C);
  mat_int64_LLL(C, A);
  printf("Output of mat_int64:\n");
  mat_int64_fprintf(stdout, C);
  mat_int64_clear(A);
  mat_int64_clear(C);
  mat_Z_clear(C_Z);
}
