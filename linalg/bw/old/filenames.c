#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "params.h"

#include <assert.h>

#include "filenames.h"

/*
 * The version that uses an auxiliary C<...> file to store the dot
 * products is in the 20010526 tarball.
 */
int current_bank;
char wdir_filename[FILENAME_LENGTH];
char wip_tag_filename[FILENAME_LENGTH];
char w_matrix_filename[FILENAME_LENGTH];
char w_indexes_filename[FILENAME_LENGTH];
char w_values_filename[FILENAME_LENGTH];
char a_meta_filename[FILENAME_LENGTH];
char a_sub_meta_filename[FILENAME_LENGTH];
char s_meta_filename[FILENAME_LENGTH];
char ds_meta_filename[FILENAME_LENGTH];
char v_meta_filename[FILENAME_LENGTH];
char f_meta_filename[FILENAME_LENGTH];
char f_sub_meta_filename[FILENAME_LENGTH];
char f_base_filename[FILENAME_LENGTH];
char pi_meta_filename[FILENAME_LENGTH];
char h_meta_filename[FILENAME_LENGTH];
char w_filename[FILENAME_LENGTH];
char l_meta_filename[FILENAME_LENGTH];
char x_meta_filename[FILENAME_LENGTH];
char y_meta_filename[FILENAME_LENGTH];
char z_meta_filename[FILENAME_LENGTH];
char valu_meta_filename[FILENAME_LENGTH];
char i_matrix_filename[FILENAME_LENGTH];
char i_vector_filename[FILENAME_LENGTH];
char i_result_filename[FILENAME_LENGTH];
char pre_stats_filename[FILENAME_LENGTH];
char mvec_filename[FILENAME_LENGTH];
char x0_filename[FILENAME_LENGTH];
char certif_meta_filename[FILENAME_LENGTH];
const char modulus_info_meta_filename[] = "input/%s";

#define wdir_PRIMNAME		"output/%d"
#define wip_tag_PRIMNAME	"output/%d/WORKING"
#define w_matrix_PRIMNAME	"output/%d/matrix"
#define w_indexes_PRIMNAME	"output/%d/indexes"
#define w_values_PRIMNAME	"output/%d/values"
#define a_meta_PRIMNAME		"output/%d/A-%%02d-%%02d.%%03d"
#define a_sub_meta_PRIMNAME	"output/%d/sA-%%02d-%%02d.%%03d-%%03d"
#define s_meta_PRIMNAME		"output/%d/S-%%02d.%%03d"
#define ds_meta_PRIMNAME	"output/%d/ndS-%%02d.%%03d-%%03d"
#define v_meta_PRIMNAME		"output/%d/V-%%02d.%%03d"
#define f_meta_PRIMNAME		"output/%d/rF-%%02d-%%02d.%%04d"
#define f_sub_meta_PRIMNAME	"output/%d/srF-%%02d-%%02d.%%04d.%%d-%%d"
#define f_base_PRIMNAME		"output/%d/F_INIT"
#define pi_meta_PRIMNAME	"output/%d/P-%%d-%%d"
#define h_meta_PRIMNAME		"output/%d/H.%%02d"
#define l_meta_PRIMNAME		"output/%d/LOCAL_INFO-%%02d"
#define x_meta_PRIMNAME		"output/%d/X%%02d"
#define w_PRIMNAME		"output/%d/W"
#define y_meta_PRIMNAME		"output/%d/Y%%02d"
#define z_meta_PRIMNAME		"output/%d/Z%%02d"
#define valu_meta_PRIMNAME	"output/%d/VALUATION-%%02d.%%04d"
#define i_matrix_PRIMNAME	"input/matrix.%d"
#define i_vector_PRIMNAME	"input/vector.%d"
#define i_result_PRIMNAME	"input/result.%d"
#define pre_stats_PRIMNAME	"output/%d/FIRST.INFO"
#define mvec_PRIMNAME		"output/%d/M-vector"
#define x0_PRIMNAME		"output/%d/X0-vector"
#define certif_meta_PRIMNAME	"output/%d/CERTIF-%%d-%%d-%%d"

static int in_run_wd=0;

void set_all_filenames(int banknum)
{
	if (in_run_wd == 0) {
		if (chdir("run") < 0) {
			perror("chdir() to run wd");
			exit(1);
		}
		in_run_wd = 1;
	}
	current_bank = banknum;
#ifndef HAS_NOT_SNPRINTF
#define SINGLE_STRING(n,f)	snprintf(n,FILENAME_LENGTH,f,banknum)
#else
#define SINGLE_STRING(n,f)	sprintf(n,f,banknum)
#endif

	SINGLE_STRING(wdir_filename,wdir_PRIMNAME);
	SINGLE_STRING(wip_tag_filename,wip_tag_PRIMNAME);
	SINGLE_STRING(w_matrix_filename,w_matrix_PRIMNAME);
	SINGLE_STRING(w_indexes_filename,w_indexes_PRIMNAME);
	SINGLE_STRING(w_values_filename,w_values_PRIMNAME);
	SINGLE_STRING(a_meta_filename,a_meta_PRIMNAME);
	SINGLE_STRING(a_sub_meta_filename,a_sub_meta_PRIMNAME);
	SINGLE_STRING(s_meta_filename,s_meta_PRIMNAME);
	SINGLE_STRING(ds_meta_filename,ds_meta_PRIMNAME);
	SINGLE_STRING(v_meta_filename,v_meta_PRIMNAME);
	SINGLE_STRING(f_meta_filename,f_meta_PRIMNAME);
	SINGLE_STRING(f_sub_meta_filename,f_sub_meta_PRIMNAME);
	SINGLE_STRING(f_base_filename,f_base_PRIMNAME);
	SINGLE_STRING(pi_meta_filename,pi_meta_PRIMNAME);
	SINGLE_STRING(h_meta_filename,h_meta_PRIMNAME);
	SINGLE_STRING(w_filename,w_PRIMNAME);
	SINGLE_STRING(l_meta_filename,l_meta_PRIMNAME);
	SINGLE_STRING(x_meta_filename,x_meta_PRIMNAME);
	SINGLE_STRING(y_meta_filename,y_meta_PRIMNAME);
	SINGLE_STRING(z_meta_filename,z_meta_PRIMNAME);
	SINGLE_STRING(valu_meta_filename,valu_meta_PRIMNAME);
	SINGLE_STRING(i_matrix_filename,i_matrix_PRIMNAME);
	SINGLE_STRING(i_vector_filename,i_vector_PRIMNAME);
	SINGLE_STRING(i_result_filename,i_result_PRIMNAME);
	SINGLE_STRING(pre_stats_filename,pre_stats_PRIMNAME);
	SINGLE_STRING(mvec_filename,mvec_PRIMNAME);
	SINGLE_STRING(x0_filename,x0_PRIMNAME);
	SINGLE_STRING(certif_meta_filename,certif_meta_PRIMNAME);
}
