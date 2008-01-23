#ifndef FILENAMES_H_
#define FILENAMES_H_

#ifdef	__cplusplus
extern "C" {
#endif

/*
 * The version that uses an auxiliary C<...> file to store the dot
 * products is in the 20010526 tarball.
 */
extern int current_bank;
extern char wdir_filename[FILENAME_LENGTH];
extern char wip_tag_filename[FILENAME_LENGTH];
extern char w_matrix_filename[FILENAME_LENGTH];
extern char w_indexes_filename[FILENAME_LENGTH];
extern char w_values_filename[FILENAME_LENGTH];
extern char a_meta_filename[FILENAME_LENGTH];
extern char a_sub_meta_filename[FILENAME_LENGTH];
extern char s_meta_filename[FILENAME_LENGTH];
extern char ds_meta_filename[FILENAME_LENGTH];
extern char v_meta_filename[FILENAME_LENGTH];
extern char f_sub_meta_filename[FILENAME_LENGTH];
extern char f_meta_filename[FILENAME_LENGTH];
extern char f_base_filename[FILENAME_LENGTH];
extern char pi_meta_filename[FILENAME_LENGTH];
extern char h_meta_filename[FILENAME_LENGTH];
extern char w_filename[FILENAME_LENGTH];
extern char l_meta_filename[FILENAME_LENGTH];
extern char x_meta_filename[FILENAME_LENGTH];
extern char y_meta_filename[FILENAME_LENGTH];
extern char z_meta_filename[FILENAME_LENGTH];
extern char valu_meta_filename[FILENAME_LENGTH];
extern char i_matrix_filename[FILENAME_LENGTH];
extern char i_vector_filename[FILENAME_LENGTH];
extern char i_result_filename[FILENAME_LENGTH];
extern char pre_stats_filename[FILENAME_LENGTH];
extern char mvec_filename[FILENAME_LENGTH];
extern char x0_filename[FILENAME_LENGTH];
extern char certif_meta_filename[FILENAME_LENGTH];
extern const char modulus_info_meta_filename[];

extern void set_all_filenames(int);

#ifdef	__cplusplus
}
#endif

#endif /* FILENAMES_H_ */
