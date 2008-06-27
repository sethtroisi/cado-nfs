extern void mpi_start_proc(char *outname, int mpi_rank, int mpi_size, sparse_mat_t *mat, FILE *purgedfile, char *purgedname, int forbw, double ration, int coverNmax);
extern void mpi_send_inactive_rows(int i);
extern void mpi_add_rows(sparse_mat_t *mat, int m, INT j, INT *ind);
extern void mpi_load_rows_for_j(sparse_mat_t *mat, int m, INT j);
