extern void mpi_send_inactive_rows(int i);
extern void mpi_add_rows(sparse_mat_t *mat, int m, INT j, INT *ind);
extern void mpi_load_rows_for_j(sparse_mat_t *mat, int m, INT j);
