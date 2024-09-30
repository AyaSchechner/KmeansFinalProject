#ifndef SYMNMF_H
#define SYMNMF_H

void free_mat(double ** matrix, int rows_num);
double ** create_mat(int vec_num, int vec_len);
double ** create_mat2(int vec_num, int vec_len);
double sqr_euclidean_dis(double * vec1, double * vec2, int vec_len);
double ** matrix_multiplication(double ** mat1, double ** mat2, int n);
double ** symmat(double ** vec_mat, int n, int d);
double **ddgmat(double **vec_mat, int n, int d);
double ** normmat(double ** vec_mat, int n, int d);
double ** matrix_multiplication2(double ** mat1, double ** mat2, int n, int k);
double ** matrix_multiplication3(double ** mat1, double ** mat2, int n, int k);
double ** transpose_matrix(double ** mat, int n, int k);
double f_norm(double ** mat, int n, int k);
double ** matrix_sub(double ** mat1, double ** mat2, int n, int k);
double ** finalmat(double ** init_mat, double ** norm_mat, int n, int k);
void print_mat(double ** vectors_mat, int vec_len, int vec_num);


#endif
