#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


struct cord{
    double value;
    struct cord *next;
};


struct vector{
    struct vector *next;
    struct cord *cords;
};


double eps = 0.0001;
int max_iter = 300;
double beta = 0.5;


void free_mat(double **, int);
void free_LinkedList(struct vector *, int);
double ** create_mat(int, int);
double ** create_mat2(int, int);
double sqr_euclidean_dis(double *, double *, int);
double ** matrix_multiplication(double ** , double ** , int);
double ** symmat(double **, int, int);
double **ddgmat(double **, int, int);
double ** normmat(double **, int, int);
double ** matrix_multiplication2(double **, double **, int, int);
double ** matrix_multiplication3(double **, double **, int, int);
double ** transpose_matrix(double **, int, int);
double f_norm(double **, int, int);
double ** matrix_sub(double **, double **, int, int);
double ** finalmat(double **, double **, int, int);
int vectors_num(struct vector *);
int vectors_len(struct vector *);
double ** insert_vectors(struct vector *, int, int);
double ** create_init_mat(struct vector *, int, int);
void print_mat(double **, int, int);


void free_mat(double ** matrix, int rows_num) {
    int i;
    for (i = 0; i < rows_num; i++) {
        free(matrix[i]);
    }
    free(matrix);
}


void free_LinkedList(struct vector *head_vec, int d) {
    struct vector *curr_vec = head_vec;
    struct vector *next_vec;
    int a = 0;

    while (a < d) {
        struct cord *curr_cord = curr_vec->cords;
        struct cord *next_cord;

        while (curr_cord != NULL) {
            next_cord = curr_cord->next;
            free(curr_cord);
            curr_cord = next_cord;
            a++;
        }

        next_vec = curr_vec->next;
        free(curr_vec);
        curr_vec = next_vec;
    }
    free(curr_vec);
}


/* create a vec_num*vec_len matrix with malloc */
double ** create_mat(int vec_num, int vec_len){
    double ** vectors_mat;
    int i = 0;
    vectors_mat = (double **) malloc(vec_num * sizeof(double *));
    if(vectors_mat == NULL) {
        printf("An Error Has Occurred");
        return NULL;
    }

    while (i < vec_num){
        vectors_mat[i] = (double *) malloc(vec_len * sizeof(double));
        if(vectors_mat[i] == NULL){
            printf("An Error Has Occurred");
            free_mat(vectors_mat, i);
            return NULL;
        }
        i++;
    }
    return vectors_mat;
}


/* create a vec_num*vec_len matrix with calloc */
double ** create_mat2(int vec_num, int vec_len){
    double ** vectors_mat;
    int i = 0;
    vectors_mat = (double **) malloc(vec_num * sizeof(double *));
    if(vectors_mat == NULL) {
        printf("An Error Has Occurred");
        return NULL;
    }

    while (i < vec_num){
        vectors_mat[i] = (double *) calloc(vec_len, sizeof(double));
        if(vectors_mat[i] == NULL){
            printf("An Error Has Occurred");
            free_mat(vectors_mat, i);
            return NULL;
        }
        i++;
    }
    return vectors_mat;
}


/* compute  the square Euclidean distance of two vectors */
double sqr_euclidean_dis(double * vec1, double * vec2, int vec_len) {
    int i;
    double dis = 0.0;
    for (i = 0; i < vec_len; i++) {
        dis += (((vec1[i]) - (vec2[i])) * ((vec1[i]) - (vec2[i])));
    }
    return dis;
}


/* multiplying n*n matrices */
double ** matrix_multiplication(double ** mat1, double ** mat2, int n){
    double ** mat;
    int i;
    int j;
    int k;

    mat = create_mat2(n, n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                mat[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return mat;
}


/* computing the Similarity Matrix */
double ** symmat(double ** vec_mat, int n, int d) {
    double **sym_mat;
    int i;
    int j;
    sym_mat = create_mat(n, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                sym_mat[i][j] = 0;
            } else {
                sym_mat[i][j] = exp(-(sqr_euclidean_dis(vec_mat[i], vec_mat[j], d) / 2));
            }
        }
    }
    return sym_mat;
}


/* computing the diagonal degree Matrix */
double **ddgmat(double **vec_mat, int n, int d) {
    int i;
    int j;
    double d_i;
    double **ddg_mat;
    double **sym_mat = symmat(vec_mat, n, d);

    ddg_mat = create_mat2(n, n);

    for (i = 0; i < n; i++) {
        d_i = 0;
        for (j = 0; j < n; j++) {
            d_i += sym_mat[i][j];
        }
        ddg_mat[i][i] = d_i;
    }
    return ddg_mat;
}


/* computing the normalized similarity Matrix */
double ** normmat(double ** vec_mat, int n, int d){
    int i;
    double **norm_mat;
    double **sym_mat = symmat(vec_mat, n, d);
    double **ddg_mat = ddgmat(vec_mat, n, d);

    for (i = 0; i < n; i++) {
        ddg_mat[i][i] = 1/(sqrt(ddg_mat[i][i])); /* DDGMAT^(-1/2) */
    }
    norm_mat = matrix_multiplication(ddg_mat, sym_mat, n);
    norm_mat = matrix_multiplication(norm_mat, ddg_mat, n);

    free_mat(sym_mat, n);
    free_mat(ddg_mat, n);
    return norm_mat;
}


/* code for finding H: */


/* multiplying n*k matrix by k*n matrix */
double ** matrix_multiplication2(double ** mat1, double ** mat2, int n, int k) {
    int i;
    int j;
    int x;
    double **result = create_mat(n, n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            result[i][j] = 0;
            for (x = 0; x < k; x++) {
                result[i][j] += mat1[i][x] * mat2[x][j];
            }
        }
    }
    return result;
}


/* multiplying n*n matrix by n*k matrix */
double ** matrix_multiplication3(double ** mat1, double ** mat2, int n, int k) {
    int i;
    int j;
    int x;
    double **result = create_mat2(n, k);

    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            for (x = 0; x < n; x++) {
                result[i][j] += mat1[i][x]*mat2[x][j];
            }
        }
    }
    return result;
}


/* transposing a n*k matrix */
double ** transpose_matrix(double ** mat, int n, int k) {
    int i;
    int j;
    double **transposed = create_mat(k, n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            transposed[j][i] = mat[i][j];
        }
    }

    return transposed;
}


/* computing frobenius norm of a matrix */
double f_norm(double ** mat, int n, int k){
    int i;
    int j;
    double norm = 0;
    double val;
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            val = mat[i][j];
            if(val < 0){
                val = (-1)*val;
            }
            norm += pow(val, 2);

        }
    }
    return norm;
}


/* subtracting between n*k matrices */
double ** matrix_sub(double ** mat1, double ** mat2, int n, int k) {
    int i;
    int j;

    double **sub_mat = create_mat(n, k);

    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            sub_mat[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return sub_mat;
}


/* computing H */
/* the function requires the initialized H, W, n, k (H dim) */
double ** finalmat(double ** init_mat, double ** norm_mat, int n, int k){
    int i;
    int j;
    double ** WHmat;
    double ** tmp;
    double ** trans_mat;
    double ** triple_mat;
    double ** prev_mat = create_mat(n, k);
    double ** curr_mat = create_mat(n, k);
    int cnt = 0;
    double ** sub_mat;
    double fnorm;

    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            prev_mat[i][j] = init_mat[i][j];
        }
    }

    WHmat = matrix_multiplication3(norm_mat, prev_mat, n, k);
    trans_mat = transpose_matrix(prev_mat, n, k);
    tmp = matrix_multiplication2(prev_mat, trans_mat, n, k);
    triple_mat = matrix_multiplication3(tmp, prev_mat, n ,k);
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            curr_mat[i][j] = prev_mat[i][j](1 - beta + beta(WHmat[i][j]/triple_mat[i][j]));
        }
    }
    cnt++;
    sub_mat = matrix_sub(curr_mat, prev_mat, n, k);
    fnorm = f_norm(sub_mat, n, k);

    while (cnt < max_iter && fnorm >= eps) {
        /update prev_mat to curr_mat/
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                prev_mat[i][j] = curr_mat[i][j];
            }
        }

        WHmat = matrix_multiplication3(norm_mat, prev_mat, n, k);
        trans_mat = transpose_matrix(prev_mat, n, k);
        tmp = matrix_multiplication2(prev_mat, trans_mat, n, k);
        triple_mat = matrix_multiplication3(tmp, prev_mat, n ,k);
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                curr_mat[i][j] = prev_mat[i][j](1 - beta + beta(WHmat[i][j]/triple_mat[i][j]));
            }
        }
        sub_mat = matrix_sub(curr_mat, prev_mat, n, k);
        fnorm = f_norm(sub_mat, n, k);
        cnt ++;
    }

    free_mat(WHmat, n);
    free_mat(tmp, n);
    free_mat(trans_mat, k);
    free_mat(triple_mat, n);
    free_mat(prev_mat, n);
    free_mat(sub_mat, n);

    return curr_mat;
}


/* computes number of vectors in the input */
int vectors_num(struct vector * vectors){
    int N = 0;
    struct vector * curr = vectors;
    while (curr->next != NULL){
        N++;
        curr = curr->next;
    }
    return N;
}


/* computes the length of the vectors in the input */
int vectors_len(struct vector * vectors){
    int d = 0;
    struct vector * vec = vectors;
    struct cord * curr = vec->cords;
    while (curr != NULL){
        d++;
        curr = curr->next;
    }
    return d;
}


/* inserting vectors into the matrix from the linked list */
double ** insert_vectors(struct vector * vectors_ll, int vec_num, int vec_len){
    double ** matrix = create_mat(vec_num, vec_len);
    int i;
    int j;
    struct vector *curr = vectors_ll;
    struct cord *curr1;
    for (i = 0; i < vec_num; ++i) {
        curr1 = curr->cords;
        for (j = 0; j < vec_len ; ++j) {
            matrix[i][j] = curr1->value;
            curr1 = curr1->next;
        }
        curr = curr->next;
    }
    return matrix;
}


/* creating matrix from the given data */
double ** create_init_mat(struct vector * input_data, int vec_num, int vec_len){
    double ** vectors_mat;
    vectors_mat = insert_vectors(input_data, vec_num, vec_len);
    return vectors_mat;
}


/* print the matrix */
void print_mat(double ** vectors_mat, int vec_len, int vec_num){
    int i;
    int j;
    int d;
    for (i = 0; i < vec_num; i++) {
        d = vec_len;
        for (j = 0; j < vec_len; j++) {
            printf("%.4f", vectors_mat[i][j]);
            if (d > 1){
                printf(",");
                d--;
            }
            else{
                printf("\n");
            }
        }
    }
}


int main(int argc, char **argv) {
    char *goal;
    char *file1;
    double **init_mat;
    int vec_num;
    int vec_len;
    FILE *file = NULL;

    struct vector *head_vec;
    struct vector *curr_vec;
    struct cord *head_cord;
    struct cord *curr_cord;
    double n;
    char c;

    int total_cord_num = 0;
    double **mat = NULL;


    if (argc != 3){
        printf("An Error Has Occurred\n");
    }

    goal = (char *) argv[1];
    file1 = argv[2];

    file = fopen(file1, "r");
    if (file == NULL) {
        printf("Failed to open the file\n");
        return 1;
    }

    head_cord = malloc(sizeof(struct cord));
    if(head_cord == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = malloc(sizeof(struct vector));
    if(head_vec == NULL){
        printf("An Error Has Occurred");
        free(head_cord);
        return 1;
    }
    curr_vec = head_vec;
    curr_vec->next = NULL;


    while (fscanf(file, "%lf%c", &n, &c) == 2){
        if (c == '\n'){
            curr_cord->value = n;
            curr_vec->cords = head_cord;
            curr_vec->next = malloc(sizeof(struct vector));
            if(curr_vec->next == NULL){
                printf("An Error Has Occurred");
                free_LinkedList(head_vec,total_cord_num);
                return 1;
            }
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_cord = malloc(sizeof(struct cord));
            total_cord_num++;
            if(head_cord == NULL){
                printf("An Error Has Occurred");
                free_LinkedList(head_vec,total_cord_num);
                return 1;
            }
            curr_cord = head_cord;
            curr_cord->next = NULL;
            continue;
        }
        curr_cord->value = n;
        curr_cord->next = malloc(sizeof(struct cord));
        if(curr_cord->next == NULL){
            printf("An Error Has Occurred");
            free_LinkedList(head_vec,total_cord_num);
            return 1;
        }
        curr_cord = curr_cord->next;
        curr_cord->next = NULL;
    }

    fclose(file);

    vec_num = vectors_num(head_vec);
    vec_len = vectors_len(head_vec);

    init_mat = create_init_mat(head_vec, vec_num, vec_len);
    free_LinkedList(head_vec,total_cord_num);
    free(head_cord);


    if(strcmp(goal, "sym") == 0){
        mat = symmat(init_mat , vec_num, vec_len);
    }

    else if(strcmp(goal, "ddg") == 0){
        mat = ddgmat(init_mat, vec_num, vec_len);
    }

    else if(strcmp(goal, "norm") == 0){
        mat = normmat(init_mat, vec_num, vec_len);
    }

    print_mat(mat, vec_num, vec_num);

    free_mat(mat, vec_num);
    free_mat(init_mat, vec_num);

    return 0;
}
