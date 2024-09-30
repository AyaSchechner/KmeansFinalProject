import pandas as pd
import numpy as np
import sys
import math
import random
import mysymnmfsp
np.random.seed(0)


# create a sorted vectors table from the input files via pandas
def create_table(input_file1):
    df1 = pd.read_csv(input_file1, header=None)
    return df1


# create the vectors matrix from the pandas table
def create_mat(vectors):
    vectors_mat = []
    for i in range(vectors.shape[0]):
        curr_vec = vectors.iloc[i, 0:].tolist()
        vectors_mat.append(curr_vec)
    return vectors_mat


# print matrix
def print_mat(mat, vec_len, vec_num):
    for i in range(vec_num):
        d = vec_len
        vec = mat[i]
        for j in range(vec_len):
            if d > 1:
                print('{:.4f}'.format(vec[j]) + ",", end="")
                d -= 1
            else:
                print('{:.4f}'.format(vec[j]) + '\n', end="")


# computing the goal matrix
def main_func(k, goal, file):
    help_table = create_table(file)
    vec_num = help_table.shape[0]
    vec_len = help_table.shape[1]
    vec_mat = create_mat(help_table)
    if goal == "sym":
        answer = mysymnmfsp.sym(vec_len, vec_num, vec_mat)
        print_mat(answer, vec_num, vec_num)
        return answer
    
    if goal == "ddg":
        answer = mysymnmfsp.ddg(vec_len, vec_num, vec_mat)
        print_mat(answer, vec_num, vec_num)
        return answer
    
    if goal == "norm":
        answer = mysymnmfsp.norm(vec_len, vec_num, vec_mat)
        print_mat(answer, vec_num, vec_num)
        return answer
    
    if goal == "symnmf":
        similarity_mat = mysymnmfsp.sym(vec_len, vec_num, vec_mat)
        diagonal_mat = mysymnmfsp.ddg(vec_len, vec_num, vec_mat)
        norm_mat = mysymnmfsp.norm(vec_len, vec_num, vec_mat)
        m = 0
        for i in range(vec_num):
            for j in range(vec_num):
                m += norm_mat[i][j]
        m = m/(vec_num*vec_num)
        H = [[0 for j in range(k)] for i in range(vec_num)]
        for i in range(vec_num):
            for j in range(k):
                H[i][j] = np.random.uniform(0, 2*(math.sqrt(m/k)))
        answer = mysymnmfsp.symnmf(vec_num, k, H, norm_mat)
        print_mat(answer, k, vec_num)
        return answer


args = sys.argv
args_len = len(args)
k1 = 0
goal1 = None
file1 = None
if args_len != 4:
    print("An Error Has Occurred")
    exit()
try:
    k1 = int(args[1])
    goal1 = str(args[2])
    file1 = args[3]
except Exception:
    print("An Error Has Occurred kkk")
    exit()

main_func(k1, goal1, file1)
