import math
import sys
import pandas as pd
import numpy as np
import random
import mysymnmfsp
from sklearn.metrics import silhouette_score
np.random.seed(0)


## kmeans implementation from HW1 ##
def extract_vectors(inpt):
    with open(inpt, "r") as file:
        vectors = []
        for line in file:
            if len(line) > 1:
                vectors.append([float(coordinate) for coordinate in line[:len(line)].split(',')])
        return vectors


def euclidean_dis(vec1, vec2):
    d = len(vec1)
    dis = 0
    for i in range(d):
        a = vec1[i]
        b = vec2[i]
        dis += ((vec1[i]-vec2[i])**2)
    dis = (math.sqrt(dis))
    return dis


def assign_centroids(clust, vec, k):
    for i in range(k):
        clust.append([vec[i]])


def make_centroids(vec, k):
    centroids = []
    for i in range(k):
        centroids.append(vec[i])
    return centroids


def assign_vectors(clust, vec, centroids, index):
    for i in range(len(vec)):
        dis_from_first_cent = euclidean_dis(vec[i], centroids[0])
        closest_cent = 0
        for j in range(1, len(centroids)):
            dis_from_curr_cent = euclidean_dis(vec[i], centroids[j])
            if dis_from_curr_cent < dis_from_first_cent:
                dis_from_first_cent = dis_from_curr_cent
                closest_cent = j
        clust[closest_cent].append(vec[i])
        index[i] = closest_cent
            

def assign_vectors2(vectors, centroids, k, index):
    new_cluster = [[] for i in range(k)]
    cnt = 0
    for i in vectors:
        dis_from_first_cent = euclidean_dis(i, centroids[0])
        closest_cent = 0
        for j in range(1, len(centroids)):
            dis_from_curr_cent = euclidean_dis(i, centroids[j])
            if dis_from_curr_cent < dis_from_first_cent:
                dis_from_first_cent = dis_from_curr_cent
                closest_cent = j
        new_cluster[closest_cent].append(i)
        index[cnt] = closest_cent
        cnt += 1
    return new_cluster


def calc_new_cent(cluster):
    new_cent = [0 for i in range(len(cluster[0]))]
    for i in cluster:
        for j in range(len(cluster[0])):
            new_cent[j] += i[j]
    for p in range(len(new_cent)):
        c = new_cent[p]
        d = len(cluster)
        new_cent[p] = ((new_cent[p]) / len(cluster))
        t = new_cent[p]
        f = new_cent[p]
    return new_cent


def adjust_centroids(clusters, centroids, k):
    delta = [0 for i in range(k)]
    for i in range(k):
        new_cent = calc_new_cent(clusters[i])
        delta[i] = euclidean_dis(centroids[i], new_cent)
        centroids[i] = new_cent
    return delta


def check_centroids(centroids):
    for i in centroids:
        if i >= 0.0001:
            return False
    return True


def print_centroids(centroids, d):
    for vector in centroids:
        centroid_len = d
        for i in range(centroid_len):
            if centroid_len > 1:
                print('{:.4f}'.format(vector[i]) + ",", end="")
                centroid_len -= 1
            else:
                print('{:.4f}'.format(vector[i]) + '\n', end="")


def k_means(K, iters, input_data):
    vectors = extract_vectors(input_data)
    index = [0 for i in range(len(vectors))]
    iter_num = 0
    k1 = int(K)  # k->k1
    centroids = make_centroids(vectors, k1)
    clusters = [[] for i in range(k1)]
    assign_vectors(clusters, vectors, centroids, index)
    deltas = adjust_centroids(clusters, centroids, k1)
    check = check_centroids(deltas)
    if check:
        for i in centroids:
            return centroids, index, vectors
    while iter_num < iters:
        clusters = assign_vectors2(vectors, centroids, k1, index)
        deltas = adjust_centroids(clusters, centroids, k1)
        check = check_centroids(deltas)
        if check:
            return centroids, index, vectors
        iter_num += 1
    return centroids, index, vectors


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


# main function
def main_func(k, file2):
    epsilon = 0.0001
    max_iter = 300
    help_table = create_table(file2)
    vec_num = help_table.shape[0]
    vec_len = help_table.shape[1]
    vec_mat = create_mat(help_table)
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
    symnmf_answer = mysymnmfsp.symnmf(vec_num, k, H, norm_mat)
    
    
    # computing Silhouette coefficient
    npH = np.array(symnmf_answer)
    symnmf_array = np.argmax(npH, axis=1)
    kmeans_table, kmeans_cluster, vectors = k_means(k, 300, file2)
    symnmf_score = "{:.4f}".format(silhouette_score(vectors, symnmf_array))
    kmeans_score = "{:.4f}".format(silhouette_score(vectors, kmeans_cluster))
    print("nmf:", symnmf_score)
    print("kmeans:", kmeans_score)
    
    
args = sys.argv
args_len = len(args)
k1 = 0
file1 = None
if args_len != 3:
    print("An Error Has Occurred")
    exit()
try:
    k1 = int(args[1])
    file1 = args[2]
except Exception:
    print("An Error Has Occurred kkk")
    exit()

main_func(k1, file1)
