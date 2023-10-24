import numpy as np
import numba
import pandas as pd

@numba.jit(nopython=True)
def cum_sum(mat):
    # cum_mat[i, j] = sum(mat[:(i+1), :(j+1)])
    cum_mat = mat.copy()
    for i in range(mat.shape[0]):
        for j in range(1, mat.shape[0]):
            cum_mat[i, j] += cum_mat[i, j-1]
        if i > 0:
            for j in range(mat.shape[0]):
                cum_mat[i, j] += cum_mat[i-1, j]
    return cum_mat

@numba.jit(nopython=True)
def normalization(mat, max_sz):
    # Normalize the diagonals [0, 2*max_sz-1]
    for i in range(2*max_sz):
        diag = np.diag(mat, i)
        diag_min = np.min(diag)
        diag_max = np.percentile(diag, 99.9)
        if diag_max == 0 or diag_max == diag_min:
            continue
        for j in range(mat.shape[0]-i):
            mat[j, j+i] = mat[j+i, j] = min(1, (diag[j] - diag_min) / \
                                        (diag_max - diag_min))
    return mat

def _load_mat(file_name, max_sz):
    # load the contact matrix
    f = open(file_name, 'r')
    line = f.readline()
    if len(line.split()) == 1:
        with open(file_name, 'r') as f:
            for i, line in enumerate(f):
                line = line.split(',')
                if i == 0:
                    mat_size = len(line)
                    mat = np.zeros((mat_size, mat_size))
                for j in range(max(0, i-2*max_sz+1), min(mat_size, i+2*max_sz)):
                    mat[i, j] = float(line[j])
    else:
        with open(file_name, 'r') as f:
            for i, line in enumerate(f):
                line = line.split()
                if i == 0:
                    mat_size = len(line)
                    mat = np.zeros((mat_size, mat_size))
                for j in range(max(0, i-2*max_sz+1), min(mat_size, i+2*max_sz)):
                    mat[i, j] = float(line[j])
    return mat

def load_mat(file_name, max_sz):
    # load and cum_sum matrix of file_name
    mat = _load_mat(file_name, max_sz)
    mat = normalization(mat, max_sz)
    cum_mat = cum_sum(mat)
    return mat, cum_mat

def shuffle(mat, max_sz):
    # shuffle mat at each geometric distance within 2*max_sz
    mat_shuffle = mat.copy()
    l = mat_shuffle.shape[0]
    for i in range(2*max_sz):
        diag = np.diag(mat_shuffle, i)
        diag.setflags(write=True)
        np.random.shuffle(diag)
        #np.random.shuffle(np.diag(mat_shuffle, i))
    for i in range(l):
        for j in range(i+1, min(l, i+2*max_sz)):
            mat_shuffle[j, i] = mat_shuffle[i, j]
    return mat_shuffle

def get_shuffle_mat(mat, max_sz):
    # return shuffled mat and its cum_sum matrix
    mat_shuffle = shuffle(mat, max_sz)
    cum_mat_shuffle = cum_sum(mat_shuffle)
    return mat_shuffle, cum_mat_shuffle
