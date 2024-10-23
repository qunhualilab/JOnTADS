import numba
import pandas as pd
import numpy as np

@numba.njit
def get_dp_mat(bds, mat, cum_mat, max_sz, min_sz, penalty):
    # 0: intermediate point position
    # 1: outer TAD whether included
    # 2: left boundary whether included
    # 3: right boundary whether included
    # 4: sum of entry values
    # 5: num of points used
    # 6: sum of TAD scores
    # 7: TAD score
    # 8: the maximum TAD size of current hierarchy
    l = len(bds)
    dp_mat = []
    for i in range(l):
        dp_mat.append([[-1, 0, 0, 0, 0., 0, 0., 0., 0]
                       for _ in range(l)])
    for i in range(1, int(np.ceil((max_sz-1)/np.ceil((min_sz-1)/2))+1)+1):
        for j in range(len(bds)-i):
            if bds[i + j] - bds[j] + 1 < min_sz:
                # TAD size >= min_sz
                continue
            if bds[i + j] - bds[j] + 1 <= max_sz:
                # TAD size <= max_sz
                tad_size = bds[i + j] - bds[j] + 1
                sub_mat_sum, sub_mat_num_points = sum_sub_mat(bds[j], bds[j + i], cum_mat)
                
                outer_include = 1
                corner_sum, corner_num = corner(bds[j], bds[j + i], max(int((min_sz+1)/2), int((bds[j+i]-bds[j]+1)/10)), cum_mat)
                if (corner_sum * sub_mat_num_points) < 0.6 * (sub_mat_sum * corner_num):
                    outer_include = 0
                
                sub_mat_surr_width = int(round(sub_mat_num_points /
                                       (2. * (tad_size-1))))
                sub_mat_surr_mean = cal_surr_mean(bds[j], bds[i + j],
                                          cum_mat, sub_mat_surr_width)
                # outer TAD score
                size = int(bds[j+i] - bds[j] + 1)
                score = 0
                score_max = max(0, sub_mat_sum / sub_mat_num_points - sub_mat_surr_mean - penalty[size-min_sz])
                
                inter = -1
                for k in range(j+1, j+i):
                    # nested TADs should be at least min_sz-1 smaller
                    if dp_mat[j][k][8] > tad_size - min_sz + 1 or dp_mat[k][j+i][8] > tad_size - min_sz + 1 or \
                    tad_size < 1.1 * max(dp_mat[j][k][8], dp_mat[k][j+i][8]):
                    #if dp_mat[j][k][8] > tad_size - min_sz + 1 or dp_mat[k][j+i][8] > tad_size - min_sz + 1:
                        if dp_mat[j][k][6] + dp_mat[k][j+i][6] > score_max:
                            score_max = dp_mat[j][k][6] + dp_mat[k][j+i][6]
                            score = 0
                            inter = k
                        continue

                    sub_mat_sum_left = sub_mat_sum
                    sub_mat_num_points_left = sub_mat_num_points
                    if dp_mat[j][k][4] != 0:
                        sub_mat_sum_left -= dp_mat[j][k][4]
                        sub_mat_num_points_left -= dp_mat[j][k][5]
                    if dp_mat[k][j+i][4] != 0:
                        sub_mat_sum_left -= dp_mat[k][j+i][4]
                        sub_mat_num_points_left -= dp_mat[k][j+i][5]
                    if dp_mat[j][k][3] != 0 and dp_mat[k][j+i][2] != 0:
                        sub_mat_sum_left += mat[bds[k], bds[k]]
                        sub_mat_num_points_left += 1
                    sub_mat_mean_left = sub_mat_sum_left / \
                                        sub_mat_num_points_left

                    sub_mat_surr_width_left = int(round(
                        sub_mat_num_points_left / (2. * (tad_size-1))))
                    #sub_mat_surr_width_left = max(int(round(
                    #    sub_mat_num_points_left / (2. * (tad_size-1)))), min_sz-1)
                    #sub_mat_surr_width_left = tad_size - 1
                    sub_mat_surr_mean_left = cal_surr_mean(
                        bds[j], bds[j+i], cum_mat, sub_mat_surr_width_left)
                    if outer_include == 1 and sub_mat_mean_left - sub_mat_surr_mean_left - penalty[size-min_sz] > 0:
                        if sub_mat_mean_left - sub_mat_surr_mean_left + \
                                dp_mat[k][j+i][6] + dp_mat[j][k][6] - \
                                penalty[size-min_sz] > score_max:
                            score_max = sub_mat_mean_left - \
                                        sub_mat_surr_mean_left + \
                                        dp_mat[j][k][6] + dp_mat[k][j+i][6] - penalty[size-min_sz]
                            score = sub_mat_mean_left - \
                                    sub_mat_surr_mean_left - penalty[size-min_sz]
                            inter = k
                    else:
                        if dp_mat[j][k][6] + dp_mat[k][j+i][6] > score_max:
                            score_max = dp_mat[j][k][6] + dp_mat[k][j+i][6]
                            score = 0
                            inter = k
                            
                if score > 0:
                    dp_mat[j][j+i][0] = inter
                    dp_mat[j][j+i][1] = 1
                    dp_mat[j][j+i][2] = 1
                    dp_mat[j][j+i][3] = 1
                    dp_mat[j][j+i][4] = sub_mat_sum
                    dp_mat[j][j+i][5] = sub_mat_num_points
                    dp_mat[j][j+i][6] = score_max
                    dp_mat[j][j+i][7] = score
                    dp_mat[j][j+i][8] = bds[j+i] - bds[j] + 1
                elif inter != -1:
                    dp_mat[j][j+i][0] = inter
                    dp_mat[j][j+i][1] = 0
                    if dp_mat[j][inter][2] == 1:
                        dp_mat[j][j+i][2] = 1
                    if dp_mat[inter][j+i][3] == 1:
                        dp_mat[j][j+i][3] = 1
                    if dp_mat[j][inter][3] == 1 and \
                            dp_mat[inter][j+i][2] == 1:
                        dp_mat[j][j+i][4] = dp_mat[j][inter][4] + \
                                            dp_mat[inter][j+i][4] - \
                                            mat[bds[inter], bds[inter]]
                        dp_mat[j][j+i][5] = dp_mat[j][inter][5] + \
                                            dp_mat[inter][j+i][5] - 1
                    else:
                        dp_mat[j][j+i][4] = dp_mat[j][inter][4] + \
                                            dp_mat[inter][j+i][4]
                        dp_mat[j][j+i][5] = dp_mat[j][inter][5] + \
                                            dp_mat[inter][j+i][5]
                    dp_mat[j][j+i][6] = dp_mat[j][inter][6] + \
                                        dp_mat[inter][j+i][6]
                    dp_mat[j][j+i][8] = max(dp_mat[j][inter][8],
                                            dp_mat[inter][j+i][8])
                elif outer_include == 1 and score_max > 0:
                    dp_mat[j][j+i][0] = inter
                    dp_mat[j][j+i][1] = 1
                    dp_mat[j][j+i][2] = 1
                    dp_mat[j][j+i][3] = 1
                    dp_mat[j][j+i][4] = sub_mat_sum
                    dp_mat[j][j+i][5] = sub_mat_num_points
                    dp_mat[j][j+i][6] = score_max
                    dp_mat[j][j+i][7] = score_max
                    dp_mat[j][j+i][8] = bds[j+i] - bds[j] + 1
                else:
                    continue
            else:
                continue
    
    return dp_mat

@numba.njit
def corner(bd1, bd2, interval, cum_mat):
    sub_mat_sum = cum_mat[bd1+interval-1, bd2] - cum_mat[bd1+interval-1, bd2-interval] - cum_mat[bd1-1, bd2] + cum_mat[bd1-1, bd2-interval]
    return sub_mat_sum, interval**2

@numba.njit
def sum_sub_mat(l_b, r_b, cum_mat):
    # l_b: left boundary
    # r_b: right boundary
    sub_mat_size = r_b - l_b + 1

    if l_b == 0:
        return cum_mat[r_b, r_b], sub_mat_size**2

    sub_mat_sum = cum_mat[r_b, r_b] - cum_mat[r_b, l_b-1] - \
                  cum_mat[l_b-1, r_b] + cum_mat[l_b-1, l_b-1]
    return sub_mat_sum, sub_mat_size**2

@numba.njit
def cal_surr_mean(l_b, r_b, cum_mat, surr_width):
    # up_mean: [max(0, l_b-surr_width) : l_b-1, l_b+1 : r_b]
    # ri_mean: [l_b : r_b-1, r_b+1 : min(r_b, size(cum_mat)-1)]
    if l_b == 0:
        up_mean = 0
    else:
        up_b = max(0, l_b-surr_width)
        if up_b == 0:
            up_sum = cum_mat[l_b-1, r_b] - cum_mat[l_b-1, l_b]
        else:
            up_sum = cum_mat[l_b-1, r_b] - cum_mat[l_b-1, l_b] - \
                  cum_mat[up_b-1, r_b] + cum_mat[up_b-1, l_b]
        up_mean = up_sum / ((r_b-l_b)*(l_b-up_b))

    if r_b == cum_mat.shape[0]-1:
        ri_mean = 0
    else:
        ri_b = min(cum_mat.shape[0]-1, r_b+surr_width)
        if l_b == 0:
            ri_sum = cum_mat[r_b-1, ri_b] - cum_mat[r_b-1, r_b]
        else:
            ri_sum = cum_mat[r_b-1, ri_b] - cum_mat[r_b-1, r_b] - \
                     cum_mat[l_b-1, ri_b] + cum_mat[l_b-1, r_b]
        ri_mean = ri_sum / ((r_b-l_b)*(ri_b-r_b))
    return max(up_mean, ri_mean)

def get_tads(bds, dp_mat, max_sz, min_sz):
    # tads: [[position, [[tad], [tad]], max_score]]
    # ensure the 'inter' is an integer
    for i in range(len(bds)):
        for j in range(len(bds)):
            dp_mat[i][j][0] = int(dp_mat[i][j][0])
            
    tads = [[0, [], 0]]
    for i, bd, in enumerate(bds[1:]):
        # TAD should be smaller than (or equal to) max_sz
        if bd - bds[tads[-1][0]] + 1 > max_sz:
            tads = [tads[-1]]
            tads[-1][0] = i + 1
        else:
            tads_last = tads[-1][1]  # [[tad], [tad]]
            max_score_last = tads[-1][2]  # float

            # tads_cur: [position, [[tad], [tad]], max_score]
            for j, tads_cur in enumerate(tads[::-1]):
                if bd - bds[tads_cur[0]] + 1 > max_sz:
                    tads = tads[-j:]
                    break

                # TAD should be larger than (or equal to) mim_sz
                if bd - bds[tads_cur[0]] + 1 < min_sz:
                    continue
                else:
                    max_score = tads_cur[2] + dp_mat[tads_cur[0]][i+1][6]
                    if max_score > max_score_last:
                        tads_last = tads_cur[1] + sub_get_tads(
                            dp_mat, tads_cur[0], i+1, bds, [])
                        max_score_last = max_score
            tads.append([i+1, tads_last, max_score_last])

    return np.array(tads[-1][1])

def sub_get_tads(dp_mat, l_b, r_b, bds, tads):
    # tads: [[l_b, r_b, size, score]]
    # in order: left to right, outer to inner
    if dp_mat[l_b][r_b][1] == 1:
        tads.append([bds[l_b], bds[r_b], dp_mat[l_b][r_b][7], bds[r_b]-bds[l_b]+1])
    if dp_mat[l_b][r_b][0] != -1:
        tads = sub_get_tads(dp_mat, l_b, dp_mat[l_b][r_b][0], bds, tads)
        tads = sub_get_tads(dp_mat, dp_mat[l_b][r_b][0], r_b, bds, tads)
    return tads
