from data_process import load_mat, get_shuffle_mat
from get_bds import get_bds, select_bds
from dp import get_dp_mat, get_tads
from quant_reg import quant_reg
import numpy as np
import time

def ks_h(mat, bd, min_sz=7, max_sz=200, p=1, threshold=0.1, mean_threshold = 0.1):
    if bd == 0 or bd > len(mat) - min_sz - 2:
        return []
    ks_stat = []
    ans1 = (mat[bd, (bd+1):min(len(mat)-1, bd+max_sz)] - mat[bd-1, (bd-1+1):min(len(mat)-2, bd-1+max_sz)]) 
    ans2 = (mat[bd, (bd+1):min(len(mat)-1, bd+max_sz)] - mat[bd+1, (bd+1+1):min(len(mat), bd+1+max_sz)]) 
    ans = np.array([ans1[i]+ans2[i] if ans1[i]>threshold and ans2[i]>threshold else 0 for i in range(len(ans1))])
    if np.sum(ans>0) < 20:
        return []
    NR, Nmiss = 0, 0
    for i in range(len(ans)):
        NR += ans[i] ** p if ans[i]>0 else 0
        Nmiss += 1 if ans[i]<=0 else 0
    if NR == 0:
        return []
    hit, miss = 0, 0
    for i in range(len(ans)):
        if i < min_sz-3:
            hit += ans[i] ** p if ans[i]>0 else 0
            miss += 1 if ans[i]<=0 else 0
        else:
            hit += ans[i] ** p if ans[i]>0 else 0
            miss += 1 if ans[i]<=0 else 0
            try:
                ks_stat.append(hit/NR-miss/Nmiss)
            except:
                ks_stat.append(hit/NR)
    max_idx = np.argmax(ks_stat)
    if max_idx == 0 or np.mean(ks_stat[:(max_idx+1)]) < mean_threshold:
        return []
    else:
        return [[bd], [bd, bd+min_sz+max_idx-2]]
def ks_v(mat, bd, min_sz=7, max_sz=200, p=1, threshold=0.1, mean_threshold = 0.1):
    if bd == len(mat)-1 or bd < min_sz + 1:
        return []
    ks_stat = []
    ans1 = (mat[max(1, bd-max_sz+1):bd, bd] - mat[max(0, bd-1-max_sz+1):bd-1, bd-1]) 
    ans2 = (mat[max(1, bd-max_sz+1):bd, bd] - mat[max(2, bd+1-max_sz+1):bd+1, bd+1]) 
    ans = np.array([ans1[i]+ans2[i] if ans1[i]>threshold and ans2[i]>threshold else 0 for i in range(len(ans1))])
    if np.sum(ans>0) < 20:
        return []
    NR, Nmiss = 0, 0
    for i in range(len(ans)):
        NR += ans[i] ** p if ans[i]>0 else 0
        Nmiss += 1 if ans[i]<=0 else 0
    if NR == 0:
        return []
    hit, miss = 0, 0
    for i in range(len(ans)):
        if i < min_sz-3:
            hit += ans[-(i+1)] ** p if ans[-(i+1)]>0 else 0
            miss += 1 if ans[-(i+1)]<=0 else 0
        else:
            hit += ans[-(i+1)] ** p if ans[-(i+1)]>0 else 0
            miss += 1 if ans[-(i+1)]<=0 else 0
            try:
                ks_stat.append(hit/NR-miss/Nmiss)
            except:
                ks_stat.append(hit/NR)
    max_idx = np.argmax(ks_stat)
    if max_idx == 0 or np.mean(ks_stat[:(max_idx+1)]) < mean_threshold:
        return []
    else:
        return [[bd-min_sz-max_idx+2, bd], [bd]]
def _get_stripe(mat, bds, min_sz=7, max_sz=200, p=1, threshold=0.1, mean_threshold=0.1):
    stripe = []
    for bd in bds:
        stripe_h = ks_h(mat, int(bd), min_sz, max_sz, p, threshold, mean_threshold)
        stripe_v = ks_v(mat, int(bd), min_sz, max_sz, p, threshold, mean_threshold)
        if len(stripe_h) == 2:
            stripe.append(stripe_h)
        if len(stripe_v) == 2:
            stripe.append(stripe_v)
    return stripe

def get_stripe(file_name, min_sz=7, max_sz=200, p=1, threshold=0.1, mean_threshold=0.1):
    mat, cum_mat = load_mat(file_name, max_sz)
    
    bds = np.arange(len(mat))
    
    stripes = _get_stripe(mat, bds, min_sz, max_sz, p, threshold, mean_threshold)
    
    return stripes

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Jointly identify nested TADs.')
    parser.add_argument('-F', '--file_name', type=str)
    parser.add_argument('-O', '--output', type=str)
    parser.add_argument('-C', '--chr', type=str)
    parser.add_argument('-MAXSZ', '--max_sz', type=int, default=200,
                        help='Maximum size of TADs')
    parser.add_argument('-MINSZ', '--min_sz', type=int, default=7,
                        help='Minimum size of TADs')
    parser.add_argument('-P', '--p', type=float, default=1)
    parser.add_argument('-T', '--threshold', type=float, default=0.1)
    args = parser.parse_args()
    stripes = get_stripe(args.file_name, max_sz=args.max_sz, min_sz=args.min_sz, p=args.p, threshold=args.threshold, mean_threshold=args.threshold)

    res = []
    for stripe in stripes:
        if len(stripe[0]) == 1:
            res.append(['chr'+args.chr, stripe[0][0], stripe[0][0], 'chr'+args.chr, stripe[1][0], stripe[1][1]])
        else:
            res.append(['chr'+args.chr, stripe[0][0], stripe[0][1], 'chr'+args.chr, stripe[1][0], stripe[1][0]])
    np.savetxt(args.output, res, fmt='%s', delimiter='\t')
