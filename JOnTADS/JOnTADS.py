import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from data_process import load_mat, get_shuffle_mat
from get_bds import get_bds, select_bds
from dp import get_dp_mat, get_tads
from quant_reg import quant_reg
import numpy as np
import time
import argparse
from get_stripe import get_stripe


def JOnTADS_one(file_name, output, max_sz, min_sz, significance, percentile, num_ft_cuts, lmda):
    print('Load contact matrix and shuffle it, spend time: ')
    start = time.time()
    mat, cum_mat = load_mat(file_name, max_sz)
    mat_shuffle, cum_mat_shuffle = get_shuffle_mat(mat, max_sz)
    end = time.time()
    print('%.2fs' % (end-start))
    
    # identify boundaries
    print('Identify boundaries, spend time: ')
    start = time.time()
    bds = get_bds(mat, mat_shuffle, max_sz, min_sz, percentile, num_ft_cuts)
    bds_shuffle = get_bds(mat_shuffle, mat_shuffle, max_sz, min_sz, percentile, num_ft_cuts)
    end = time.time()
    print('%.2fs' % (end - start))

    # get tads
    print('Get tads, spend time: ')
    start = time.time()
    dp_mat_shuffle = get_dp_mat(np.array(bds_shuffle), mat_shuffle, cum_mat_shuffle, max_sz, min_sz, np.zeros(max_sz-min_sz+1))
    tads_shuffle = get_tads(bds_shuffle, dp_mat_shuffle, max_sz, min_sz)
    dp_mat = get_dp_mat(np.array(bds), mat, cum_mat, max_sz, min_sz, np.zeros(max_sz-min_sz+1))
    tads = get_tads(bds, dp_mat, max_sz, min_sz)
    end = time.time()
    print('%.2fs' % (end - start))
    
    # filter tads and boundaries
    print('Filter tads and boundaries, spend time: ')
    start = time.time()
    idx, curve = quant_reg(tads, tads_shuffle, max_sz, min_sz, significance, lmda)
    bds = np.unique(tads[idx][:, :2])
    end = time.time()
    print('%.2fs' % (end - start))

    # use filtered boundaries to identify tads
    print('Use filtered boundaries to identify tads, spend time: ')
    start = time.time()
    dp_mat = get_dp_mat(np.array(bds, dtype=int), mat, cum_mat, max_sz, min_sz, np.zeros(max_sz-min_sz+1))
    tads = get_tads(bds, dp_mat, max_sz, min_sz)
    end = time.time()
    print('%.2fs' % (end - start))

    # save tads
    print('Save tads in the output file.')
    np.savetxt(output, tads, fmt='%d, %d, %.3g, %d')

def JOnTADS(file_names, outputs, max_sz, min_sz, significance, percentile, num_ft_cuts, lmda):
    if len(file_names) == 1:
        JOnTADS_one(file_names[0], outputs[0], max_sz, min_sz, significance, percentile, num_ft_cuts, lmda)
        return 1
    
    def declutter_bds(bds_all, mats, max_sz, min_sz):
        max_num = 0
        idx_max_num = -1
        bds_num = []
        for i in range(len(bds_all)):
            bds_num.append(len(bds_all[i]))
            if len(bds_all[i]) > max_num:
                max_num = len(bds_all[i])
                idx_max_num = i
        ratio = []
        for i in range(len(bds_all)):
            ratio.append(len(bds_all[i])/max_num)
        ft_all = []
        for i in range(len(mats)):
            ft = []
            for j in range(mats[i].shape[0]):
                ft.append(np.mean(np.sort(mats[i][j, max(0, j-max_sz+1):min(mats[i].shape[0], j+max_sz)])[-min_sz:]))
            ft_all.append(ft)
        ft_all = np.array(ft_all)
        ft_pooled = np.mean(ft_all, axis=0)
        bds_all = np.unique([j for i in bds_all for j in i])
        bds_all = [np.array(select_bds(bds_all, ft_pooled, int((min_sz+1)/2))) for i in range(len(mats))]
        bds_all = [np.sort(np.array(bds_all[i])[np.argsort(ft_all[i][bds_all[i]])[int(-len(bds_all[idx_max_num])*ratio[i]):]]) for i in range(len(mats))]
        return bds_all
        
    print('Load contact matrices and shuffle them, spend time: ')
    start = time.time()
    mats, cum_mats, mats_shuffle, cum_mats_shuffle = [], [], [], []   
    for file_name in file_names:
        _mat, _cum_mat = load_mat(file_name, max_sz)
        _mat_shuffle, _cum_mat_shuffle = get_shuffle_mat(_mat, max_sz)
        mats.append(_mat)
        cum_mats.append(_cum_mat)
        mats_shuffle.append(_mat_shuffle)
        cum_mats_shuffle.append(_cum_mat_shuffle)
    end = time.time()
    print('%.2fs' % (end - start))

    print('Get boundaries, spend time: ')
    start = time.time()
    bds_each = []
    bds_all, bds_all_shuffle = [], []  
    bds_vas = []
    for mat, mat_shuffle in zip(mats, mats_shuffle):
        _bds = get_bds(mat, mat_shuffle, max_sz, min_sz, percentile, num_ft_cuts)
        bds_all += _bds
        bds_each.append(_bds)
        bds_all_shuffle.append(get_bds(mat_shuffle, mat_shuffle, max_sz, min_sz, percentile, num_ft_cuts))
    bds_all = declutter_bds(bds_each, mats, max_sz, min_sz)
    end = time.time()
    print('%.2fs' % (end - start))
    
    print('Get tads, spend time: ')
    start = time.time()
    tads_all = []
    for i in range(len(file_names)):
        _dp_mat_shuffle = get_dp_mat(np.array(bds_all_shuffle[i]), mats_shuffle[i], cum_mats_shuffle[i], max_sz, min_sz, np.zeros(max_sz-min_sz+1))
        _tads_shuffle = get_tads(bds_all_shuffle[i], _dp_mat_shuffle, max_sz, min_sz)
        _dp_mat = get_dp_mat(np.array(bds_all[i]), mats[i], cum_mats[i], max_sz, min_sz, np.zeros(max_sz-min_sz+1))
        _tads = get_tads(bds_all[i], _dp_mat, max_sz, min_sz)
        _idx, curve = quant_reg(_tads, _tads_shuffle, max_sz, min_sz, significance, lmda)
        bds_all[i] = np.unique(_tads[_idx][:, :2])

        _dp_mat = get_dp_mat(np.array(bds_all[i], dtype=int), mats[i], cum_mats[i], max_sz, min_sz, np.zeros(max_sz-min_sz+1))
        _tads = get_tads(bds_all[i], _dp_mat, max_sz, min_sz)
        tads_all.append(_tads)
    end = time.time()
    print('%.2fs' % (end - start))
        
    for i in range(len(outputs)):
        np.savetxt(outputs[i], tads_all[i], fmt='%d, %d, %.3g, %d')
    
def get_args():
    parser = argparse.ArgumentParser(description='Jointly identify nested TADs.')
    parser.add_argument('-F', '--file_name', type=str, nargs='+',
                        help='Input file path')
    parser.add_argument('-O', '--output', type=str, nargs='+',
                        help='Output file path')
    parser.add_argument('-MAXSZ', '--max_sz', type=int, default=200,
                        help='Maximum size of TADs')
    parser.add_argument('-MINSZ', '--min_sz', type=int, default=7,
                        help='Minimum size of TADs')
    parser.add_argument('--significance', type=float, default=0.1,
                        help='Significance for filtering TADs and boundaries')
    parser.add_argument('--percentile', type=float, default=90,
                        help='Percentile of boundaries considered')
    parser.add_argument('--num_ft_cuts', type=int, default=5,
                        help='Number of feature cut points')
    parser.add_argument('--num_rg_cuts', type=int, default=5,
                        help='Number of regression cut points')
    parser.add_argument('--seed', type=int, default=0,
                        help='Initialize the random number generator')
    parser.add_argument('--lmda', type=float, default=0.001,
                        help='Lambda of L1 norm in nonparametric quantile regression')

    parser.add_argument('--stripe', type=bool, default=False,
                        help='Call stripes or not')
    parser.add_argument('--stripe_output', type=str, default='',
                        help='Output file to save stripe results')
    parser.add_argument('--chr', type=str, default='',
                        help='Chromosome information for calling stripes')
    parser.add_argument('--p', type=float, default=1,
                        help='Exponent for measuring similarity')
    parser.add_argument('--threshold', type=float, default=0.1,
                        help='Threshold to call stripes')
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    np.random.seed(args.seed)
    if args.stripe == True:
        stripes = get_stripe(args.file_name[0], chr_num=args.chr, max_sz=args.max_sz, min_sz=args.min_sz, p=args.p, threshold=args.threshold, mean_threshold=args.threshold)
        np.savetxt(args.stripe_output, stripes, fmt='%s', delimiter='\t')
    print('Finished stripe identification.')
    np.random.seed(args.seed)
    JOnTADS(args.file_name, args.output,
               args.max_sz, args.min_sz,
               args.significance, args.percentile,
               args.num_ft_cuts, args.lmda)
    return

if __name__ == '__main__':
    main()
