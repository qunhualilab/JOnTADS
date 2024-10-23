import numpy as np
import bisect
import multiprocessing

def maximize_sum_ft(bds, ft, interval, pos):
    # pos: points that must be included
    p = 0
    res = [[[-interval], -1]] # [[[bds], sum]]
    for i in range(len(bds)):
        if p < len(pos) and bds[i] == pos[p]:
            max_sum = -1
            for j in range(len(res)):
                res_cur = res[-(j+1)]
                bds_cur = res_cur[0]
                sum_cur = res_cur[1]
                if bds[i] - bds_cur[-1] + 1 < interval:
                    if sum_cur + ft[bds[i]] - ft[bds_cur[-1]] > max_sum:
                        max_sum = sum_cur + ft[bds[i]] - ft[bds_cur[-1]]
                        bds_next = bds_cur[:-1] + [bds[i]]
                else:
                    if ft[bds[i]] + sum_cur > max_sum:
                        bds_next = bds_cur + [bds[i]]
                        max_sum = ft[bds[i]] + sum_cur
            p += 1
            res.append([bds_next, max_sum])
            res = [res[-1]]
        else:
            bds_next = res[-1][0]
            max_sum = res[-1][1]
            for j in range(len(res)):
                res_cur = res[-(j+1)]
                bds_cur = res_cur[0]
                sum_cur = res_cur[1]
                if bds[i] - bds_cur[-1] + 1 < interval:
                    if ft[bds[i]] > ft[bds_cur[-1]]:
                        if sum_cur + ft[bds[i]] - ft[bds_cur[-1]] > max_sum:
                            max_sum = sum_cur + ft[bds[i]] - ft[bds_cur[-1]]
                            bds_next = bds_cur[:-1] + [bds[i]]
                else:
                    if ft[bds[i]] + sum_cur > max_sum:
                        bds_next = bds_cur + [bds[i]]
                        max_sum = ft[bds[i]] + sum_cur
            res.append([bds_next, max_sum])
            res = res[-interval+1:]
    return res[-1][0][1:]
    
def _select_bds(bds, ft, interval):
    pos = bds[np.where(ft[bds] == np.max(ft[bds]))[0]]
    pos = maximize_sum_ft(pos, ft, interval, [-1])
    bds = maximize_sum_ft(bds, ft, interval, pos)
    return bds    

def select_bds(bds, ft, interval):
    bds_selected = []
    n = len(bds)
    i = 0
    while i <= n-1:
        bds_part = [bds[i]]
        for j in range(1, n-i+1):
            if i+j <= n-1 and bds[i+j] - bds_part[-1] + 1 < interval:
                bds_part.append(bds[i+j])
            else:
                i += j
                break
        bds_selected += _select_bds(np.array(bds_part), ft, interval)
    return bds_selected

def get_diag_ft_diff(mat, max_sz, min_sz, num):
    # get feature difference of diagonal elements
    cut_p = np.linspace(min_sz, max_sz, num, dtype=int)
    diag_h, diag_v = [[] for i in range(num)], [[] for i in range(num)]
    diag_h_diff, diag_v_diff = [[] for i in range(num)], [[] for i in range(num)]                    
    for i in range(num):
        for j in range(mat.shape[0]):
            if j <= mat.shape[0] - cut_p[i]:
                diag_h[i].append(np.mean(np.sort(mat[j, j+2 : j+cut_p[i]])[-min_sz+2:]))
            if j >= cut_p[i]-1:
                diag_v[i].append(np.mean(np.sort(mat[j-cut_p[i]+1 : j-1, j])[-min_sz+2:]))
    
    for i in range(num):
        diag_h[i] = np.array(diag_h[i])
        diag_v[i] = np.array(diag_v[i][::-1])
    
    for i in range(num):
        diag_h_diff[i] = np.diff(diag_h[i])
        diag_v_diff[i] = np.diff(diag_v[i])
        
    return diag_h_diff, diag_v_diff

from scipy.stats import norm
def norm_ref(c):
    return 1-norm.cdf(norm.ppf(0.9)+c)
def cur_prop(diff, c):
    return np.sum(diff > np.percentile(diff, 90) + c * np.std(diff))/len(diff)
def bisect_trsd(diff):
    left = 1
    right = (np.max(diff) - np.percentile(diff, 90))/np.std(diff)
    mid = (left+right)/2
    for _ in range(10):
        mid = (left+right)/2
        if norm_ref(mid)/cur_prop(diff,mid) < 0.2:
            right = mid
        elif norm_ref(mid)/cur_prop(diff,mid) > 0.2:
            left = mid
        else:
            return mid
    return mid
from scipy.stats import norm
def gmm(data):
    tau, sigma1, sigma2 = 0.5, np.std(data), np.std(data)/100
    Qold = -np.inf
    T = np.ones(len(data))*0.5
    Qnew = 0
    for i in range(len(data)):
        Qnew += T[i] * (np.log(tau)-np.log(sigma1)-1/2*data[i]**2/sigma1**2-1/2*np.log(2*np.pi)) + (1-T[i]) * (np.log(1-tau)-np.log(sigma2)-1/2*data[i]**2/sigma2**2-1/2*np.log(2*np.pi))
    while Qnew-Qold > 0.001:
        Qold = Qnew
        T = np.zeros(len(data))
        for i in range(len(data)):
            T[i] = tau*norm.pdf(data[i], loc=0, scale = sigma1) / (tau*norm.pdf(data[i], loc=0, scale = sigma1) + (1-tau)*norm.pdf(data[i], loc=0, scale = sigma2))
        Qnew = 0
        for i in range(len(data)):
            Qnew += T[i] * (np.log(tau)-np.log(sigma1)-1/2*data[i]**2/sigma1**2-1/2*np.log(2*np.pi)) + (1-T[i]) * (np.log(1-tau)-np.log(sigma2)-1/2*data[i]**2/sigma2**2-1/2*np.log(2*np.pi))
        tau = np.mean(T)
        sigma1 = np.sqrt(np.sum(T*data**2)/np.sum(T))
        sigma2 = np.max([np.sqrt(np.sum((1-T)*data**2)/np.sum(1-T)), sigma1/1000])
    return tau, sigma1, sigma2

from scipy.optimize import minimize
from scipy.stats import norm
def optim(data):
    bnds = ((0, 1), (1e-10, None), (1e-10, None))
    def llh(pars, data):
        tau, sigma1, sigma2 = pars[0], pars[1], pars[2]
        return -np.sum(np.log(tau*norm.pdf(data, loc=0, scale=sigma1)+(1-tau)*norm.pdf(data, loc=0, scale=sigma2)+1e-10))
    pars = np.array([0.5, np.std(data), np.std(data)/10000])
    res = minimize(llh, pars, args=data, bounds=bnds, tol=1e-6, method='L-BFGS-B')
    return res.x

def _get_bds(diag_ft_diff, diag_ft_diff_shuffle, percentile):
    tau, sigma1, sigma2 = optim(diag_ft_diff)
    if sigma2 > np.max(diag_ft_diff):
        tau, sigma2 = 1, 0
    if sigma1 > np.max(diag_ft_diff):
        tau, sigma1, sigma2 = 1, sigma2, 0
    if sigma2 > sigma1:
        tau, sigma1, sigma2 = 1-tau, sigma2, sigma1
    return np.where(diag_ft_diff >= (2*sigma1+(1-tau)*sigma2))[0] + 1

def para_bds(diag_diff, ft, interval, percentile):
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    tasks = [(diag_diff[i], diag_diff[i], percentile) for i in range(len(diag_diff))]
    res = pool.starmap(_get_bds, tasks)
    pool.close()
    pool.join()
    return res
    
def get_bds(mat, mat_shuffle, max_sz, min_sz, percentile, num):
    ft = []
    for i in range(mat.shape[0]):
        ft.append(np.mean(np.sort(mat[i, max(0, i-max_sz+1):min(mat.shape[0], i+max_sz)])[-min_sz:]))
    ft = np.array(ft)
    
    n = mat.shape[0]
    diag_h_diff, diag_v_diff = get_diag_ft_diff(mat, max_sz, min_sz, num)
    diag_h_diff_shuffle, diag_v_diff_shuffle = get_diag_ft_diff(mat_shuffle, max_sz, min_sz, num)

    interval = int((min_sz+1)/2)
    bds_h, bds_v = [], []
    resh = para_bds(diag_h_diff, ft, interval, percentile)
    for _bds_h in resh:
        bds_h += select_bds(_bds_h, ft, interval)
    resv = para_bds(diag_v_diff, ft, interval, percentile)
    for _bds_v in resv:
        bds_v += (n-1- np.array(select_bds(_bds_v, ft[::-1], interval))).tolist()
    #for i in range(len(diag_h_diff)):
        #_bds_h = _get_bds(diag_h_diff[i], diag_h_diff_shuffle[i], percentile)
        #_bds_h = select_bds(_bds_h, ft, interval)
        #bds_h += _bds_h
        #_bds_v = _get_bds(diag_v_diff[i], diag_v_diff_shuffle[i], percentile)
        #_bds_v = (n-1- np.array(select_bds(_bds_v, ft[::-1], interval))).tolist()
        #bds_v += _bds_v
        
    bds = np.unique(bds_h+bds_v)
    return select_bds(bds, ft, interval)

if __name__ == '__main__':
    pass
