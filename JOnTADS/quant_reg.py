import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from qpsolvers import solve_qp
import scipy
from scipy.spatial.distance import pdist, squareform

def kernel(x1, x2, sigma):
    return np.exp((-((x1-x2)**2))/(2*sigma**2))
def d1kernel(x1, x2, sigma):
    return (x2-x1)/sigma**2*kernel(x1, x2, sigma)
def nqr(x, y, tau, sigma, penalty):
    K = np.exp(-squareform(pdist(np.expand_dims(x, axis=1), 'euclidean')) ** 2 / (2 * sigma ** 2))
    D2K = (np.tile(x, (len(x), 1)).T-x)/sigma**2*K
    D1K = D2K.T
    D1D2K = (1/sigma**2-(np.tile(x, (len(x), 1)).T-x)**2/sigma**4)*K
    P = np.vstack((np.hstack((K,D1K)), np.hstack((D2K,D1D2K))))
    P += 1e-6*np.eye(2*len(x))
    q = np.hstack((-y, np.zeros(len(x))))
    C = 1/(penalty*len(x))
    lb = np.hstack((np.ones(len(x))*C*(tau-1), np.zeros(len(x))))
    ub = np.hstack((np.ones(len(x))*C*tau, np.ones(len(x))*np.Inf))
    A = np.hstack((np.ones(len(x)), np.zeros(len(x))))
    b = np.zeros(1)
    res = solve_qp(P=P, q=q, A=A, b=b, lb=lb, ub=ub, solver="quadprog")
    alpha = res[:len(x)]
    beta = res[len(x):]
    for i in range(len(x)):
        if np.abs(alpha[i]-C*(tau-1))>1e-2*C*tau and np.abs(alpha[i]-C*tau)>1e-2*C*tau:
            b = y[i] - sum(alpha*K[i, :] + beta*D1K[i, :])
    return alpha, beta, b
def force_decreasing(curve):
    curve[curve<0] = 0
    for i in range(1, len(curve)):
        if curve[i] > curve[i-1]:
            curve[i] = curve[i-1]
    return curve
def _quant_reg(x, y, max_sz, min_sz, tau, lmda):
    fitted = []
    sigma = (max_sz-min_sz)/10
    alpha, beta, b = nqr(x, y, tau, sigma, lmda)
    for i in range(min_sz, max_sz+1):
        fitted.append(np.sum(alpha*kernel(i, x, sigma)+beta*d1kernel(i, x, sigma))+b)
    return force_decreasing(np.array(fitted))

def quant_reg(tads, tads_shuffle, max_sz, min_sz, significance, lmda):
    def helper(tads_shuffle, min_sz, max_sz):
        ps = np.percentile(tads_shuffle[:, 3], [60, 90])
        idx1 = np.random.choice(np.arange(len(tads_shuffle))[tads_shuffle[:, 3]<=ps[0]], 100)
        idx2 = np.random.choice(np.arange(len(tads_shuffle))[np.array(tads_shuffle[:, 3]>ps[0]) & np.array(tads_shuffle[:, 3]<ps[1])], 100)
        idx3 = np.random.choice(np.arange(len(tads_shuffle))[tads_shuffle[:, 3]>=ps[1]], 100)
        idx4 = []
        #p = (1-50/len(tads_shuffle))*100.
        for i in range(min_sz, max_sz+1):
            scores = tads_shuffle[:, 2][tads_shuffle[:, 3]==i]
            if len(scores) != 0:
                # idx4 += (np.arange(len(tads_shuffle))[np.array(tads_shuffle[:, 2]>np.percentile(scores, p)) & np.array(tads_shuffle[:, 3]==i)]).tolist()
                idx4 += (np.arange(len(tads_shuffle))[np.array(tads_shuffle[:, 2]==np.max(scores)) & np.array(tads_shuffle[:, 3]==i)]).tolist()
        idx4 = np.array(idx4)
        idx = np.array(idx1.tolist()+idx2.tolist()+idx3.tolist()+idx4.tolist())
        return tads_shuffle[:, 3][idx], tads_shuffle[:, 2][idx]
    
    size_p, score_p = helper(tads_shuffle, min_sz, max_sz)
    
    left = 0.5
    right = 1
    diff = 1
    idx_out, curve_out = [], []
    for _ in range(10):
        mid = (left+right)/2
        curve = _quant_reg(size_p, score_p, max_sz, min_sz, mid, lmda)
        if (np.max(curve)-np.min(curve)) > (np.std(tads_shuffle[:, 2]) * 0.15):
            num, idx = above_curve(tads, curve, min_sz)
            num_shuffle, _ = above_curve(tads_shuffle, curve, min_sz)
        else:
            while True:
                print('!')
                size_p, score_p = helper(tads_shuffle, min_sz, max_sz)
                curve = _quant_reg(size_p, score_p, max_sz, min_sz, mid+np.random.uniform(-0.0005, 0.0005, 1)[0], lmda)
                if (np.max(curve)-np.min(curve)) > (np.std(tads_shuffle[:, 2]) * 0.15):
                    num, idx = above_curve(tads, curve, min_sz)
                    num_shuffle, _ = above_curve(tads_shuffle, curve, min_sz)
                    break
        if num_shuffle / (num_shuffle+num) > significance:
            left = mid
        elif num_shuffle / (num+num_shuffle) < significance:
            right = mid
        else:
            #plt.plot(tads[:, 3], tads[:, 2], '.')
            #plt.plot(np.arange(min_sz, max_sz+1), curve)
            #import datetime
            #time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
            #plt.show()
            #plt.close()
            return idx, curve
        if np.abs(num_shuffle / (num_shuffle+num) - significance) < diff:
            diff = np.abs(num_shuffle / (num_shuffle+num) - significance)
            idx_out = idx
            curve_out = curve

    #plt.plot(tads[:, 3], tads[:, 2], '.')
    #plt.plot(np.arange(min_sz, max_sz+1), curve_out)
    #import datetime
    #time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    #plt.show()
    #plt.close()
    if left == 0.1 or right == 1:
        print('TADs might not exist.')
    return idx_out, curve_out
    
def above_curve(tads, curve, min_sz):
    size = tads[:, 3]
    score = tads[:, 2]
    ans = 0
    idx = []
    for i, sz in enumerate(np.array(size, dtype=int)):
        if score[i] >= curve[sz-min_sz]:
            ans += 1
            idx.append(i)
    ans = len(np.unique(tads[idx, :2]))
    return ans, idx

if __name__ == '__main__':
    pass
