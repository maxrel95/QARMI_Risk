# ==========================================================================
# Quantitative Asset and Risk Management
# Risk Management Project : Comparison of VaR estimation techniques
#
# Author: Maxime Borel
# Date: August 2020
# ==========================================================================

import numpy as np
import pandas as pd
from scipy import stats
from scipy import optimize
import matplotlib.pyplot as plt
from statsmodels.stats.diagnostic import acorr_ljungbox as ljungbox
from statsmodels.tsa.stattools import pacf
import seaborn as sns
from arch import arch_model


def failurerate(rt, var, conflvlvar, conflvltest):
    ret=-rt
    p=1 - conflvlvar
    T=ret.shape[0]
    ret[rt > var[:, 0]]=1
    ret[ret != 1]=0
    N=ret.sum()
    failrate=N / T

    if failrate != 0:
        LR=2 * ((N * np.log(failrate) + np.log(1 - failrate) * (T - N))
                - (N * np.log(p)) + np.log(1 - p) * (T - N))
    else:
        LR=0

    criticalvalue=stats.chi2.ppf(conflvltest, df=1)

    if LR > criticalvalue:
        resultat=1
    else:
        resultat=0

    ff=0
    fs=0
    sf=0
    ss=0

    for i in range(1, T):
        if ret[i] == 1 and ret[i - 1] == 1:
            ff=ff + 1
        elif ret[i] == 1 and ret[i - 1] == 0:
            fs=fs + 1
        elif ret[i] == 0 and ret[i - 1] == 1:
            sf=sf + 1
        elif ret[i] == 0 and ret[i - 1] == 0:
            ss=ss + 1

    if failrate == 0:
        LRCCI=0
    else:
        pi0=sf / (ss + sf)
        pi1=ff / (fs + ff)
        pi2=(sf + ff) / (ss + ff + sf + fs)
        upper=((1 - pi2) ** (ss + fs)) * (pi2 ** (sf + ff))
        lower=((1 - pi0) ** ss) * (pi0 ** sf) * ((1 - pi1) ** fs) * (pi1 ** ff)
        LRCCI=-2 * np.log(upper / lower)

    CC=LR + LRCCI
    result=np.array([failrate, LR, resultat, LRCCI, CC])

    return result


def test(rt, var):
    ret=-rt
    T=ret.shape[0]
    ret[rt > var[:, 0]]=1
    ret[ret != 1]=0
    N=ret.sum()
    failrate=N / T

    return failrate


df = pd.read_excel('classeur5.xlsx', index_col=0)
price = df['Adj Close']
returns = np.log(price).diff().dropna()
T = returns.shape[0]
evalwindow = np.arange(T-3344+1, T)
estimwindow = np.arange(0, T-3344+1)

rt_eval = returns.iloc[evalwindow]
rt_estim = returns.iloc[estimwindow]
t = rt_eval.shape[0]
varlvl = np.array([0.05, 0.01, 0.005, 0.001, 0.0001])
cv = np.ones((1, t))*stats.chi2.ppf(0.95, 1)
cv2 = np.ones((1, t))*stats.chi2.ppf(0.95, 2)
VarNames1 = ['VaR estimate', 'Failure rate', 'LR tstat', 'LRCCI', 'CV5', 'CC', 'CV25']
VarNames2 = ['Average VaR estimate', 'Failure rate', 'LR', 'LRCCI', 'CV5', 'CC', 'CV225']

# descriptive statistics
sumstat_rt = pd.DataFrame([returns.mean(), returns.median(), returns.std(), returns.min(),
                           returns.max(), stats.skew(returns), stats.kurtosis(returns)],
                          columns=['Descriptive statistic of returns'],
                          index=['Mean', 'Median', 'Standard Deviation', 'Minimum',
                                 'Maximum', 'Skewness', 'Excess Kurtosis'])
print(sumstat_rt)
sumstat_rt.to_excel('sumstat_rt.xlsx')

# check the autocorrelation in the squared eps for arch effect
centered_rt = (returns-returns.mean())**2
centered_rt.autocorr(lag=10)
tempotable = pd.DataFrame(np.zeros((10, 4)), index=['1', '2', '3', '4', '5', '6',
                                                    '7', '8', '9', '10'],
                          columns=['ACF', 'PACF', 'Ljungbox', 'pvalue'])
for i in range(1, 10):
    tempotable.values[i, 0] = centered_rt.autocorr(lag=i)
tempotable.values[:, 1] = pacf(centered_rt, nlags=10)[1:]
tempotable.values[:, 2:] = ljungbox(centered_rt, lags=10, return_df=True)
print(tempotable)
tempotable.to_excel('tempotable.xlsx')

fig1 = plt.figure()
plt.plot(returns)
plt.ylabel('daily Returns')
plt.xlabel('Time')
plt.title('Daily return of the SMI')
plt.show()
fig1.savefig('dailyreturn.png')

fig2 = plt.figure()
plt.hist(returns, bins=30, density=True)
sns.kdeplot(returns, color='red')
plt.ylabel('# number of occurence')
plt.xlabel('return')
plt.title('Histogram of the return')
plt.show()
fig2.savefig('histreturn.png')

## Static models for VaR estimation
# Fitting Normal distribution and testing parameters


def loglik_norm(params, rt):

    if params[1] <= 0:
        return 11111110

    mu = params[0]
    sigma = params[1]
    loglik = stats.norm.logpdf(rt, loc=mu, scale=sigma).sum()

    return - loglik


norm_optim = optimize.minimize(loglik_norm, x0=np.array([0, 1]),
                      args=(rt_estim,), method='BFGS', options={'disp': True})

tstat = norm_optim.x/np.sqrt(np.diag(norm_optim.hess_inv))

static_norm = pd.DataFrame(np.array([[norm_optim.x[0], tstat[0]],
                                    [norm_optim.x[1], tstat[1]]]), index=['Mu', 'sigma'],
                           columns=['Estimates', 'tStat'])
print(static_norm)
static_norm.to_excel('static_norm.xlsx')

# Fitting the student t distribution and testing parameters


def loglik_student(params, rt):
    mu = params[0]
    sigma = params[1]
    nu = params[2]

    if sigma <= 0:
        return 10000000

    if nu <= 0:
        return 10000000

    loglik = stats.t.logpdf(rt, df=nu, loc=mu, scale=sigma).sum()

    return - loglik


student_optim = optimize.minimize(loglik_student, x0=np.array([0, 1, 3]),
                                  args=(rt_estim,), method='BFGS', options={'disp': True})

tstat = student_optim.x/np.sqrt(np.diag(student_optim.hess_inv))

static_t = pd.DataFrame(np.array([[student_optim.x[0], tstat[0]],
                                  [student_optim.x[1], tstat[1]],
                                 [student_optim.x[2], tstat[2]]]), index=['Mu', 'sigma', 'nu'],
                        columns=['Estimates', 'tStat'])
print(static_t)
static_t.to_excel('static_t.xlsx')

# Fitting the skewed t distribution and testing parameters


def loglik_skewt(params, rt):
    mu = params[0]
    sigma = params[1]
    nu = params[2]
    lamb = params[3]

    if sigma <= 0:
        return 10000000
    if nu <= 0:
        return 10000000
    if lamb >= 1 or lamb <= -1:
        return 1000000

    loglik = stats.nct.logpdf(rt, df=nu, nc=lamb, loc=mu, scale=sigma).sum()

    return -loglik


skewed_optim = optimize.minimize(loglik_skewt, x0=np.array([0, 1, 3, 0]),
                                  args=(rt_estim,), method='BFGS', options={'disp': True})

tstat = skewed_optim.x/np.sqrt(np.diag(skewed_optim.hess_inv))

static_skewed = pd.DataFrame(np.array([[skewed_optim.x[0], tstat[0]],
                                       [skewed_optim.x[1], tstat[1]],
                                       [skewed_optim.x[2], tstat[2]],
                                       [skewed_optim.x[3], tstat[3]]]),
                             index=['Mu', 'sigma', 'nu', 'lambda'], columns=['Estimates', 'tStat'])
print(static_skewed)
static_skewed.to_excel('static_skewed.xlsx')

# computing the unconditional VaR for all distribution
# Normal
n = varlvl.shape[0]
FR_s_n = np.zeros((n, 1))
FR_s_t = np.zeros((n, 1))
FR_s_sk = np.zeros((n, 1))
LR_s_n = np.zeros((n, 1))
LR_s_t = np.zeros((n, 1))
LR_s_sk = np.zeros((n, 1))
rlt_s_n = np.zeros((n, 1))
rlt_s_t = np.zeros((n, 1))
rlt_s_sk = np.zeros((n, 1))
LRCCI_s_n = np.zeros((n, 1))
LRCCI_s_t = np.zeros((n, 1))
LRCCI_s_sk = np.zeros((n, 1))
CC_s_n = np.zeros((n, 1))
CC_s_t = np.zeros((n, 1))
CC_s_sk = np.zeros((n, 1))
var_s_n = np.zeros((n, 1))
var_s_t = np.zeros((n, 1))
var_s_sk = np.zeros((n, 1))

for i in varlvl:
    var_s_n[varlvl == i] = -(np.exp(static_norm.values[0, 0] + static_norm.values[1, 0]*
                                    stats.norm.ppf(i))-1)
    vecvar = np.ones((t, 1))*var_s_n[varlvl == i]
    [FR_s_n[varlvl == i], LR_s_n[varlvl == i], rlt_s_n[varlvl == i], LRCCI_s_n[varlvl == i],
     CC_s_n[varlvl == i]] = failurerate(rt_eval.values, vecvar, 1-i, 0.95)


    var_s_t[varlvl == i] = -(np.exp(static_t.values[0, 0] + static_t.values[1, 0] *
                                  stats.t.ppf(i, df=static_t.values[-1, 0])) - 1)
    vecvar1 = np.ones((t, 1))*var_s_t[varlvl == i]
    [FR_s_t[varlvl == i], LR_s_t[varlvl == i], rlt_s_t[varlvl == i], LRCCI_s_t[varlvl == i],
     CC_s_t[varlvl == i]] = failurerate(rt_eval.values, vecvar1, 1-i, 0.95)


    var_s_sk[varlvl == i] = -(np.exp(static_skewed.values[0, 0] + static_skewed.values[1, 0] *
                                  stats.nct.ppf(i, df=static_skewed.values[-2, 0],
                                                nc=static_skewed.values[-1, 0])) - 1)
    vecvar2 = np.ones((t, 1))*var_s_sk[varlvl == i]
    [FR_s_sk[varlvl == i], LR_s_sk[varlvl == i], rlt_s_sk[varlvl == i], LRCCI_s_sk[varlvl == i],
     CC_s_sk[varlvl == i]] = failurerate(rt_eval.values, vecvar2, 1-i, 0.95)

# graph several var static
fig3 = plt.figure()
plt.plot(rt_eval)
plt.plot(rt_eval.index, -np.ones((t, 1))*var_s_n[0], 'r')
plt.title('Static VaR at 5% Normal')
plt.xlabel('Time')
plt.ylabel('Log returns')
plt.legend(['Return', 'VaR 5%'])
plt.show()
fig3.savefig('varnormalstatic.png')

fig4 = plt.figure()
plt.plot(rt_eval)
plt.plot(rt_eval.index, -np.ones((t, 1))*var_s_t[0], 'r')
plt.title('Static VaR at 5% Student')
plt.xlabel('Time')
plt.ylabel('Log returns')
plt.legend(['Return', 'VaR 5%'])
plt.show()
fig4.savefig('vartstatic.png')

fig5 = plt.figure()
plt.plot(rt_eval)
plt.plot(rt_eval.index, -np.ones((t, 1))*var_s_sk[0], 'r')
plt.title('Static VaR at 5% Skewed T')
plt.xlabel('Time')
plt.ylabel('Log returns')
plt.legend(['Return', 'VaR 5%'])
plt.show()
fig5.savefig('varSKstatic.png')

# Table of results static
tablename = ['VaR Estimate', 'Failure rate', 'LR test', 'LRCCI', 'CC']
tablename1 = ['VaR Average', 'Failure rate', 'LR test', 'LRCCI', 'CC']

result_s_n = pd.DataFrame(np.array(
    [var_s_n*100, FR_s_n*t, LR_s_n, LRCCI_s_n, CC_s_n]).reshape((5, 5)).T,
                          index=varlvl,
                          columns=tablename)
result_s_n.to_excel('result_s_n.xlsx')
result_s_t = pd.DataFrame(np.array(
    [var_s_t*100, FR_s_t*t, LR_s_t, LRCCI_s_t, CC_s_t]).reshape((5, 5)).T,
                          index=varlvl,
                          columns=tablename)
result_s_t.to_excel('result_s_t.xlsx')
result_s_sk = pd.DataFrame(np.array([var_s_sk*100, FR_s_sk*t, LR_s_sk,
                                     LRCCI_s_sk, CC_s_sk]).reshape((5, 5)).T,
                          index=varlvl,
                          columns=tablename)
result_s_sk.to_excel('result_s_sk.xlsx')

## GARCH Model
FR_g_n = np.zeros((n, 1))
FR_g_t = np.zeros((n, 1))
FR_g_sk = np.zeros((n, 1))
LR_g_n = np.zeros((n, 1))
LR_g_t = np.zeros((n, 1))
LR_g_sk = np.zeros((n, 1))
rlt_g_n = np.zeros((n, 1))
rlt_g_t = np.zeros((n, 1))
rlt_g_sk = np.zeros((n, 1))
LRCCI_g_n = np.zeros((n, 1))
LRCCI_g_t = np.zeros((n, 1))
LRCCI_g_sk = np.zeros((n, 1))
CC_g_n = np.zeros((n, 1))
CC_g_t = np.zeros((n, 1))
CC_g_sk = np.zeros((n, 1))
var_g_n = np.zeros((t, n))
var_g_t = np.zeros((t, n))
var_g_sk = np.zeros((t, n))

eps_garch_eval = rt_eval-rt_estim.mean()
eps_garch_estim = rt_estim-rt_estim.mean()
garch_n = arch_model(100*eps_garch_estim, mean='zero', vol='GARCH', p=1, q=1, dist='Normal')
garch_n_result = garch_n.fit()
print(garch_n_result)
garch_t = arch_model(100*eps_garch_estim, mean='zero', vol='GARCH', p=1, q=1, dist='t')
garch_t_result = garch_t.fit()
print(garch_t_result)

plt.plot(garch_n_result.resid**2)
plt.title('GARCH eps squared')
plt.xlabel('Date')
plt.ylabel('eps squared')
plt.show()

garch_n_param = pd.concat([garch_n_result.params, garch_n_result.tvalues], axis=1)
garch_n_param.to_excel('para_n_garch.xlsx')

sigma_normal = np.zeros((t, 1))
sigma_normal[0, 0] = garch_n_param.values[1, 0]+\
                     garch_n_param.values[2, 0]*garch_n_result.resid[-1]**2+\
                     garch_n_param.values[-1, 0]*garch_n_result.conditional_volatility[-1]**2
for i in range(1, t):
    sigma_normal[i, 0] = garch_n_param.values[1, 0]+\
                       garch_n_param.values[2, 0]*eps_garch[i-1]**2+\
                       garch_n_param.values[-1, 0]*sigma_normal[i-1, 0]
sigma_normal = np.vstack([garch_n_result.conditional_volatility[-1]**2, sigma_normal[:-1]])

for i in varlvl:
    var_g_n[:, varlvl == i] = -(np.exp(garch_n_param.values[0, 0] + np.sqrt(sigma_normal)*
                                       stats.norm.ppf(i))-1)
    vecvar = np.ones((t, 1))*var_g_n[:, varlvl == i]
    [FR_g_n[varlvl == i], LR_g_n[varlvl == i], rlt_g_n[varlvl == i], LRCCI_g_n[varlvl == i],
     CC_g_n[varlvl == i]] = failurerate(rt_eval.values, vecvar, 1-i, 0.95)















