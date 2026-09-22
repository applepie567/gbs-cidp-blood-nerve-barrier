import numpy as np
from scipy import stats


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ix = np.flatnonzero(np.isfinite(p))
    ix = ix[np.argsort(p[ix])]
    if len(ix):
        q[ix] = np.minimum(1, np.minimum.accumulate((p[ix] * len(ix) / np.arange(1, len(ix) + 1))[::-1])[::-1])
    return q


def hc3(y, meta, covariates=()):
    y = np.asarray(y, float)
    x = [np.ones(len(y)), meta.disease.eq('CIDP').to_numpy(float)]
    for c in covariates:
        if c == 'center':
            x.extend(meta.center.eq(l).to_numpy(float) for l in sorted(meta.center.unique())[1:])
        else:
            x.append(meta[c].to_numpy(float))
    X = np.column_stack(x)
    rank = np.linalg.matrix_rank(X)
    if rank != X.shape[1] or len(y) <= rank:
        return dict(delta=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan,
                    p_hc3=np.nan, residual_df=len(y)-rank, estimable=False)
    bread = np.linalg.inv(X.T @ X)
    beta = bread @ X.T @ y
    h = np.sum((X @ bread) * X, axis=1)
    if np.max(h) >= 1-1e-10:
        raise ValueError('HC3 leverage is one')
    residual = y - X @ beta
    covariance = bread @ (X.T @ np.diag((residual/(1-h))**2) @ X) @ bread
    se = np.sqrt(max(0, covariance[1, 1]))
    df = len(y) - rank
    ci = stats.t.ppf(.975, df) * se
    return dict(delta=float(beta[1]), se=float(se), ci_low=float(beta[1]-ci),
                ci_high=float(beta[1]+ci), p_hc3=float(2*stats.t.sf(abs(beta[1]/se), df)),
                residual_df=int(df), estimable=True)


def meta_two(beta, se):
    beta, se = np.asarray(beta, float), np.asarray(se, float)
    k = len(beta)
    v = se**2
    w = 1/v
    fixed = np.average(beta, weights=w)
    fixed_se = np.sqrt(1/w.sum())
    Q = np.sum(w*(beta-fixed)**2)
    tau2 = max(0, (Q-(k-1))/(w.sum()-(w*w).sum()/w.sum()))
    wr = 1/(v+tau2)
    random = np.average(beta, weights=wr)
    random_se = np.sqrt(1/wr.sum())
    hk_se = np.sqrt(max(1, np.sum(wr*(beta-random)**2)/(k-1))/wr.sum())
    z = stats.norm.ppf(.975)
    t = stats.t.ppf(.975, k-1)
    return dict(fixed_beta=fixed, fixed_se=fixed_se, fixed_low=fixed-z*fixed_se,
                fixed_high=fixed+z*fixed_se, fixed_p=2*stats.norm.sf(abs(fixed/fixed_se)),
                random_beta=random, random_se=random_se, random_low=random-z*random_se,
                random_high=random+z*random_se, random_p=2*stats.norm.sf(abs(random/random_se)),
                mhk_se=hk_se, mhk_low=random-t*hk_se, mhk_high=random+t*hk_se,
                mhk_p=2*stats.t.sf(abs(random/hk_se), k-1), Q=Q, tau2=tau2,
                I2_percent=max(0,(Q-k+1)/Q)*100 if Q > 0 else 0)
