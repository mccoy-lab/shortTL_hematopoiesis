import numpy as np
from scipy.stats import binom, chi2
from scipy.special import gamma
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.integrate import quad


def freq_age_density(s, mu=1e-9, N=1e5, tau=1):
    """Density of VAF for individual at a specific age."""
    assert s >= 0
    assert mu > 0
    assert N > 0
    assert tau > 0
    theta = N * tau * mu
    assert (theta > 0) and (theta < 1)
    ntilde = lambda t: (np.exp(s * t) - 1) / (s * tau) if s > 0 else t / tau
    rho = (
        lambda f, t: (2 * N) ** theta
        / (f ** (1 - theta) * (1 - 2 * f) ** (1 + theta))
        * np.exp((-2 * N * f) / (ntilde(t) * (1 - 2 * f)))
        / (gamma(theta) * ntilde(t) ** theta)
    )
    return rho


def downsample(n, x, max_n=100):
    """Downsample super high clinical read depths"""
    p = x / n
    x = binom.rvs(n=max_n, p=p)
    return max_n, x


def macgregor_single_likelihood(
    r, d, t=30, s=0.01, mu=1e-9, read_threshold=100, **kwargs
):
    """Likelihood for a sample with actual read information."""
    assert r >= 0
    assert d >= r
    assert t >= 0
    assert read_threshold > 0
    x = freq_age_density(s=s, mu=mu, **kwargs)
    if d > read_threshold:
        d2, r2 = downsample(n=d, x=r)
        integrand = quad(lambda f: binom.pmf(r2, n=d2, p=f) * x(f, t), 0, 0.5)
    else:
        integrand = quad(lambda f: binom.pmf(r, n=d, p=f) * x(f, t), 0, 0.5)
    return integrand


def macgregor_unobserved_likelihood(t=30, s=1e-2, mu=1e-9, f0=1e-2, **kwargs):
    """Likelihood for an unobserved mutation."""
    assert f0 > 0
    integrand = quad(lambda f: freq_age_density(s=s, mu=mu, **kwargs)(f, t), 0.0, f0)
    return integrand


def macgregor_total_loglik(rs, ds, ts, f0=1e-2, s=1e-4, mu=1e-9, seed=42, **kwargs):
    """The total log-likelihood for observed and unobserved mutations."""
    assert rs.ndim == 1
    assert ds.ndim == 1
    assert rs.size == ds.size
    assert ts.size == rs.size
    assert seed > 0
    # Seed is only really relevant to downsampling ...
    np.random.seed(seed)
    loglik = 0.0
    for r, d, t in zip(rs, ds, ts):
        if (r > 0) & (d > 0):
            loglik += np.log(macgregor_single_likelihood(r, d, t=t, s=s, mu=mu)[0])
        else:
            loglik += np.log(macgregor_unobserved_likelihood(t=t, s=s, mu=mu, f0=f0)[0])
    return loglik


def reverse_cumulative(r, d, **kwargs):
    """Reverse cumulative distribution at a given age threshold."""
    assert (r != 0) & (d != 0)
    assert d >= r
    cum_density = np.sum(
        [macgregor_single_likelihood(rp, ds, **kwargs) for rp in range(r, d)]
    )
    return cum_density


def mle_inference_sonly(rs, ds, ts, mu=1e-9, **kwargs):
    """Single-variate MLE inference for only selection; mutation rate is provided."""
    f = lambda x: -macgregor_total_loglik(rs, ds, ts, mu=mu, s=x, **kwargs)
    opt_res = minimize_scalar(f, bounds=(1e-10, 0.5), options={"disp": True})
    return opt_res


def mle_inference(rs, ds, ts, **kwargs):
    """Numerical estimation of the MLE for joint mutation and selection parameters."""
    f = lambda x: -macgregor_total_loglik(rs, ds, ts, mu=x[0], s=x[1], **kwargs)
    opt_res = minimize(
        f,
        x0=[1e-9, 1e-2],
        bounds=[(1e-20, 1e-6), (1e-10, 0.5)],
        method="Powell",
        options={"disp": True},
    )
    return opt_res


def confint_s(rs, ds, ts, s_hat, mu=1e-9, q=0.95, df=1, **kwargs):
    """Estimate the confidence interval for s using a profile likelihood."""
    wilks = lambda x: 2 * (
        macgregor_total_loglik(rs, ds, ts, mu=mu, s=s_hat, **kwargs)
        - macgregor_total_loglik(rs, ds, ts, mu=mu, s=x, **kwargs)
    )
    qval = chi2.ppf(q, df=df)
    # NOTE: these are typically use to avoid roundoff errors in integration
    lower = brentq(lambda x: wilks(x) - qval, 1e-9, s_hat)
    upper = brentq(lambda x: wilks(x) - qval, s_hat, 0.4 - 1e-9)
    return (lower, upper)


def confint_mu(rs, ds, ts, mu_hat, s=1e-2, q=0.95, df=1, **kwargs):
    """Estimate the confidence interval for mu using a profile likelihood."""
    wilks = lambda x: 2 * (
        macgregor_total_loglik(rs, ds, ts, mu=mu_hat, s=s, **kwargs)
        - macgregor_total_loglik(rs, ds, ts, mu=x, s=s, **kwargs)
    )
    qval = chi2.ppf(q, df=df)
    # NOTE: these are typically use to avoid roundoff errors in integration
    lower = brentq(lambda x: wilks(x) - qval, 1e-12, mu_hat)
    upper = brentq(lambda x: wilks(x) - qval, mu_hat, 1e-6)
    return (lower, upper)


def clonesize_density(s=1e-2, mu=1e-9, N=1e5, tau=1):
    """Obtain the density of clone sizes as a function of n."""
    pass
