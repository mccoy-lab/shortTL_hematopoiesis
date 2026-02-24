import numpy as np
from scipy.stats import truncnorm
from scipy.optimize import minimize


def fit_age_dist_quantiles(
    xs=np.array([19.83847, 24.07392, 37.93292, 54.02875, 65.57153, 77.43737, 84.799]),
    qs=np.array([0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]),
    min_age=18,
    max_age=106,
    **kwargs
):
    """Fit the normal age distribution based on quantiles."""
    assert min_age > 0
    assert max_age > min_age
    f = lambda x: np.sum(
        (
            truncnorm.cdf(
                xs,
                (min_age - x[0]) / x[1],
                (max_age - x[0]) / x[1],
                loc=x[0],
                scale=x[1],
            )
            - qs
        )
        ** 2
    )
    opt_res = minimize(
        f, x0=[np.mean(xs), np.std(xs)], bounds=[(min_age, max_age), (0, 1e5)], **kwargs
    )
    # Return the actual truncated normal distribution...
    a = (min_age - opt_res.x[0]) / opt_res.x[1]
    b = (max_age - opt_res.x[0]) / opt_res.x[1]
    rv = truncnorm(a, b, loc=opt_res.x[0], scale=opt_res.x[1])
    return rv
