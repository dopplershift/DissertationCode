import numpy as np

try:
    from scipy import linalg
except ImportError:
    from numpy import linalg

# A collection of computational tools
def linear_regression(x, y, weights=None, mask=None, intercept=True):
    # Include intercept term if necessary
    if intercept:
        x = np.concatenate((np.ones_like(y).reshape(-1, 1), x), 1)

    # Save x-values so we can return the value of the fit at all
    # points
    full_x = x

    # If no weights given, just use ones
    if weights is None:
        weights = np.ones_like(y).reshape(-1, 1)

    # Mask everything if necessary
    if mask is not None:
        weights = weights[mask]
        x = x[mask]
        y = y[mask]

    # Do the least squares fit
    coeffs,resid,rank,sigma = linalg.lstsq(weights * x, weights * y)

    # Calculate the value of the fit at all x values
    fit = np.dot(full_x, coeffs).squeeze()
    coeffs = coeffs.squeeze()
    return coeffs, fit

def power_law_fit(x, y, weights=None, mask=None):
    # A power law fit is just a linear fit in the log-transformed domain.
    # So just transform the data and pass to a linear fit.
    logX = np.log(x)
    logY = np.log(y)
    coeffs, fit = linear_regression(logX, logY, weights, mask, intercept=True)

    # Need to transform the intercept coefficient into the multiplicative
    # factor.
    coeffs[0] = np.exp(coeffs[0])

    # Untransform the fit values as well.
    return coeffs, np.exp(fit)
