import numpy as np
import scipy.stats as stat


def _parse_file(file_handle):
    if isinstance(file_handle, str):
        return open(file_handle, 'r')
    else:
        return file_handle


def _extract_data(file_handle, column, warmup):
    # Parse file name
    file_handle = _parse_file(file_handle)
    file_handle.seek(0)

    # Read data from file and copy specified column to 1-D Numpy array
    all_data = np.genfromtxt(file_handle, delimiter=',', comments='#')
    data = all_data[:, column]

    # Check for NaN in data
    if np.isnan(np.min(data)):
        raise ValueError("Data contains NaN.")

    # Make sure data size > warmup size
    if data.size < warmup:
        raise ValueError("Warmup length is greater than sample size.")

    # Return warmup and production arrays
    return data[:warmup], data[warmup:]


def _detect_warmup(file_handle, column):
    # Extract data
    _, data = _extract_data(file_handle, column, 0)

    # Pre-block-average in chunks of 5 steps and truncate integer counts of those steps
    m = 5
    len_data = len(data)
    if len_data % m != 0:
        data_pad = np.pad(data, (0, m - len_data % m), mode='mean', stat_length=len_data % m)
        data_block_avg = np.mean(data_pad.reshape(int(len_data/m) + 1, m), axis=1)
    else:
        data_block_avg = np.mean(data.reshape(int(len_data/m), m), axis=1)

    # Compute standard error of mean
    std_err_mean_list = []
    full_size = len(data_block_avg)
    while True:
        std_err_mean = stat.sem(data_block_avg)
        std_err_mean_list.append(std_err_mean)
        len_data_remaining = len(data_block_avg)
        if len_data_remaining <= 0.25 * full_size or len_data_remaining < 5:
            break
        data_block_avg = data_block_avg[1:]

    # Find index that minimizes the variance
    index_min_var = np.argmin(std_err_mean_list)*m

    # Divide data into warmup and production
    warmup_data = data[:index_min_var]
    production_data = data[index_min_var:]

    return warmup_data, production_data, index_min_var


def compute_stats(data):
    # Compute statistics
    nobs, minmax, mean, variance, skewness, kurtosis = stat.describe(data)

    # Compute standard error of mean
    std_err_mean = stat.sem(data)

    # Compute autocorrelation function
    data_shifted = data - mean
    correlation = np.correlate(data_shifted, data_shifted, mode='same')/variance
    autocorrelation = correlation[int(correlation.size/2):]
    autocorrelation = autocorrelation/np.arange(nobs - 1, nobs - 1 - autocorrelation.size, -1)

    # Choose where to cutoff the autocorrelation time sum
    cutoff = autocorrelation.size
    i = 0
    while i < cutoff:
        if autocorrelation[i] < np.sqrt(2.0/(nobs - i)):
            cutoff = np.minimum(cutoff, 5*i)
        i += 1

    # Compute correlation time
    kappa = 1.0 + 2.0*np.sum(autocorrelation[1:int(2.0*cutoff/5.0)])

    # Update standard error of the mean for a correlation correction
    semcc = std_err_mean*np.sqrt(kappa)

    return nobs, minmax, mean, semcc, kappa, variance, autocorrelation


def compute_stats_from_file(file_handle, column):
    _, data, _ = _detect_warmup(file_handle, column)
    return compute_stats(data)
