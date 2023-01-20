import sys
import numpy as np

def calculate_full_fisher(*args):

    if len(args) == 8:
        f, segDuration, csd_f, sigma_sq_inv_f,combined_antenna_response,t_delay, H_f, notch = args

    print('Calculating full-fisher matrix for frequency {} Hz.       '.format(f))
    sys.stdout.write("\033[F")

    gamma_star = combined_antenna_response * np.exp(-1j*(2*np.pi*f)*t_delay)
    fisher_f =  (np.dot(np.multiply(sigma_sq_inv_f, np.conjugate(gamma_star).T), gamma_star) * segDuration ** 2) * H_f ** 2

    if not(notch):
        return np.zeros_like(fisher_f)
    else:
        return fisher_f
