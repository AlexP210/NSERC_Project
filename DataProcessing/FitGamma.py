import math as m
import numpy as np

def gamma_probability(x, theta, k, mu):
    return 1/( m.factorial(1-k) * theta ** k) * (x - mu)**(k-1)*math.e**(-(x-mu)/theta)

def log_gamma_probability(x, theta, k, mu):
    return -m.log(m.factorial(1-k)) - k*m.log(theta) + (k-1)*m.log( x-mu ) - (x-mu)/theta

def fit_gamma(X)
    theta_hypotheses = list(np.linspace(0, 1, 100))
    k_hypotheses = list(np.linspace(0, 1, 100))
    mu_hypotheses = list(np.linspace(0, 1, 100))

    log_likelihood_space = np.ones((len(theta_hypotheses), len(k_hypotheses), len(mu_hypotheses)))
    max_log_likelihood = float("inf")
    max_log_likelihood_idx = (None, None, None)
    for theta_idx in range(len(theta_hypotheses)):
        for k_idx in range(len(k_hypotheses)):
            for mu_idx in range(len(mu_hypotheses)):
                # Set the hypothesis
                theta = theta_hypotheses[theta_idx]
                k = k_hypotheses[k_idx]
                mu = mu_hypotheses[mu_idx]
                # Get the log-likelihood
                log_likelihood = sum( [log_gamma_probability(x, theta, k, mu) for x in X] )
                if log_likelihood > max_log_likelihood: 
                    max_log_likelihood = log_likelihood
                    max_log_likelihood_idx = (theta_idx, k_idx, mu_idx)
                log_likelihood_space[theta_idx, k_idx, mu_idx] = log_likelihood
    best_theta, best_k, best_mu = *max_log_likelihood_idx

    return lambda x: gamma_probability(x, theta_hypotheses[best_theta], k_hypotheses[best_k], mu_hypotheses[best_mu])
                





