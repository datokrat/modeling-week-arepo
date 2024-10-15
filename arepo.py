import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import umbridge as um
from matplotlib.animation import FuncAnimation
from Riemann import RiemannProblem

from scipy.integrate import solve_ivp

import tinyDA as tda

class UmbridgeArepoModel(um.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [4]
    
    def get_output_sizes(self, config):
        return [4]
    
    def __call__(self, parameters, config={}):
        RESOLUTION = 200
        input = parameters[0]
        rho_L, v_L, p_L = input[0]
        rho_R, v_R, p_R = input[1]
        T_max, n_T = input[2]
        L_max, n_L = input[3]

        time = np.linspace(0, T_max, n_T + 1)
        solutions = [
            RiemannProblem(np.linspace(0, L_max, RESOLUTION), L_max / 2.0, [rho_L, v_L, p_L], [rho_R, v_R, p_R], 1.4, t)[1]
            for t in time
        ]

        def convert_space_index(input):
            return int(input * (RESOLUTION - 1) / n_L)

        return [
            [
                [], # x coordinates
                [
                    [
                        solutions[t][convert_space_index(s), 0] # 0 = density
                        for s in range(n_L + 1)
                    ]
                    for t in range(len(time))
                ],
                [
                    [
                        solutions[t][convert_space_index(s), 1] # 1 = velocity
                        for s in range(n_L + 1)
                    ]
                    for t in range(len(time))
                ],
                [
                    [
                        solutions[t][convert_space_index(s), 2] # 2 = pressure
                        for s in range(n_L + 1)
                    ]
                    for t in range(len(time))
                ]
            ]
        ]

    def supports_evaluate(self):
        return True
    
class NpArepoModel(um.Model):
    def __init__(self, inner_model, n_L=5, n_T=20, L_max=20, T_max=20, verbose=False):
        self.inner = inner_model
        self.n_L = n_L
        self.n_T = n_T
        self.L_max = L_max
        self.T_max = T_max
        self.verbose = verbose

    def __call__(self, parameters):
        left = [0, 0, 0]
        right = [0, 0, 0]
        left[0] = np.clip(parameters[0], 1e-6, 1)
        left[1] = 0.0
        left[2] = np.clip(parameters[1], 1e-6, 1)
        right[0] = np.clip(parameters[2], 1e-6, 1)
        right[1] = 0.0
        right[2] = np.clip(parameters[3], 1e-6, 1)
        inner_args = [[left, right, [self.T_max, self.n_T], [self.L_max, self.n_L]]]
        if self.verbose:
            print("Inner input:")
            print(inner_args)
        inner_result = np.array(self.inner(inner_args)[0][3])
        if self.verbose:
            print("Inner output:")
            print(inner_result)
        if inner_result.shape[0] != self.n_T + 1:
            raise ValueError(f"Expected n_T + 1 = {self.n_T + 1} time steps, but got {inner_result.shape[0]}")
        if inner_result.shape[1] != self.n_L + 1:
            raise ValueError(f"Expected n_L + 1 = {self.n_L + 1} spatial steps, but got {inner_result.shape[1]}")
        return np.array(inner_result).flatten(), True
    
class MultivariatePrior:
    def __init__(self, low, high):
        self.low = np.array(low)
        self.high = np.array(high)
        self.dim = len(low)

    def logpdf(self, x):
        x = np.array(x)
        in_range = np.all((x >= self.low) & (x <= self.high))
        return 0 if in_range else -np.inf  # np.where(in_range, 0, -np.inf)

    def rvs(self, size=1):
        if isinstance(size, int):
            size = (size, self.dim)
        elif isinstance(size, tuple) and len(size) == 1:
            size = size + (self.dim,)
        samples = np.random.uniform(self.low, self.high, size=size)
        if size == 1 or (isinstance(size, tuple) and size[0] == 1):
            return samples[0]
        return samples
    

iterations = 12000
burnin = 2000

# Set up the prior
low = np.array([0, 0, 0, 0])
high = np.array([1, 1, 1, 1])
my_prior = MultivariatePrior(low, high)

# draw some samples and plot them
prior_samples = my_prior.rvs(size=(1000, 4))

def inverse_uq(sigma, L_max, n_L, T_max, n_T, rho_L, p_L, rho_R, p_R):
    true_parameters = np.array([rho_L, p_L, rho_R, p_R])

    my_exact_model = NpArepoModel(UmbridgeArepoModel(), n_L=n_L, n_T=n_T, L_max=L_max, T_max=T_max)

    # initalise the true model and solve it
    # http_model = um.HTTPModel("http://localhost:4242", "shockwave1D")
    # my_model = NpArepoModel(http_model, n_L=n_L, n_T=n_T, L_max=L_max, T_max=T_max, verbose=False)
    my_model = NpArepoModel(UmbridgeArepoModel(), n_L=n_L, n_T=n_T, L_max=L_max, T_max=T_max, verbose=False)
    y_true = my_model(true_parameters)[0]
    y_true

    # set the noise level
    experimental_measurements = np.array(my_exact_model(true_parameters)[0])

    noise = np.random.normal(scale=sigma, size=experimental_measurements.size) # fine noise
    data = experimental_measurements + noise # noisy fine data.
    data[data < 0] = 0 # make sure all the data is positive.

    # define the likelihood
    cov_likelihood = sigma**2*np.eye(data.size)

    my_loglike = tda.GaussianLogLike(data, cov_likelihood)

    # set up the link factories
    my_posterior = tda.Posterior(my_prior, my_loglike, my_model)

    my_posteriors = [my_posterior]

    MAP = tda.get_MAP(my_posterior)

    am_cov = np.eye(true_parameters.size)
    am_t0 = 100
    am_sd = None
    am_epsilon = 1e-6
    am_adaptive = True
    my_proposal = tda.AdaptiveMetropolis(C0=am_cov, t0=am_t0, sd=am_sd, epsilon=am_epsilon)

    my_chain = tda.sample(my_posteriors, my_proposal, iterations=iterations, n_chains=2, initial_parameters=MAP, subchain_length=5, adaptive_error_model='state-independent', force_sequential=True)

    return tda.to_inference_data(my_chain, level='fine', burnin=burnin)