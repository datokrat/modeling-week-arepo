from Riemann import RiemannProblem
import numpy as np

def sample_riemann(xx, x0, W_L, W_R, gamma, time):
    return RiemannProblem(xx, x0, W_L, W_R, gamma, time)


