## Ornstein Uhlerbeck
# dX(t) = θ * [μ – X(t)] * dt + σ * dW(t),
# where dX(t) is an increment of the process X between t and dt, and σ > 0 is the instantaneous diffusion term, used to measure the volatility of the process, which is assumed
# to be constant. On the other side, μ is the process long-term expected value, and θ > 0 is
# the speed or reversion of X(t) toward μ, both also assumed to be constant. Finally, dW(t)
# is an increment during the interval (t,t + dt) of a standard Brownian motion under the
# real probability measure P, which follows a normal distribution with expected value 0 and
# variance t.
import numpy as np
#import matplotlib.pyplot as plt
import random


def ou_proc(dt, tau, T, sigma, mu):
    np.random.seed(10)
    n = int(T / dt)  # Number of time steps.
    t = np.linspace(0., T, n)  # Vector of times.
    
    sigma_bis = sigma * np.sqrt(2. / tau)
    sqrtdt = np.sqrt(dt)
    x = np.zeros(n)
    for i in range(n - 1):
        x[i + 1] = x[i] + dt * (-(x[i] - mu) / tau) + sigma_bis * sqrtdt * np.random.normal(0,1)#,1) #mu should be positive, and different than DC_bias
    # fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    # plt.plot(t, x, lw=2)
    # plt.xlim(0,500)
    # plt.title('The OU Process with tau {} ms'.format(tau))
    # plt.ylabel('OU-value (mV)',fontsize=15)
    # plt.xlabel('time (ms)',fontsize=15)
    return x