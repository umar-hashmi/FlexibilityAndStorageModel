# -*- coding: utf-8 -*-
"""
Created on Wed May 15 14:49:16 2024

@author: uhashmi
"""

import numpy as np
from scipy.optimize import linprog
import time

# Load data from a .mat file
from scipy.io import loadmat
data = loadmat('Load_price_data.mat')
real = data['real'].flatten()
load = data['load'].flatten()

N = len(real)

# Battery parameters
e_ch = 0.95
e_dis = 0.95
del_max = 500
del_min = -del_max
b_0 = 200
b_max = 1000
b_min = 200
h = 0.25
kappa = 1

start_time = time.time()

t_a = 24
t_d = 72

K1 = np.zeros(N)
K2 = np.zeros(N)

real_buy = real / e_ch
real_sell = kappa * real * e_dis
real_interme = real * e_dis
real_interme_negload = kappa * real / e_ch

A_buy = np.diag(real_buy)
A_sell = np.diag(real_sell)
A_interm = np.diag(real_interme)
A_intermnegload = np.diag(real_interme_negload)
A_minus = -1 * np.eye(N)
A_zero = np.zeros((N, N))

K_mat = np.diag(np.ones(N))
K_mat += np.diag(-np.ones(N-1), -1)

Scenario_max = 1000

results = []

for id in range(1, Scenario_max + 1):
    del_max_red = id * del_max / Scenario_max
    del_min_red = -del_max_red

    # LP matrix formulation
    A = np.vstack([
        np.hstack([A_buy, A_minus]),
        np.hstack([A_interm, A_minus]),
        np.hstack([A_intermnegload, A_minus]),
        np.hstack([A_sell, A_minus]),
        np.hstack([np.tril(np.ones((N, N)), -1) + np.eye(N), A_zero]),
        np.hstack([-np.tril(np.ones((N, N)), -1) - np.eye(N), A_zero]),
        np.hstack([K_mat, A_zero]),
        np.hstack([-K_mat, A_zero])
    ])

    b = np.hstack([
        K1, K1, K2, K2,
        (b_max - b_0) * np.ones(N),
        (b_0 - b_min) * np.ones(N),
        del_max * h,
        del_max_red * h * np.ones(N-1),
        -del_min * h,
        -del_min_red * h * np.ones(N-1)
    ])

    lb = np.hstack([del_min * h * np.ones(N), -1000000 * np.ones(N)])
    ub = np.hstack([del_max * h * np.ones(N), 1000000 * np.ones(N)])

    Aeq = None
    beq = None
    f = np.hstack([np.zeros(N), np.ones(N)])

    result = linprog(f, A_ub=A, b_ub=b, A_eq=Aeq, b_eq=beq, bounds=list(zip(lb, ub)), method='highs')
    x_state = result.x

    x = x_state[:N]

    cost_of_consumption_nominal = (np.sum(real * np.maximum(0, load)) - kappa * np.sum(real * np.maximum(0, -load))) / 1000

    x_ch = np.maximum(0, x)
    x_ds = -np.minimum(0, x)
    lhouse = load + x_ch / e_ch - x_ds * e_dis

    bat_out = x_ch / e_ch - x_ds * e_dis

    profit_only_arbitrage = (np.sum(real * np.maximum(0, bat_out)) - kappa * np.sum(real * np.maximum(0, -bat_out))) / 1000

    x_adj = x / b_max

    results.append([id, del_max_red, -profit_only_arbitrage, np.max(np.abs(np.diff(x_adj)))])

end_time = time.time()
print(f"Execution time: {end_time - start_time} seconds")

# Convert results to a numpy array for further analysis if needed
results = np.array(results)
