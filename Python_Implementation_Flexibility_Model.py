# -*- coding: utf-8 -*-
"""
Created on Wed May 15 15:02:21 2024

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

# Flexibility parameters
del_flex_max = 4
del_flex_min = 0
h = 0.25
kappa = 1

start_time = time.time()

t_a = 24
t_d = 72

K1 = np.zeros(N)
K2 = np.zeros(N)

real_buy = real
real_sell = kappa * real

A_buy = np.diag(real_buy)
A_sell = np.diag(real_sell)

A_minus = -1 * np.eye(N)
A_zero = np.zeros((N, N))

K_mat = np.diag(np.ones(N))
K_mat += np.diag(-np.ones(N-1), -1)

K_mat[:t_a-1, :] = 0
K_mat[:, :t_a-1] = 0
K_mat[t_d:, :] = 0
K_mat[:, t_d:] = 0

Tag_mat = np.zeros(N)
Tag_mat[t_a:t_d] = 1

flex_profile_nominal = np.zeros(N)
flex_profile_nominal[t_a:48] = 1
t = np.arange(0.25, 24.25, 0.25)

gh_mat = np.zeros(N)
gh_mat[t_a:t_d] = 1

K_agg = np.sum(flex_profile_nominal)
epsilon_k = 1e-3

i_st = 30
i_end = 1000
i_max = 1000

results = []

for id in range(i_st, i_end + 1):
    flex_LB = np.zeros(N)
    flex_LB[t_a:t_d] = del_flex_min * h
    flex_UB = np.zeros(N)
    flex_UB[t_a:t_d] = del_flex_max * h

    del_max_red = id * del_flex_max / i_max
    del_min_red = -del_max_red

    A = np.vstack([
        np.hstack([A_buy, A_minus]),
        np.hstack([A_sell, A_minus]),
        np.hstack([Tag_mat, np.zeros(N)]),
        np.hstack([-Tag_mat, np.zeros(N)]),
        np.hstack([K_mat, A_zero]),
        np.hstack([-K_mat, A_zero])
    ])

    b = np.hstack([
        K1, K2,
        K_agg + epsilon_k, -K_agg + epsilon_k,
        del_max_red * h * np.ones(N), -del_min_red * h * np.ones(N)
    ])

    lb = np.hstack([del_flex_min * h * gh_mat, -1000000 * np.ones(N)])
    ub = np.hstack([del_flex_max * h * gh_mat, 1000000 * np.ones(N)])

    Aeq = None
    beq = None
    f = np.hstack([np.zeros(N), np.ones(N)])

    result = linprog(f, A_ub=A, b_ub=b, A_eq=Aeq, b_eq=beq, bounds=list(zip(lb, ub)), method='highs')
    x_state = result.x

    x_flex = x_state[:N]

    cost_of_consumption_nominal = np.sum(real * np.maximum(0, flex_profile_nominal)) - kappa * np.sum(real * np.maximum(0, -flex_profile_nominal))

    x_ch = np.maximum(0, x_flex)
    x_ds = -np.minimum(0, x_flex)

    flex_out = x_ch - x_ds

    profit_only_arbitrage = cost_of_consumption_nominal - (np.sum(real * np.maximum(0, flex_out)) - kappa * np.sum(real * np.maximum(0, -flex_out)))

    results.append([id, del_max_red, profit_only_arbitrage, np.max(np.abs(np.diff(x_flex)))])

end_time = time.time()
print(f"Execution time: {end_time - start_time} seconds")

# Convert results to a numpy array for further analysis if needed
results = np.array(results)
