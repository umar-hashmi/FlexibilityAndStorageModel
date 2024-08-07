# Walking through code
## Initialization and Setup:

The script starts with clearing the workspace, closing all figures, and clearing the command window to ensure a fresh environment.
Data is loaded from a .mat file, and relevant parameters are initialized.

## Precomputations:

Several constants and matrices are precomputed to avoid redundancy within the loop.

## Matrix Construction:

Matrices for the linear programming problem (A_buy, A_sell, A_minus, A_zero, and K_mat) are constructed and modified according to the constraints.

## Main Loop:

For each iteration, the flexibility bounds are adjusted, and the linear programming problem is formulated and solved using linprog.
The cost and profit are calculated based on the results of the optimization.

## Results Collection and Display:

Results are collected into an array and displayed after the loop execution.



# The proposed battery and flexibility model
Energy storage and flexibility model with ramp rate constraint: we improve the flexibility model by considering the ramp rate constraint. The flexibility model also includes a desired time window of operation and a deadline constraint ensures that the energy consumed is
temporally optimized with this time window. The battery models also consider ramp and capacity constraints. The proposed model is linear and can be efficiently solved using off-the-shelf solvers. The run-time for 1000 Monte Carlo (MC) simulation days with 96 time steps for each MC
scenario takes less than 30 seconds. Due to its computational efficiency, the proposed models are suitable for (near) real-time operation

# GAMS Implementation
To accommodate linear programming constraints in GAMS, two pseudo variables are introduced - x_ch and x_dis - to represent the positive and negative values of the charging variable x. 

This results in comparable results between the Matlab, Python, and GAMS code. One thing to note is that in Matlab in the initial time step, the ramp rate constraint has been relaxed so that any level of charging within the charging limits is accepted. In the GAMs version, the first time step has been constrained by the ramp rate assuming starting at a 0 charging rate. Thus the results slightly differ.

# Citation
Paper title: Linear storage and flexibility model with ramp rate, ramping, deadline and capacity constraints

Conference: IEEE PES ISGT Europe 2024 (ISGT Europe 2024)
