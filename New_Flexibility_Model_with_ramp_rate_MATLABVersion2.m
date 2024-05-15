%% Energy Arbitrage with Linear Flexibility Model
% Including Ramp Constraint & Deadline Constraint
% Author: Md Umar Hashmi
% Date: 13 May 2024
% Location: Leuven, Belgium

% Clear workspace, close all figures, and clear command window
clear;
close all;
clc;

% Load electricity price and (inelastic load - renewable generation) data
load('Load_price_data.mat'); 

% Initialize parameters
N = length(real);  % Length of horizon in number of samples (unit: kWh)
del_flex_max = 4;  % Maximum charging power (kW)
del_flex_min = 0;  % Minimum discharging power (kW)
h = 0.25;          % Sampling time (hours)
kappa = 1;         % Ratio of selling price to buying price

% Flexibility arrival and departure time
t_a = 24;
t_d = 72;

% Precompute constants
K1 = zeros(N, 1);
K2 = zeros(N, 1);
real_buy = real;
real_sell = kappa * real;
flex_profile_nominal = zeros(N,1);
flex_profile_nominal(t_a:48) = 1;
t=0.25:0.25:24;
K_agg = sum(flex_profile_nominal);
epsilon_k = 1e-3;
Tag_mat = zeros(1, N);
Tag_mat(t_a:t_d) = 1;
gh_mat = zeros(N, 1);
gh_mat(t_a:t_d) = 1;


% ramp rate variables
i_st = 30;
i_end = 1000;
i_max = 1000;

% Construct matrix operators
A_buy = diag(real_buy);
A_sell = diag(real_sell);
A_minus = -eye(N);
A_zero = zeros(N);

K_mat = diag(ones(1, N)) + diag(-ones(1, N-1), -1);
K_mat(1:t_a-1, :) = 0;
K_mat(:, 1:t_a-1) = 0;
K_mat(t_d+1:end, :) = 0;
K_mat(:, t_d+1:end) = 0;

% Initialize results matrix
results = zeros(i_end - i_st + 1, 4);

% Start timer
tic;

% Main loop for energy arbitrage calculation
for id = i_st:i_end
    % Adjust flexibility bounds
    del_max_red = id * del_flex_max / i_max;
    del_min_red = -del_max_red;
    flex_LB = zeros(N, 1);
    flex_LB(t_a:t_d) = del_flex_min * h;
    flex_UB = zeros(N, 1);
    flex_UB(t_a:t_d) = del_flex_max * h;
    
    % Formulate linear programming problem
    A = [A_buy A_minus; A_sell A_minus; Tag_mat zeros(1, N); -Tag_mat zeros(1, N); K_mat A_zero; -K_mat A_zero];
    b = [K1; K2; K_agg + epsilon_k; -K_agg + epsilon_k; del_max_red * h * ones(N, 1); -del_min_red * h * ones(N, 1)];
    lb = [del_flex_min * h * ones(N, 1) .* gh_mat; -1000000 * ones(N, 1)];
    ub = [del_flex_max * h * ones(N, 1) .* gh_mat; 1000000 * ones(N, 1)];
    f = [zeros(N, 1); ones(N, 1)];
    
    % Solve linear programming problem
    x_state = linprog(f, A, b, [], [], lb, ub);
    x_flex = x_state(1:N);
    
    % Calculate cost and profit
    cost_of_consumption_nominal = sum(real' * max(0, flex_profile_nominal) - kappa * real' * max(0, -flex_profile_nominal));
    x_ch = max(0, x_flex);
    x_ds = -min(0, x_flex);
    flex_out = x_ch - x_ds;
    profit_only_arbitrage = cost_of_consumption_nominal - (sum(real' * max(0, flex_out) - kappa * real' * max(0, -flex_out)));
    
    % Record results
    results(id - i_st + 1, :) = [id, del_max_red, profit_only_arbitrage, max(abs(diff(x_flex)))];
end

% End timer and display execution time
toc;

% Display results (optional)
disp('Results:');
disp(results);
%%
figure
plot(results(:,1)/1000, results(:,3))

figure
plot(results(:,1)/1000, results(:,3)*100/results(end,3))

figure
plot(results(1:end-1,1)/1000, diff(results(:,3)*100/results(end,3)))

figure
plot(t,flex_profile_nominal)
hold on
plot(t,x_flex)