%% Energy Arbitrage with linear storage model including the ramp constraint
% Code Author: Md Umar Hashmi
% Date: 13 May 2024
% Location: Leuven, Belgium


clear
close all
clc

load('Load_price_data.mat');        % Load electricity price and (inelastic load - renewable generation) data
N=length(real);                     % Length of horizon in number of samples
% Unit of electricity price cents per kWh
% Unit of load is kWh

% battery parameters
e_ch=0.95;                          % Charging efficiency
e_dis =0.95;                        % Discharging efficiency
del_max = 500;                     % Maximum charging rate

del_min = -del_max;                 % Minimum discharging rate

b_0 = 200;                         % Initial battery capacity
b_max = 1000;                       % Maximum battery capacity
b_min = 200;                        % Minimum permissible battery capacity
h=0.25;                             % Sampling time

kappa = 1;                          % the ratio of selling price and buying price

tic


% K1 = -real.*load;
% K2 = -kappa*real.*load;

t_a = 24;
t_d = 72;

K1 = zeros(N,1);
K2 = zeros(N,1);

real_buy = real/e_ch;
real_sell = kappa*real*e_dis;
real_interme = real*e_dis;
real_interme_negload = kappa*real/e_ch;

A_buy = diag(real_buy);
A_sell = diag(real_sell);
A_interm = diag(real_interme);
A_intermnegload = diag(real_interme_negload);
A_minus = -1*eye(N);
A_zero = zeros(N);

K_mat = diag(ones(1, N)); % Create a diagonal matrix with 1s on the diagonal
K_mat = K_mat + diag(-ones(1, N-1), -1); % Add lower diagonal elements as -1

Scenario_max = 1000;

for id = 1:Scenario_max
    
    del_max_red = id*del_max/Scenario_max;
    del_min_red = -del_max_red;
    
    % LP matrix formulation
    A = [A_buy A_minus; A_interm A_minus; A_intermnegload A_minus; A_sell A_minus; tril(ones(N,N),-1) + eye(N)  A_zero; -tril(ones(N,N),-1) - eye(N)  A_zero;...
        K_mat A_zero; -K_mat A_zero];
    b=[K1; K1; K2; K2;(b_max-b_0)*ones(N,1); (b_0-b_min)*ones(N,1); del_max*h; del_max_red*h*ones(N-1,1); -del_min*h;-del_min_red*h*ones(N-1,1)];
    % b=[K1; K1; K2; K2;(b_max-b_0)*ones(N,1); (b_0-b_min)*ones(N,1); del_max_red*h*ones(N,1); -del_min_red*h*ones(N,1)];
    
    % A = [A_buy A_minus; A_interm A_minus; A_intermnegload A_minus; A_sell A_minus; tril(ones(N,N),-1) + eye(N)  A_zero; -tril(ones(N,N),-1) - eye(N)  A_zero];
    % b=[K1; K1; K2; K2;(b_max-b_0)*ones(N,1); (b_0-b_min)*ones(N,1)];
    
    lb=[del_min*h*ones(N,1); -1000000*ones(N,1)];
    ub=[del_max*h*ones(N,1); 1000000*ones(N,1)];
    
    Aeq=[];
    beq=[];
    f=[zeros(N,1); ones(N,1)];
    
    x_state = linprog(f,A,b,Aeq,beq,lb,ub);
    
    x= x_state(1:N);
    
    cost_of_consumption_nominal = sum(real'*subplus(load)-kappa*real'*subplus(-load))/1000;
    
    x_ch = max(0,x);
    x_ds = -min(0,x);
    lhouse = load+x_ch/e_ch - x_ds*e_dis;
    
    bat_out = x_ch/e_ch - x_ds*e_dis;
    
    profit_only_arbitrage =  sum(real'*subplus(bat_out)-kappa*real'*subplus(-bat_out))/1000;
    toc
    
    x_adj= x/b_max;
    
    res(id,:) = [id  del_max_red  -profit_only_arbitrage   max(abs(diff(x_adj)))];
    
end

t = 0.25:0.25:24;
figure
plot(t,x_adj)


%%
des = 0.1:0.1:1;
idx = [1:100:Scenario_max];
figure
% hold on
plot(res(:,1)/Scenario_max, res(:,3)*100/res(end,3))
hold on
plot(res(idx,1)/Scenario_max, res(idx,3)*100/res(end,3), '*r')
hold on
xline(res(idx,1)/Scenario_max,'--','color', [.5 .5 .5])
hold on
yline(res(idx,3)*100/res(end,3),'--','color', [.5 .5 .5])

%%
figure
plot(res(1:end-1,1)/Scenario_max, diff(res(:,3)*100/res(end,3)))

%%
figure
plot(res(:,1)/Scenario_max, res(:,4))