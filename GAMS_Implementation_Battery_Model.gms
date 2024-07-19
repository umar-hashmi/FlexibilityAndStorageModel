* Energy Arbitrage with linear storage model including the ramp constraint
* Code Author: Aleen van der Meer and Md Umar Hashmi
* Date: 1 July 2024
* Location: Leuven, Belgium

* Clear previous data
$clear


*Define scalars
Scalar
    e_ch    charging efficiency /0.95/
    e_dis   discharging efficiency /0.95/
    
    del_max maximum charging rate /500.0/
    del_min minimum discharging rate /-500.0/
    
    ramp_rate  maximum ramp rate as percent of battery charging rate
    
    b_0     initial battery capacity /200.0/
    b_max   maximum battery capacity /1000.0/
    b_min   minimum permissible battery capacity /200.0/
    h       sampling time   /0.25/
    
    kappa   ratio of selling price and buying price /1.0/   ;


*Define sets and parameters
sets
    i time intervals
    / 1*96 /
    ifirst(i) first period
    ilast(i) last period

    m data_types
    / real, load /
    
    id set to drive iterations
    / 1*1000 /;
    
ifirst(i) = yes$(ord(i) eq 1);
ilast(i) = yes$(ord(i) eq card(i) );


* Load price and load data from .csv file
Table data_input(i,m) 'price and load data' 
$ondelim
$include load_price_data.csv
$offdelim
; 

* Create parameters for the price and load data
parameter real(i), load(i);
real(i) = data_input(i, "real");
load(i) = data_input(i, "load");

* Create parameter for iterations and number of scenarios 
parameter
    profit(id)  used to hold the profit of each iteration
    val(id)     numeric value of the iterations
    scen_max    maximum number of scenarios ;

val(id) = ord(id);
scen_max = card(id);

* Define variables
positive variables x_ch(i);
negative variables x_dis(i);
variables x(i), b(i), profit_only_arbitrage;

* Define equations
equations obj,
    batt_soc_constraint(i),
    ramp_rate_constraint_up(i),
    ramp_rate_constraint_down(i),
    sum_constraint(i);

obj.. profit_only_arbitrage =e= (sum(i, real(i)*x_ch(i)/e_ch) + sum(i, real(i)*x_dis(i)*e_dis))/1000;

batt_soc_constraint(i).. b(i) =e= b(i-1) + b_0$ifirst(i) + x(i);

ramp_rate_constraint_up(i).. x(i) - x(i-1) =l= ramp_rate*h*del_max;

ramp_rate_constraint_down(i).. x(i-1) - x(i) =l= ramp_rate*h*del_max;

sum_constraint(i).. x_ch(i) + x_dis(i) =g= x(i);

* Bounds
x.up(i) = del_max * h;
x.lo(i) = del_min * h;

b.up(i) = b_max;
b.lo(i) = b_min;


* Solve
model arbitrage /all/;
loop(id,
    ramp_rate = val(id)/scen_max;
    
    solve arbitrage using lp minimizing profit_only_arbitrage;
    
    profit(id) = profit_only_arbitrage.l
) ;

* Results
display profit