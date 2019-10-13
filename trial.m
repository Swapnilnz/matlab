clear, clc, close all


% Fixed parameters
rho = 998; % [kg/m^3]
g = 9.81; % [m/s^2)
K = 1.55; % [m^1.5/s]
L = 300; % [m]
A_NI = 620e6; % [km^2]
A_SI = 350e6; % [km^2]
h_NI_min = 355.85; % [m]
h_SI_min = 402; % [m]
h_NI_max = 357.25; % [m]
h_SI_max = 410; % [m]
h_NI_w = h_NI_max; % [m]
h_SI_w = h_SI_max; % [m]
h_NI_gen = 80; % [m]
h_SI_gen = 0; % [m]



% Support functions for differential equations
cf = @(t) 0.41 + 0.12.*sin(2.*pi.*(t - 5660) ./ 8760);

f_NI_in = @(t) 345 + 73.*sin(2.*pi.*(t - 3624) ./ 8760);
f_NI_gen = @(t) f_NI_in(t);


% Power functions, (pf = power function)
pf_SI_demand = 1940; % [MW]
pf_NI_demand = @(t) 4065 + 1.4e6.*normpdf(t, 5000, 1000); % [MW]
pf_geo = 1525; % [MW]
pf_NI_wind = @(t) 4500.*cf(t); % [MW] (Maximum capacity is a guess)
pf_SI_wind = @(t) 4500.*cf(t); % [MW] (Maximum capacity is a guess)
pf_NI_hydro = @(t, h) (0.9.*f_NI_gen(t).*rho.*g.*(h - h_NI_gen))./ 1e6; % [MW]
pf_HVDC = @(t, h) pf_NI_demand(t) - pf_geo - pf_NI_wind(t) - pf_NI_hydro(t, h); % [MW]
pf_SI_hydro = @(t, h) pf_SI_demand + pf_HVDC(t, h) - pf_SI_wind(t); % [MW]

% SI in and gen
f_SI_in = @(t) 593 - 183.*sin(2.*pi.*(t - 2320) ./ 8760); % [m^3/s]
f_SI_gen = @(t, h) (pf_SI_hydro(t, h).*10^6) ./ (0.9.*rho.*g.*(h - h_SI_gen)); % [m^3/s]

% Defiine a set of 4 variables for differential equations

% NI hydroelectric lake differential equation
dh_NI = @(t, h)...
    (f_NI_in(t) - f_NI_gen(t) - f_spill(K, L, h, h_NI_w)) ./ A_NI;

% SI hydroelectric lake differential equation
dh_SI = @(t, h)...
    (f_SI_in(t) - f_SI_gen(t, h) - f_spill(K, L, h, h_SI_w)) ./ A_SI;

% V North Island differential equation
dv_NI = @(h) f_spill(K, L, h, h_NI_w) .* 3600;

% V South  Island differential equation
dv_SI = @(h) f_spill(K, L, h, h_SI_w) .* 3600;  

% Euler's 
t_array = [0:1:8760]; %#ok<*NBRAK> % [h]
h_NI_array = 356.55;
h_SI_array = 406;
v_NI_array = 0;
v_SI_array = 0;


for t = t_array(1:end-1) 
    
    h_NI = h_NI_array(end);
    h_SI = h_SI_array(end);
    h_NI_array(end + 1) = h_NI + dh_NI(t, h_NI); %#ok<*SAGROW>
    h_SI_array(end + 1) = h_SI + dh_SI(t, h_SI);
    
    
    v_NI = v_NI_array(end);
    v_SI = v_SI_array(end);
    v_NI_array(end + 1) = v_NI + dv_NI(h_NI);
    v_SI_array(end + 1) = v_SI + dv_SI(h_SI);
    
end

% Task 1?
disp(h_NI_array)
disp(h_SI_array(end))
disp(v_NI_array(end))
disp(v_SI_array(end))

% Plots
figure(1)
plot(t_array, h_SI_array)
figure(2)
plot(t_array,h_NI_array)
