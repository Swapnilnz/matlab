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

% Initialisation of Euler's method arrays
t0 = 0;
tf = 8760;
n = 8760;
t_array = linspace(t0, tf, n+1); % [h]
h_NI_array = 356.55;
h_SI_array = 406;
v_NI_array = 0;
v_SI_array = 0;


for t = t_array(1:end-1) 
    
    
    % Initialise key variables for Euler's method
    h_NI = h_NI_array(end);
    h_SI = h_SI_array(end);
    v_NI = v_NI_array(end);
    v_SI = v_SI_array(end);
    
    % 1. Calculate the North Island and South Island demand in MW (14, 15).
    p_NI_demand = 4065 + 1.4e6.*normpdf(t, 5000, 1000);
    p_SI_demand = 1940;
    
    % 2. Calculate the inlet flow rate into the North Island and SI Lakes
    % (5,6)
    f_NI_in = 345 + 73.*sin((2.*pi.*(t - 3624)) ./ 8760);
    f_SI_in = 593 - 183.*sin((2.*pi.*(t - 2320)) ./ 8760);
    
    % 3. Set the generating flow for the North Island (8).
    f_NI_gen = f_NI_in;
    
    % 4. Calculate the amount of hydro-electric generation in the NI (7).
    p_NI_hydro = (0.9.*f_NI_gen.*rho.*g.*(h_NI - h_NI_gen))./ (10^6);
    
    % 5. Calculate the capacity factor for wind (same in both islands)
    % (4). Also calculate NI and SI wind generation in MW (guessed).
    cf = 0.41 + 0.12.*sin(2.*pi.*(t - 5660) ./ 8760);
    p_NI_wind = 2000.*cf;
    p_SI_wind = 2130.*cf;
    
    % 6. Calculate the total of hydro, geothermal and wind generation in
    % the North Island in MW.
    p_geo = 1525;
    p_NI_total = p_NI_hydro + p_NI_wind + p_geo;
    
    % 7. Calculate power needed to be transmitted from or to the SI.
    p_HVDC = p_NI_demand - p_NI_total;
    
    % 8. From the South Island demand and wind production, determine the
    % amount of hydro power required (3).
    p_SI_hydro = p_SI_demand + p_HVDC - p_SI_wind;
    
    % 9. Determine the generating flow rate of hydro in the South Island
    % (7)
    f_SI_gen = (p_SI_hydro.*10^6) ./ (0.9.*rho.*g.*(h_SI - h_SI_gen));
    
    % 10. Determine if there is any spillway flow for each lake (9).
    if h_NI > h_NI_w
        NI_spill = K.*L.*((h_NI - h_NI_w).^1.5);
    else
        NI_spill = 0;
    end
    if h_SI > h_SI_w
        SI_spill = K.*L.*((h_SI - h_SI_w).^1.5);
    else
        SI_spill = 0;
    end
    
    % 11. Determine the four derivaties (10-13)
    dv_NI_spill = NI_spill .* 3600;
    dv_SI_spill = SI_spill .* 3600;
    dh_NI = (f_NI_in - f_NI_gen - NI_spill)./ A_NI;
    dh_SI = (f_SI_in - f_SI_gen - SI_spill)./ A_SI;
    
    % 12. Use Euler’s method to take a step in time (17).
    h_NI_array(end + 1) = h_NI + dh_NI; %#ok<*SAGROW>
    h_SI_array(end + 1) = h_SI + dh_SI;
    
    v_NI_array(end + 1) = v_NI + dv_NI_spill;
    v_SI_array(end + 1) = v_SI + dv_SI_spill;
    
    
end

% Task 1 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% Plots (b)
figure(1)
plot(t_array, h_SI_array)
title('South Island lake level vs time (2035)')
xlabel('Time in year 2035 (hours)')
ylabel('Height of South Island lake (m)')

figure(2)
plot(t_array, h_NI_array)
title('North Island lake level vs time (2035)')
xlabel('Time in year 2035 (hours)')
ylabel('Height of North Island lake (m)')

% Minimum lake levels
disp(min(h_SI_array))
disp(min(h_NI_array))


