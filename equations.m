%(1)
dh = (f_in - f_out) ./ A;

%(2)
p_HVDC = p_NI_demand - p_geo - p_NI_wind - p_NI_hydro;

%(3)
p_SI_hydro = p_SI_demand + p_HVDC - p_SI_wind;

%(4)
cf = 0.41 + 0.12.*sin(2.*pi.*(t - 5660) ./ 8760);

%(5)
f_NI_in = 345 + 73.*sin(2.*pi.*(t - 3624) ./ 8760);

%(6) 
f_SI_in = 593 - 183.*sin(2.*pi.*(t - 2320) ./ 8760);

%(7)
p_hydro = (0.9.*f_gen.*rho.*g.*(h - h_gen))./1e6;

%(8)
f_NI_gen = f_NI_in;

%(9)
f_spill = K.*L.*((h - h_w).^1.5);

%(10)
dv_NI_spill = f_NI_spill .* 3600;

%(11)
dv_SI_spill = f_SI_spill .* 3600;

%(12)
dh_NI = (f_NI_in - f_NI_gen - f_NI_spill)./ A_NI;

%(13)
dh_SI = (f_SI_in - f_SI_gen - f_SI_spill)./ A_SI;

%(14)
p_NI_demand = 4065 + 1.4e6.*normpdf(t, 5000, 1000);

%(15)
p_SI_demand = 1940;

%(16)
value = (0.9.*rho.*v_spill.*g.*(h - h_gen))./(3600e6) .* 100;






