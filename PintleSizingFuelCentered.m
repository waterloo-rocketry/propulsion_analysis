% This program generates pintle injector parameters for a given engine
% mission 

clear all;

%% Constants - Checked
KELVIN = 273.15;
P_ATM = 14.69594878; %psi

%% Extract NOS Data - Checked
liqdata = readtable('liquidnosdata.xlsx');
vapdata = readtable('vapournosdata.xlsx');

liqdata = table2array(liqdata);
vapdata = table2array(vapdata);

% Usage: 
% value = findit(lookupvalue, lookupindex, outputvalueindex, liqdata);

% Array  Indexing guide
% For Liquid data file
% [01] Temperature (K)	
% [02] Pressure (psia)	
% [03] Density (kg/m3)	
% [04] Volume (m3/kg)	
% [05] Internal Energy (kJ/kg)	
% [06] Enthalpy (kJ/kg)	
% [07] Entropy (J/g*K)	
% [08] Cv (J/g*K)	
% [09] Cp (J/g*K)	
% [10] Sound Spd. (m/s)	
% [11] Joule-Thomson (F/psia)	
% [12] Surface tension (N/m)

% For Vapour data file
% [01] Temperature (K)	
% [02] Pressure (psia)	
% [03] Density (kg/m3)	
% [04] Volume (m3/kg)	
% [05] Internal Energy (kJ/kg)	
% [06] Enthalpy (kJ/kg)	
% [07] Entropy (J/g*K)	
% [08] Cv (J/g*K)	
% [09] Cp (J/g*K)	
% [10] Sound Spd. (m/s)	
% [11] Joule-Thomson (F/psia)

%% Engine Mission - CHECKED
F_thrust = 2500; %N
P_chamber = 600; %psi
OF = 5.482; %optimal at 600 psi is 5.482
Mdot_fu = 0.17028; %kg/s
Mdot_ox = Mdot_fu * OF; %kg/s
Mdot_total = Mdot_fu + Mdot_ox; %kg/s

%% Basic Pintle Geometry - CHECKED
% Using an estimate from "Spray characteristics of a liquid-liquid" pintle
% injector by S.Ninish et al. and then getting the floor that estimate
% Mass flow rate is converted to g/s for the emperical correlation
D_pintle_estimate = sqrt(Mdot_total * 1000); %mm
D_pintle = floor(D_pintle_estimate); %mm

% Chamber diameter for a pintle injector must be between 3 and 5 times the
% pintle diameter. As long as D_cc is between these 2 values, it works
D_cc_minimum = 3 * D_pintle; %mm
D_cc_maximum = 5 * D_pintle; %mm
D_cc = 5/0.0393701; %mm, 5 inches

D_pintle = 34;
D_pintle = D_pintle / 1000; %m

%% Pintle Hole Geometry - CHECKED

% The ratio of the diameter of the second row of holes relative to the
% first row of holes as well as the ratio of the third row of holes
% relative to the first row of holes
secondary_primary_hole_diameter_ratio = 0.9; % unitless
tertiary_primary_hole_diameter_ratio = 0; % unitless
hole_diameter_ratios = [1 secondary_primary_hole_diameter_ratio tertiary_primary_hole_diameter_ratio];

% Number of each size of holes
n_holes_primary = 8; % number of holes
n_holes_secondary = 16; % number of holes
n_holes_tertiary = 0; % number of holes
n_holes = [n_holes_primary n_holes_secondary n_holes_tertiary];

% Initial guess of the discharge coefficient of holes before the program
% generates a more accurate estimate
Cd_fuel = 0.6; % unitless

%% Pintle Annulus Geometry - CHECKED
% length of oxidizer annulus
l = 10/1000; % m

% Initial guess of the discharge coefficient of the annulus before the
% program generates a more accurate estimate
Cd_oxidizer = 0.35; % unitless

%% Propellant Inlet Properties - CHECKED
% Ambient temperature at which this engine operates
T_outside = 25; % Celsius
T_outside_k = KELVIN + T_outside; % Kelvin

% Fuel Tank 
    % Fuel Tank pressure
    P_fuel = 800; %psi
    
    % Fuel Density
    Rho_fuel = 789.00; %kg/m^3
    
    % Pressure drop in the fuel tank
    Pdrop_fuel = P_fuel - P_chamber; %psi
    
    % Viscosity of fuel
    % https://www.engineeringtoolbox.com/ethanol-dynamic-kinematic-viscosity-temperature-pressure-d_2071.html
    mu_fu = 1074/1000000; % Pa*s

% Oxidizer Tank
    % Oxidizer Tank pressure
    P_oxidizer = 800; %psi
    
    % Oxidizer Vapour Pressure
    P_oxidizer_vap = findit(T_outside_k, 1, 2, liqdata);
    
    % Oxidizer Density
    Rho_oxidizer_in = findit(P_oxidizer, 2, 3, liqdata); %kg/m^3
    
    % Pressure drop across the injector
    Pdrop_oxidizer = P_oxidizer - P_chamber; %psi
    
    % Entropy at the inlet of the combustion chamber
    Entropy_liquid_oxidizer_tank = findit(P_oxidizer, 2, 7, liqdata);
    
    % Enthalpy at the inlet of the combustion chamber
    Enthalpy_oxidizer_in = findit(P_oxidizer, 2, 6, liqdata);
    Enthalpy_oxidizer_in = Enthalpy_oxidizer_in * 1000;

%% Propellant Outlet Properties - CHECKED
% Oxidizer Outlet Properties
    % Entropy at the exit of the injector
    Entropy_liquid_oxidizer_chamber = findit(P_chamber, 2, 7, liqdata);
    Entropy_gas_oxidizer_chamber = findit(P_chamber, 2, 7, vapdata);
    
    % Find exit vapour fraction assuming isentropic process
    X_vapour_fraction = (Entropy_liquid_oxidizer_tank-Entropy_liquid_oxidizer_chamber)/(Entropy_gas_oxidizer_chamber - Entropy_liquid_oxidizer_chamber);

    % Find the oxidizer's exit enthalpy
    Enthalpy_oxidizer_out_liquid = findit(P_chamber, 2, 6, liqdata);
    Enthalpy_oxidizer_out_gas = findit(P_chamber, 2, 6, vapdata);
    Enthalpy_oxidizer_out = (1-X_vapour_fraction) * Enthalpy_oxidizer_out_liquid + X_vapour_fraction * Enthalpy_oxidizer_out_gas; % kJ/kg
    Enthalpy_oxidizer_out = Enthalpy_oxidizer_out * 1000; % J/kg

    % Find the oxidizer's exit density
    Rho_oxidizer_out_gas = findit(P_chamber, 2, 3, vapdata); %kg/m^3
    Rho_oxidizer_out_liq = findit(P_chamber, 2, 3, liqdata); %kg/m^3
    Rho_oxidizer_out = (1-X_vapour_fraction) * Rho_oxidizer_out_liq + (X_vapour_fraction) * Rho_oxidizer_out_gas; %kg/m^3

    % Viscosity of Oxidizer at the given temperature
    % Correlation, may not be accurate
    [mu_ox_g] = dynamic_viscosity_nos_g(T_outside_k);
    [mu_ox_l] = dynamic_viscosity_nos_l(T_outside_k);
    mu_ox = (1-X_vapour_fraction) * mu_ox_l + X_vapour_fraction * mu_ox_g;
    
%% Calculate Required Injector Areas and output velocities - CHECKED
    % Annulus, oxidizer stream - CHECKED
    [A_oxidizer] = oxidizer_area(Mdot_ox, Rho_oxidizer_in, Rho_oxidizer_out, P_oxidizer, P_oxidizer_vap, P_chamber, Enthalpy_oxidizer_in, Enthalpy_oxidizer_out);
    [A_oxidizer, Cd_oxidizer] = oxidizer_area_corrected(Mdot_ox, Rho_oxidizer_in, Rho_oxidizer_out, P_oxidizer, P_oxidizer_vap, P_chamber, Enthalpy_oxidizer_in, Enthalpy_oxidizer_out, A_oxidizer, D_pintle, mu_ox, l); 
    
    % Calculate the annulus width of the oxidizer - CHECKED
    A_a = A_oxidizer + pi * (D_pintle/2)^2; 
    D_a = diameter_of_hole(A_a);
    t = 0.5 * (D_a - D_pintle);
    t_mm = t * 1000;
    
    % Calculate the exit velocity of the oxidizer - CHECKED
    Enthalpy_loss = Enthalpy_oxidizer_in - Enthalpy_oxidizer_out; 
    P_drop = P_oxidizer - P_chamber;
    [frac_SPI, frac_HEM] = HEM_fractions(P_oxidizer, P_oxidizer_vap, P_chamber);
    V_oxidizer = hole_mdot(frac_SPI, frac_HEM, Rho_oxidizer_in, Rho_oxidizer_out, Enthalpy_loss, Cd_oxidizer, P_drop)/Rho_oxidizer_out;
    
    % Hole, fuel stream
    % Estimate needed area for fuel stream
    A_fuel = fuel_area(Mdot_fu, Rho_fuel, Pdrop_fuel);
    [A_1, A_2, A_3] = split_area(A_fuel, n_holes, hole_diameter_ratios);
    
    
    m_dot_act_1 = A_1 * Cd_fuel * sqrt(2 * Rho_fuel * (P_fuel-P_chamber) * 6894.76);
	m_dot_act_2 = A_2 * Cd_fuel * sqrt(2 * Rho_fuel * (P_fuel-P_chamber) * 6894.76);
	m_dot_act_3 = A_3 * Cd_fuel * sqrt(2 * Rho_fuel * (P_fuel-P_chamber) * 6894.76);
    
    Cd_1 = Cd_fuel;
    Cd_2 = Cd_fuel;
    Cd_3 = Cd_fuel;
    
    m_dot_act = m_dot_act_1 + m_dot_act_2 + m_dot_act_3;
    
    while m_dot_act < Mdot_fu
        [A_1, A_2, A_3] = split_area(A_fuel, n_holes, hole_diameter_ratios);
        
        V_fuel_1 = m_dot_act_1 / (A_1 * Rho_fuel);
        V_fuel_2 = m_dot_act_2 / (A_2 * Rho_fuel);
        V_fuel_3 = m_dot_act_3 / (A_3 * Rho_fuel);
        
        D_fu = diameter_of_hole(A_fuel);
        D_primary = diameter_of_hole(A_1);
        D_secondary = diameter_of_hole(A_2);
        D_tertiary = diameter_of_hole(A_3);
        
        l = (D_pintle - D_fu)/2;
        
        Re_1 = (Rho_fuel * V_fuel_1 * D_primary)/mu_fu;
        Re_2 = (Rho_fuel * V_fuel_2 * D_secondary)/mu_fu;
        Re_3 = (Rho_fuel * V_fuel_3 * D_tertiary)/mu_fu;
        
        Cd_1 = discharge_coefficient(Re_1, l, D_primary);
        Cd_2 = discharge_coefficient(Re_2, l, D_secondary);
        Cd_3 = discharge_coefficient(Re_3, l, D_tertiary);
        
        m_dot_act_1 = A_1 * Cd_1 * sqrt(2 * Rho_fuel * (P_fuel-P_chamber) * 6894.76);
        m_dot_act_2 = A_2 * Cd_2 * sqrt(2 * Rho_fuel * (P_fuel-P_chamber) * 6894.76);
        m_dot_act_3 = A_3 * Cd_3 * sqrt(2 * Rho_fuel * (P_fuel-P_chamber) * 6894.76);
        
        
        if n_holes(2) == 0
            m_dot_act_2 = 0;
            V_fuel_2 = 0;
        end
        if n_holes(3) == 0
            m_dot_act_3 = 0;
            V_fuel_3 = 0;
        end
        
        m_dot_act = m_dot_act_1 + m_dot_act_2 + m_dot_act_3;
        
        if m_dot_act < Mdot_fu
            A_fuel = 1.00001 * A_fuel;
        end
    end
    
    % Calculate fuel hole diameters in mm
    D_primary_mm = 1000 * D_primary;
    D_secondary_mm = 1000 * D_secondary;
    D_tertiary_mm = 1000 * D_tertiary;
    
    
%% Calculate TMR
V_fuel_average = (n_holes(1) * A_1 * V_fuel_1 + n_holes(2) * A_2 * V_fuel_2 + n_holes(3) * A_3 * V_fuel_3)/A_fuel;

TMR = (Mdot_fu * V_fuel_average)/(Mdot_ox * V_oxidizer);

% Calculate angle of exit based on cheng's liquid-liquid pintle model
theta_cheng = acosd(1/(1+TMR));

% Calculate angle of exit based on heister's gas-liquid pintle model
theta_heister = (180/pi)*TMR^0.5;

% Estimate actual angle (very bad, but idk what else to do)
theta_est = theta_cheng*(1 - X_vapour_fraction) + theta_heister * X_vapour_fraction;


%% Calculate LMR

LMR_primary = (m_dot_act_1 * V_fuel_1)/((D_primary/(D_pintle*pi/n_holes(1))) * Mdot_ox * V_oxidizer);
LMR_secondary = (m_dot_act_2 * V_fuel_2)/((D_secondary/(D_pintle*pi/n_holes(2))) * Mdot_ox * V_oxidizer);
LMR_tertiary = (m_dot_act_3 * V_fuel_3)/((D_tertiary/(D_pintle*pi/n_holes(3))) * Mdot_ox * V_oxidizer);

theta_primary = acosd(1/(1+LMR_primary));
theta_secondary = acosd(1/(1+LMR_secondary));
theta_tertiary = acosd(1/(1+LMR_tertiary));

%% HEM Model, find SPI, HEM fractions - CHECKED
function [frac_SPI, frac_HEM] = HEM_fractions(P, P_vap, P_out)
    k = sqrt((P-P_out)/(P_vap-P_out));
    % Equation from model: K = sqrt((P_supply-P_exit)/(P_vapour-P_exit))
    frac_SPI = k/(k+1);
    frac_HEM = 1/(k+1);
end

%% Calculate Required Area for the Oxidizer - CHECKED
function [A_ox] = oxidizer_area(mdot_ox, Rho_oxidizer_in, Rho_oxidizer_out, P_oxidizer, P_oxidizer_vap, P_chamber, Enthalpy_oxidizer_in, Enthalpy_oxidizer_out)
    % Calculate fractions of single phase incompressible flow and HEM model
    [frac_SPI, frac_HEM] = HEM_fractions(P_oxidizer, P_oxidizer_vap, P_chamber);
    
    m_SPI = frac_SPI * sqrt(2 * Rho_oxidizer_in * (P_oxidizer - P_chamber) * 6894.76);
    m_HEM = frac_HEM * Rho_oxidizer_out * sqrt(2 * (Enthalpy_oxidizer_in - Enthalpy_oxidizer_out));
    
    A_ox = mdot_ox/(m_SPI + m_HEM);
end

%% Calculate the mass flow rate through a single hole given entrance and exit properties - CHECKED
function [m_dot_hole] = hole_mdot(frac_SPI, frac_HEM, Rho_in, Rho_out, Enthalpy_loss, Cd, P_drop)
    
    m_SPI_hole = Cd * sqrt(2 * Rho_in * P_drop * 6894.76);
    m_HEM_hole = Cd * Rho_out * sqrt(2 * Enthalpy_loss);
    
    m_dot_hole = frac_SPI * m_SPI_hole + frac_HEM * m_HEM_hole;
    
    if isnan(m_dot_hole)
        m_dot_hole = 0;
    end
end

%% Split Area Between 3 sets of holes - CHECKED
function [A_primary, A_secondary, A_tertiary] = split_area(A, n_holes, hole_diameter_ratios)
    A_primary = A/(n_holes(3) * hole_diameter_ratios(3)^2 + n_holes(2) * hole_diameter_ratios(2)^2 + n_holes(1));
	A_secondary = A_primary * hole_diameter_ratios(2)^2;
	A_tertiary = A_primary * hole_diameter_ratios(3)^2;
end

%% Correct Area of Oxidizer Annulus to account for frictional effects - CHECKED
function [A_ox_corrected, Cd] = oxidizer_area_corrected(mdot_ox, Rho_oxidizer_in, Rho_oxidizer_out, P_oxidizer, P_oxidizer_vap, P_chamber, Enthalpy_oxidizer_in, Enthalpy_oxidizer_out, A_ox, D_pintle, mu_ox, l)
    
    m_dot_act = 0;
    Enthalpy_loss = Enthalpy_oxidizer_in - Enthalpy_oxidizer_out; 
    P_drop = P_oxidizer - P_chamber;
    [frac_SPI, frac_HEM] = HEM_fractions(P_oxidizer, P_oxidizer_vap, P_chamber);
    A_pintle = pi * (D_pintle/2)^2;
    
    while m_dot_act < mdot_ox
        
        A_a = A_ox + A_pintle; 
        D_a = diameter_of_hole(A_a);
        D_h = D_a - D_pintle; % m, hydraulic diameter of annulus
        
        % Reynolds number through annulus
        V_ox = mdot_ox/Rho_oxidizer_out/A_ox;
        Re = (Rho_oxidizer_in * V_ox * D_h)/mu_ox;
        
        % Discharge coefficients through annulus
        Cd = discharge_coefficient(Re, l, D_h);
        
        % Find actual mdot
        m_dot_act = A_ox * hole_mdot(frac_SPI, frac_HEM, Rho_oxidizer_in, Rho_oxidizer_out, Enthalpy_loss, Cd, P_drop);
        
        if m_dot_act < mdot_ox
            A_ox = 1.0001*A_ox;
        end
    end
    
    A_ox_corrected = A_ox;
end

%% Calculate Hole Diameter from a given hole area - CHECKED
function [D_hole] = diameter_of_hole(A_hole)
    D_hole = sqrt((4 * A_hole)/pi);
end

%% Calculate Required Area for the Fuel - CHECKED

function [A_fu] = fuel_area(mdot_fu, Rho_fuel, Pdrop_fuel)
    A_fu = (mdot_fu/sqrt(2 * Rho_fuel * Pdrop_fuel * 6894.76));
end

%% Lookup function with interpolation
function [ytar] = findit(xtar, nx, ny, data)
    [x1, x2, y1, y2] = lookup(xtar, nx, ny, data);
    [ytar] = linint(xtar, x1, x2, y1, y2);
end

%% Lookup function for the arrays
function [x1, x2, y1, y2] = lookup(xmax, nx, ny, data)
    x1 = 0;
    x2 = 1;
    y1 = 0;
    y2 = 0;
    i = 1;
    while x2 < xmax
        try
            x1 = data(i, nx);
            x2 = data(i+1, nx);
            y1 = data(i, ny);
            y2 = data(i+1, ny);
            i = i + 1;
        catch
            x1 = 0;
            x2 = xmax;
            y1 = 0;
            y2 = 0;
        end
    end
end

%% Function for linear interpolation between 2 points
function [ytar] = linint(xtar, x1, x2, y1, y2)
    m = (y2-y1)/(x2-x1);
    delx = xtar - x1;
    ytar = m * delx + y1;
end

%% Dynamic Viscosity Liquid NOS
function [mu_NOS_l] = dynamic_viscosity_nos_l(T_NOS)
    b1 = 1.6088;
    b2 = 2.0439;
    b3 = 5.24;
    b4 = 0.0293423;
    
    T_c = 309.57;
    
    theta = (T_c - b3)/(T_NOS - b3);
    
    mu_NOS_l = b4*exp(b1*(theta-1)^(1/3) + b2*(theta-1)^(4/3));
    mu_NOS_l = mu_NOS_l / 1000;
end

%% Dynamic Viscosity Gaseous NOS
function [mu_NOS_g] = dynamic_viscosity_nos_g(T_NOS)
    b1 = 3.3281;
    b2 = -1.18237;
    b3 = -0.055155;
    
    T_c = 309.57;
    
    T_r = T_NOS/T_c;
    
    mu_NOS_g = exp(b1 + b2*(1/T_r - 1)^(1/3) + b3*(1/T_r - 1)^(4/3));
    mu_NOS_g = mu_NOS_g / 1000000;
end

%% Discharge Coefficient Calculator - Checked
function [Cd0] = discharge_coefficient(Re, l, d)
    Cf = 0.0791*(Re)^(-0.25);
    if Re < 2300
        Cf = 16/Re;
    end
    K = 2.28; % Assumed minor loss coefficient for sharp edged entrance
    Cd0 = sqrt(1/(4*Cf*(l/d)+K));
end
