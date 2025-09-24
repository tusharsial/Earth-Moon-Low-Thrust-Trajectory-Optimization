%%% This script is to reproduce the results from the paper titled 
% "Three-stage approach to optimal low-thrust Earth-moon trajectories."

clc;
clear;
close all;

%% Input
mu_e = 3.986004418*1e14; % in m^3/s^-2
g = 9.80665; % standard gravity acceleration at sea level (in m/s^2)
Re = earthRadius; % Radius of earth (in m)

m_leo = 1e5; % initial mass at LEO (in kg) 
r_leo = Re + 315*1e3; % initial orbit radius (in m)
v_theta_leo = sqrt(mu_e/r_leo);

% Define time vector
t_c = 24*60*60; % constant to convert from days to seconds
t0 = 0*t_c; % in seconds
t1_vec = linspace(2.1, 2.37, 6)*t_c; % Final escape time vector in seconds
dt = 100; % in seconds

% Initial conditions (from Equations of Motion)
x0 = zeros(4,1);

x0(1) = r_leo;
x0(3) = v_theta_leo;

%% Curve Fitting (for Earth-Escape Maximum Spiral Trajectory
vr_vec = zeros(1, length(t1_vec)); % Final radial velocity vector
vtheta_vec = zeros(1, length(t1_vec)); % Final tangential velocity vector 
rf_vec = zeros(1, length(t1_vec)); % Final radius vector

for i = 1:length(t1_vec)
    t1 = t1_vec(i);
    %% Solve (Using fmincon)
    % Initial guess for missing initial conditions or parameters
    s0 = [-mu_e/r_leo^2; 0; -v_theta_leo];
    
    % Optimization options for fmincon (SQP)
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'TolFun', 1e-8);

    % Optimize using fmincon
    s_solution = fmincon(@(s) Escape_Cost_Function(s, x0, t0, t1), s0, [], [], [], [], [], [], @(s) Escape_Constraint_Function(s, x0, t0, t1), options);

    %% Solution
    % Solve ODEs using the correct initial conditions
    options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
    [t, y] = ode45(@(t,y) Escape_ODE(t, y), [0, t1], initial_conditions(s_solution, x0), options_ode);
    
    %% Update vectors (for curve fitting)
    rf_vec(i) = y(end, 1);
    vr_vec(i) = y(end, 2);
    vtheta_vec(i) = y(end, 3);
end

%% Curve fit the following functions
% Choose a degree n
n = 5;

% Fit polynomials
p_vr = polyfit(rf_vec, vr_vec, n);       % vr = f1(rf)
p_vtheta = polyfit(rf_vec, vtheta_vec, n); % vtheta = f2(rf)
p_tf = polyfit(rf_vec, t1_vec, n);       % tf = f3(rf)

% Evaluate at a new rf (e.g., rf_new)
rf_new = 65562185.8252543;  % example value
vr_new = polyval(p_vr, rf_new);
vtheta_new = polyval(p_vtheta, rf_new);
tf_new = polyval(p_tf, rf_new);

% Save the polynomial coefficients to a .mat file
save('Escape_Curve_Fitting_Function.mat', 'p_vr', 'p_vtheta', 'p_tf');