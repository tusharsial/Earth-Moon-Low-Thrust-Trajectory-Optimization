%% This script is for computing the trajectory for the coasting phase

clc;
clear;
close all;

% Load the polynomial coefficients from the .mat file
load('Escape_Curve_Fitting_Function.mat');
load('Capture_Spline_Function.mat');

%% Input
mu_m = 4.9048695*1e12; % in m^3/s^-2
mu_e = 3.986004418*1e14; % in m^3/s^-2
Re = earthRadius; % Radius of earth (in m)
Rm = 1737.1*1e3; % Radius of moon (in m)
D = 3.84399*1e8; % Distance between Earth and Moon (in m)
mu_e = 3.986004418*1e14; % in m^3/s^-2
mu_m = 4.9048695*1e12; % in m^3/s^-2
mu = mu_e + mu_m; % in m^3/s^-2
omg = sqrt(mu/D^3); % Constant angular rate of Earth-Moon system
m_leo = 1e5; % Initial mass at start of spiral
mdot = 107.5/3600; % time rate of change of mass due to thrust (in kg/s)
t_c = 24*60*60; % constant to convert from days to seconds
t0 = 0*t_c; % in seconds

%% Optimisation
% escape and capture curve fitting function
e_cf = [p_vr; p_vtheta; p_tf];

% Initial guess
%s0 = [9.69*Re; 50*pi/180; 0.93*m_leo; 3*t_c];
%s0 = [9.47*Re; 136*pi/180; 0.93*m_leo; 4.5*t_c];
%s0 = [11.29*Re; 10*pi/180; 0.87*m_leo; 1*t_c];
s0 = [8.29*Re; 105*pi/180; 0.96*m_leo; 5.2*t_c];


%lb = [6*Re, 0, 0.06*m_leo, 0.001*t_c]; % Lower bounds
%ub = [20*Re, 2*pi, 0.96*m_leo, 10*t_c];

% Optimization options for fmincon (SQP)
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'MaxFunctionEvaluations', 10000, 'StepTolerance', 1e-6,'OptimalityTolerance', 1e-6,'ConstraintTolerance', 1e-6);

% Optimize using fmincon
s_solution = fmincon(@(s) Coasting_Cost_Function(s, e_cf, vr_spline, vtheta_spline, tf_spline), s0, [], [], [], [], [], [], @(s) Coasting_Constraint_Function(s, e_cf, vr_spline, vtheta_spline, tf_spline, m_leo, mdot), options);

% options = optimoptions('patternsearch', ...
%     'Display', 'iter', ...          % Show iteration details
%     'MaxFunctionEvaluations', 10000, ... % Max function evaluations
%     'PlotFcn', {@psplotbestf, @psplotfuncount}); % Optional: Plot progress
% 
% s_solution = patternsearch(@(s) Coasting_Cost_Function(s, e_cf, vr_spline, vtheta_spline, tf_spline), s0, [], [], [], [], [], [], @(s) Coasting_Constraint_Function(s, e_cf, vr_spline, vtheta_spline, tf_spline, m_leo, mdot), options);

%% Simulation of Optimal Coast trajectory
r0 = s_solution(1); % intial radius at start of coast
theta0 = s_solution(2); % initial angle at start of coast
m_llo = s_solution(3); % Final mass in LLO (in kg)
t1 = s_solution(4); % duration of coast

vr_0 = polyval(p_vr, r0);
vtheta_0 = polyval(p_vtheta, r0) - omg*r0;

x0 = [r0; vr_0; vtheta_0; theta0];

% Solve ODEs using the guessed parameters
options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
[~, y] = ode45(@(t,y) Coasting_ODE(t, y), [0, t1], x0, options_ode);

% Extract final values at t1
y_final = y(end, :);

rf = y_final(1);
vr_f = y_final(2);
vtheta_f = y_final(3);
theta_f = y_final(4);

[r2_f, theta2_f, vr2_f, vtheta2_f] = ECRF_to_MCRF(rf, theta_f, vr_f, vtheta_f, D, omg);

%% Optimal Solution Analysis
t_esc = polyval(p_tf, r0); 
t_capt = tf_spline(r2_f, m_llo);

%% Plot Trajectory
r = y(:,1);
theta = y(:,4);

% Position in Cartesian Coordinates
x_vec = r.*cos(theta);
y_vec = r.*sin(theta);

figure;
plot(x_vec/Re, y_vec/Re)

axis('equal')
xlabel(' X (normalised with Earth radii)')
ylabel(' Y (normalised with Earth radii)')
title('Coasting Trajectory')
