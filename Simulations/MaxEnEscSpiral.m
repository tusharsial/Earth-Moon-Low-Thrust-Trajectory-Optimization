%%% This script is to reproduce the results from the paper titled 
% "Three-stage approach to optimal low-thrust Earth-moon trajectories."

clc;
clear; 
close all;

%% Input
mu_e = 3.986004418*1e14; % in m^3/s^-2
%Re = 6371*1e3; % in m
Re = earthRadius;
T = 2942; % Constant Thrust magnitude (in N)
m_leo = 1e5; % initial mass at LEO (in kg) 
r_leo = Re + 315*1e3; % initial orbit radius (in m)
v_theta_leo = sqrt(mu_e/r_leo);
mdot = 107.5/3600; % time rate of change of mass due to thrust (in kg/s)

% Define time vector
c = 24*60*60; % constant to convert from days to seconds
t0 = 0*c; % in seconds
t1 = 2.37*c; % in seconds
dt = 10; % in seconds
t = t0:dt:t1;

% Compute Thrust acceleration magnitude
m = m_leo - mdot*t;
a_T = T./m;   

% Initial conditions (from Equations of Motion)
x0 = zeros(4,1);
x0(1) = r_leo;
x0(3) = v_theta_leo;

%% Solve (Using fmincon)
tic;
% Initial guess for missing initial conditions or parameters
s0 = [-mu_e/r_leo^2; 0; -v_theta_leo];

% Optimization options for fmincon (SQP)
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'TolFun', 1e-10);

% Optimize using fmincon
s_solution = fmincon(@(s) cost_function(s, x0, t0, t1), s0, [], [], [], [], [], [], @(s) constraint_function(s, x0, t0, t1), options);

%% Solution
% Solve ODEs using the correct initial conditions
options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
[t, y] = ode45(@(t,y) myODEs(t, y), [0, t1], initial_conditions(s_solution, x0), options_ode);
toc;

%% Plot Trajectory and Optimal Thrust angle Profile 
r = y(:,1);
theta = y(:,4);

% Position in Cartesian Coordinates
x_vec = r.*cos(theta);
y_vec = r.*sin(theta);

% Optimal Control

u = atan2(-y(:,6),-y(:,7)); % Optimal control 
gam = atan2(y(:,2), y(:,3)); % Flight path angle

figure(1);
plot(x_vec/Re,y_vec/Re)
axis('equal')
xlabel(' X (normalised with Earth radii)')
ylabel(' Y (normalised with Earth radii)')
title('Maximum-energy Earth-escape spiral for t_f = 2.37 days')

figure(2);
plot(t/c, u)
xlabel('time (in days)')
ylabel('Steering and Flight path angles (in deg.)')
hold on
plot(t/c, gam, '--')
legend('Optimal Thrust Angle', 'Flight Path Angle')
title('Optimal thrust direction and flight path angles: Earth Escape Phase.')

