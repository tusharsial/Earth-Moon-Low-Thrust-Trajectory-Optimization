%% This script is for computing the thrust angle to maximize the spacecraft 
% energy at the end of start of a fixed-time capture spiral in LLO. 

clc;
clear;
close all;

global m_c;

%% Input
mu_m = 4.9048695*1e12; % in m^3/s^-2
g = 9.80665; % standard gravity acceleration at sea level (in m/s^2)
Rm = 1737.1*1e3; % Radius of moon (in m)
mdot = -107.5/3600; % time rate of change of mass due to thrust (in kg/s)
%m_llo = 0.89*1e5; % mass left at LLO (in kg) 
m_llo = 92990.1164250887; % mass left at LLO (in kg) 

%Isp = 10046.5287022250; % Constant Specific Impulse value 
%T_max = 2942;
%c = Isp*g; % Exhaust velocity (in m/s)

r_llo = Rm + 100*1e3; % initial orbit radius (in m)
v_theta_llo = -sqrt(mu_m/r_llo);

% Define time vector
t_c = 24*60*60; % constant to convert from days to seconds
t0 = 0*t_c; % in seconds
%t1 = 15.6*60*60; % Final escape time vector in seconds
t1 = 40502.1515870877;

m_c = m_llo; % mass left at LL0 (in kg) 

% Initial conditions (from Equations of Motion)
x0 = zeros(4,1);
%x0 = zeros(5,1);

x0(1) = r_llo;
x0(3) = v_theta_llo;
%x0(5) = m_leo;

%% Solve (Using fmincon)
% Initial guess for missing initial conditions or parameters
s0 = [-mu_m/r_llo^2; 0; -v_theta_llo];
%s0 = [-mu_e/r_leo^2; 10; -v_theta_leo; (T_max/m_leo)*v_theta_leo];

% Optimization options for fmincon (SQP)
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'TolFun', 1e-8);

% Optimize using fmincon
s_solution = fmincon(@(s) Capture_Cost_Function(s, x0, t0, t1), s0, [], [], [], [], [], [], @(s) Capture_Constraint_Function(s, x0, t0, t1), options);

%% Solution
% Solve ODEs using the correct initial conditions
options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
[t, y] = ode45(@(t,y) Capture_ODE(t, y), [t0, t1], initial_conditions(s_solution, x0), options_ode);

%[t, y] = ode45(@(t,y) VariableThrustODE(t, y), [0, t1], initial_conditions(s_solution, x0), options_ode);

%% Data for testing curve fits
rf = y(end, 1);
vr_f = y(end, 2);
vtheta_f = y(end, 3);
mf = m_llo - mdot*t1;
J = (vr_f^2 + vtheta_f^2)/2 - mu_m/rf;
%% Plot Trajectory and Optimal Thrust angle Profile 
r = y(:,1);
theta = y(:,4);

% Position in Cartesian Coordinates
x_vec = r.*cos(theta);
y_vec = r.*sin(theta);

% Optimal Control
u = atan2(-y(:,6),-y(:,7))+pi; % Optimal control 
%u = atan2(-y(:,7),-y(:,8)); % Optimal control 
gam = atan2(y(:,2), y(:,3)); % Flight path angle

figure;
plot(x_vec/Rm,y_vec/Rm)

axis('equal')
xlabel(' X (normalised with Earth radii)')
ylabel(' Y (normalised with Earth radii)')
title('Maximum-energy Earth-escape spiral for t_f = 2.37 days')

figure;
plot(t/t_c, u*180/pi)
xlabel('time (in days)')
ylabel('Steering and Flight path angles (in deg.)')
hold on
plot(t/t_c, gam*180/pi, '--')
legend('Optimal Thrust Angle', 'Flight Path Angle')
title('Optimal thrust direction and flight path angles: Earth Escape Phase.')

