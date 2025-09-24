%% This script is for implementing the hybrid direct/indirect optimizer for 
% generating optimal low-thrust earth-moon trajectory. 

clc;
clear;
close all;

%% Input
Re = earthRadius; % Radius of earth (in m)
Rm = 1737.1*1e3; % Radius of moon (in m)
soi_m = 66100*1e3; % Moon's SOI radius (km)
D = 3.84399*1e8; % Distance between Earth and Moon (in m)
mu_e = 3.986004418*1e14; % in m^3/s^-2
mu_m = 4.9048695*1e12; % in m^3/s^-2
mu = mu_e + mu_m; % in m^3/s^-2
omg = sqrt(mu/D^3); % Constant angular rate of Earth-Moon system
m_leo = 1e5; % Initial mass at start of spiral
mdot = 107.5/3600; % time rate of change of mass due to thrust (in kg/s)
t_c = 24*60*60; % constant to convert from days to seconds
t0 = 0*t_c; % in seconds

% LEO Data
r_leo = Re + 315*1e3; % initial orbit radius (in m)
v_theta_leo = sqrt(mu_e/r_leo);

% LLO Data
r_llo = Rm + 100*1e3; % initial orbit radius (in m)
v_theta_llo =sqrt(mu_m/r_llo); 

%% Optimisation
x0_1 = [r_leo; 0; v_theta_leo - omg*r_leo]; % intial conditions for escape + partial coast phase
x0_2 = [r_llo; 0; v_theta_llo - omg*r_llo]; % intial conditions for escape + partial coast phase

theta_0 = 0.*pi/180;
theta_f = 0.2*pi/180;

% theta_0 = 2.2532;
% theta_f = -5.3835;

% theta_0 = -2.2532;
% theta_f = 0.1;

lam_0 = [-mu_e/r_leo^2; 0; -v_theta_leo];
lam_f = [-mu_m/r_llo^2; 0; -v_theta_llo];

t_esc0 = 197293.978080361;
t_capt0 = 49789.9037231618;
t_coast0 = 218388.470468064;

% t_esc0 = 192915.837887292;
% t_capt0 = 38657.0176538703;
% t_coast0 = 409697.479401245;

% t_esc0 = 193347.111945247;
% t_capt0 = 39753.4248502906;
% t_coast0 = 509369.174009751;

s0 = [theta_0; theta_f; lam_0; lam_f; t_esc0; t_capt0; t_coast0]; % parameters to optimize: theta_0, theta_f, lam_r0, lam_vr0, lam_vtheta0, lam_r_tf, lam_vr_tf, lam_vtheta_tf, t_esc, t_capt, t_coast 

tmax = 3.5*t_c; % Initial maximum value for computing t_coast2

% Optimization options for fmincon (SQP) %'MaxFunctionEvaluations', 10000
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'StepTolerance', 1e-6,'OptimalityTolerance', 1e-6,'ConstraintTolerance', 1e-6);

% Optimize using fmincon
% options = optimoptions('patternsearch', ...
%      'Display', 'iter', ...          % Show iteration details
%      'MaxFunctionEvaluations', 10000, ... % Max function evaluations
%      'PlotFcn', {@psplotbestf, @psplotfuncount}); % Optional: Plot progress

tic;
s_opt = fmincon(@(s) Hybrid_Cost_Function(s), s0, [], [], [], [], [], [], @(s) Hybrid_Constraint_Function(s, x0_1, x0_2, t0, m_leo, mdot, soi_m, tmax), options);
toc;
%s_opt = patternsearch(@(s) Hybrid_Cost_Function(s), s0, [], [], [], [], [], [], @(s) Hybrid_Constraint_Function(s, x0_1, x0_2, t0, m_leo, mdot, soi_m, tmax), options);

%% Simulations
theta0 = s_opt(1);
thetaf = s_opt(2);
lam_0 = [s_opt(3); s_opt(4); s_opt(5)];
t_esc = s_opt(9);
t_capt = s_opt(10);
t_coast = s_opt(11);
tf = t_esc+t_coast+t_capt;

%% ODE 2: Capture + partial coast (t = 0 (LLO) till t = t_capt + t_coast2) (Backward prograpagate)
m_llo = m_leo - mdot*(t_esc+t_coast+t_capt);
y0_2 = [x0_2; thetaf; lam_f; 0]; % intial 8 by 1 state

tspan2 = t0:100:tmax;

options_ode1 = odeset('Events', @soi_event, 'RelTol', 1e-6, 'AbsTol', 1e-8); % Event function for detecting crossing of SOI
[t2, y2] = ode45(@(t,y) Hybrid_ODE_2(t, y, t_capt, m_llo), tspan2, y0_2, options_ode1); 

t_coast2 = t2(end) - t_capt;

%% ODE 1: Escape + partial coast (t = 0 (LEO) till t = t_esc + t_coast1) (Forward propagate)
y0_1 = [x0_1; theta0; lam_0; 0]; % intial 8 by 1 state
t_coast1 = t_coast - t_coast2;
t1 = t_esc + t_coast1;

tspan1 = t0:100:t1;

%options_ode2 = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
[t1_vec, y1] = ode45(@(t,y) Hybrid_ODE_1(t, y, t_esc), tspan1, y0_1); 

%% Plot Trajectory and Optimal Thrust angle Profile 
r1 = y1(:,1);
vr1 = y1(:,2);

vtheta1 = y1(:,3);
theta1 = y1(:,4);

r2 = y2(:,1);
vr2 = y2(:,2);
vtheta2 = y2(:,3);
theta2 = y2(:,4);

%[rm_1, thetam_1, vrm_1, vthetam_1] = ECRF_to_MCRF(r1, theta1, v_r1, v_theta1, D, omg); 

% Position in Cartesian Coordinates
x_vec1 = r1.*cos(theta1);
y_vec1 = r1.*sin(theta1);

x_vec2 = r2.*cos(theta2) - D;
y_vec2 = r2.*sin(theta2);

%% for sphere of influence
% Create angle values from 0 to 2*pi
theta_soi = linspace(0, 2*pi, 100);

% Compute x and y coordinates of the circle
x_soi = soi_m*cos(theta_soi)-D;
y_soi = soi_m*sin(theta_soi);

figure;
plot(x_vec1(t1_vec<=t_esc)/Re,y_vec1(t1_vec<=t_esc)/Re);
hold on;
plot(x_vec1(t1_vec(1:end-200)>t_esc)/Re,y_vec1(t1_vec(1:end-200)>t_esc)/Re, 'k--');
plot(x_vec2(t2>t_capt)/Re,y_vec2(t2>t_capt)/Re, 'k--');
plot(x_vec2(t2<=t_capt)/Re,y_vec2(t2<=t_capt)/Re);
plot(x_soi/Re, y_soi/Re, '--')
axis('equal')
xlabel(' X (normalised with Earth radii)')
ylabel(' Y (normalised with Earth radii)')
legend('Thrust arc', 'Coast arc')
title('Optimal retrograde trajectory: rotating coordinates')
xlim([-75, 10])

%% Flight path angle and Thrust (Escape)
lam_vr1 = y1(:,6);
lam_vtheta2 = y1(:,7);

thrust_idx1 = t1_vec <= t_esc;

u1_esc = atan2(-lam_vr1(thrust_idx1),-lam_vtheta2(thrust_idx1)); % Optimal control 
u_coast1 = zeros(sum(~thrust_idx1), 1);
u1 = [u1_esc; u_coast1];

gam1 = atan2(y1(:,2), y1(:,3)); % Flight path angle

figure;
plot(t1_vec/t_c, u1*180/pi)
xlabel('time (in days)')
ylabel('Steering and Flight path angles (in deg.)')
hold on
plot(t1_vec/t_c, gam1*180/pi, '--')
legend('Optimal Thrust Angle', 'Flight Path Angle')
title('Optimal thrust direction and flight path angles: Earth Escape Phase.')
xlim([0, 3])

%% Flight path angle and thrust (Capture)
lam_vr2 = y2(:,6);
lam_vtheta2 = y2(:,7);

thrust_idx2 = t2 <= t_capt;

u1_capt = atan2(-lam_vr2(thrust_idx2),-lam_vtheta2(thrust_idx2))+pi; % Optimal control 
u_coast2 = zeros(sum(~thrust_idx2), 1);
u2 = [u1_capt; u_coast2];

gam2 = atan2(y2(:,2), y2(:,3)); % Flight path angle

figure;
plot((tf-t2)/t_c +2, u2*180/pi)
xlabel('time (in days)')
ylabel('Steering and Flight path angles (in deg.)')
hold on
plot((tf-t2)/t_c +2, gam2*180/pi, '--')
legend('Optimal Thrust Angle', 'Flight Path Angle')
title('Optimal thrust direction and flight path angles: Earth Escape Phase.')
%xlim([0, 3])