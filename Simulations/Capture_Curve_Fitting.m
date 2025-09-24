%%% This script is to reproduce the results from the paper titled 
% "Three-stage approach to optimal low-thrust Earth-moon trajectories."

clc;
clear;
close all;

global m_c;

%% Input
mu_m = 4.9048695*1e12; % in m^3/s^-2
g = 9.80665; % standard gravity acceleration at sea level (in m/s^2)
Rm = 1737.1*1e3; % Radius of moon (in m)
mdot = -107.5/3600; % time rate of change of mass due to thrust (in kg/s)

%Isp = 10046.5287022250; % Constant Specific Impulse value 
%T_max = 2942;
%c = Isp*g; % Exhaust velocity (in m/s)

r_llo = Rm + 100*1e3; % initial orbit radius (in m)
v_theta_llo = sqrt(mu_m/r_llo);

% Define time vector
t_c = 24*60*60; % constant to convert from days to seconds
t0 = 0*t_c; % in seconds

% Initial conditions (from Equations of Motion)
x0 = zeros(4,1);
x0(1) = r_llo;
x0(3) = v_theta_llo;

%% Curve Fitting (for Earth-Escape Maximum Spiral Trajectory
m_llo_vec = linspace(0.86, 0.93, 4)*1e5; % vector of mass left at LLO (in kg) 
rf_vec = linspace(7, 30, 5)*Rm; % Final radius vector vector at start of capture (in m)

%m_llo_vec = 0.86*1e5;
%rf_vec = 7*Rm;

vr_vec = zeros(length(m_llo_vec), length(rf_vec)); % Final radial velocity vector
vtheta_vec = zeros(length(m_llo_vec), length(rf_vec)); % Final tangential velocity vector 
t1_vec = zeros(length(m_llo_vec), length(rf_vec)); % Capture time vector 

itr_max = 1000; % max number of iterations for convergence

for i = 1:length(m_llo_vec)
    m_c = m_llo_vec(i);
    for j = 1:length(rf_vec)
        t1 = 15.6*3600; % intial guess for t1 (in s)
        for k = 1:itr_max
            disp(k)
            %% Solve (Using fmincon)
            % Initial guess for missing initial conditions or parameters
            s0 = [-mu_m/r_llo^2; 0; -v_theta_llo];
            
            % Optimization options for fmincon (SQP)
            options = optimoptions('fmincon', 'Algorithm', 'sqp', 'TolFun', 1e-8);

            % Optimize using fmincon
            s_solution = fmincon(@(s) Capture_Cost_Function(s, x0, t0, t1), s0, [], [], [], [], [], [], @(s) Capture_Constraint_Function(s, x0, t0, t1), options);

            %% Solution
            % Solve ODEs using the correct initial conditions
            options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
            [t, y] = ode45(@(t,y) Capture_ODE(t, y), [t0, t1], initial_conditions(s_solution, x0), options_ode);

            %[t, y] = ode45(@(t,y) VariableThrustODE(t, y), [0, t1], initial_conditions(s_solution, x0), options_ode);

            %% Newton's Method
            rf = y(end, 1);
            rf_dot = y(end, 2);
            vtheta_f = y(end, 3);

            del_t1 = (rf_vec(j) - rf)/rf_dot;
            if abs(del_t1) < 0.02
                t1_vec(i,j) = t1;
                vr_vec(i,j) = rf_dot;
                vtheta_vec(i,j) = vtheta_f;
                break;
            end

            % update step
            t1 = t1-del_t1;
        end

    end
end


%% Curve fit the following functions
% Evaluate at a specific point
rf_query = 680016782.639722;
m_llo_query = 0.94*1e5;

% Create the spline interpolants
[vr_spline, vtheta_spline, tf_spline] = Capture_Splines(rf_vec, m_llo_vec, t1_vec, vr_vec, vtheta_vec);

vr_val = vr_spline(rf_query, m_llo_query);
vtheta_val = vtheta_spline(rf_query, m_llo_query);
tf_val = tf_spline(rf_query, m_llo_query);

% Save the polynomial coefficients to a .mat file
save('Capture_Spline_Function.mat', 'vr_spline', 'vtheta_spline', 'tf_spline');