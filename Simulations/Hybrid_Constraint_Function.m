%% Constraint Function (Ensures Boundary Conditions are Satisfied)
function [c, ceq] = Hybrid_Constraint_Function(s, x0_1, x0_2, t0, m_leo, mdot, soi_m, tmax)

% Transform from Earth-Centered to Moon-centered
D = 3.84399*1e8; % Distance between Earth and Moon (in m)
mu_e = 3.986004418*1e14; % in m^3/s^-2
mu_m = 4.9048695*1e12; % in m^3/s^-2
mu = mu_e + mu_m; % in m^3/s^-2
omg = sqrt(mu/D^3); % Constant angular rate of Earth-Moon system

%% Extract parameters
theta0 = s(1);
thetaf = s(2);

lam_0 = [s(3); s(4); s(5)];
lam_f = [s(6); s(7); s(8)];

t_esc = s(9);
t_capt = s(10);
t_coast = s(11);

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
[~, y1] = ode45(@(t,y) Hybrid_ODE_1(t, y, t_esc), tspan1, y0_1); 

%% Extract;
y_final1 = y1(end, :);
r1 = y_final1(1);
v_r1 = y_final1(2);
v_theta1 = y_final1(3); 
theta1 = y_final1(4);

y_final2 = y2(end, :);
r2 = y_final2(1);
theta2 = y_final2(4);
v_r2 = y_final2(2);
v_theta2 = y_final2(3);

% Convert from Earth-centered rotating to Moon-Centered Rotating
[rm_1, thetam_1, vrm_1, vthetam_1] = ECRF_to_MCRF(r1, theta1, v_r1, v_theta1, D, omg); 

%% Constraint
c = []; % No inequality constraints
ceq = [rm_1 - r2; mod(thetam_1, 2*pi) - mod(theta2, 2*pi); vrm_1 - v_r2; vthetam_1 - v_theta2]; % 4 Equality Constraints