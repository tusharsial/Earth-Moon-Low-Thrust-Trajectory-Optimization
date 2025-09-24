%% Constraint Function (Ensures Boundary Conditions are Satisfied)
function [c, ceq] = Coasting_Constraint_Function(s, e_cf, vr_spline, vtheta_spline, tf_spline, m_leo, mdot)

% Transform from Earth-Centered to Moon-centered
D = 3.84399*1e8; % Distance between Earth and Moon (in m)
mu_e = 3.986004418*1e14; % in m^3/s^-2
mu_m = 4.9048695*1e12; % in m^3/s^-2
mu = mu_e + mu_m; % in m^3/s^-2
omg = sqrt(mu/D^3); % Constant angular rate of Earth-Moon system

%% Escape curve fitting functions
p_vr = e_cf(1, :);
p_vtheta = e_cf(2, :);
p_tf = e_cf(3, :);

%% Initial value for the EoM 
r0 = s(1); % intial radius at start of coast
theta0 = s(2); % initial angle at start of coast
m_llo = s(3); % Final mass in LLO (in kg)
t1 = s(4); % duration of coast

vr_0 = polyval(p_vr, r0);
vtheta_0 = polyval(p_vtheta, r0) - omg*r0;

x0 = [r0; vr_0; vtheta_0; theta0];

%% Solve ODEs using the guessed parameters
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
[~, y] = ode45(@(t,y) Coasting_ODE(t, y), [0, t1], x0, options); 

% Extract final values at t1
y_final = y(end, :);
rf = y_final(1);
vr_f = y_final(2);
vtheta_f = y_final(3);
theta_f = y_final(4);

[r2_f, theta2_f, vr2_f, vtheta2_f] = ECRF_to_MCRF(rf, theta_f, vr_f, vtheta_f, D, omg);
% r2_f = sqrt(rf^2 + D^2 - 2*rf*D*cos(theta_f)); % radius at tf in moon-centered frame
% vr2_f = vr_f * cos(theta_f) + vtheta_f * sin(theta_f);
% vtheta2_f = -vr_f * sin(theta_f) + vtheta_f * cos(theta_f) - omg*r2_f;

% No inequality constraints (c = [])
c = [-r0; -t1];

% Compute escape and capture times
t_esc = polyval(p_tf, r0);
t_capt = tf_spline(r2_f, m_llo);

% Equality constraints (terminal conditions)
ceq = [vr2_f - vr_spline(r2_f, m_llo); abs(vtheta2_f) - abs(vtheta_spline(r2_f, m_llo) - omg*r2_f); m_llo - m_leo + mdot*(t_esc + t_capt)];
