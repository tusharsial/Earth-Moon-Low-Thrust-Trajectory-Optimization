function J = Coasting_Cost_Function(s, e_cf, vr_spline, vtheta_spline, tf_spline)

%% Constants
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
vtheta_0 = polyval(p_vtheta, r0) - omg*r0; % velocity in earth-centered rotating frame

x0 = [r0; vr_0; vtheta_0; theta0];

% Solve ODEs using the guessed parameters
options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); 
[~, y] = ode45(@(t,y) Coasting_ODE(t, y), [0, t1], x0, options_ode); 

% Extract final values at t1
y_final = y(end, :);
rf = y_final(1);
theta_f = y_final(4);

r2_f = sqrt(rf^2 + D^2 - 2*rf*D*cos(theta_f)); % radius at tf in moon-centered frame

% Compute escape and capture times
t_esc = polyval(p_tf, r0);
t_capt = tf_spline(r2_f, m_llo);

% Cost function
J = t_esc + t_capt;
